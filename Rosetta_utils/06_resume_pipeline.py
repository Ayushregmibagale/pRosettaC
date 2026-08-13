#!/usr/bin/env python3
"""
Resume pRosettaC pipeline from the constraint_generation step.

Assumes local docking is already done (pd.*_docking_????.pdb exist in Patchdock_Results/).
Steps performed:
  1. Score existing docking PDBs → local.fasc  (skipped if already complete)
  2. Clean up failed constraint_generation outputs (empty/stale files)
  3. Run constraint_generation for each docking solution
  4. Run clustering

Usage:
  cd <run_dir_containing_Protac_params.txt>
  python /path/to/resume_pipeline.py
"""
import sys
import os
import glob
import subprocess
import concurrent.futures

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS_DIR)

# Apply all patches before doing anything else
from pyrosetta_adapter import apply_patches
apply_patches()

sys.path.insert(0, str(os.path.join(os.path.expanduser('~'), 'PRosettaC')))
import utils

RUN_DIR = os.getcwd()
PATCHDOCK = os.path.join(RUN_DIR, 'Patchdock_Results')


def read_params(params_file='Protac_params.txt'):
    params = {}
    with open(params_file) as f:
        for line in f:
            if ':' in line:
                k, v = line.split(':', 1)
                params[k.strip()] = v.strip()
    return params


def score_docking_pdbs(pt0_params, pt1_params):
    """Score all docking PDBs and write local.fasc; skip if already complete."""
    fasc = os.path.join(PATCHDOCK, 'local.fasc')
    pdbs = sorted(glob.glob(os.path.join(PATCHDOCK, 'pd.*_docking_????.pdb')))
    if os.path.exists(fasc):
        with open(fasc) as f:
            n = sum(1 for line in f if line.startswith('SCORE:'))
        if n >= len(pdbs) * 0.9:
            print(f'local.fasc already has {n}/{len(pdbs)} entries — skipping scoring.')
            return
    print(f'Scoring {len(pdbs)} docking PDBs for local.fasc...', flush=True)
    # Rewrite paths through ~/protac symlink so PyRosetta init (which splits on
    # whitespace) doesn't break on the space in "PROTAC degradation Prediction".
    def _sl(p):
        return p.replace(str(RUN_DIR).split('pRosettaC')[0].rstrip('/'),
                         os.path.expanduser('~/protac'))
    import pyrosetta
    pyrosetta.init(
        f'-extra_res_fa {_sl(pt0_params)} -extra_res_fa {_sl(pt1_params)} '
        '-ignore_unrecognized_res true -load_PDB_components false -out:level 100'
    )
    scorefxn = pyrosetta.create_score_function('ref2015')
    with open(fasc, 'w') as lf:
        for k, pdb in enumerate(pdbs):
            try:
                pose = pyrosetta.pose_from_file(pdb)
                score = scorefxn(pose)
            except Exception as e:
                print(f'  SKIP {os.path.basename(pdb)}: {e}')
                score = 0.0
            desc = os.path.basename(pdb).replace('.pdb', '')
            lf.write(f'SCORE:  {score:.3f}  0  0  0  {score:.3f}  {desc}\n')
            if (k + 1) % 200 == 0:
                print(f'  scored {k+1}/{len(pdbs)}', flush=True)
    print(f'local.fasc written ({len(pdbs)} entries)')


def cleanup_failed_constraint_gen():
    """Remove stale files from the previous failed constraint_generation run."""
    removed = 0
    stale_patterns = ['confs_*.sdf', 'v_*.sdf', 'docked_*.sdf',
                      'PT_*.pdb', 'PT_*.params']
    # Only delete combined PDBs and score.sc if docking poses exist to regenerate them;
    # if docking PDBs are gone (deleted to save space), combined PDBs + score.sc are
    # the only results and must be preserved for clustering.
    if glob.glob(os.path.join(PATCHDOCK, 'pd.*_docking_????.pdb')):
        stale_patterns += ['combined_*.pdb', 'score.sc']
    for pattern in stale_patterns:
        for f in glob.glob(os.path.join(PATCHDOCK, pattern)):
            os.remove(f)
            removed += 1
    # Keep docked_*.pdb (they have HETATM content that is reused)
    for f in glob.glob(os.path.join(PATCHDOCK, 'docked_*.pdb')):
        if os.path.getsize(f) == 0:
            os.remove(f)
            removed += 1
    print(f'Cleaned up {removed} stale files.')


def run_constraint_generation(heads, linkers, chains):
    """Run constraint_generation wrapper for each docking solution in parallel."""
    os.chdir(PATCHDOCK)
    wrapper = os.path.join(SCRIPTS_DIR, 'constraint_generation.py')
    docking_solutions = sorted(glob.glob('*_docking_????.pdb'))
    suffix_list = []
    for s in docking_solutions:
        parts = s.split('.')[1].split('_')
        sfx = parts[0] + '_' + str(int(parts[2]))
        suffix_list.append((s, sfx))

    chains_str = ''.join(chains)
    # Use absolute paths directly; the original "../" prefix breaks when heads/linkers
    # are absolute paths (f"../{abs_path}" → "..//home/..." = wrong relative path).
    def _abspath(p):
        return p if os.path.isabs(p) else os.path.abspath(os.path.join(RUN_DIR, p))

    h0 = _abspath(heads[0])
    h1 = _abspath(heads[1])
    lk = _abspath(linkers)
    commands = [
        f'python "{wrapper}" "{h0}" "{h1}" "{lk}" '
        f'{sfx} {pdb} {chains_str}'
        for pdb, sfx in suffix_list
    ]

    n_total = len(commands)
    print(f'Running constraint_generation for {n_total} poses '
          f'(~{n_total // max(1, os.cpu_count()-2)} min expected)...', flush=True)

    # Conservative default: 4 workers leaves ~92 of 96 cores free for other users
    # on the shared ernie node. Override with N_WORKERS_RESUME env var if you have
    # the node to yourself (e.g. N_WORKERS_RESUME=16 for a dedicated run).
    # Previous default of cpu_count()-2 (~94) caused load avg 550 (see Bug 11).
    n_workers = int(os.environ.get('N_WORKERS_RESUME', 4))

    def _run(cmd):
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            # Only print stderr excerpt if GenConstConf didn't just produce no conformers
            if 'No conformations were generated' not in r.stderr and r.returncode != 0:
                print(f'[CG] fail rc={r.returncode}: {cmd[:70]}')
        return r.returncode

    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_run, c): c for c in commands}
        for f in concurrent.futures.as_completed(futs):
            done += 1
            if done % 200 == 0:
                print(f'  {done}/{n_total} done', flush=True)

    os.chdir(RUN_DIR)

    n_combined = len(glob.glob(os.path.join(PATCHDOCK, 'combined_*_0001.pdb')))
    print(f'Constraint generation complete: {n_combined} combined PDBs produced.')


def shift_scores_if_needed():
    """
    Shift score.sc values to be negative if needed.
    clustering.py rejects all entries where score >= 0 (hardcoded `< 0` filter).
    PackRotamersMover produces positive Rosetta scores, so we shift:
      shifted = original - (max_score + 100)
    Backs up the original as score.sc.orig and records the shift value in
    score.sc.shift for downstream analysis scripts to recover original values.
    """
    sc_file = os.path.join(PATCHDOCK, 'score.sc')
    if not os.path.exists(sc_file):
        print('No score.sc found — skipping score shift check.')
        return

    with open(sc_file) as f:
        lines = f.readlines()

    score_lines = [l for l in lines if l.startswith('SCORE:') and len(l.split()) > 2]
    if not score_lines:
        print('score.sc is empty — skipping score shift.')
        return

    scores = []
    for l in score_lines:
        parts = l.split()
        try:
            scores.append(float(parts[1]))
        except (IndexError, ValueError):
            continue

    max_score = max(scores) if scores else 0.0
    if max_score < 0:
        print(f'score.sc already negative (max={max_score:.1f}) — no shift needed.')
        return

    shift = max_score + 100.0
    print(f'Shifting score.sc by {shift:.3f} (max={max_score:.1f} → all scores < 0)')

    backup = sc_file + '.orig'
    if not os.path.exists(backup):
        import shutil
        shutil.copy2(sc_file, backup)

    new_lines = []
    for l in lines:
        if l.startswith('SCORE:') and len(l.split()) > 2:
            parts = l.split()
            try:
                parts[1] = f'{float(parts[1]) - shift:.3f}'
                new_lines.append('  '.join(parts) + '\n')
            except (IndexError, ValueError):
                new_lines.append(l)
        else:
            new_lines.append(l)

    with open(sc_file, 'w') as f:
        f.writelines(new_lines)

    with open(sc_file + '.shift', 'w') as f:
        f.write(f'{shift:.6f}\n')

    print(f'score.sc shifted; original backed up to score.sc.orig')


def run_clustering(chains):
    """Run clustering (matches main.py call)."""
    os.system('cat Init0.pdb Init1.pdb > Init.pdb')
    sys.path.insert(0, str(os.path.join(os.path.expanduser('~'), 'PRosettaC')))
    import importlib
    cl = importlib.import_module('clustering')
    cl.main('clustering.py', ['1000', '200', '4', chains[1]])
    if os.path.isdir('Results/'):
        print('SUCCESS — Results/ directory created.')
    else:
        print('Clustering done but no Results/ — check score.sc/local.fasc/combined_*.pdb')


def main():
    params = read_params()
    raw_heads = params.get('Heads', '').split()
    # main.py adds _H.sdf suffix during preprocessing
    heads = [h.rsplit('.', 1)[0] + '_H.sdf' for h in raw_heads]
    linkers = params.get('Protac', '')
    chains = params.get('Chains', 'ABC A').split()

    pt0 = os.path.join(RUN_DIR, 'PT0.params')
    pt1 = os.path.join(RUN_DIR, 'PT1.params')

    print(f'[resume] RUN_DIR: {RUN_DIR}')
    print(f'[resume] Heads: {heads}')
    print(f'[resume] Linker: {linkers}')
    print(f'[resume] Chains: {chains}')

    score_docking_pdbs(pt0, pt1)
    cleanup_failed_constraint_gen()
    run_constraint_generation(heads, linkers, chains)
    shift_scores_if_needed()
    run_clustering(chains)


if __name__ == '__main__':
    main()
