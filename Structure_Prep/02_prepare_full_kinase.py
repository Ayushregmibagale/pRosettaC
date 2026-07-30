#!/usr/bin/env python3
"""
Prepare full-kinase-domain structures for RIPK2 and MAPK14 for pRosettaC re-run.

Problem: pRosettaC concatenates VHL (chain A, res 1-103) and POI (chain A, original
numbering). For RIPK2 (4C8B, res 8-316) and MAPK14 (3GCP, res 4-352), residues 8-103
and 4-103 respectively overlap numerically with VHL, causing pRosettaC to silently drop
the kinase N-lobe. The N-lobe contains the ATP pocket where SJF8240 binds.

Fix: renumber POI residues to start at 200 (well above VHL's 1-103), preserving the
full kinase domain. Atomic coordinates are unchanged; only residue numbers shift.

Output:
  data/structures/POI/high_affinity_degraders/RIPK2/RIPK2_protein_clean_full.pdb
  data/structures/POI/high_affinity_degraders/MAPK14/MAPK14_protein_clean_full.pdb
  pRosettaC/runs/high_affinity_degraders/RIPK2_full/Protac_params.txt
  pRosettaC/runs/high_affinity_degraders/MAPK14_full/Protac_params.txt
"""

import shutil
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
RUNS_ROOT = PROJECT / 'pRosettaC' / 'runs' / 'high_affinity_degraders'
POI_ROOT = PROJECT / 'data' / 'structures' / 'POI' / 'high_affinity_degraders'
VHL_ROOT = PROJECT / 'data' / 'structures' / 'VHL'


def renumber_pdb(src_pdb, dst_pdb, offset):
    """Shift all residue numbers in src by +offset and write to dst."""
    lines_out = []
    with open(src_pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                try:
                    orig_resnum = int(line[22:26])
                    new_resnum = orig_resnum + offset
                    line = line[:22] + f'{new_resnum:4d}' + line[26:]
                except ValueError:
                    pass
            lines_out.append(line)
    with open(dst_pdb, 'w') as f:
        f.writelines(lines_out)
    return dst_pdb


def write_protac_params(run_dir, poi_pdb, poi_warhead_sdf, poi_warhead_sdf_H, anchor_atom):
    """Write Protac_params.txt for a pRosettaC run."""
    vhl_pdb = VHL_ROOT / '4W9H_VBC_clean.pdb'
    vhl_warhead_sdf = VHL_ROOT / 'SJF8240_vhll_placed.sdf'
    protac_smi = PROJECT / 'data' / 'structures' / 'PROTAC' / 'SJF8240' / 'SJF8240.smi'

    content = f"""Structures: {vhl_pdb} {poi_pdb}
Chains: ABC A
Heads: {vhl_warhead_sdf} {poi_warhead_sdf}
Anchor atoms: 6 {anchor_atom}
Protac: {protac_smi}
Full: True
ClusterName: LOCAL
"""
    params_file = run_dir / 'Protac_params.txt'
    params_file.write_text(content)
    print(f'  [wrote] {params_file}')

    # Also write the _H version used by resume_pipeline
    content_H = content.replace(str(poi_warhead_sdf), str(poi_warhead_sdf_H))
    content_H = content_H.replace(str(vhl_warhead_sdf),
                                  str(VHL_ROOT / 'SJF8240_vhll_placed_H.sdf'))
    (run_dir / 'Protac_params_H.txt').write_text(content_H)


def setup_poi(poi_name, src_clean_pdb, offset, anchor_atom):
    poi_dir = POI_ROOT / poi_name
    src_warhead_sdf = poi_dir / f'{poi_name}_warhead_placed.sdf'
    src_warhead_sdf_H = poi_dir / f'{poi_name}_warhead_placed_H.sdf'

    # 1. Renumber PDB
    dst_pdb = poi_dir / f'{poi_name}_protein_clean_full.pdb'
    renumber_pdb(src_clean_pdb, dst_pdb, offset)

    # Verify residue range
    resnums = set()
    with open(dst_pdb) as f:
        for line in f:
            if line.startswith('ATOM'):
                resnums.add(int(line[22:26]))
    print(f'{poi_name}: renumbered to res {min(resnums)}-{max(resnums)} '
          f'({len(resnums)} unique residues, offset +{offset})')
    print(f'  Full PDB: {dst_pdb}')

    # 2. Set up run directory (copy from existing run, clear Patchdock_Results)
    old_run = RUNS_ROOT / poi_name
    new_run = RUNS_ROOT / f'{poi_name}_full'
    new_run.mkdir(exist_ok=True)

    # Copy any reference files needed
    for fname in ['SJF8240.params', 'SJF8240.pdb']:
        src = old_run / fname
        if src.exists():
            shutil.copy2(src, new_run / fname)

    # 3. Write Protac_params.txt
    write_protac_params(new_run, dst_pdb, src_warhead_sdf, src_warhead_sdf_H, anchor_atom)

    print(f'  Run dir: {new_run}')
    print()
    return new_run


def main():
    print('Setting up full-kinase-domain runs for RIPK2 and MAPK14')
    print('=' * 60)

    # RIPK2: 4C8B res 8-316 → renumber to 200-508 (offset +192)
    # Warhead anchor atom: 1 (same as existing run)
    setup_poi(
        poi_name='RIPK2',
        src_clean_pdb=POI_ROOT / 'RIPK2' / 'RIPK2_protein_clean.pdb',
        offset=192,
        anchor_atom=1,
    )

    # MAPK14: 3GCP res 4-352 → renumber to 200-548 (offset +196)
    setup_poi(
        poi_name='MAPK14',
        src_clean_pdb=POI_ROOT / 'MAPK14' / 'MAPK14_protein_clean.pdb',
        offset=196,
        anchor_atom=1,
    )

    print('Next steps:')
    print('  1. Verify warhead alignment (warhead SDF coordinates unchanged — still valid)')
    print('  2. Rsync new run dirs to ernie:')
    print('       rsync -av "RIPK2_full/" ernie:~/protac/pRosettaC/runs/high_affinity_degraders/RIPK2_full/')
    print('       rsync -av "MAPK14_full/" ernie:~/protac/pRosettaC/runs/high_affinity_degraders/MAPK14_full/')
    print('  3. On ernie (after MET finishes), inside tmux:')
    print('       cd ~/protac/pRosettaC/runs/high_affinity_degraders/RIPK2_full')
    print('       python ~/protac/scripts/rosetta_utils/run_local_prosettac.py Protac_params.txt --n_dock 200')
    print('       # then MAPK14_full after RIPK2_full finishes')


if __name__ == '__main__':
    main()
