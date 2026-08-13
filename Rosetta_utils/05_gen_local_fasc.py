#!/usr/bin/env python3
"""
Score existing docking PDB files and write local.fasc.
Run from the Patchdock_Results/ directory.

Usage:
  cd Patchdock_Results
  python /path/to/gen_local_fasc.py <PT0.params> <PT1.params>
"""
import sys
import os
import glob

def main():
    if len(sys.argv) < 3:
        print('Usage: gen_local_fasc.py <PT0.params> <PT1.params>')
        sys.exit(1)

    pt0, pt1 = sys.argv[1], sys.argv[2]

    import pyrosetta
    pyrosetta.init(
        f'-extra_res_fa {pt0} -extra_res_fa {pt1} '
        '-ignore_unrecognized_res true -load_PDB_components false -out:level 100'
    )
    scorefxn = pyrosetta.create_score_function('ref2015')

    pdbs = sorted(glob.glob('pd.*_docking_????.pdb'))
    print(f'Scoring {len(pdbs)} docking PDBs...', flush=True)

    with open('local.fasc', 'w') as lf:
        for k, pdb in enumerate(pdbs):
            try:
                pose = pyrosetta.pose_from_file(pdb)
                score = scorefxn(pose)
            except Exception as e:
                print(f'  SKIP {pdb}: {e}')
                score = 0.0
            desc = pdb.replace('.pdb', '')
            lf.write(f'SCORE:  {score:.3f}  0  0  0  {score:.3f}  {desc}\n')
            if (k + 1) % 100 == 0:
                print(f'  {k+1}/{len(pdbs)}', flush=True)

    print('Done → local.fasc')

if __name__ == '__main__':
    main()
