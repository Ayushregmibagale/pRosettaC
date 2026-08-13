#!/usr/bin/env python3
"""
Find all PASS/MARGINAL PROTAC bridging poses in Patchdock_Results for given POIs.

PROTAC (chain X) must be within 5 Å of both VHL (chain A res 1-103) and POI
to be PASS, or within 10 Å of both to be MARGINAL.

Run this ON ERNIE from ~/protac/:
    python find_feasible_poses.py DDR2 EPHB2

Writes ~/ddr2_feasible.txt, ~/ephb2_feasible.txt with lines:
    PASS<tab>relative/path/to/combined_X_Y_0001.pdb
    MARGINAL<tab>relative/path/...
"""

import sys
import glob
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

RUNS_ROOT = Path.home() / 'protac' / 'pRosettaC' / 'runs'
PROTAC_ROOT = Path.home() / 'protac'

POI_CONFIGS = {
    # Phase 2 degrader set (n_dock=200, 10000 poses each)
    'MET':         {'poi_min': 1053, 'poi_max': 1345, 'category': 'high_affinity_degraders'},
    'DDR2':        {'poi_min': 553,  'poi_max': 858,  'category': 'high_affinity_degraders'},
    'EPHB2':       {'poi_min': 615,  'poi_max': 894,  'category': 'high_affinity_degraders'},
    # RIPK2_full / MAPK14_full: renumbered to res 200+ to fix N-lobe artifact (Bug 21)
    # Supersede the old RIPK2/MAPK14 runs which had VHL cleft overwritten by the N-lobe
    'RIPK2_full':  {'poi_min': 200,  'poi_max': 508,  'category': 'high_affinity_degraders'},
    'MAPK14_full': {'poi_min': 200,  'poi_max': 548,  'category': 'high_affinity_degraders'},
}


def check_pdb(args):
    pdb, poi_min, poi_max = args
    vhl, poi, protac = [], [], []
    try:
        with open(pdb) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                chain = line[21]
                resnum = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                if chain == 'X':
                    protac.append([x, y, z])
                elif chain == 'A' and resnum <= 103:
                    vhl.append([x, y, z])
                elif chain == 'A' and poi_min <= resnum <= poi_max:
                    poi.append([x, y, z])
    except Exception:
        return pdb, 'FAIL'
    if not protac or not vhl or not poi:
        return pdb, 'FAIL'
    pa = np.array(protac)
    va = np.array(vhl)
    da = np.array(poi)
    mv = np.min(np.linalg.norm(va[:, None] - pa[None, :], axis=2))
    mp = np.min(np.linalg.norm(da[:, None] - pa[None, :], axis=2))
    if mv < 5 and mp < 5:
        return pdb, 'PASS'
    if mv < 10 and mp < 10:
        return pdb, 'MARGINAL'
    return pdb, 'FAIL'


def run_poi(poi_name):
    cfg = POI_CONFIGS[poi_name]
    patchdock = RUNS_ROOT / cfg['category'] / poi_name / 'Patchdock_Results'
    pdbs = sorted(patchdock.glob('combined_*_0001.pdb'))
    if not pdbs:
        print(f'{poi_name}: no combined PDBs found in {patchdock}')
        return

    print(f'{poi_name}: checking {len(pdbs)} poses with 8 workers...', flush=True)
    args = [(str(p), cfg['poi_min'], cfg['poi_max']) for p in pdbs]

    feasible = []
    done = 0
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(check_pdb, a): a for a in args}
        for f in as_completed(futs):
            pdb_path, verdict = f.result()
            done += 1
            if verdict in ('PASS', 'MARGINAL'):
                rel = str(Path(pdb_path).relative_to(PROTAC_ROOT))
                feasible.append(f"{verdict}\t{rel}")
            if done % 1000 == 0:
                print(f'  {done}/{len(pdbs)} checked, {len(feasible)} feasible', flush=True)

    outfile = Path.home() / f'{poi_name.lower()}_feasible.txt'
    with open(outfile, 'w') as fh:
        fh.write('\n'.join(feasible) + '\n')
    print(f'{poi_name}: {len(feasible)} PASS/MARGINAL → {outfile}')


def main():
    pois = sys.argv[1:] if len(sys.argv) > 1 else list(POI_CONFIGS)
    for poi in pois:
        if poi not in POI_CONFIGS:
            print(f'Unknown POI: {poi}')
            continue
        run_poi(poi)


if __name__ == '__main__':
    main()
