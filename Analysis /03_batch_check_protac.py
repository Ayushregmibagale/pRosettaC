#!/usr/bin/env python3
"""
Batch PROTAC placement check across all combined PDBs for one or more POIs.

Usage:
    python batch_check_protac.py RIPK2 EPHB2 MAPK14
    python batch_check_protac.py --all

Checks whether the PROTAC (chain X) is positioned between VHL (chain A res 1-103)
and the POI (chain A res poi_min-poi_max) in each combined_*_0001.pdb.

Verdict:
    PASS     — closest PROTAC atom < 5 A to both VHL and POI
    MARGINAL — closest PROTAC atom < 10 A to both
    FAIL     — PROTAC too far from one or both proteins
"""

import sys
import glob
import argparse
from pathlib import Path
import numpy as np

RUNS_ROOT = Path(__file__).resolve().parents[2] / 'pRosettaC' / 'runs'

POI_CONFIGS = {
    'MET':    {'poi_min': 1053, 'poi_max': 1345, 'category': 'high_affinity_degraders'},
    'DDR2':   {'poi_min': 553,  'poi_max': 858,  'category': 'high_affinity_degraders'},
    # RIPK2: N-lobe (8-103) absent; residue numbering immediately follows VHL (1-103).
    # Using range-based check is invalid — VHL and RIPK2 are directly adjacent in chain A.
    # Instead, use specific C-lobe ATP-site residues (derived from RIPK2_warhead_placed_H.sdf, 5A cutoff).
    'RIPK2':  {'poi_residues': {105,109,112,113,122,144,153,163,164,165,
                                261,264,265,267,268,269,271,281,286,289,
                                290,292,293,294,296,313},
               'category': 'high_affinity_degraders'},
    'EPHB2':  {'poi_min': 615,  'poi_max': 894,  'category': 'high_affinity_degraders'},
    'MAPK14': {'poi_min': 104,  'poi_max': 352,  'category': 'high_affinity_degraders'},
    'AXL':    {'poi_min': 473,  'poi_max': 712,  'category': 'high_affinity_no_degradation'},
    'SLK':    {'poi_min': 1,    'poi_max': 340,  'category': 'high_affinity_no_degradation'},
    'ABL1':   {'poi_min': 229,  'poi_max': 500,  'category': 'high_affinity_no_degradation'},
    'EPHA2':  {'poi_min': 596,  'poi_max': 896,  'category': 'high_affinity_no_degradation'},
    'MAP4K5': {'poi_min': 1,    'poi_max': 290,  'category': 'high_affinity_no_degradation'},
}


def check_pdb(pdb, cfg):
    poi_residues = cfg.get('poi_residues')
    poi_min = cfg.get('poi_min')
    poi_max = cfg.get('poi_max')

    vhl, poi, protac = [], [], []
    with open(pdb) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain  = line[21]
            resnum = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if chain == 'X':
                protac.append([x, y, z])
            elif chain == 'A' and resnum <= 103:
                vhl.append([x, y, z])
            elif chain == 'A':
                if poi_residues is not None:
                    if resnum in poi_residues:
                        poi.append([x, y, z])
                elif poi_min <= resnum <= poi_max:
                    poi.append([x, y, z])
    if not protac or not vhl or not poi:
        return 'FAIL'
    pa = np.array(protac)
    va = np.array(vhl)
    da = np.array(poi)
    min_vhl = np.min(np.linalg.norm(va[:, None] - pa[None, :], axis=2))
    min_poi = np.min(np.linalg.norm(da[:, None] - pa[None, :], axis=2))
    if min_vhl < 5 and min_poi < 5:
        return 'PASS'
    elif min_vhl < 10 and min_poi < 10:
        return 'MARGINAL'
    return 'FAIL'


def run_poi(poi_name):
    cfg = POI_CONFIGS[poi_name]
    pattern = str(RUNS_ROOT / cfg['category'] / poi_name /
                  'Patchdock_Results' / 'combined_*_0001.pdb')
    pdbs = sorted(glob.glob(pattern))
    if not pdbs:
        print(f'{poi_name}: no combined PDBs found at {pattern}')
        return
    counts = {'PASS': 0, 'MARGINAL': 0, 'FAIL': 0}
    for i, pdb in enumerate(pdbs):
        counts[check_pdb(pdb, cfg)] += 1
        if (i + 1) % 500 == 0:
            print(f'  {poi_name}: {i+1}/{len(pdbs)} checked...', flush=True)
    total = len(pdbs)
    pct = 100 * (counts['PASS'] + counts['MARGINAL']) / total if total else 0
    print(f"{poi_name} ({total} poses) — "
          f"PASS: {counts['PASS']}  MARGINAL: {counts['MARGINAL']}  FAIL: {counts['FAIL']}  "
          f"({pct:.1f}% feasible)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pois', nargs='*', choices=list(POI_CONFIGS) + [[]], metavar='POI',
                        help='POI name(s) to check')
    parser.add_argument('--all', action='store_true', help='Check all configured POIs')
    args = parser.parse_args()

    targets = list(POI_CONFIGS) if args.all else args.pois
    if not targets:
        parser.print_help()
        sys.exit(1)

    for poi in targets:
        if poi not in POI_CONFIGS:
            print(f'Unknown POI: {poi}. Choose from: {list(POI_CONFIGS)}')
            continue
        run_poi(poi)


if __name__ == '__main__':
    main()
#!/usr/bin/env python3
"""
Batch PROTAC placement check across all combined PDBs for one or more POIs.

Usage:
    python batch_check_protac.py RIPK2 EPHB2 MAPK14
    python batch_check_protac.py --all

Checks whether the PROTAC (chain X) is positioned between VHL (chain A res 1-103)
and the POI (chain A res poi_min-poi_max) in each combined_*_0001.pdb.

Verdict:
    PASS     — closest PROTAC atom < 5 A to both VHL and POI
    MARGINAL — closest PROTAC atom < 10 A to both
    FAIL     — PROTAC too far from one or both proteins
"""

import sys
import glob
import argparse
from pathlib import Path
import numpy as np

RUNS_ROOT = Path(__file__).resolve().parents[2] / 'pRosettaC' / 'runs'

POI_CONFIGS = {
    'MET':    {'poi_min': 1053, 'poi_max': 1345, 'category': 'high_affinity_degraders'},
    'DDR2':   {'poi_min': 553,  'poi_max': 858,  'category': 'high_affinity_degraders'},
    # RIPK2: N-lobe (8-103) absent; residue numbering immediately follows VHL (1-103).
    # Using range-based check is invalid — VHL and RIPK2 are directly adjacent in chain A.
    # Instead, use specific C-lobe ATP-site residues (derived from RIPK2_warhead_placed_H.sdf, 5A cutoff).
    'RIPK2':  {'poi_residues': {105,109,112,113,122,144,153,163,164,165,
                                261,264,265,267,268,269,271,281,286,289,
                                290,292,293,294,296,313},
               'category': 'high_affinity_degraders'},
    'EPHB2':  {'poi_min': 615,  'poi_max': 894,  'category': 'high_affinity_degraders'},
    'MAPK14': {'poi_min': 104,  'poi_max': 352,  'category': 'high_affinity_degraders'},
    'AXL':    {'poi_min': 473,  'poi_max': 712,  'category': 'high_affinity_no_degradation'},
    'SLK':    {'poi_min': 1,    'poi_max': 340,  'category': 'high_affinity_no_degradation'},
    'ABL1':   {'poi_min': 229,  'poi_max': 500,  'category': 'high_affinity_no_degradation'},
    'EPHA2':  {'poi_min': 596,  'poi_max': 896,  'category': 'high_affinity_no_degradation'},
    'MAP4K5': {'poi_min': 1,    'poi_max': 290,  'category': 'high_affinity_no_degradation'},
}


def check_pdb(pdb, cfg):
    poi_residues = cfg.get('poi_residues')
    poi_min = cfg.get('poi_min')
    poi_max = cfg.get('poi_max')

    vhl, poi, protac = [], [], []
    with open(pdb) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain  = line[21]
            resnum = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if chain == 'X':
                protac.append([x, y, z])
            elif chain == 'A' and resnum <= 103:
                vhl.append([x, y, z])
            elif chain == 'A':
                if poi_residues is not None:
                    if resnum in poi_residues:
                        poi.append([x, y, z])
                elif poi_min <= resnum <= poi_max:
                    poi.append([x, y, z])
    if not protac or not vhl or not poi:
        return 'FAIL'
    pa = np.array(protac)
    va = np.array(vhl)
    da = np.array(poi)
    min_vhl = np.min(np.linalg.norm(va[:, None] - pa[None, :], axis=2))
    min_poi = np.min(np.linalg.norm(da[:, None] - pa[None, :], axis=2))
    if min_vhl < 5 and min_poi < 5:
        return 'PASS'
    elif min_vhl < 10 and min_poi < 10:
        return 'MARGINAL'
    return 'FAIL'


def run_poi(poi_name):
    cfg = POI_CONFIGS[poi_name]
    pattern = str(RUNS_ROOT / cfg['category'] / poi_name /
                  'Patchdock_Results' / 'combined_*_0001.pdb')
    pdbs = sorted(glob.glob(pattern))
    if not pdbs:
        print(f'{poi_name}: no combined PDBs found at {pattern}')
        return
    counts = {'PASS': 0, 'MARGINAL': 0, 'FAIL': 0}
    for i, pdb in enumerate(pdbs):
        counts[check_pdb(pdb, cfg)] += 1
        if (i + 1) % 500 == 0:
            print(f'  {poi_name}: {i+1}/{len(pdbs)} checked...', flush=True)
    total = len(pdbs)
    pct = 100 * (counts['PASS'] + counts['MARGINAL']) / total if total else 0
    print(f"{poi_name} ({total} poses) — "
          f"PASS: {counts['PASS']}  MARGINAL: {counts['MARGINAL']}  FAIL: {counts['FAIL']}  "
          f"({pct:.1f}% feasible)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pois', nargs='*', choices=list(POI_CONFIGS) + [[]], metavar='POI',
                        help='POI name(s) to check')
    parser.add_argument('--all', action='store_true', help='Check all configured POIs')
    args = parser.parse_args()

    targets = list(POI_CONFIGS) if args.all else args.pois
    if not targets:
        parser.print_help()
        sys.exit(1)

    for poi in targets:
        if poi not in POI_CONFIGS:
            print(f'Unknown POI: {poi}. Choose from: {list(POI_CONFIGS)}')
            continue
        run_poi(poi)


if __name__ == '__main__':
    main()
