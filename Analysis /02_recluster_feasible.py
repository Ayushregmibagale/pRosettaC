#!/usr/bin/env python3
"""
Re-cluster only the PASS/MARGINAL (bridging) poses, following the pRosettaC paper:
  1. Filter combined PDBs to only those where PROTAC bridges VHL and POI
  2. Cluster by protein Cα RMSD (4 Å threshold, average linkage)
  3. Rank clusters by size (most populated = most sampled geometry)
  4. Take top 3 clusters, picking best PASS > best MARGINAL > best score representative

Usage:
    python recluster_feasible.py DDR2
    python recluster_feasible.py MET RIPK2_full MAPK14_full
"""

import sys
import glob
import json
import csv
import numpy as np
from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

RUNS_ROOT = Path(__file__).resolve().parents[2] / 'pRosettaC' / 'runs'

POI_CONFIGS = {
    # poi_min/poi_max: full residue range of the POI in chain A
    # poi_check_min: stricter lower bound for PROTAC bridging check
    #   (for N-lobe truncated kinases where POI starts at 104, adjacent to VHL res 103,
    #    requiring PROTAC to reach residues deep in the C-lobe avoids the artifact
    #    where 92%+ of poses trivially pass because PROTAC is near the VHL-POI junction)
    'MET':    {'poi_min': 1053, 'poi_max': 1345, 'poi_check_min': 1053,
               'category': 'high_affinity_degraders', 'poi_full': 293, 'degrader': True},
    'DDR2':   {'poi_min': 553,  'poi_max': 858,  'poi_check_min': 553,
               'category': 'high_affinity_degraders', 'poi_full': 306, 'degrader': True},
    # RIPK2_full / MAPK14_full: renumbered to res 200+ to fix N-lobe artifact (Bug 21)
    # These supersede the old RIPK2/MAPK14 runs (N-lobe overwrite → invalid combined PDBs)
    'RIPK2_full': {'poi_min': 200, 'poi_max': 508, 'poi_check_min': 290,
                   'category': 'high_affinity_degraders', 'poi_full': 309, 'degrader': True},
    'EPHB2':  {'poi_min': 615,  'poi_max': 894,  'poi_check_min': 615,
               'category': 'high_affinity_degraders', 'poi_full': 280, 'degrader': True},
    'MAPK14_full': {'poi_min': 200, 'poi_max': 548, 'poi_check_min': 290,
                    'category': 'high_affinity_degraders', 'poi_full': 349, 'degrader': True},
    'AXL':    {'poi_min': 473,  'poi_max': 712,  'poi_check_min': 473,
               'category': 'high_affinity_no_degradation', 'poi_full': 270, 'degrader': False},
    'SLK':    {'poi_min': 1,    'poi_max': 340,  'poi_check_min': 1,
               'category': 'high_affinity_no_degradation', 'poi_full': 288, 'degrader': False},
    'ABL1':   {'poi_min': 229,  'poi_max': 500,  'poi_check_min': 229,
               'category': 'high_affinity_no_degradation', 'poi_full': 263, 'degrader': False},
    'EPHA2':  {'poi_min': 596,  'poi_max': 896,  'poi_check_min': 596,
               'category': 'high_affinity_no_degradation', 'poi_full': 300, 'degrader': False},
    'MAP4K5': {'poi_min': 1,    'poi_max': 290,  'poi_check_min': 1,
               'category': 'high_affinity_no_degradation', 'poi_full': 290, 'degrader': False},
}

FEASIBLE_DIR = Path.home()
RMSD_CUTOFF = 4.0
MIN_CLUSTER_SIZE = 1
# For POIs where the kinase construct is large relative to VHL, nearly all poses can be
# "feasible" by loose geometry. Re-clustering thousands of poses is expensive and gives the
# same result as the original clustering — cap at MAX_POSES_FULL_CLUSTER.
# Fall back to the existing pRosettaC Results/ cluster structure when feasible set is too large.
MAX_POSES_FOR_RECLUSTERING = 500


def load_ca_coords(pdb_path):
    """Extract Cα coordinates keyed by residue number (chain A only)."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            if line[21] != 'A':
                continue
            resnum = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords[resnum] = np.array([x, y, z])
    return coords


def kabsch_rmsd(a, b):
    """RMSD after optimal superposition (Kabsch algorithm)."""
    ca = a.mean(0)
    cb = b.mean(0)
    a = a - ca
    b = b - cb
    H = a.T @ b
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    a_rot = a @ R.T
    return float(np.sqrt(((a_rot - b) ** 2).sum() / len(a)))


def pairwise_rmsd(coords_list):
    n = len(coords_list)
    dist = np.zeros((n, n))
    # Find common residues across all structures for a consistent comparison
    common = set(coords_list[0].keys())
    for c in coords_list[1:]:
        common &= set(c.keys())
    common = sorted(common)
    if not common:
        return dist
    mats = [np.array([c[r] for r in common]) for c in coords_list]
    for i in range(n):
        for j in range(i + 1, n):
            rmsd = kabsch_rmsd(mats[i].copy(), mats[j].copy())
            dist[i, j] = dist[j, i] = rmsd
    return dist


def load_scores(score_sc):
    scores = {}
    if not score_sc.exists():
        return scores
    for line in score_sc.read_text().splitlines():
        if line.startswith('SCORE:') and 'total_score' not in line:
            parts = line.split()
            try:
                scores[parts[-1]] = float(parts[1])
            except (ValueError, IndexError):
                pass
    return scores


def protac_quality(pdb, poi_check_min, poi_max):
    """PROTAC quality check using the strict poi_check_min residue lower bound."""
    vhl, poi, protac = [], [], []
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
            elif chain == 'A' and poi_check_min <= resnum <= poi_max:
                poi.append([x, y, z])
    if not protac or not vhl or not poi:
        return 2
    pa, va, da = np.array(protac), np.array(vhl), np.array(poi)
    mv = np.min(np.linalg.norm(va[:, None] - pa[None, :], axis=2))
    mp = np.min(np.linalg.norm(da[:, None] - pa[None, :], axis=2))
    if mv < 5 and mp < 5:
        return 0
    if mv < 10 and mp < 10:
        return 1
    return 2


def run_poi(poi_name):
    cfg = POI_CONFIGS[poi_name]
    run_dir = RUNS_ROOT / cfg['category'] / poi_name
    patchdock = run_dir / 'Patchdock_Results'
    score_sc = patchdock / 'score.sc'
    out_dir = run_dir / 'Results_feasible'

    # --- 1. Load feasible poses ---
    feasible_file = FEASIBLE_DIR / f'{poi_name.lower()}_feasible.txt'
    if not feasible_file.exists():
        print(f'{poi_name}: no feasible file at {feasible_file} — run find_feasible_poses first')
        return

    project_root = RUNS_ROOT.parents[1]   # ~/PROTAC degradation Prediction
    poi_check_min = cfg['poi_check_min']
    poi_max = cfg['poi_max']
    needs_recheck = poi_check_min != cfg['poi_min']

    feasible_paths = []
    verdicts = {}
    n_recheck_fail = 0
    for line in feasible_file.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) != 2:
            continue
        verdict, path = parts
        p = Path(path)
        if not p.is_absolute():
            p = project_root / path
        if not p.exists():
            continue
        # Re-check using poi_check_min (stricter lower bound) to require PROTAC near
        # C-lobe residues, avoiding the VHL-adjacency artifact where nearly all poses pass.
        if needs_recheck:
            q = protac_quality(p, poi_check_min, poi_max)
            if q == 2:
                n_recheck_fail += 1
                continue
            verdicts[str(p)] = q
        else:
            verdicts[str(p)] = 0 if verdict == 'PASS' else 1
        feasible_paths.append(p)

    if needs_recheck and n_recheck_fail > 0:
        print(f'  [recheck] {n_recheck_fail} poses failed strict criterion (poi_check_min={poi_check_min}), kept {len(feasible_paths)}')

    if not feasible_paths:
        # No bridging poses in Patchdock_Results — try scanning Results/ cluster dirs directly
        print(f'{poi_name}: no valid paths in feasible file; scanning Results/ cluster dirs for PASS/MARGINAL...')
        results_dir = run_dir / 'Results'
        poi_check_min = cfg['poi_check_min']
        poi_max = cfg['poi_max']
        for cluster_dir in sorted(results_dir.glob('cluster*'),
                                  key=lambda d: int(d.name.replace('cluster', ''))):
            for p in cluster_dir.glob('combined_*_0001.pdb'):
                q = protac_quality(p, poi_check_min, poi_max)
                if q < 2:
                    feasible_paths.append(p)
                    verdicts[str(p)] = q
        if feasible_paths:
            print(f'  Found {len(feasible_paths)} PASS/MARGINAL in Results/ cluster dirs')
        else:
            print(f'{poi_name}: no bridging poses found anywhere — skipping')
            return

    n = len(feasible_paths)
    print(f'{poi_name}: {n} feasible poses → clustering (RMSD cutoff {RMSD_CUTOFF} Å)...', flush=True)

    if n < 3:
        print(f'{poi_name}: fewer than 3 feasible poses — skipping clustering, taking all as candidates')
        scores = load_scores(score_sc)
        candidates = sorted(feasible_paths,
                            key=lambda p: (verdicts.get(str(p), 2), scores.get(p.stem, float('inf'))))
        _write_results(poi_name, candidates[:3], [1] * len(candidates), run_dir, cfg, scores, verdicts)
        return

    if n > MAX_POSES_FOR_RECLUSTERING:
        # Too many feasible poses to efficiently re-cluster (typically N-lobe truncated kinases
        # where the bridging check doesn't filter meaningfully). Fall back to the existing
        # pRosettaC Results/ cluster directories, picking the best PASS/MARGINAL representative.
        print(f'  {n} poses exceeds {MAX_POSES_FOR_RECLUSTERING} limit — using existing pRosettaC clusters')
        scores = load_scores(score_sc)
        _fallback_to_existing_clusters(poi_name, run_dir, cfg, scores, verdicts, feasible_paths)
        return

    # --- 2. Compute pairwise RMSD ---
    print(f'  Loading Cα coordinates...', flush=True)
    coords_list = []
    valid_paths = []
    for p in feasible_paths:
        try:
            c = load_ca_coords(p)
            if len(c) > 10:
                coords_list.append(c)
                valid_paths.append(p)
        except Exception as e:
            print(f'  WARNING: could not load {p.name}: {e}')

    print(f'  Computing {len(valid_paths)}×{len(valid_paths)} RMSD matrix...', flush=True)
    dist_matrix = pairwise_rmsd(coords_list)

    # --- 3. Hierarchical clustering ---
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')
    labels = fcluster(Z, t=RMSD_CUTOFF, criterion='distance')

    # --- 4. Rank clusters by size ---
    from collections import Counter
    cluster_counts = Counter(labels)
    ranked = [cl for cl, _ in cluster_counts.most_common()]

    scores = load_scores(score_sc)
    _write_results(poi_name, valid_paths, labels, run_dir, cfg, scores, verdicts, ranked)


def _fallback_to_existing_clusters(poi_name, run_dir, cfg, scores, verdicts, feasible_paths):
    """Use existing pRosettaC Results/ cluster dirs; pick best PASS/MARGINAL rep per cluster."""
    results_dir = run_dir / 'Results'
    # Build stem → quality map so cluster dir copies are recognised regardless of path
    stem_quality = {p.stem: verdicts[str(p)] for p in feasible_paths}

    cluster_dirs = sorted(results_dir.glob('cluster*'),
                          key=lambda d: int(d.name.replace('cluster', '')))
    if not cluster_dirs:
        print(f'{poi_name}: no cluster dirs found in {results_dir}')
        return

    candidates = []
    for cluster_dir in cluster_dirs:
        pdbs = sorted(cluster_dir.glob('combined_*_0001.pdb'))
        if not pdbs:
            continue
        # Prefer PASS > MARGINAL > score within this cluster
        best, best_q, best_sc = None, 2, float('inf')
        for p in pdbs:
            q = stem_quality.get(p.stem, 2)
            sc = scores.get(p.stem, float('inf'))
            if (q, sc) < (best_q, best_sc):
                best_q, best_sc, best = q, sc, p
        if best is None:
            continue
        q_str = ['PASS', 'MARGINAL', 'FAIL'][best_q]
        candidates.append({'cluster': cluster_dir.name, 'n_members': len(pdbs),
                            'representative': str(best), 'name': best.name,
                            'score': best_sc if best_sc != float('inf') else float('nan'),
                            'quality': q_str})

    # Rank by cluster size (n_members), take top 3
    candidates.sort(key=lambda c: -c['n_members'])
    top3 = candidates[:3]

    print(f'\n{"="*60}')
    print(f'  {poi_name} — Existing pRosettaC Clusters (PASS-aware reps)')
    print(f'{"="*60}')
    print(f'  {len(candidates)} clusters | using top 3 by size')
    for c in top3:
        print(f"  {c['cluster']} (n={c['n_members']}): {c['name']}  score={c['score']:.1f}  {c['quality']}")

    out_csv = run_dir / 'feasible_md_candidates.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['rank', 'cluster_size', 'representative', 'name', 'score', 'quality'])
        w.writeheader()
        for i, c in enumerate(top3, 1):
            w.writerow({'rank': i, 'cluster_size': c['n_members'],
                        'representative': c['representative'], 'name': c['name'],
                        'score': c['score'], 'quality': c['quality']})
    print(f'\n[saved] {out_csv}')


def _write_results(poi_name, paths, labels, run_dir, cfg, scores, verdicts, ranked=None):
    from collections import defaultdict
    clusters = defaultdict(list)
    for p, lbl in zip(paths, labels):
        clusters[lbl].append(p)

    if ranked is None:
        ranked = sorted(clusters, key=lambda cl: len(clusters[cl]), reverse=True)

    out_dir = run_dir / 'Results_feasible'
    out_dir.mkdir(exist_ok=True)

    print(f'\n{"="*60}')
    print(f'  {poi_name} — Feasible-only Clusters')
    print(f'{"="*60}')
    print(f'  {len(clusters)} clusters from {len(paths)} bridging poses')

    candidates = []
    for rank, cl in enumerate(ranked, 1):
        members = clusters[cl]
        # pick best: PASS first, then MARGINAL, then score
        best = min(members, key=lambda p: (verdicts.get(str(p), 2),
                                            scores.get(p.stem, float('inf'))))
        q = verdicts.get(str(best), 2)
        sc = scores.get(best.stem, float('nan'))
        q_str = ['PASS', 'MARGINAL', 'FAIL'][q]
        print(f'  Cluster {rank} (n={len(members)}): {best.name}  score={sc:.1f}  {q_str}')
        candidates.append({'rank': rank, 'cluster_size': len(members),
                            'representative': str(best), 'name': best.name,
                            'score': sc, 'quality': q_str})
        if rank >= 3:
            break

    # Save candidate list
    out_csv = run_dir / 'feasible_md_candidates.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['rank', 'cluster_size', 'representative', 'name', 'score', 'quality'])
        w.writeheader()
        w.writerows(candidates)
    print(f'\n[saved] {out_csv}')
    return candidates


def main():
    pois = sys.argv[1:] if len(sys.argv) > 1 else list(POI_CONFIGS)
    for poi in pois:
        if poi not in POI_CONFIGS:
            print(f'Unknown POI: {poi}')
            continue
        run_poi(poi)


if __name__ == '__main__':
    main()
