#!/usr/bin/env python3
"""
Generalised cluster feature extraction + pre-MD check for any pRosettaC POI.

Usage:
    python run_poi_full_analysis.py <POI>          # e.g. DDR2
    python run_poi_full_analysis.py --all          # run all configured POIs

Outputs per POI (written to the POI's run directory):
    cluster_features.csv      — flat ML feature table (all clusters)
    cluster_features.json     — nested cluster data
    premd_descriptors.csv     — pre-MD feature table (top 3 MD candidates)
    premd_descriptors.json    — full nested pre-MD data
    premd_report.txt          — human-readable narrative
"""

import os, sys, glob, math, csv, json, warnings, argparse
import numpy as np
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
#  POI configuration table
# ══════════════════════════════════════════════════════════════════════════════
_RUNS_ROOT = Path(__file__).resolve().parents[2] / 'pRosettaC' / 'runs'

# Legacy alias for backwards compat within this file
RUNS_BASE = _RUNS_ROOT / 'high_affinity_degraders'

POI_CONFIGS = {
    'MET': {
        'poi_min': 1053, 'poi_max': 1345, 'poi_full': 293,
        'score_shift': 9576.0,
        'degrader': True, 'affinity': 'high',
        'category': 'high_affinity_degraders',
        'struct_note': None,
    },
    'DDR2': {
        'poi_min': 553,  'poi_max': 858,  'poi_full': 306,
        'score_shift': 0.0,
        'degrader': True, 'affinity': 'high',
        'category': 'high_affinity_degraders',
        'struct_note': 'Structure: 7AYM (human DDR2, 2.12Å); replaces DDR1/3ZOS',
    },
    'RIPK2': {
        'poi_min': 104,  'poi_max': 316,  'poi_full': 213,
        'score_shift': 1250.0,
        'degrader': True, 'affinity': 'high',
        'category': 'high_affinity_degraders',
        'struct_note': ('N-lobe residues 8-103 absent — numbering overlaps VHL 1-103; '
                        'only C-lobe (104-316) is modelled'),
    },
    'EPHB2': {
        'poi_min': 615,  'poi_max': 894,  'poi_full': 280,
        'score_shift': 0.0,
        'degrader': True, 'affinity': 'high',
        'category': 'high_affinity_degraders',
        'struct_note': 'Structure: 3ZFM (human EphB2 apo, 2.27Å); replaces ADP-bound 2HEN',
    },
    'MAPK14': {
        'poi_min': 104,  'poi_max': 352,  'poi_full': 237,
        'score_shift': 2660.647,
        'degrader': True, 'affinity': 'high',
        'category': 'high_affinity_degraders',
        'struct_note': ('N-lobe residues 4-103 absent — numbering overlaps VHL 1-103; '
                        'only C-lobe (104-352) is modelled'),
    },
    # ── Priority 2: High affinity, no degradation ─────────────────────────────
    'AXL': {
        'poi_min': 521,  'poi_max': 813,  'poi_full': 270,
        'score_shift': 0.0,
        'degrader': False, 'affinity': 'high',
        'category': 'high_affinity_no_degradation',
        'struct_note': None,
    },
    'SLK': {
        'poi_min': 104,  'poi_max': 308,  'poi_full': 288,
        'score_shift': 0.0,
        'degrader': False, 'affinity': 'high',
        'category': 'high_affinity_no_degradation',
        'struct_note': ('Residues 21-103 overlap VHL chain A numbering; '
                        'interface analysis uses only residues 104-308'),
    },
    'ABL1': {
        'poi_min': 235,  'poi_max': 498,  'poi_full': 263,
        'score_shift': 0.0,
        'degrader': False, 'affinity': 'high',
        'category': 'high_affinity_no_degradation',
        'struct_note': None,
    },
    'EPHA2': {
        'poi_min': 600,  'poi_max': 896,  'poi_full': 301,
        'score_shift': 0.0,
        'degrader': False, 'affinity': 'high',
        'category': 'high_affinity_no_degradation',
        'struct_note': None,
    },
    'MAP4K5': {
        'poi_min': 104,  'poi_max': 320,  'poi_full': 320,
        'score_shift': 0.0,
        'degrader': False, 'affinity': 'high',
        'category': 'high_affinity_no_degradation',
        'struct_note': ('Residues 1-103 overlap VHL chain A numbering; '
                        'interface analysis uses only residues 104-320'),
    },
    # ── Priority 3: Low affinity, degraders ───────────────────────────────────
    'TNIK': {
        'poi_min': 104,  'poi_max': 310,  'poi_full': 289,
        'score_shift': 0.0,
        'degrader': True, 'affinity': 'low',
        'category': 'low_affinity_degraders',
        'struct_note': ('Residues 11-103 overlap VHL chain A numbering; '
                        'interface analysis uses only residues 104-310'),
    },
    'PIP4K2C': {
        'poi_min': 104,  'poi_max': 417,  'poi_full': 243,
        'score_shift': 0.0,
        'degrader': True, 'affinity': 'low',
        'category': 'low_affinity_degraders',
        'struct_note': ('Non-kinase phosphatase; residues 45-103 overlap VHL; '
                        'interface analysis uses residues 104-417'),
    },
    'CDK17': {
        'poi_min': 104,  'poi_max': 350,  'poi_full': 350,
        'score_shift': 0.0,
        'degrader': True, 'affinity': 'low',
        'category': 'low_affinity_degraders',
        'struct_note': ('Residues 1-103 overlap VHL chain A numbering; '
                        'interface analysis uses only residues 104-350'),
    },
}

# Fixed across all POIs
VHL_CHAIN, VHL_MIN, VHL_MAX     = 'A', 1, 103
VHL_FULL                        = 103
ELONGB_CHAIN, ELONGB_MIN, ELONGB_MAX = 'B', 17, 111
ELONGB_FULL                     = 95
ELONGC_CHAIN, ELONGC_MIN, ELONGC_MAX = 'C', 62, 202
ELONGC_FULL                     = 141
PROTAC_CHAIN                    = 'X'
PROTAC_N_ATOMS_EXPECTED         = 79


# ══════════════════════════════════════════════════════════════════════════════
#  PDB parsing helpers
# ══════════════════════════════════════════════════════════════════════════════
def parse_pdb(path):
    atoms = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        try:
            atoms.append({
                'record': line[:6].strip(),
                'name':   line[12:16].strip(),
                'resn':   line[17:20].strip(),
                'chain':  line[21],
                'resi':   int(line[22:26]),
                'xyz':    np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                'elem':   line[76:78].strip() if len(line) > 76 else '',
            })
        except (ValueError, IndexError):
            continue
    return atoms


def heavy(atoms):
    """Heavy atoms only (excludes H and records with empty element field)."""
    return [a for a in atoms if a['elem'] not in ('H', '')]


def coords(atoms):
    return np.array([a['xyz'] for a in atoms])


def select(atoms, chain=None, resi_range=None, name=None, resn=None):
    out = atoms
    if chain:
        out = [a for a in out if a['chain'] == chain]
    if resi_range:
        lo, hi = resi_range
        out = [a for a in out if lo <= a['resi'] <= hi]
    if name:
        out = [a for a in out if a['name'] == name]
    if resn:
        out = [a for a in out if a['resn'] == resn]
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  Score loading
# ══════════════════════════════════════════════════════════════════════════════
def load_scores(run_dir, score_shift):
    sc_file = run_dir / 'Patchdock_Results' / 'score.sc'
    scores = {}
    if not sc_file.exists():
        return scores
    for line in sc_file.read_text().splitlines():
        if line.startswith('SCORE:') and 'total_score' not in line:
            parts = line.split()
            try:
                scores[parts[-1]] = float(parts[1]) + score_shift
            except (ValueError, IndexError):
                pass
    return scores


# ══════════════════════════════════════════════════════════════════════════════
#  Feature extraction (per cluster representative)
# ══════════════════════════════════════════════════════════════════════════════
def _protac_quality_rank(pdb_path, cfg):
    """Return 0=PASS, 1=MARGINAL, 2=FAIL based on PROTAC bridging geometry."""
    poi_residues = cfg.get('poi_residues')
    poi_min = cfg.get('poi_min', 0)
    poi_max = cfg.get('poi_max', 99999)
    vhl, poi, protac = [], [], []
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                chain = line[21]
                resnum = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                if chain == 'X':
                    protac.append((x, y, z))
                elif chain == 'A' and resnum <= 103:
                    vhl.append((x, y, z))
                elif chain == 'A':
                    if poi_residues is not None:
                        if resnum in poi_residues:
                            poi.append((x, y, z))
                    elif poi_min <= resnum <= poi_max:
                        poi.append((x, y, z))
    except Exception:
        return 2
    if not protac or not vhl or not poi:
        return 2
    pa = np.array(protac)
    va = np.array(vhl)
    da = np.array(poi)
    mv = np.min(np.linalg.norm(va[:, None] - pa[None, :], axis=2))
    mp = np.min(np.linalg.norm(da[:, None] - pa[None, :], axis=2))
    if mv < 5 and mp < 5:
        return 0
    if mv < 10 and mp < 10:
        return 1
    return 2


def pick_representative(cluster_dir, scores, cfg=None):
    """Pick cluster representative preferring PASS > MARGINAL > score-based."""
    pdbs = sorted(cluster_dir.glob('combined_*_0001.pdb'))
    if not pdbs:
        return None, float('inf')
    best, best_sc, best_q = pdbs[0], float('inf'), 2
    for p in pdbs:
        sc = scores.get(p.stem, float('inf'))
        q = _protac_quality_rank(p, cfg) if cfg is not None else 2
        if (q, sc) < (best_q, best_sc):
            best_q, best_sc, best = q, sc, p
    return best, best_sc


def interface_contacts(atoms, cfg, cutoff=5.0):
    poi_ha  = heavy(select(atoms, chain='A', resi_range=(cfg['poi_min'], cfg['poi_max'])))
    vhl_ha  = heavy(select(atoms, chain='A', resi_range=(VHL_MIN, VHL_MAX)))
    elb_ha  = heavy([a for a in atoms if a['chain'] == ELONGB_CHAIN])
    elc_ha  = heavy([a for a in atoms if a['chain'] == ELONGC_CHAIN])
    vhl_all = vhl_ha + elb_ha + elc_ha
    if not poi_ha or not vhl_all:
        return 0, 0.0, 0, 0
    pc, vc = coords(poi_ha), coords(vhl_all)
    dists = np.sqrt(((pc[:, None, :] - vc[None, :, :]) ** 2).sum(axis=2))
    n_contacts = int((dists < cutoff).sum())
    area = n_contacts * 15.0
    n_poi_iface = int((dists < cutoff).any(axis=1).sum())
    n_vhl_iface = int((dists < cutoff).any(axis=0).sum())
    return n_contacts, area, n_poi_iface, n_vhl_iface


def hbond_count(atoms, cfg, d_cutoff=3.5):
    poi_donors = [a for a in atoms
                  if a['chain'] == 'A' and cfg['poi_min'] <= a['resi'] <= cfg['poi_max']
                  and a['elem'] in ('N', 'O')]
    vhl_acc = [a for a in atoms
               if ((a['chain'] == 'A' and VHL_MIN <= a['resi'] <= VHL_MAX)
                   or a['chain'] in (ELONGB_CHAIN, ELONGC_CHAIN))
               and a['elem'] in ('N', 'O')]
    if not poi_donors or not vhl_acc:
        return 0
    dc, vc = coords(poi_donors), coords(vhl_acc)
    dists = np.sqrt(((dc[:, None, :] - vc[None, :, :]) ** 2).sum(axis=2))
    return int((dists < d_cutoff).sum())


def protac_geometry(atoms):
    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))
    if not ptc:
        return 0, 0.0, 0, 0
    pts = coords(ptc)
    d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(d, 999)
    nn = d.min(axis=1)
    disconn = int((nn > 3.0).sum())
    overlaps = int((nn < 1.3).sum())
    np.fill_diagonal(d, 0)
    e2e = float(d.max())
    return len(ptc), e2e, disconn, overlaps


def exposed_lysines(atoms, cfg):
    poi_lys_nz = [a for a in atoms
                  if a['chain'] == 'A' and cfg['poi_min'] <= a['resi'] <= cfg['poi_max']
                  and a['resn'] == 'LYS' and a['name'] == 'NZ']
    vhl_all_ha = heavy([a for a in atoms
                        if (a['chain'] == 'A' and VHL_MIN <= a['resi'] <= VHL_MAX)
                        or a['chain'] in (ELONGB_CHAIN, ELONGC_CHAIN)])
    poi_ha = heavy(select(atoms, chain='A', resi_range=(cfg['poi_min'], cfg['poi_max'])))
    if not poi_lys_nz:
        return []
    vc = coords(vhl_all_ha)
    pc = coords(poi_ha)
    results = []
    for lys in poi_lys_nz:
        nz = lys['xyz']
        d_vhl = float(np.sqrt(((vc - nz) ** 2).sum(axis=1)).min()) if len(vc) else 999
        d_poi = np.sqrt(((pc - nz) ** 2).sum(axis=1))
        n_buried = int((d_poi < 5.0).sum()) - 1
        results.append((lys['resi'], float(d_vhl), n_buried))
    results.sort(key=lambda x: x[1])
    return results


def e2_distances(atoms, cfg):
    elb_ca = [a for a in atoms if a['chain'] == ELONGB_CHAIN and a['name'] == 'CA']
    poi_lys_nz = [a for a in atoms
                  if a['chain'] == 'A' and cfg['poi_min'] <= a['resi'] <= cfg['poi_max']
                  and a['resn'] == 'LYS' and a['name'] == 'NZ']
    if not elb_ca or not poi_lys_nz:
        return []
    eb_centroid = coords(elb_ca).mean(0)
    out = [(lys['resi'], round(float(np.linalg.norm(lys['xyz'] - eb_centroid)), 1))
           for lys in poi_lys_nz]
    out.sort(key=lambda x: x[1])
    return out


def clash_check(atoms, cfg, cutoff=2.0):
    poi = heavy(select(atoms, chain='A', resi_range=(cfg['poi_min'], cfg['poi_max'])))
    vhl = heavy(select(atoms, chain='A', resi_range=(VHL_MIN, VHL_MAX)))
    elb = heavy([a for a in atoms if a['chain'] == ELONGB_CHAIN])
    elc = heavy([a for a in atoms if a['chain'] == ELONGC_CHAIN])
    vhl_all = vhl + elb + elc
    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))
    prot = heavy([a for a in atoms if a['chain'] != PROTAC_CHAIN])

    def _clashes(A, B):
        if not A or not B:
            return 0
        d = np.sqrt(((coords(A)[:, None, :] - coords(B)[None, :, :]) ** 2).sum(axis=2))
        return int((d < cutoff).sum())

    return _clashes(poi, vhl_all), _clashes(ptc, prot)


def multi_clash(atoms, cfg):
    out = {}
    for label, thr in [('severe_15A', 1.5), ('bad_20A', 2.0), ('warn_25A', 2.5)]:
        ci, cp = clash_check(atoms, cfg, cutoff=thr)
        out[f'clash_iface_{label}'] = ci
        out[f'clash_protac_{label}'] = cp
    return out


def analyse_cluster(cluster_id, cluster_dir, scores, cfg):
    rep, rep_score = pick_representative(cluster_dir, scores, cfg)
    if rep is None:
        return None
    atoms = parse_pdb(rep)
    n_poses = len(list(cluster_dir.glob('combined_*_0001.pdb')))

    n_cont, area, n_poi_i, n_vhl_i = interface_contacts(atoms, cfg)
    n_hb = hbond_count(atoms, cfg)
    n_ptc, e2e, disconn, overlaps = protac_geometry(atoms)
    lys_list = exposed_lysines(atoms, cfg)
    e2_list  = e2_distances(atoms, cfg)
    clashes  = clash_check(atoms, cfg, cutoff=2.0)

    close_lys = [(r, d) for r, d, nb in lys_list if d < 25.0 and nb < 8]
    e2_opt    = [(r, d) for r, d in e2_list if 35.0 < d < 75.0]

    return {
        'cluster': cluster_id,
        'n_poses': n_poses,
        'representative': rep.name,
        'score_orig': rep_score,
        'interface_contacts_5A': n_cont,
        'interface_area_est_A2': round(area),
        'poi_iface_atoms': n_poi_i,
        'vhl_iface_atoms': n_vhl_i,
        'polar_hbond_pairs': n_hb,
        'protac_n_atoms': n_ptc,
        'protac_end_to_end_A': round(e2e, 1),
        'protac_disconnected': disconn,
        'protac_overlaps': overlaps,
        'clash_interface': clashes[0],
        'clash_protac_protein': clashes[1],
        'n_lys_within_25A_vhl': len(close_lys),
        'closest_lys_to_vhl': close_lys[:3],
        'n_lys_e2_optimal': len(e2_opt),
        'best_e2_lys': e2_opt[:3],
    }


# ══════════════════════════════════════════════════════════════════════════════
#  Chain completeness check
# ══════════════════════════════════════════════════════════════════════════════
def chain_completeness(atoms, cfg):
    def _res(chain, lo, hi):
        return {a['resi'] for a in atoms
                if a['chain'] == chain and lo <= a['resi'] <= hi and a['record'] == 'ATOM'}

    vhl  = _res('A',           VHL_MIN,    VHL_MAX)
    poi  = _res('A',           cfg['poi_min'], cfg['poi_max'])
    elb  = _res(ELONGB_CHAIN,  ELONGB_MIN, ELONGB_MAX)
    elc  = _res(ELONGC_CHAIN,  ELONGC_MIN, ELONGC_MAX)

    def _gaps(present, lo, hi):
        return sorted(r for r in range(lo, hi + 1) if r not in present)

    vhl_gaps = _gaps(vhl, VHL_MIN, VHL_MAX)
    poi_gaps = _gaps(poi, cfg['poi_min'], cfg['poi_max'])
    elb_gaps = _gaps(elb, ELONGB_MIN, ELONGB_MAX)
    elc_gaps = _gaps(elc, ELONGC_MIN, ELONGC_MAX)

    # Interface centroid
    vhl_ha = heavy(select(atoms, chain='A', resi_range=(VHL_MIN, VHL_MAX)))
    poi_ha = heavy(select(atoms, chain='A', resi_range=(cfg['poi_min'], cfg['poi_max'])))
    if_cen = (coords(vhl_ha).mean(0) + coords(poi_ha).mean(0)) / 2 if vhl_ha and poi_ha else None

    def _near_iface(gap_list, chain, cutoff=12.0):
        if if_cen is None:
            return []
        near = []
        for r in gap_list:
            flanks = [a for a in atoms if a['chain'] == chain
                      and abs(a['resi'] - r) <= 2 and a['name'] == 'CA']
            if not flanks:
                continue
            min_d = min(np.linalg.norm(a['xyz'] - if_cen) for a in flanks)
            if min_d <= cutoff:
                near.append((r, round(float(min_d), 1)))
        return near

    return {
        'vhl_present': len(vhl), 'vhl_full': VHL_FULL,
        'vhl_pct': round(100 * len(vhl) / VHL_FULL, 1),
        'vhl_gaps': vhl_gaps,
        'vhl_gaps_near_iface': _near_iface(vhl_gaps, 'A'),
        'poi_present': len(poi), 'poi_full': cfg['poi_full'],
        'poi_pct': round(100 * len(poi) / cfg['poi_full'], 1),
        'poi_gaps': poi_gaps,
        'poi_gaps_near_iface': _near_iface(poi_gaps, 'A'),
        'elongB_present': len(elb) > 0, 'elongB_res': len(elb), 'elongB_gaps': elb_gaps,
        'elongC_present': len(elc) > 0, 'elongC_res': len(elc), 'elongC_gaps': elc_gaps,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  Params check
# ══════════════════════════════════════════════════════════════════════════════
def check_params(run_dir):
    info = {}
    for tag in ('PT0', 'PT1'):
        p = run_dir / f'{tag}.params'
        if not p.exists():
            info[tag] = {'exists': False}
            continue
        lines = p.read_text().splitlines()
        n_atom  = sum(1 for l in lines if l.startswith('ATOM '))
        n_bond  = sum(1 for l in lines if l.startswith('BOND'))
        n_icoor = sum(1 for l in lines if l.startswith('ICOOR'))
        charges = []
        charges_ok = True
        for l in lines:
            if l.startswith('ATOM '):
                try:
                    charges.append(float(l.split()[3]))
                except (ValueError, IndexError):
                    charges_ok = False
        info[tag] = {
            'exists': True, 'n_atoms': n_atom, 'n_bonds': n_bond, 'n_icoors': n_icoor,
            'icoor_ok': n_icoor == n_atom,
            'bonds_ok': n_bond >= n_atom - 1,
            'charges_ok': charges_ok,
            'charge_sum': round(sum(charges), 3),
        }
    return info


# ══════════════════════════════════════════════════════════════════════════════
#  Key starting distances
# ══════════════════════════════════════════════════════════════════════════════
def key_distances(atoms, cfg):
    poi_lys_nz = [a for a in atoms
                  if a['chain'] == 'A' and cfg['poi_min'] <= a['resi'] <= cfg['poi_max']
                  and a['resn'] == 'LYS' and a['name'] == 'NZ']
    vhl_ha = heavy([a for a in atoms
                    if (a['chain'] == 'A' and VHL_MIN <= a['resi'] <= VHL_MAX)
                    or a['chain'] in (ELONGB_CHAIN, ELONGC_CHAIN)])
    elb_ca = [a for a in atoms if a['chain'] == ELONGB_CHAIN and a['name'] == 'CA']

    lys_to_vhl, lys_to_e2 = [], []
    if poi_lys_nz and vhl_ha:
        vc = coords(vhl_ha)
        eb_cen = coords(elb_ca).mean(0) if elb_ca else None
        for lys in poi_lys_nz:
            d_vhl = float(np.sqrt(((vc - lys['xyz']) ** 2).sum(axis=1)).min())
            lys_to_vhl.append((lys['resi'], round(d_vhl, 1)))
            if eb_cen is not None:
                lys_to_e2.append((lys['resi'], round(float(np.linalg.norm(lys['xyz'] - eb_cen)), 1)))

    lys_to_vhl.sort(key=lambda x: x[1])
    lys_to_e2.sort(key=lambda x: x[1])
    e2_opt = [(r, d) for r, d in lys_to_e2 if 35.0 < d < 75.0]

    poi_ha = heavy(select(atoms, chain='A', resi_range=(cfg['poi_min'], cfg['poi_max'])))
    vhl_ha2 = heavy([a for a in atoms
                     if (a['chain'] == 'A' and VHL_MIN <= a['resi'] <= VHL_MAX)
                     or a['chain'] in (ELONGB_CHAIN, ELONGC_CHAIN)])
    n_contacts = 0
    if poi_ha and vhl_ha2:
        d = np.sqrt(((coords(poi_ha)[:, None, :] - coords(vhl_ha2)[None, :, :]) ** 2).sum(axis=2))
        n_contacts = int((d < 5.0).sum())

    ptc = heavy(select(atoms, chain=PROTAC_CHAIN))
    e2e = 0.0
    if ptc:
        d = np.sqrt(((coords(ptc)[:, None, :] - coords(ptc)[None, :, :]) ** 2).sum(axis=2))
        e2e = float(d.max())

    return {
        'closest_lys_resid': lys_to_vhl[0][0] if lys_to_vhl else None,
        'closest_lys_A':     lys_to_vhl[0][1] if lys_to_vhl else None,
        'top3_lys_vhl':      lys_to_vhl[:3],
        'n_lys_e2_optimal':  len(e2_opt),
        'best_e2_resid':     e2_opt[0][0] if e2_opt else None,
        'best_e2_A':         e2_opt[0][1] if e2_opt else None,
        'top3_e2_lys':       e2_opt[:3],
        'interface_contacts_5A': n_contacts,
        'protac_e2e_A':      round(e2e, 1),
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MD candidate selection
# ══════════════════════════════════════════════════════════════════════════════
def select_md_candidates(all_features):
    valid = [f for f in all_features
             if f['n_poses'] >= 5
             and f['clash_interface'] <= 20
             and f['protac_disconnected'] == 0]

    if not valid:
        valid = [f for f in all_features if f['n_poses'] >= 1]

    cand1 = max(valid, key=lambda x: x['n_poses'])
    others = [f for f in valid if f['cluster'] != cand1['cluster']]
    cand2 = min(others, key=lambda x: x['score_orig']) if others else None
    already = {cand1['cluster']} | ({cand2['cluster']} if cand2 else set())
    rest = [f for f in valid if f['cluster'] not in already]
    cand3 = max(rest, key=lambda x: x['n_lys_e2_optimal']) if rest else None

    return [c for c in [cand1, cand2, cand3] if c]


# ══════════════════════════════════════════════════════════════════════════════
#  MD verdict
# ══════════════════════════════════════════════════════════════════════════════
def md_verdict(chain_info, clashes, params_info, topo, cfg, cluster_id):
    issues, warnings, actions = [], [], []

    if chain_info['poi_pct'] < 85:
        issues.append(f"POI only {chain_info['poi_pct']}% complete")
    elif chain_info['poi_pct'] < 95:
        warnings.append(f"POI {chain_info['poi_pct']}% complete")

    if chain_info['poi_gaps_near_iface']:
        nr = [f"{r}({d}Å)" for r, d in chain_info['poi_gaps_near_iface']]
        warnings.append(f"POI gaps near interface: {', '.join(nr)}")
    if chain_info['vhl_gaps_near_iface']:
        nr = [f"{r}({d}Å)" for r, d in chain_info['vhl_gaps_near_iface']]
        warnings.append(f"VHL gaps near interface: {', '.join(nr)}")

    if cfg.get('struct_note'):
        warnings.append(f"Structural limitation: {cfg['struct_note']}")

    if not chain_info['elongB_present']:
        issues.append("ElonginB missing — E2-facing geometry undefined")

    for tag in ('PT0', 'PT1'):
        p = params_info.get(tag, {})
        if not p.get('exists'):
            issues.append(f"{tag}.params missing")
        elif not p.get('icoor_ok'):
            issues.append(f"{tag}.params ICOOR incomplete")

    if topo['disconnected'] > 0:
        issues.append(f"PROTAC {topo['disconnected']} disconnected heavy atoms")
    if topo['n_atoms'] != PROTAC_N_ATOMS_EXPECTED:
        issues.append(f"PROTAC atom count {topo['n_atoms']} ≠ expected {PROTAC_N_ATOMS_EXPECTED}")

    if clashes.get('clash_iface_severe_15A', 0) > 3:
        issues.append(f"Severe interface clashes ({clashes['clash_iface_severe_15A']}) → minimisation required")
        actions.append('minimisation_required')
    elif clashes.get('clash_iface_bad_20A', 0) > 5:
        warnings.append(f"Interface bad clashes ({clashes['clash_iface_bad_20A']}) → minimisation recommended")
        actions.append('minimisation_recommended')

    ready = len(issues) == 0
    return {'ready': ready, 'verdict': 'READY' if ready else 'NEEDS_ACTION',
            'issues': issues, 'warnings': warnings, 'actions': actions}


# ══════════════════════════════════════════════════════════════════════════════
#  Report writer
# ══════════════════════════════════════════════════════════════════════════════
def _b(ok): return '✓' if ok else '✗'

def write_premd_report(poi, candidates, run_dir, params_info, cfg, fh):
    def pr(msg=''):
        print(msg)
        print(msg, file=fh)

    pr(f'\n{"="*80}')
    pr(f'Pre-MD Structural Check  —  {poi}–VHL–SJF8240')
    pr(f'{"="*80}')
    if cfg.get('struct_note'):
        pr(f'  ⚠  STRUCTURAL NOTE: {cfg["struct_note"]}')
        pr()

    for i, (feat, chain_info, clashes, dist, verdict, _) in enumerate(candidates):
        pr(f'\n{"─"*80}')
        pr(f'CANDIDATE {i+1}: {feat["cluster"].upper()}  —  {feat.get("_reason", "")}')
        pr(f'  PDB     : {feat["representative"]}')
        pr(f'  Score   : {feat["score_orig"]:.1f} REU')
        pr(f'  Verdict : {verdict["verdict"]}')
        if verdict['issues']:
            pr(f'  ISSUES  : ' + ' | '.join(verdict['issues']))
        if verdict['warnings']:
            pr(f'  Warnings: ' + ' | '.join(verdict['warnings']))
        if verdict['actions']:
            pr(f'  Actions : ' + ' | '.join(verdict['actions']))

        pr(f'\n  Chain completeness:')
        pr(f'    VHL  (A 1-103)               : {chain_info["vhl_present"]}/{chain_info["vhl_full"]} ({chain_info["vhl_pct"]}%)'
           + (f'  gaps={chain_info["vhl_gaps"]}' if chain_info['vhl_gaps'] else ''))
        pr(f'    {poi} (A {cfg["poi_min"]}-{cfg["poi_max"]}) : '
           f'{chain_info["poi_present"]}/{chain_info["poi_full"]} ({chain_info["poi_pct"]}%)'
           + (f'  ⚠ gaps near iface={chain_info["poi_gaps_near_iface"]}' if chain_info['poi_gaps_near_iface'] else ''))
        pr(f'    ElonginB (B)                 : {_b(chain_info["elongB_present"])}  '
           f'{chain_info["elongB_res"]}/{ELONGB_FULL} res')
        pr(f'    ElonginC (C)                 : {_b(chain_info["elongC_present"])}  '
           f'{chain_info["elongC_res"]}/{ELONGC_FULL} res')
        pr(f'    VHL complex                  : '
           f'{"COMPLETE (E2 geometry defined)" if chain_info["elongB_present"] else "INCOMPLETE — ElonginB missing"}')

        pr(f'\n  PROTAC:')
        pr(f'    Heavy atoms  : {feat["protac_n_atoms"]} {_b(feat["protac_n_atoms"] == PROTAC_N_ATOMS_EXPECTED)}')
        pr(f'    Disconnected : {feat["protac_disconnected"]} {_b(feat["protac_disconnected"] == 0)}')
        pr(f'    End-to-end   : {feat["protac_end_to_end_A"]:.1f} Å')
        for tag in ('PT0', 'PT1'):
            p = params_info.get(tag, {})
            if p.get('exists'):
                pr(f'    {tag}.params  : {p["n_atoms"]}at {p["n_bonds"]}bd {p["n_icoors"]}ic  '
                   f'charge_sum={p["charge_sum"]:+.3f}  ICOOR {_b(p["icoor_ok"])}  bonds {_b(p["bonds_ok"])}')

        pr(f'\n  Clashes (heavy atoms):')
        pr(f'    Threshold  │ Interface │ PROTAC–prot')
        for label, tk in [('1.5Å severe','severe_15A'),('2.0Å bad','bad_20A'),('2.5Å warn','warn_25A')]:
            ci = clashes.get(f'clash_iface_{tk}', 0)
            cp = clashes.get(f'clash_protac_{tk}', 0)
            flag = ' ⚠' if (tk=='severe_15A' and ci>3) or (tk=='bad_20A' and ci>5) else ''
            pr(f'    {label:12s} │ {ci:5d}{flag:2s}    │ {cp:5d}')

        pr(f'\n  Key starting distances:')
        pr(f'    Interface contacts (≤5Å)  : {dist["interface_contacts_5A"]}')
        pr(f'    PROTAC end-to-end         : {dist["protac_e2e_A"]:.1f} Å')
        pr(f'    Closest {poi} Lys → VHL:')
        for r, d in dist['top3_lys_vhl']:
            pr(f'      K{r}  {d:.1f} Å' + (' ⚠ very close' if d < 15 else ''))
        pr(f'    E2-optimal Lys (35-75Å from ElonginB): {dist["n_lys_e2_optimal"]}')
        for r, d in dist['top3_e2_lys']:
            pr(f'      K{r}  {d:.1f} Å  ✓')


# ══════════════════════════════════════════════════════════════════════════════
#  Save outputs
# ══════════════════════════════════════════════════════════════════════════════
def save_cluster_features(poi, all_features, md_candidates, run_dir):
    cand_rank = {c['cluster']: i + 1 for i, c in enumerate(md_candidates)}
    csv_path = run_dir / 'cluster_features.csv'
    cols = ['poi','cluster','n_poses','representative','score_orig_REU',
            'interface_contacts_5A','interface_area_est_A2','poi_iface_atoms','vhl_iface_atoms',
            'polar_hbond_pairs','protac_disconnected','protac_overlaps',
            'clash_interface','clash_protac_protein',
            'n_lys_within_25A_vhl','n_lys_e2_optimal','best_e2_lys_top3','md_candidate']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for feat in all_features:
            w.writerow({
                'poi': poi,
                'cluster': feat['cluster'],
                'n_poses': feat['n_poses'],
                'representative': feat['representative'],
                'score_orig_REU': round(feat['score_orig'], 1),
                'interface_contacts_5A': feat['interface_contacts_5A'],
                'interface_area_est_A2': feat['interface_area_est_A2'],
                'poi_iface_atoms': feat['poi_iface_atoms'],
                'vhl_iface_atoms': feat['vhl_iface_atoms'],
                'polar_hbond_pairs': feat['polar_hbond_pairs'],
                'protac_disconnected': feat['protac_disconnected'],
                'protac_overlaps': feat['protac_overlaps'],
                'clash_interface': feat['clash_interface'],
                'clash_protac_protein': feat['clash_protac_protein'],
                'n_lys_within_25A_vhl': feat['n_lys_within_25A_vhl'],
                'n_lys_e2_optimal': feat['n_lys_e2_optimal'],
                'best_e2_lys_top3': ';'.join(f"K{r}({d}A)" for r, d in feat['best_e2_lys']),
                'md_candidate': cand_rank.get(feat['cluster'], ''),
            })
    with open(run_dir / 'cluster_features.json', 'w') as f:
        json.dump({'poi': poi, 'clusters': all_features,
                   'md_candidates': [{'rank': i+1, 'cluster': c['cluster'],
                                       'representative': c['representative']}
                                      for i, c in enumerate(md_candidates)]},
                  f, indent=2)
    print(f'[saved] {csv_path}')


def save_premd_descriptors(poi, candidates_full, run_dir, cfg):
    cols = [
        'poi', 'cluster', 'n_poses', 'representative', 'score_REU', 'md_verdict',
        'degrader', 'affinity',
        'vhl_pct', 'vhl_n_gaps', 'vhl_gaps_near_iface',
        'poi_pct', 'poi_n_gaps', 'poi_gaps_near_iface',
        'struct_limitation',
        'elongB_present', 'elongC_present',
        'protac_n_atoms', 'protac_disconnected', 'protac_e2e_A',
        'pt0_atoms', 'pt0_bonds', 'pt0_charge_sum',
        'pt1_atoms', 'pt1_bonds', 'pt1_charge_sum',
        'clash_iface_severe_15A', 'clash_iface_bad_20A', 'clash_iface_warn_25A',
        'clash_protac_severe_15A', 'clash_protac_bad_20A',
        'interface_contacts_5A',
        'closest_lys_resid', 'closest_lys_A',
        'n_lys_e2_optimal', 'best_e2_resid', 'best_e2_A', 'top3_e2_lys',
        'needs_minimisation', 'issues', 'warnings',
    ]
    records = []
    for feat, chain_info, clashes, dist, verdict, params_info in candidates_full:
        pt0 = params_info.get('PT0', {})
        pt1 = params_info.get('PT1', {})
        records.append({
            'poi': poi, 'cluster': feat['cluster'], 'n_poses': feat['n_poses'],
            'representative': feat['representative'],
            'score_REU': round(feat['score_orig'], 1),
            'md_verdict': verdict['verdict'],
            'degrader': cfg['degrader'], 'affinity': cfg['affinity'],
            'vhl_pct': chain_info['vhl_pct'],
            'vhl_n_gaps': len(chain_info['vhl_gaps']),
            'vhl_gaps_near_iface': ';'.join(f"{r}({d}A)" for r, d in chain_info['vhl_gaps_near_iface']),
            'poi_pct': chain_info['poi_pct'],
            'poi_n_gaps': len(chain_info['poi_gaps']),
            'poi_gaps_near_iface': ';'.join(f"{r}({d}A)" for r, d in chain_info['poi_gaps_near_iface']),
            'struct_limitation': cfg.get('struct_note') or '',
            'elongB_present': chain_info['elongB_present'],
            'elongC_present': chain_info['elongC_present'],
            'protac_n_atoms': feat['protac_n_atoms'],
            'protac_disconnected': feat['protac_disconnected'],
            'protac_e2e_A': feat['protac_end_to_end_A'],
            'pt0_atoms': pt0.get('n_atoms'), 'pt0_bonds': pt0.get('n_bonds'),
            'pt0_charge_sum': pt0.get('charge_sum'),
            'pt1_atoms': pt1.get('n_atoms'), 'pt1_bonds': pt1.get('n_bonds'),
            'pt1_charge_sum': pt1.get('charge_sum'),
            'clash_iface_severe_15A': clashes.get('clash_iface_severe_15A'),
            'clash_iface_bad_20A':    clashes.get('clash_iface_bad_20A'),
            'clash_iface_warn_25A':   clashes.get('clash_iface_warn_25A'),
            'clash_protac_severe_15A': clashes.get('clash_protac_severe_15A'),
            'clash_protac_bad_20A':   clashes.get('clash_protac_bad_20A'),
            'interface_contacts_5A': dist['interface_contacts_5A'],
            'closest_lys_resid': dist['closest_lys_resid'],
            'closest_lys_A': dist['closest_lys_A'],
            'n_lys_e2_optimal': dist['n_lys_e2_optimal'],
            'best_e2_resid': dist['best_e2_resid'],
            'best_e2_A': dist['best_e2_A'],
            'top3_e2_lys': ';'.join(f"K{r}({d}A)" for r, d in dist['top3_e2_lys']),
            'needs_minimisation': 'minimisation' in ' '.join(verdict['actions']),
            'issues': ' | '.join(verdict['issues']),
            'warnings': ' | '.join(verdict['warnings']),
        })

    csv_path = run_dir / 'premd_descriptors.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in records:
            w.writerow(r)

    with open(run_dir / 'premd_descriptors.json', 'w') as f:
        json.dump({'poi': poi, 'candidates': records}, f, indent=2, default=str)

    print(f'[saved] {csv_path}')
    return records


# ══════════════════════════════════════════════════════════════════════════════
#  Per-POI runner
# ══════════════════════════════════════════════════════════════════════════════
def run_poi(poi_name):
    cfg = POI_CONFIGS[poi_name]
    category = cfg.get('category', 'high_affinity_degraders')
    run_dir = _RUNS_ROOT / category / poi_name
    results_dir = run_dir / 'Results'

    if not results_dir.exists():
        print(f'[SKIP] {poi_name}: Results/ not found at {results_dir}')
        return

    print(f'\n{"#"*80}')
    print(f'  {poi_name}')
    print(f'{"#"*80}')

    scores = load_scores(run_dir, cfg['score_shift'])
    params_info = check_params(run_dir)

    # ── Extract features for ALL clusters ───────────────────────────────────
    all_features = []
    cluster_dirs = sorted(results_dir.glob('cluster*'),
                          key=lambda p: int(p.name.replace('cluster', '')))
    print(f'  Extracting features from {len(cluster_dirs)} clusters...', flush=True)
    for cl_dir in cluster_dirs:
        feat = analyse_cluster(cl_dir.name, cl_dir, scores, cfg)
        if feat:
            all_features.append(feat)
        if len(all_features) % 10 == 0:
            print(f'    {len(all_features)}/{len(cluster_dirs)} done', flush=True)

    print(f'  {len(all_features)} clusters extracted')

    # ── Select MD candidates ─────────────────────────────────────────────────
    md_candidates = select_md_candidates(all_features)
    reasons = ['largest cluster', 'best Rosetta score', 'most E2-optimal Lys']
    for i, c in enumerate(md_candidates):
        c['_reason'] = reasons[i] if i < len(reasons) else ''

    # ── Save cluster features ────────────────────────────────────────────────
    save_cluster_features(poi_name, all_features, md_candidates, run_dir)

    # ── Pre-MD checks on candidates ──────────────────────────────────────────
    print(f'  Running pre-MD checks on {len(md_candidates)} candidates...', flush=True)
    candidates_full = []
    for feat in md_candidates:
        pdb_path = results_dir / feat['cluster'] / feat['representative']
        atoms = parse_pdb(pdb_path)
        chain_info = chain_completeness(atoms, cfg)
        clashes = multi_clash(atoms, cfg)
        dist = key_distances(atoms, cfg)
        n_ptc = len(heavy(select(atoms, chain=PROTAC_CHAIN)))
        e2e_val = feat['protac_end_to_end_A']
        disconn = feat['protac_disconnected']
        topo = {'n_atoms': n_ptc, 'disconnected': disconn, 'e2e': e2e_val}
        verdict = md_verdict(chain_info, clashes, params_info, topo, cfg, feat['cluster'])
        candidates_full.append((feat, chain_info, clashes, dist, verdict, params_info))

    # ── Write report ─────────────────────────────────────────────────────────
    report_path = run_dir / 'premd_report.txt'
    with open(report_path, 'w') as fh:
        write_premd_report(poi_name, candidates_full, run_dir, params_info, cfg, fh)
    print(f'[saved] {report_path}')

    # ── Save pre-MD descriptors ───────────────────────────────────────────────
    save_premd_descriptors(poi_name, candidates_full, run_dir, cfg)

    # ── Summary ──────────────────────────────────────────────────────────────
    print(f'\n  MD candidates:')
    for i, (feat, _, _, dist, verdict, _) in enumerate(candidates_full):
        e2_str = ', '.join(f"K{r}({d}A)" for r, d in dist['top3_e2_lys'])
        print(f'    [{i+1}] {feat["cluster"]:10s}  n={feat["n_poses"]:3d}  '
              f'score={feat["score_orig"]:8.1f}  contacts={dist["interface_contacts_5A"]:4d}  '
              f'e2_lys={dist["n_lys_e2_optimal"]:2d}  verdict={verdict["verdict"]}')
        print(f'         top E2 Lys: {e2_str}')
        if verdict['warnings']:
            print(f'         warnings: {"; ".join(verdict["warnings"][:2])}')


# ══════════════════════════════════════════════════════════════════════════════
#  Entry point
# ══════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('poi', nargs='?', help='POI name, e.g. DDR2')
    parser.add_argument('--all', action='store_true', help='Run all configured POIs')
    args = parser.parse_args()

    if args.all:
        for poi in POI_CONFIGS:
            run_poi(poi)
    elif args.poi:
        if args.poi not in POI_CONFIGS:
            sys.exit(f'Unknown POI "{args.poi}". Options: {list(POI_CONFIGS.keys())}')
        run_poi(args.poi)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
