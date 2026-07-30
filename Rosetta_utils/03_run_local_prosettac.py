#!/usr/bin/env python3
"""
Local pRosettaC runner — no cluster scheduler, no Rosetta binaries, no PatchDock.

Usage (run from the target's run directory):
  cd ~/PROTAC\ degradation\ Prediction/pRosettaC/runs/<category>/<target>
  python ~/PROTAC\ degradation\ Prediction/scripts/02_rosetta_utils/03_run_local_prosettac.py \
      Protac_params.txt [--n_dock 200] [--n_workers 4]

What it does:
  1. Applies pyrosetta_adapter patches (PyRosetta + Python OpenBabel replace executables)
  2. Delegates to pRosettaC's main.py with ClusterName=LOCAL (runs locally via
     Python's ProcessPoolExecutor instead of PBS/SGE/SLURM)

Prerequisites (all already installed):
  - PyRosetta       (~/miniconda3, Python 3.13)
  - openbabel-wheel (Python OpenBabel)
  - RDKit, scipy, numpy
  - PRosettaC repo  (~/PRosettaC)
"""

import sys
import os
import argparse
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────

# Path to the folder containing this script.
_SCRIPT_DIR   = Path(__file__).parent

# Project root for the surrounding repository.
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent

# Path to the installed PRosettaC codebase in the user's home directory.
_PROSETTAC    = Path.home() / 'PRosettaC'

# ── inject PRosettaC into path so its modules are importable ──────────────────

# Add the PRosettaC repository to Python's import path.
# This allows lines like `import main as prosettac_main` to work.
if str(_PROSETTAC) not in sys.path:
    sys.path.insert(0, str(_PROSETTAC))

# ── apply monkey-patches BEFORE importing pRosettaC ──────────────────────────

# Add this script directory to the import path so Python can find
# `pyrosetta_adapter.py`, which contains local compatibility patches.
sys.path.insert(0, str(_SCRIPT_DIR))
from pyrosetta_adapter import apply_patches

# Apply the patches before importing PRosettaC itself.
# This is important because the imported PRosettaC code will use the patched behavior.
apply_patches()

# ── now import pRosettaC main ─────────────────────────────────────────────────

# After patching, import the PRosettaC main entry point.
import main as prosettac_main


# Parse command-line arguments for the local runner.
def parse_args():
    p = argparse.ArgumentParser(description='Local pRosettaC runner')
    p.add_argument('params_file', help='Protac_params.txt path')
    p.add_argument('--n_dock',    type=int, default=200,
                   help='Max global docking poses to generate (default 200)')
    p.add_argument('--n_workers', type=int, default=None,
                   help='Parallel workers for local docking (default: CPU count - 1)')
    return p.parse_args()


# Main workflow for running a PRosettaC job locally.
def main():
    args = parse_args()

    # Convert the provided params file path into a Path object and make sure it exists.
    params_path = Path(args.params_file)
    if not params_path.exists():
        sys.exit(f'ERROR: {params_path} not found.\n'
                 f'Run this script from the target run directory.')

    # If the user requested a specific number of workers, patch the LOCAL cluster
    # object so PRosettaC uses that worker count when parallelizing locally.
    if args.n_workers:
        from cluster.LOCAL.LOCAL import LOCAL as LocalCluster
        LocalCluster.__init__ = lambda self, _=None: setattr(self, 'n_workers', args.n_workers)

    # Optionally scale the number of global docking poses.
    # PRosettaC normally asks for 1000 poses when Full=True (or 500 otherwise).
    # This wrapper replaces that value with a smaller user-chosen cap, which is
    # often more practical on a local workstation.
    import pyrosetta_adapter as _adp
    _orig_global_dock = _adp.global_dock

    def _scaled_global_dock(structs, anchors, min_dist, max_dist,
                             num_results=1000, threshold=2.0):
        return _orig_global_dock(structs, anchors, min_dist, max_dist,
                                 num_results=min(num_results, args.n_dock),
                                 threshold=threshold)

    # Replace PRosettaC's patchdock-like global docking function with the capped version.
    import utils
    utils.patchdock = _scaled_global_dock

    print(f'[run_local] Starting pRosettaC for: {params_path.resolve()}')
    print(f'[run_local] Max global poses: {args.n_dock}  |  Workers: {args.n_workers or "auto"}')

    # Hand off the main run to PRosettaC's own main function.
    prosettac_main.main(sys.argv[0], [str(params_path)])

    # PRosettaC's main.py builds some constraint-generation paths incorrectly when
    # absolute paths are used for heads/linker inputs. To fix that, run the follow-up
    # resume script, which redoes constraint generation and clustering with correct paths.
    import subprocess
    print('\n[run_local] Re-running constraint_gen + clustering via resume_pipeline...')
    subprocess.run(
        [sys.executable, str(_SCRIPT_DIR / 'resume_pipeline.py')],
        check=True,
        cwd=os.getcwd(),
    )


# Standard Python entry point: only run main() when this file is executed directly.
if __name__ == '__main__':
    main()
