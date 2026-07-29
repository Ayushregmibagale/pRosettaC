"""
Download and clean PDB structures for all POIs and VHL.
Outputs cleaned PDBs into data/structures/POI/<category>/<TARGET>/
"""

import os
import subprocess
import urllib.request
from pathlib import Path

# Find the project root based on where this script lives.
# `parents[2]` means "go up two folders from this file".
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# All downloaded and cleaned structure files will be stored under this folder.
DATA_DIR = PROJECT_ROOT / "data" / "structures"

# A table describing each protein of interest (POI).
# Each row contains:
# (target name, category folder, PDB ID to download, chain to keep, priority)
POI_TABLE = [
    # (target, category, pdb_id, chain, priority)
    ("MET",     "high_affinity_degraders",     "3DKC", "A", 1),
    ("DDR2",    "high_affinity_degraders",     "3ZOS", "A", 1),
    ("RIPK2",   "high_affinity_degraders",     "5AR3", "A", 1),
    ("EPHB2",   "high_affinity_degraders",     "2HEN", "A", 1),
    ("MAPK14",  "high_affinity_degraders",     "3GCP", "A", 1),
    ("AXL",     "high_affinity_no_degradation","5U6B", "A", 2),
    ("SLK",     "high_affinity_no_degradation","2J51", "A", 2),
    ("ABL1",    "high_affinity_no_degradation","2HYY", "A", 2),
    ("EPHA2",   "high_affinity_no_degradation","3MBW", "A", 2),
    ("MAP4K5",  "high_affinity_no_degradation","6YSS", "A", 2),
    ("TNIK",    "low_affinity_degraders",      "4UXT", "A", 3),
    # MAPK14 reuses same PDB as above — symlink rather than re-download
    ("PIP4K2C", "low_affinity_degraders",      "2GK9", "A", 3),
    # CDK17 has no crystal structure — placeholder for AlphaFold model
]

# VHL is handled separately from the POI table.
# This tuple stores: (name, PDB ID, chain to keep).
VHL_PDB = ("VHL", "4W9H", "B")   # chain B is VHL in 4W9H


def fetch_pdb(pdb_id: str, out_path: Path) -> None:
    """Download a PDB file from the RCSB website."""
    # Convert the PDB ID to uppercase and build the download URL.
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"  Fetching {pdb_id} → {out_path.name}")

    # Save the downloaded file to the requested output path.
    urllib.request.urlretrieve(url, out_path)



def extract_chain(in_pdb: Path, chain: str, out_pdb: Path) -> None:
    """Keep only ATOM/HETATM records for the requested chain."""
    # Collect the lines we want to keep from the original PDB file.
    lines = []
    with open(in_pdb) as fh:
        for line in fh:
            # Keep atom records only if they belong to the requested chain.
            # In the PDB text format, the chain ID is stored at character index 21.
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
                lines.append(line)
            # Also preserve TER/END records so the output remains a valid PDB-style file.
            elif line.startswith(("TER", "END")):
                lines.append(line)

    # Write the filtered lines into a new chain-specific PDB file.
    with open(out_pdb, "w") as fh:
        fh.writelines(lines)
    print(f"    Chain {chain} extracted → {out_pdb.name}")



def main():
    # First handle VHL, which uses its own PDB ID and chain.
    vhl_dir = DATA_DIR / "VHL"
    raw_vhl = vhl_dir / f"{VHL_PDB[1]}.pdb"

    # Download the full VHL PDB only if it is not already present.
    if not raw_vhl.exists():
        fetch_pdb(VHL_PDB[1], raw_vhl)

    # Create a second file that contains only the VHL chain.
    chain_vhl = vhl_dir / f"{VHL_PDB[1]}_chain{VHL_PDB[2]}.pdb"
    if not chain_vhl.exists():
        extract_chain(raw_vhl, VHL_PDB[2], chain_vhl)

    # Now process every POI listed in the table above.
    for target, category, pdb_id, chain, priority in POI_TABLE:
        # Build the folder where this target's files should live.
        out_dir = DATA_DIR / "POI" / category / target

        # Create the folder (and any missing parent folders) if needed.
        out_dir.mkdir(parents=True, exist_ok=True)

        # Path for the full downloaded PDB file.
        raw_pdb = out_dir / f"{pdb_id}.pdb"
        if not raw_pdb.exists():
            fetch_pdb(pdb_id, raw_pdb)
        else:
            print(f"  {target} ({pdb_id}) already present, skipping download.")

        # Path for the cleaned file containing only the requested chain.
        chain_pdb = out_dir / f"{target}_{pdb_id}_chain{chain}.pdb"
        if not chain_pdb.exists():
            extract_chain(raw_pdb, chain, chain_pdb)

    # Special case: MAPK14 in the low-affinity set reuses the same files
    # as the high-affinity set, so create symlinks instead of copying/downloading again.
    src = DATA_DIR / "POI" / "high_affinity_degraders" / "MAPK14"
    dst = DATA_DIR / "POI" / "low_affinity_degraders" / "MAPK14"
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        link = dst / f.name
        if not link.exists():
            link.symlink_to(f)
    print("  MAPK14 low_affinity_degraders symlinked to high_affinity entry.")

    # Final status message and pointer to the next step in the workflow.
    print("\nDone. Next step: place warhead coordinates in each POI structure.")
    print("See scripts/structure_prep/place_warhead.py")


# Run `main()` only when this file is executed directly.
# If this file is imported into another Python script, `main()` will not run automatically.
if __name__ == "__main__":
    main()
