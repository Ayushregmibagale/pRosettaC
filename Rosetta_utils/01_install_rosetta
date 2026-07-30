#!/bin/bash
# Rosetta + pRosettaC installation guide for SJF8240 campaign
#
# Step 1: Get Rosetta (free academic license)
#   https://www.rosettacommons.org/software/license-and-download
#   Download: rosetta_src_<version>_bundle.tgz
#
# Step 2: Compile Rosetta (or download pre-compiled binaries)
#   tar -xzf rosetta_src_<version>_bundle.tgz
#   cd rosetta_src_<version>/main/source
#   ./scons.py -j$(nproc) mode=release bin/rosetta_scripts
#   ./scons.py -j$(nproc) mode=release bin/docking_protocol
#
# Step 3: Set environment variable
#   export ROSETTA_HOME=/path/to/rosetta_src_<version>
#   export ROSETTA_BIN=$ROSETTA_HOME/main/source/bin
#   Add these to ~/.bashrc
#
# Step 4: Get pRosettaC
#   git clone https://github.com/LiorZ/pRosettaC.git ~/pRosettaC
#   export PROSETTAC_HOME=~/pRosettaC
#
# Step 5: Regenerate SJF8240 params with proper CHI dihedrals
# (run this once ROSETTA_HOME is set)

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SDF="${PROJECT_ROOT}/data/structures/PROTAC/SJF8240/SJF8240.sdf"
PARAMS_DIR="${PROJECT_ROOT}/data/ligand_params"

if [[ -z "${ROSETTA_HOME:-}" ]]; then
    echo "ERROR: ROSETTA_HOME not set."
    echo "Download from https://www.rosettacommons.org/software/license-and-download"
    exit 1
fi

PARAMS_SCRIPT="${ROSETTA_HOME}/main/source/scripts/python/public/molfile_to_params.py"
if [[ ! -f "$PARAMS_SCRIPT" ]]; then
    echo "ERROR: molfile_to_params.py not found at ${PARAMS_SCRIPT}"
    exit 1
fi

echo "Regenerating SJF8240.params with Rosetta molfile_to_params.py..."
python "$PARAMS_SCRIPT" \
    "$SDF" \
    -n SJF \
    -p "${PARAMS_DIR}/SJF8240.params" \
    --conformers-in-one-file \
    --no-pdb

echo ""
echo "Regenerating fragment params..."
for frag in warhead vhll; do
    python "$PARAMS_SCRIPT" \
        "${PROJECT_ROOT}/data/structures/PROTAC/SJF8240/SJF8240_${frag}.sdf" \
        -n $(echo $frag | tr '[:lower:]' '[:upper:]' | cut -c1-3) \
        -p "${PARAMS_DIR}/SJF8240_${frag}.params" \
        --no-pdb
done

echo ""
echo "All params files regenerated in ${PARAMS_DIR}/"
echo ""
echo "Next: run batch_run.sh to launch all pRosettaC jobs."
echo "  cd ${PROJECT_ROOT}/scripts/rosetta_utils"
echo "  PROSETTAC_HOME=~/pRosettaC ROSETTA_BIN=\$ROSETTA_BIN ./batch_run.sh"
