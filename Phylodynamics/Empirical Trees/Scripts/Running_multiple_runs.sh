#!/usr/bin/env bash

set -euo pipefail

XML_IN="emp_combined_compact_subsampled_data1_AL.xml"
SBATCH_TEMPLATE="sub.sh"
NRUNS=8

[[ -f "$XML_IN"        ]] || { echo "❌  XML file $XML_IN not found!"; exit 1; }
[[ -f "$SBATCH_TEMPLATE" ]] || { echo "❌  SBATCH template $SBATCH_TEMPLATE not found!"; exit 1; }

# ---------------------------------------------------------------------------
# 2.  Create r1 … r7, drop files inside, and submit
# ---------------------------------------------------------------------------
for i in $(seq 1 "$NRUNS"); do
    RUNDIR="r${i}"
    echo "▶  Setting up $RUNDIR"

    mkdir -p   "$RUNDIR"

    ln -sf "../$XML_IN"           "$RUNDIR/$XML_IN"
    cp      "$SBATCH_TEMPLATE"    "$RUNDIR/run.sh"

    sed -i "s|emp_combined_compact_subsampled_data1_AL.xml|$XML_IN|" "$RUNDIR/run.sh"

    chmod +x "$RUNDIR/run.sh"

    # -----------------------------------------------------------------------
    # 3.  Submit
    # -----------------------------------------------------------------------
    ( cd "$RUNDIR" && sbatch run.sh )
done

echo "✅  All $NRUNS jobs queued."
