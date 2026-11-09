#!/bin/bash
#
# quick check of the conjugate conduction operator of the CHT scalars (see README.md):
# the operator the harness assembles is the one CaNS applies (one RK3 step, both geometry paths),
# and the steady order on a reduced set of cases
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
WORKDIR=$CANSDIR/run/cht_order_ci
rm -rf "${WORKDIR:?}"; mkdir -p "$WORKDIR"
echo "INFO: operator fidelity"
python3 cht_order.py fidelity --cans "$CANSDIR/build/cans" --work "$WORKDIR"
echo "INFO: steady order, reduced set"
python3 cht_order.py run --quick --cans "$CANSDIR/build/cans" --work "$WORKDIR" --paths sdf \
        --families flat-on,flat-off,corr --kappas 2 --aspects 1 --jobs 4 --out results.json
echo "INFO: assertions"
CHT_ORDER_WORK="$WORKDIR" pytest -q test.py
rm -rf "${WORKDIR:?}"
