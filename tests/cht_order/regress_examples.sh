#!/bin/bash
#
# regress_examples.sh <cans> <refdir> <rundir> [case ...] [-- key=value ...]
#
# Bit-for-bit regression of the CHT examples (examples/cht-levelset, default cases: slab wavy).
# Each case runs with <cans> in <rundir>/<case> (one MPI task, from the example's input.nml and
# level set). If <refdir>/<case>.md5 does not exist, the checksums of the outputs (data/fld_*.bin,
# data/forcing*.out, data/ibm_force.out) are SAVED there as the reference; otherwise they are
# CHECKED against it. check.py is run in both cases. Entries after `--` are set in input.nml
# before the run (e.g. `-- is_cht_tang(:)=T`), a key that is absent is appended to &scalar.
# Exit status 1 if any checksum differs. The reference is only meaningful for the same compiler,
# flags and machine (runs are bit-reproducible there).
#
set -e
CANS=$(realpath "$1"); REF=$(realpath -m "$2"); RUN=$(realpath -m "$3"); shift 3
CASES=(); SETS=()
while [ $# -gt 0 ]; do
  if [ "$1" == "--" ]; then shift; SETS=("$@"); break; fi
  CASES+=("$1"); shift
done
[ ${#CASES[@]} -eq 0 ] && CASES=(slab wavy)
EX=$(cd "$(dirname "$0")"/../../examples/cht-levelset && pwd)
mkdir -p "$REF" "$RUN"
status=0
for c in "${CASES[@]}"; do
  D="$RUN/$c"
  rm -rf "${D:?}"; mkdir -p "$D/data"
  cp "$EX/$c/input.nml" "$D/"
  (cd "$EX" && python3 setup.py "$c" > /dev/null) && cp "$EX/$c/data/sdfp.bin" "$D/data/"
  for kv in "${SETS[@]}"; do
    k=${kv%%=*}; v=${kv#*=}
    kre=$(printf '%s' "$k" | sed 's/[][()*.^$]/\\&/g')
    if grep -q "^[[:space:]]*$kre[[:space:]]*=" "$D/input.nml"; then
      sed -i "s|^\([[:space:]]*\)$kre[[:space:]]*=.*|\1$k = $v|" "$D/input.nml"
    else
      sed -i "/^&scalar/a $k = $v" "$D/input.nml"
    fi
  done
  echo "INFO: running $c in $D"
  (cd "$D" && mpirun -np 1 "$CANS" > run.log 2>&1) || { echo "ERROR: CaNS failed in $D"; exit 1; }
  python3 "$EX/check.py" "$c" "$D"
  (cd "$D/data" && md5sum fld_*.bin forcing*.out ibm_force.out > ../outputs.md5)
  if [ -f "$REF/$c.md5" ]; then
    if diff -q "$REF/$c.md5" "$D/outputs.md5" > /dev/null; then
      echo "PASS: $c bit-identical to $REF/$c.md5"
    else
      echo "FAIL: $c differs from $REF/$c.md5:"; diff "$REF/$c.md5" "$D/outputs.md5" || true; status=1
    fi
  else
    cp "$D/outputs.md5" "$REF/$c.md5"; echo "SAVED: reference $REF/$c.md5"
  fi
done
exit $status
