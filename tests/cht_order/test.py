#!/usr/bin/env python
#
# assertions on the output of testit.sh (CHT_ORDER_WORK: its work directory)
#
import json, os
from pathlib import Path

WORK = Path(os.environ.get('CHT_ORDER_WORK', Path(__file__).resolve().parents[2]/'run'/'cht_order_ci'))


def test_fidelity():
    # one RK3 step of CaNS and of the assembled operator agree to round-off
    for r in json.load(open(WORK/'fidelity.json')):
        assert r['rel_to_change'] < 1e-12, r
        assert r['max_u_after'] == 0.0, r


def test_geometry_and_conservation():
    res = json.load(open(WORK/'results.json'))
    assert all(r['link_mismatch'] == 0 for r in res)
    assert max(r['colsum'] for r in res) < 1e-13


def test_order_on_face():
    # Step 1 is second order when the interface lies on a row of cell faces
    res = sorted([r for r in json.load(open(WORK/'results.json')) if r['family'] == 'flat-on'], key=lambda r: r['N'])
    assert res[-1]['order_L2'] > 1.9 and res[-1]['order_Linf'] > 1.9, res[-1]


if __name__ == '__main__':
    test_fidelity()
    test_geometry_and_conservation()
    test_order_on_face()
    print('Passed!')
