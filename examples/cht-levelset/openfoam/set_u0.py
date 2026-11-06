#!/usr/bin/env python3
"""
Initial velocity of the OpenFOAM cylinder cases, as in the CaNS example (../setup.py cylinder): Poiseuille profile with
U_b = 1, u = 1.5 (1 - (2 z/lz - 1)^2), plus the vertical bump w = 0.05 exp(-((x - xc - 1.5)^2 + (z - zc - 0.3)^2)/0.25)
behind the cylinder that triggers the shedding; written into <case>/0/fluid/U from the cell centres <case>/0/fluid/C
  python3 set_u0.py <case>
"""
import re, sys
from pathlib import Path
import numpy as np
LZ, XC, ZC = 4.0, 3.0, 2.0
case = Path(sys.argv[1])
txt = (case/'0'/'fluid'/'C').read_text()
m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(', txt)
n = int(m.group(1)); body = txt[m.end():txt.index('\n)', m.end())]
c = np.array(body.replace('(', ' ').replace(')', ' ').split(), dtype=float).reshape(n, 3)
x, z = c[:, 0], c[:, 2]
u = 1.5*(1 - (2*z/LZ - 1)**2)
w = 0.05*np.exp(-((x - XC - 1.5)**2 + (z - ZC - 0.3)**2)/0.25)
U = case/'0'/'fluid'/'U'; t = U.read_text()
lst = f'internalField   nonuniform List<vector>\n{n}\n(\n' + ''.join(f'({a:.12g} 0 {b:.12g})\n' for a, b in zip(u, w)) + ')\n;'
t2 = re.sub(r'internalField\s+uniform\s+\(0 0 0\);', lambda _: lst, t)
assert t2 != t, 'internalField not replaced'
U.write_text(t2)
print(f'{case}: {n} fluid cells, U_b(cells, unweighted) = {u.mean():.4f}, max w = {w.max():.4f}')
