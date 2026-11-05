# Level-set (signed-distance) input for `&ibm geometry = 'sdf'`

The level-set IBM reads the signed distance `phi` at the cell centres from `sdf_file` (default `data/sdfp.bin`):
raw float64, `ng(1)*ng(2)*ng(3)` values, x fastest (C order `(nz, ny, nx)`), `phi > 0` in the fluid.
The solver needs the correct sign everywhere and a smooth, accurate zero level within about two cells of the
wall; `|grad phi| = 1` is needed only for the per-point wall outputs and for `is_ibm_damp`, so far from the wall
`phi` may be clipped.

* `hmap_sdf.py <case> <map.npy|map.txt> --zm <mean wall height> [--upper mirror|shift|none] [--lmap lx ly] [--np N]`:
  rough walls from a periodic height map; exact distance within a band of the wall (Newton projection on a
  bicubic spline of the map), slope-corrected vertical distance beyond it.
* `stl_sdf.py <case> wall.stl [...] [options]`: any geometry from triangulated surfaces (closed solids or
  periodic sheets); exact distance in a band (libigl), sign from pseudonormals or the winding number, flood
  fill outside the band. Needs `igl` (libigl Python bindings) and SciPy.

Both read `ng` and `l` from `<case>/input.nml` and the z grid from `<case>/grid_zf.in` (`gtype = 9`),
`<case>/data/grid.bin` or a uniform spacing, and print the quantiles of `|grad phi|` near the wall as a check.

Simple shapes can be written directly, e.g. two flat walls at z = 0.5 and 2.5 in a domain of height 3:
```python
phi = np.minimum(zc - 0.5, 2.5 - zc)[:, None, None]*np.ones((nz, ny, nx)); phi.tofile('data/sdfp.bin')
```
and combined with `min` (union of solids) and `max` (intersection), e.g. a cylinder of radius R on that wall:
`phi = min(phi_walls, sqrt((x - xc)**2 + (z - zc)**2) - R)` (with the nearest periodic image in x).
