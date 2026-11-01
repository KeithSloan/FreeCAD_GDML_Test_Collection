# FreeCAD vs Geant4 centre of mass

Run: `/Users/ksloan/Workbenches/FreeCAD_GDML_Test_Collection/autotest_results/overlap-main-fc114`

Tolerances: volume 1 %, centre of mass 1 % of radius of gyration, inertia 2 % of largest term.  Geant4 values are Monte Carlo estimates.

## fc-1.1.4

DIFF 24, MATCH 49, NO_GEANT4 7

| file | status | volume | cm off / tol (mm) | inertia | G4 load |
|---|---|---|---|---|---|
| COG_Expanded_simple_array.FCStd | NO_GEANT4 |  |  |  | ABORT |
| COG_TestExtrudeSketchBox.FCStd | DIFF | 63.76 % | 0.755 / 0.056 | 70.4 % | OK |
| MultiCommon.FCStd | DIFF | 60.17 % | 5.57 / 0.053 | 9.5 % | OK |
| PathArray.FCStd | DIFF | 0.49 % | 0.139 / 0.11 | 1.6 % | OK |
| PathArray2.FCStd | DIFF | 0.00 % | 0.01 / 0.64 | 99.6 % | OK |
| Sweep.FCStd | DIFF | 1.06 % | 0.131 / 0.16 | 2.0 % | OK |
| Target-Unmeshed.FCStd | NO_GEANT4 |  |  |  | FAILED |
| UnionTest.FCStd | DIFF | 0.68 % | 0.124 / 0.075 | 1.0 % | OK |
| arrayTest-5.FCStd | DIFF | 0.00 % | 0.0141 / 0.6 | 2.4 % | OK |
| arrayTest-6_scaled.FCStd | DIFF | 1.04 % | 77.2 / 0.98 | 183.3 % | OK |
| arrayTest-7.FCStd | DIFF | 0.17 % | 0.119 / 0.48 | 3.2 % | OK |
| arrayTest-7_scaled.FCStd | DIFF | 0.56 % | 0.454 / 0.55 | 38.0 % | OK |
| extrudeLimitation-1.FCStd | DIFF | 1.37 % | 0.0332 / 0.044 | 1.9 % | OK |
| extrudedParabola-2.FCStd | DIFF | 1.61 % | 0.0447 / 0.14 | 2.0 % | OK |
| importObjTest-MMH.FCStd | DIFF | 0.00 % | 1.15 / 0.1 | 4.0 % | OK |
| importObjTest-edited.FCStd | DIFF | 0.00 % | 1.15 / 0.1 | 4.0 % | OK |
| mirrorTest-1b.FCStd | NO_GEANT4 |  |  |  | FAILED |
| mirrorTest-1c.FCStd | NO_GEANT4 |  |  |  | FAILED |
| mirrorTest-3.FCStd | NO_GEANT4 |  |  |  | FAILED |
| mirrorTest-3b.FCStd | NO_GEANT4 |  |  |  | FAILED |
| mirrorTest-3c.FCStd | NO_GEANT4 |  |  |  | FAILED |
| orthoArray-of-Container.FCStd | DIFF | 4.06 % | 0.1 / 0.72 | 4.0 % | OK |
| orthoArray-of-polarArray-of-Part.FCStd | DIFF | 13.20 % | 0.462 / 1.1 | 13.1 % | OK |
| polarArray-of-Part-array-xformed.FCStd | DIFF | 0.00 % | 0.502 / 0.18 | 5.0 % | OK |
| polarArray-of-Part-base-xformed-offcenter.FCStd | DIFF | 0.00 % | 0.502 / 0.34 | 3.7 % | OK |
| polarArray-of-Part-base-xformed.FCStd | DIFF | 0.00 % | 0.506 / 0.18 | 2.4 % | OK |
| polarArray-of-Part.FCStd | DIFF | 0.00 % | 0.64 / 0.18 | 3.1 % | OK |
| revolveTester-5-a.FCStd | DIFF | 0.00 % | 30 / 18 | 1.4 % | OK |
| revolveTester-5.FCStd | DIFF | 0.00 % | 30 / 18 | 1.4 % | OK |
| revolveTester-7-a.FCStd | DIFF | 0.53 % | 30 / 18 | 1.8 % | OK |
| revolveTester-7.FCStd | DIFF | 0.54 % | 30 / 18 | 1.8 % | OK |
| AssemblyTest.FCStd | MATCH | 0.01 % | 0.0574 / 0.39 | 0.8 % | OK |
| COG_TestExtrudeSketch-Circle.FCStd | MATCH | 0.00 % | 0.01 / 0.032 | 0.2 % | OK |
| COG_TestExtrudeSketch-Outer.FCStd | MATCH | 0.00 % | 0.0141 / 0.057 | 0.3 % | OK |
| Chamfer_and_fillet.FCStd | MATCH | 0.00 % | 0.01 / 0.12 | 0.1 % | OK |
| InvoluteGear.FCStd | MATCH | 0.39 % | 0.01 / 0.056 | 0.8 % | OK |
| Loft.FCStd | MATCH | 0.28 % | 0.102 / 0.12 | 0.6 % | OK |
| PointArray.FCStd | MATCH | 0.76 % | 0 / 0.031 | 0.3 % | OK |
| array-of-Part.FCStd | MATCH | 0.00 % | 0.224 / 0.35 | 0.3 % | OK |
| array-of-array-of-solid-offcenter.FCStd | MATCH | 0.06 % | 0.499 / 0.92 | 1.1 % | OK |
| array-of-array-of-solid.FCStd | MATCH | 0.22 % | 0.341 / 0.89 | 0.9 % | OK |
| array-of-mother-daughter-vols-correct.FCStd | MATCH | 0.00 % | 0 / 0.12 | 0.2 % | OK |
| arrayTest-1-a.FCStd | MATCH | 0.00 % | 0 / 0.31 | 0.0 % | OK |
| arrayTest-1.FCStd | MATCH | 0.00 % | 0 / 0.31 | 0.0 % | OK |
| arrayTest-2.FCStd | MATCH | 0.35 % | 0.195 / 0.73 | 0.8 % | OK |
| arrayTest-3.FCStd | MATCH | 0.08 % | 0.0616 / 0.29 | 0.1 % | OK |
| arrayTest-4.FCStd | MATCH | 0.00 % | 0.0173 / 0.21 | 0.5 % | OK |
| assemblyTest-A.FCStd | MATCH | 0.01 % | 0.0574 / 0.39 | 0.8 % | OK |
| assemblyTest-B.FCStd | MATCH | 0.02 % | 0.07 / 0.34 | 1.3 % | OK |
| daughter_in_array.FCStd | MATCH | 0.00 % | 0.01 / 0.15 | 0.0 % | OK |
| duplicated_bolts.FCStd | MATCH | 0.81 % | 0.01 / 0.29 | 0.7 % | OK |
| duplicated_nuts.FCStd | MATCH | 0.25 % | 0.185 / 2.1 | 0.2 % | OK |
| example1b.FCStd | MATCH | 0.00 % | 0 / 0.38 | 0.0 % | OK |
| extrudeAllShapes.FCStd | MATCH | 0.06 % | 0.0245 / 0.23 | 0.2 % | OK |
| extrudeCircleTest-a.FCStd | MATCH | 0.00 % | 0.0141 / 0.077 | 0.1 % | OK |
| extrudeCircleTest.FCStd | MATCH | 0.00 % | 0.0141 / 0.077 | 0.1 % | OK |
| extrudeTester.FCStd | MATCH | 0.59 % | 0.228 / 1.2 | 0.5 % | OK |
| extrudeTesterA.FCStd | MATCH | 0.59 % | 0.228 / 1.2 | 0.5 % | OK |
| extrudedHyperbola.FCStd | MATCH | 0.60 % | 0.0332 / 0.14 | 0.7 % | OK |
| fillet.FCStd | MATCH | 0.08 % | 0 / 0.049 | 0.2 % | OK |
| fillet_Gmshed.FCtxt.FCStd | MATCH | 0.00 % | 0 / 0.049 | 0.1 % | OK |
| orthoArray_of_Cut.FCStd | MATCH | 0.29 % | 0.01 / 0.15 | 0.1 % | OK |
| pathArray-of-part.FCStd | MATCH | 0.00 % | 0.01 / 0.24 | 0.1 % | OK |
| polarArray-7-Base-rotated-translated.FCStd | MATCH | 0.12 % | 0.0566 / 0.25 | 0.1 % | OK |
| polarArray-7-array-rotated-translated.FCStd | MATCH | 0.12 % | 0.0173 / 0.27 | 0.8 % | OK |
| polarArray-7-array-rotated.FCStd | MATCH | 0.12 % | 0.0173 / 0.27 | 0.8 % | OK |
| polarArray-7-array-translated.FCStd | MATCH | 0.12 % | 0.0245 / 0.27 | 0.4 % | OK |
| polarArray-7.FCStd | MATCH | 0.12 % | 0.0245 / 0.27 | 0.8 % | OK |
| polarArray.FCStd | MATCH | 0.00 % | 0 / 0.05 | 0.0 % | OK |
| polarArrray-of-orthoArray-of-Assembly.FCStd | MATCH | 0.00 % | 0 / 0.97 | 1.8 % | OK |
| revolve.FCStd | MATCH | 0.38 % | 0.051 / 0.51 | 0.5 % | OK |
| revolveSingleEdge-Bspline.FCStd | MATCH | 0.33 % | 0.0374 / 0.28 | 0.3 % | OK |
| revolveSingleEdge.FCStd | MATCH | 0.40 % | 0.0224 / 0.34 | 0.7 % | OK |
| revolveTest-torus.FCStd | MATCH | 0.00 % | 0.0141 / 0.13 | 0.1 % | OK |
| revolveTest-xaxis.FCStd | MATCH | 0.00 % | 0.01 / 0.19 | 0.1 % | OK |
| revolveTester-6-a.FCStd | MATCH | 0.25 % | 10 / 18 | 0.4 % | OK |
| revolveTester-6.FCStd | MATCH | 0.25 % | 10 / 18 | 0.4 % | OK |
| revolveTester-7_scaled-a.FCStd | MATCH | 0.22 % | 14.1 / 18 | 0.3 % | OK |
| revolveTester-7_scaled.FCStd | MATCH | 0.22 % | 14.1 / 18 | 0.3 % | OK |
| solids.FCStd | MATCH | 0.05 % | 0.412 / 6.8 | 0.1 % | OK |

### Differences

- **COG_TestExtrudeSketchBox.FCStd** (`COG_TestExtrudeSketchBox-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - volume 2570 vs 931.3 mm^3 (63.76 %)
  - cm (0.03, -0.3, 5.32) vs (-0.47, 0.16, 4.99) mm, off by 0.755 mm (tol 0.056)
  - inertia (1.298e+05, 1.213e+05, 5.682e+04) vs (3.847e+04, 3.819e+04, 1.469e+04) mm^5 (70.4 %)
- **MultiCommon.FCStd** (`MultiCommon-worldVOL.gdml`, 2 FC / 1 G4 volumes)
  - volume 183.9 vs 73.26 mm^3 (60.17 %)
  - cm (4.33, 1.14, 3.68) vs (9.36, 1.36, 6.06) mm, off by 5.57 mm (tol 0.053)
  - inertia (4449, 1.074e+04, 7393) vs (3433, 9722, 6678) mm^5 (9.5 %)
- **PathArray.FCStd** (`PathArray-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - cm (-0, 1.22, 11.5) vs (-0.08, 1.11, 11.5) mm, off by 0.139 mm (tol 0.11)
- **PathArray2.FCStd** (`PathArray2-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - inertia (3.828e+06, 2.028e+06, 2.236e+06) vs (1.665e+04, 1.667e+04, 1.667e+04) mm^5 (99.6 %)
- **Sweep.FCStd** (`Sweep-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - volume 2.838e+04 vs 2.808e+04 mm^3 (1.06 %)
  - inertia (1.043e+07, 8.871e+06, 5.888e+06) vs (1.022e+07, 8.688e+06, 5.829e+06) mm^5 (2.0 %)
- **UnionTest.FCStd** (`UnionTest-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - cm (18.1, 0.3, 7.62) vs (18, 0.23, 7.72) mm, off by 0.124 mm (tol 0.075)
- **arrayTest-5.FCStd** (`arrayTest-5-worldVOL.gdml`, 1 FC / 8 G4 volumes)
  - inertia (9.201e+07, 9.797e+07, 1.862e+08) vs (9.411e+07, 1.005e+08, 1.907e+08) mm^5 (2.4 %)
- **arrayTest-6_scaled.FCStd** (`arrayTest-6_scaled-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - volume 1.17e+05 vs 1.182e+05 mm^3 (1.04 %)
  - cm (11.8, 30.4, 3.57) vs (86.4, 10.7, 2.66) mm, off by 77.2 mm (tol 0.98)
  - inertia (6.569e+08, 5.971e+08, 1.237e+09) vs (3.173e+08, 2.864e+09, 3.172e+09) mm^5 (183.3 %)
- **arrayTest-7.FCStd** (`arrayTest-7-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - inertia (6.866e+07, 6.866e+07, 1.326e+08) vs (7.1e+07, 7.066e+07, 1.369e+08) mm^5 (3.2 %)
- **arrayTest-7_scaled.FCStd** (`arrayTest-7_scaled-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - inertia (1.349e+08, 1.35e+08, 2.607e+08) vs (1.065e+08, 2.34e+08, 3.333e+08) mm^5 (38.0 %)
- **extrudeLimitation-1.FCStd** (`extrudeLimitation-1-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - volume 488.9 vs 482.2 mm^3 (1.37 %)
- **extrudedParabola-2.FCStd** (`extrudedParabola-2-worldVOL.gdml`, 1 FC / 1 G4 volumes)
  - volume 3110 vs 3060 mm^3 (1.61 %)
  - inertia (1.815e+06, 1.673e+06, 3.281e+06) vs (1.779e+06, 1.641e+06, 3.215e+06) mm^5 (2.0 %)
- **importObjTest-MMH.FCStd** (`importObjTest-MMH-worldVOL.gdml`, 1 FC / 8 G4 volumes)
  - cm (0.05, 1.22, 47.7) vs (0.01, 1.83, 48.6) mm, off by 1.15 mm (tol 0.1)
  - inertia (6.988e+06, 7.096e+06, 2.138e+05) vs (7.27e+06, 7.383e+06, 2.289e+05) mm^5 (4.0 %)
- **importObjTest-edited.FCStd** (`importObjTest-edited-worldVOL.gdml`, 1 FC / 8 G4 volumes)
  - cm (0.05, 1.22, 47.7) vs (0.01, 1.83, 48.6) mm, off by 1.15 mm (tol 0.1)
  - inertia (6.988e+06, 7.096e+06, 2.138e+05) vs (7.27e+06, 7.383e+06, 2.289e+05) mm^5 (4.0 %)
- **orthoArray-of-Container.FCStd** (`orthoArray-of-Container-worldVOL.gdml`, 2 FC / 13 G4 volumes)
  - volume 1.652e+04 vs 1.719e+04 mm^3 (4.06 %)
  - inertia (2.116e+08, 1.347e+08, 3.453e+08) vs (2.2e+08, 1.401e+08, 3.59e+08) mm^5 (4.0 %)
- **orthoArray-of-polarArray-of-Part.FCStd** (`orthoArray-of-polarArray-of-Part-worldVOL.gdml`, 1 FC / 40 G4 volumes)
  - volume 2e+04 vs 2.264e+04 mm^3 (13.20 %)
  - inertia (2.253e+08, 1.802e+08, 4.052e+08) vs (2.547e+08, 2.037e+08, 4.581e+08) mm^5 (13.1 %)
- **polarArray-of-Part-array-xformed.FCStd** (`polarArray-of-Part-array-xformed-worldVOL.gdml`, 2 FC / 11 G4 volumes)
  - cm (-15.1, -11.6, -0) vs (-15.6, -11.8, -0) mm, off by 0.502 mm (tol 0.18)
  - inertia (1.481e+06, 2.509e+06, 3.901e+06) vs (1.517e+06, 2.668e+06, 4.096e+06) mm^5 (5.0 %)
- **polarArray-of-Part-base-xformed-offcenter.FCStd** (`polarArray-of-Part-base-xformed-offcenter-worldVOL.gdml`, 2 FC / 11 G4 volumes)
  - cm (-24.8, -1.96, -0) vs (-25.2, -2.16, -0) mm, off by 0.502 mm (tol 0.34)
  - inertia (1.766e+06, 8.405e+06, 1.008e+07) vs (1.832e+06, 8.702e+06, 1.045e+07) mm^5 (3.7 %)
- **polarArray-of-Part-base-xformed.FCStd** (`polarArray-of-Part-base-xformed-worldVOL.gdml`, 1 FC / 10 G4 volumes)
  - cm (4.86, 8.4, -0) vs (4.4, 8.19, -0) mm, off by 0.506 mm (tol 0.18)
  - inertia (1.117e+06, 1.343e+06, 2.372e+06) vs (1.106e+06, 1.399e+06, 2.416e+06) mm^5 (2.4 %)
- **polarArray-of-Part.FCStd** (`polarArray-of-Part-worldVOL.gdml`, 1 FC / 10 G4 volumes)
  - cm (0.01, 9.7, -0) vs (-0.56, 9.41, -0) mm, off by 0.64 mm (tol 0.18)
  - inertia (1.002e+06, 1.458e+06, 2.372e+06) vs (9.946e+05, 1.531e+06, 2.437e+06) mm^5 (3.1 %)
- **revolveTester-5-a.FCStd** (`revolveTester-5-a-worldVOL.gdml`, 2 FC / 2 G4 volumes)
  - cm (-0, -0, 3.39e+03) vs (-0, 0, 3.36e+03) mm, off by 30 mm (tol 18)
- **revolveTester-5.FCStd** (`revolveTester-5-worldVOL.gdml`, 2 FC / 2 G4 volumes)
  - cm (-0, -0, 3.39e+03) vs (-0, 0, 3.36e+03) mm, off by 30 mm (tol 18)
- **revolveTester-7-a.FCStd** (`revolveTester-7-a-worldVOL.gdml`, 1 FC / 2 G4 volumes)
  - cm (0, -0, 3.38e+03) vs (0, 0, 3.35e+03) mm, off by 30 mm (tol 18)
- **revolveTester-7.FCStd** (`revolveTester-7-worldVOL.gdml`, 1 FC / 2 G4 volumes)
  - cm (0, -0, 3.38e+03) vs (0, 0, 3.35e+03) mm, off by 30 mm (tol 18)

