# GDML 2.2.13 on FreeCAD 1.1.4 — 25 November 2026

Regression run of the GDML workbench changes since 1 November.
Write-up: GDML repository `TESTING.md`, section *GDML 2.2.13 on FreeCAD 1.1.4*.

| | baseline | this run |
|---|---|---|
| run (`autotest_results/`) | `overlap-main-fc114` | `gdml-2.2.13` |
| date | 2026-11-01 | 2026-11-25 |
| FreeCAD | 1.1.4 (Python 3.11.14, OCC 7.8.1) | same |
| GDML | Main `ff23b2bf` | Main `cc28775f` (2.2.13) |
| test repo | `test_FC_1.1.4@7d44b8c` | `test_GDML_2.2.13` (from `test_FC_1.1.4@04a7905`) |

## FreeCAD side

| | result vs baseline |
|---|---|
| statuses | 98/98 OK (after the recheck below) |
| COG files (`cogs/`) | 80/80 byte-identical |
| exported GDML | 79/80 identical; `extrudedParabola-2` differs by 6e-14 mm |
| round-trip GDML | 18/18 identical |
| export warnings | unchanged: `PathArray`, `array-of-array-of-solid`, `assembly.gdml` (overlaps) |

22 files were first `OK_WITH_DIALOGS`, all *Document Recovery*: a FreeCAD was
already running when the run started. Run `gdml-2.2.13-recheck` (those 22
files, no other FreeCAD open, GDML Main `7dfca778`, 2.2.14) had no dialogs.

## Geant4 (`autotest/run_geant4.sh`)

| | baseline | this run |
|---|---|---|
| MATCH / DIFF / not compared | 49 / 24 / 7 | 49 / 24 / 7 |
| not loaded | `COG_Expanded_simple_array` (ABORT), `Target-Unmeshed`, `mirrorTest-1b`, `-1c`, `-3`, `-3b`, `-3c` | same 7 |

No file changed status.

## Scripting API self-test

GDML `freecad/gdml/tests/ai_functions/smoke_test.py`, run in the same
FreeCAD 1.1.4 through the FreeCAD MCP server; its export
`ai_functions_smoke.gdml` was added to the run's `exported_gdml/` before the
Geant4 stage.

- Found two bugs in the new API (fixed in GDML 2.2.14):
  `makePolyhedra` / `makeGenericPolyhedra` (group extension) and
  `makeTessellated` (unhashable vertices). Then 57/57 checks passed.
- The export defines `C_element` / `H_element` for the formula material, and
  Geant4 loads it (`OK`). Its 22 `GeomVol1002` overlap warnings came from
  the self-test placing most solids at the origin (changed in GDML 2.2.15).

## Reproduce

```sh
cd ~/Workbenches/FreeCAD_GDML_Test_Collection
AT_RUN_ID=gdml-2.2.13 autotest/run_versions.sh
autotest/run_geant4.sh autotest_results/gdml-2.2.13
# side-by-side with the baseline
mkdir -p autotest_results/gdml-2.2.13-vs-nov01 && cd autotest_results/gdml-2.2.13-vs-nov01
ln -s ../overlap-main-fc114/fc-1.1.4 fc-1.1.4-nov01
ln -s ../gdml-2.2.13/fc-1.1.4 fc-1.1.4-today
cd ../.. && python3 autotest/compare.py autotest_results/gdml-2.2.13-vs-nov01
```
