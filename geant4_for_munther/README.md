# Geant4 check results for Munther

Exported GDML that Munther's `load_gdml_color -printCM` rejects, with the
Geant4 output for each, from the full test-collection run on
FreeCAD 1.1.4 with GDML `Main` (`ff23b2bf`), November 2026.

Totals: 80 GDML files, 73 load OK, 6 FAILED, 1 ABORT.
`fc-1.1.4/geant4_cog.md` has the full volume / centre-of-mass / inertia
comparison with FreeCAD (49 MATCH, 24 DIFF, 7 not compared).

## Checker crashes (FAILED)

`Target-Unmeshed`, `mirrorTest-1b`, `-1c`, `-3`, `-3b`, `-3c`.
The geometry loads, then the checker terminates with
`ZMxpvInfiniteVector` (vector divided by zero). The mirror tests use
reflected placements.

## Geant4 abort (ABORT): COG_Expanded_simple_array

To discuss before fixing the exporter:

1. Duplicated array element: `av_2_impr_2_Obj001_pv_1` is placed exactly on
   `pv_0` (Geant4 GeomVol1002, overlap with volume already placed).
2. Border surfaces refer to physvols that are not exported, or are named
   differently, inside expanded assemblies.

Source test files are in `test_FC_files/`.
