# GDML autotest — run the collection against several FreeCAD builds

Runs every file in `FC_test_files/` and `test_gdml_files/` through each FreeCAD
listed in `versions.conf`, unattended, and writes a comparison report.

Per FreeCAD file: **open → recompute → export GDML → centre-of-mass / inertia**
Per GDML file: **import → recompute → re-export (round trip)**

The report compares every version against the others **and** against the
committed baselines in `FC_COGs/` and `FC_exported_gdml/`.

## Quick start (macOS)

```sh
cd ~/Workbenches/FreeCAD_GDML_Test_Collection

# 1. Check each FreeCAD sees the GDML checkout you want to test
autotest/check_setup.sh            # report
autotest/check_setup.sh --link     # symlink <Mod>/GDML -> ~/Workbenches/GDML where missing

# 2. Run (FreeCAD windows will open and close; leave them alone)
autotest/run_versions.sh                     # every version in versions.conf
autotest/run_versions.sh fc-26.3rc2          # one version
AT_ONLY='Loft*' autotest/run_versions.sh     # a subset of files

# 3. Read autotest_results/<timestamp>/report.md
```

Before the first run, open each new FreeCAD once by hand and dismiss any
first-run / Gatekeeper dialogs. Quit other FreeCAD instances while it runs.

## Why the GUI, not freecadcmd

`exportGDML.export()` takes the world volume from `FreeCADGui.Selection`, and
`calcCenterOfMass.py` walks the GUI document tree. So each version is launched
as `FreeCAD autotest/gdml_autotest_runner.py`. The runner:

- defers its work until the GUI event loop is running;
- dismisses any modal dialog after 2 s and **logs it against the file** (these
  show up in the report — useful in themselves);
- turns off auto-recovery files, so a relaunch after a crash isn't blocked by
  the Document Recovery dialog;
- writes one JSON line per file to `results.jsonl` and skips files already
  recorded, so it resumes after a crash.

`run_versions.sh` watches a heartbeat file. If a stage takes longer than
`AT_STAGE_TIMEOUT` (900 s) FreeCAD is killed; if FreeCAD crashes, the file it
was on is recorded as `CRASH`/`TIMEOUT` and FreeCAD is relaunched for the rest.

## Options (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `AT_RESULTS` | `autotest_results/` | results root (git-ignore it) |
| `AT_RUN_ID` | timestamp | reuse an id to resume an interrupted run |
| `AT_STAGE_TIMEOUT` | 900 | seconds before a stuck stage is killed |
| `AT_STARTUP_TIMEOUT` | 300 | seconds for FreeCAD to start the runner |
| `AT_MAX_RESTARTS` | 25 | relaunches per version after crashes |
| `AT_STAGES` | `open,export,cog,import` | stages to run |
| `AT_ONLY` | | glob(s) restricting file names, comma-separated |
| `AT_GEANT4` | 0 | `1` = also load every export with `load_gdml_color` |
| `AT_VERSIONS_CONF` | `autotest/versions.conf` | alternative version list |
| `AT_HEADLESS_X` | 0 | Linux only: run under `xvfb-run` |
| `GDML_REPO` | `~/Workbenches/GDML` | used by `check_setup.sh --link` |

## Output

```
autotest_results/<run>/
  report.md  report.json  run.log
  <version>/env.json          FreeCAD, Python, OCC, GDML commit, XML backend
  <version>/results.jsonl     one record per file, with tracebacks
  <version>/exported_gdml/    <stem>-<worldVol>.gdml  (same names as FC_exported_gdml)
  <version>/cogs/             <stem>.txt              (same format as FC_COGs)
  <version>/roundtrip_gdml/
  <version>/dialogs.log  freecad.log  geant4.jsonl
```

Statuses: `OK`, `OK_WITH_DIALOGS`, `FAIL_OPEN`, `FAIL_EXPORT`, `FAIL_COG`,
`FAIL_IMPORT`, `FAIL_ROUNDTRIP`, `CRASH`, `TIMEOUT`.

COG tolerances (in `compare.py`): volume 0.1 %, centre of mass 0.1 % or
0.02 mm, inertia 1 % (printed to 4 significant figures, may be Monte Carlo).
GDML exports are compared structurally — solid types and counts, facets,
materials, volumes, assemblies, physvols — not as text.

## Geant4 stage

`AT_GEANT4=1` (or `autotest/run_geant4.sh <run_dir>` afterwards) loads every
exported GDML with `load_gdml_color -printCM -m=batch.mac`, records
`OK / FAILED / ABORT / TIMEOUT` per file in `geant4.jsonl`, and the report
flags any that Geant4 rejects. Source the Geant4 environment first so
`load_gdml_color` is on `PATH`.

Comparing Geant4's centre of mass with FreeCAD's needs
`filter_load_gdml_CM.awk`, which `geant4COG.sh` uses but is not in the repo.

## Scheduling overnight

Any scheduler that can run a shell command works, e.g. `launchd`, or a
Claude scheduled task on a machine linked to the Claude desktop app:

```sh
cd ~/Workbenches/FreeCAD_GDML_Test_Collection && git -C ../GDML pull --ff-only && autotest/run_versions.sh
```

`report.json` is meant for an agent to read the next morning and triage.
