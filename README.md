# FreeCAD_GDML_Test_Collection
A set of FreeCAD files created with the GDML workbench - In Future will be used for Auto Testing.

# Macros

  * GDMLExportDirectory.FCMacro 


      - Creates gdml exports in target directory for all FC files created with GDML workbench in source directory
  * GDMLImportDirectory.FCMacro
      - For all gdml files in a directory import them into a FC document, then in a target directory create a saved FC file and an exported gdml version.
The newsly exported gdml file can be checked against the original gdml file.     
  * calcCenterOfMass.py
      - Creates CofM and Moment of Inertia Calculations for selected World Volume
  * GDMLCOGcalcDirectory.py
     - Perform CofM and Moment of Inertia Calcs for all FC files in a directory 

# Notes: 

## FreeCAD files should have been recomputed and saved

## For Geant4 COG

  load_gdml_color -m batchFile.mac

  ( batchFile.mac same as initInit.mac minus /control/execture vis.mac )

# Run a test

  * create a branch

        git branch testDate
        git checkout testDate

  * run FC Macro

        GDMLExportDirectory
             Source : test_FC_files
             Target : FC_exported_gdml

  * run FC Macro

        GDMLCOGcalcDirectory
            Source : test_FC_files
            Target : FC_COGs

  * create outputs for Geant4 COG's 

    cd  : FC_exported_gdml
```

    for f in *.gdml; do load_gdml_color $f -m=batch.mac -printCM >${f}.out; done
```
  * TODO run comparison script

            Sources : FC_COGs and Geant4_COGs
                                


