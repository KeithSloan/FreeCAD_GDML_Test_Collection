#!/bin/sh
#echo $PATH
# printCM is needed for awk script
cd FC_exported_gdml
LOGFILE="../geant4_failures.txt"
#cp ../vis.mac vis.mac
#cp ../initInit.mac initInit.mac
cp ../batch.mac batch.mac
for f in *.gdml; do
    load_gdml_color $f -printCM -m=batch.mac >${f}.out 2>&1;
    retCode=$?
    if [ $retCode -ne 0 ]; then
        echo $f "Failed " $retCode >> $LOGFILE
    fi
done    
for f in *.out; do 
    gawk -f ../filter_load_gdml_CM.awk "$f" >"${f%}.txt"
done
exit
    
