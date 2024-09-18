#!/bin/sh
#iexport Path=$PATH:.
#export Path=$PATH:$(pwd).
#echo $PATH
cd FC_exported_gdml
touch ../geant4_failures.txt
#cp ../vis.mac vis.mac
#cp ../initInit.mac initInit.mac
cp ../batch.mac batch.mac
for f in *.gdml; do load_gdml_color $f -m=batch.mac >${f}.out; done
retVal=$?
if [ $retVal -ne 0 ]; then
    echo $f "Failed " $retval >> ../geant4_failures.txt
fi    
for f in *.out; do gawk -f ../filter_load_gdml_CM.awk $f; done
exit
    
