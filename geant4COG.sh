#!/bin/sh

export Path=$PATH:.
#export Path=$PATH:$(pwd).
echo $PATH
cd FC_exported_gdml
#cp ../vis.mac vis.mac
#cp ../initInit.mac initInit.mac
#cp ../initInit.mac batch.mac
for f in *.gdml; do load_gdml_color $f -m=initInit.mac -printCM >${f}.out; done

gawk -f ../filter_gdml_CM.awk    
