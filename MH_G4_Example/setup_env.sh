#
# MacOS Homebrew Geant4 + Root6 environment for tutorial
#

#qt5 on MacOS (brew)
export PKG_CONFIG_PATH=/usr/local/opt/qt@5/lib/pkgconfig
export PATH="/usr/local/opt/qt@5/bin:$PATH"

# Root 6
if [ -x /usr/local/bin/thisroot.sh ]; then
. /usr/local/bin/thisroot.sh
fi

# Geant4
if [ -x /usr/local/geant4/11.2.2/bin/geant4.sh ]; then
    . /usr/local/geant4/11.2.2/bin/geant4.sh
    # Need this if other Geant4 installed in /usr/local
    export CMAKE_PREFIX_PATH=/usr/local/bin/geant/11.2.2/lib/Geant4-v11.2.2
fi

#export G4PIIDATA=/opt/local/share/Geant4/Data/Geant4.11.2.2/G4PII1.3
#export G4PARTICLEXSDATA=/opt/local/share/Geant4/Data/Geant4.11.2.2/G4PARTICLEXS4.0
#export G4SAIDXSDATA=/opt/local/share/Geant4/Data/Geant4.11.2.2/G4SAIDDATA2.0
#export G4DATADIR=/opt/local/share/Geant4/Data/Geant4.11.2.2
#export G4REALSURFACEDATA=/opt/local/share/Geant4/Data/Geant4.11.2.2/RealSurface2.1.1
#export G4INCLDATA=/opt/local/share/Geant4/Data/Geant4.11.2.2/G4INCL1.0
#export G4TENDL=/opt/local/share/Geant4/Data/Geant4.11.2.2/G4TENDL1.3.2
#export G4EMLOW=/opt/local/share/Geant4/Data/Geant4.11.2.2/G4EMLOW.8.5
XPC_SERVICE_NAM=E0

export GEANT4_INSTALL_DATADIR=/opt/local/share/Geant4/Data/Geant4.11.2.2

#show enveronment in prompt
PS1="(g4) ${PS1}"
