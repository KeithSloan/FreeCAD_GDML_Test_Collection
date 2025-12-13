//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef G4FINDDATADIR_HH
#define G4FINDDATADIR_HH

#define CMAKE_INSTALL_PREFIX "/Users/ksloan/geant4-v11.3.2-full-install"
#define CMAKE_INSTALL_DATADIR "share/Geant4"
#define GEANT4_INSTALL_DATADIR "share/Geant4/data"
#define GEANT4_INSTALL_FULL_DATADIR "/Users/ksloan/geant4-v11.3.2-full-install/share/Geant4/data"

struct dataset {
  const char * const env;
  const char * const dir;
};

static const struct dataset dataset_definitions[] = {
{ "G4NEUTRONHPDATA", "G4NDL4.7.1" },
{ "G4LEDATA", "G4EMLOW8.6.1" },
{ "G4LEVELGAMMADATA", "PhotonEvaporation6.1" },
{ "G4RADIOACTIVEDATA", "RadioactiveDecay6.1.2" },
{ "G4PARTICLEXSDATA", "G4PARTICLEXS4.1" },
{ "G4PIIDATA", "G4PII1.3" },
{ "G4REALSURFACEDATA", "RealSurface2.2" },
{ "G4SAIDXSDATA", "G4SAIDDATA2.0" },
{ "G4ABLADATA", "G4ABLA3.3" },
{ "G4INCLDATA", "G4INCL1.2" },
{ "G4ENSDFSTATEDATA", "G4ENSDFSTATE3.0" },
{ "G4CHANNELINGDATA", "G4CHANNELING1.0" },};

#endif
