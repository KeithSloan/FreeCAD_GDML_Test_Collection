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
/// \file persistency/gdml/G01/load_gdml.cc
/// \brief Main program of the persistency/gdml/G01 example
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - load_gdml
//
// --------------------------------------------------------------

#include <G4String.hh>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iterator>
#include <ostream>
#include <string>
#include <vector>
#include <sstream>

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"

//#include "G01PrimaryGeneratorAction.hh"
#include "G01DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"

#include "G4VSolid.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4GDMLParser.hh"
#include "G4ios.hh"

#include "GeometryMoments.hh"
#include "MomentsEstimator.hh"

using namespace std;
using namespace CLHEP;


void print_aux(const G4GDMLAuxListType* auxInfoList, G4String prepend="|")
{
  for(std::vector<G4GDMLAuxStructType>::const_iterator
      iaux = auxInfoList->begin(); iaux != auxInfoList->end(); iaux++ )
    {
      G4String str=iaux->type;
      G4String val=iaux->value;
      G4String unit=iaux->unit;

      G4cout << prepend << str << " : " << val  << " " << unit << G4endl;

      if (iaux->auxList) print_aux(iaux->auxList, prepend + "|");
    }
  return;
}

void addColor(G4String lvName, G4String color){
  std::stringstream ss;
  unsigned int r, g, b, a;

  std::string rs, gs, bs, as;
  rs = "0x"+color.substr(1,2);
  gs = "0x"+color.substr(3,2);
  bs = "0x"+color.substr(5,2);
  sscanf(rs.c_str(), "%x", &r);
  sscanf(gs.c_str(), "%x", &g);
  sscanf(bs.c_str(), "%x", &b);
  
  a = 256;
  if(color.length() == 9) { // includes the #
	as = "0x"+color.substr(7,2);
	sscanf(as.c_str(), "%x", &a);
  }

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  char cmd[1024];
  sprintf(cmd, "/vis/geometry/set/colour %s %f %f %f %f",
		  lvName.c_str(), r/256., g/256., b/256., 1.0 - a/256.);
  G4cout << cmd << G4endl;
  UImanager->ApplyCommand(cmd);
  
}

void usage()
{
   G4cout << G4endl;
   G4cout << left << "Usage: load_gdml <intput_gdml_file:mandatory> \n"
		  << setw(20) << "Optional arguments: \n"
          << setw(30) << "\t-o=<output_gdml_file>" <<  "Optional output gdml file\n"
		  << setw(30) << "\t-m=<macroFile>" << "Macro file to execute, Default initInter.mac\n"
		  << setw(30) << "\t-help:" << "print usage\n"
		  << setw(30) << "\t-printCM" << "print center of mass and moments of inertia for each placed volume.\n"
		  << setw(30) << "\t" << "The CM and a moments are culculated using Monte Carlo\n"
		  << setw(30) << "\t" << "The moments of inertia are calculated about the CM of each object\n"
		  << setw(30) << "\t" <<"and parallel to the World x, y and z axes\n"
		  << setw(30) << "\t-logNsim=[3|4|5|6|7|8]" << "log10(Nsim) of number of simulations in CM Monte Carlo\n"
		  << setw(30) << "\t" << "=3 ==> 1000 simulations, =4 ==>10000 simulations, etc.\n"
		  << setw(30) << "\t" << "default=6 (1000,000) simulations per volume.\n";
   G4cout << G4endl;
}

// return index of given argument. -1 if not found
int findArg(int argc, char **argv, G4String arg)
{
  arg = "-" + arg;
  for(int i=1; i < argc; i++) {
	G4String carg = argv[i];
	auto ieq = carg.find("=");
	if ( ieq != string::npos) {
	  carg =  carg.substr(0, ieq);
	}
	if(carg == arg) {
	  return i;
	}
  }
  return -1;
}

int getIntArg(char **argv, int iarg)
{
  string arg = argv[iarg];
  auto ieq = arg.find("=");
  int val = atoi(arg.substr(ieq+1).c_str());

  return val;
}

G4String getArg(char **argv, int iarg)
{
  string arg = argv[iarg];
  auto ieq = arg.find("=");
  G4String val = arg.substr(ieq+1);

  return val;
}

vector<G4String> validateOptions(int argc, char **argv, G4String opts[], int nopts)
{
  vector<G4String> invalidArgs;
  for(int i=1; i < argc; i++) {
	G4String  arg = argv[i];
	bool found = false;
	if (arg.find("-") == 0) {
	  for(int j=0; j < nopts; j++) {
		if(arg.find(opts[j]) != string::npos) {
		  found = true;
		}
	  }
	  if(!found) {
		invalidArgs.push_back(arg);
	  }
	}
  }
  return invalidArgs;
}

// --------------------------------------------------------------

int main(int argc,char **argv)
{
   if (argc < 2)
   {
      G4cerr << "Error! Mandatory input file is not specified!" << G4endl;
      G4cerr << G4endl;
	  usage();
      return -1;
   }

   if (argc > 7)
   {
      G4cerr << "Error! Too many arguments!" << G4endl;
      G4cerr << G4endl;
	  usage();
      return -1;
   }

   G4String options [] = {"-o=", "-m=", "-help", "-printCM", "logNsim="};
   vector<G4String> invalidOptions = validateOptions(argc, argv, options, 5);
   if (invalidOptions.size() > 0) {
	 for(size_t i=0; i < invalidOptions.size(); i++) {
	   G4cerr << "Invalid option: " << invalidOptions[i] << endl;
	 }
	 usage();
	 return -1;
   }

   if(findArg(argc, argv, "help") > 0) {
	 usage();
	 return -1;
   }

   G4GDMLParser parser;

// Uncomment the following if wish to avoid names stripping
// parser.SetStripFlag(false);

   parser.SetOverlapCheck(true);
   parser.Read(argv[1]);


   auto* runManager = G4RunManagerFactory::CreateRunManager();

   runManager->SetUserInitialization(new G01DetectorConstruction(
                                     parser.GetWorldVolume()));
   G4VModularPhysicsList* physicsList = new FTFP_BERT;
   physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
   G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
   
   physicsList->RegisterPhysics(opticalPhysics);
   runManager->SetUserInitialization(physicsList);
   runManager->SetUserInitialization(new ActionInitialization());

   runManager->Initialize();

   // Initialize visualization
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();

   // Get the pointer to the User Interface manager
   G4UImanager* UImanager = G4UImanager::GetUIpointer();

   ///////////////////////////////////////////////////////////////////////
   //
   // Example how to retrieve Auxiliary Information
   //

   G4cout << std::endl;

   const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
   std::vector<G4LogicalVolume*>::const_iterator lvciter;
   for( lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++ )
   {
     G4GDMLAuxListType auxInfo = parser.GetVolumeAuxiliaryInformation(*lvciter);

     if (auxInfo.size()>0) {
       G4cout << "Auxiliary Information is found for Logical Volume :  "
              << (*lvciter)->GetName() << G4endl;
	   for(std::vector<G4GDMLAuxStructType>::const_iterator
			 iaux = auxInfo.begin(); iaux != auxInfo.end(); iaux++ )
		 {
		   G4String str=iaux->type;
		   G4String val=iaux->value;
		   G4String unit=iaux->unit;

		   if(str == "Color") {
			 addColor((*lvciter)->GetName(), val);
		   }
		 }
	 }

     print_aux(&auxInfo);
   }

   // now the 'global' auxiliary info
   G4cout << std::endl;
   G4cout << "Global auxiliary info:" << std::endl;
   G4cout << std::endl;

   print_aux(parser.GetAuxList());

   G4cout << std::endl;

   //
   // End of Auxiliary Information block
   //
   ////////////////////////////////////////////////////////////////////////


   // runManager->BeamOn(0);

   // example of writing out

   int iOutput = findArg(argc, argv, "o");
   if (iOutput > 0) {
	 G4String outfile = getArg(argv, iOutput);
	 
/*
     G4GDMLAuxStructType mysubaux = {"mysubtype", "mysubvalue", "mysubunit", 0};
     G4GDMLAuxListType* myauxlist = new G4GDMLAuxListType();
     myauxlist->push_back(mysubaux);

     G4GDMLAuxStructType myaux = {"mytype", "myvalue", "myunit", myauxlist};
     parser.AddAuxiliary(myaux);


     // example of setting auxiliary info for world volume
     // (can be set for any volume)

     G4GDMLAuxStructType mylocalaux = {"sometype", "somevalue", "someunit", 0};

     parser.AddVolumeAuxiliary(mylocalaux,
       G4TransportationManager::GetTransportationManager()
       ->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
*/

     parser.SetRegionExport(true);
     //     parser.SetEnergyCutsExport(true);
     //     parser.SetOutputFileOverwrite(true);
     parser.Write(outfile, G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
   }

   int iprintCM = findArg(argc, argv, "printCM");
   G4cout << "** iprintCM= " << iprintCM << endl;
   int nsim = 1000000;
   if (iprintCM > 0) {
	 int insim = findArg(argc, argv, "logNsim");
	 if (insim > 0) {
	   int logn = getIntArg(argv, insim);
	   if(logn > 2 && logn < 9) {
		 nsim = 10;
		 int i = 1;
		 while(i++ < logn) {
		   nsim *= 10;
		 }
	   }
	   else {
		 G4cout << "logNsim outside allowed range ([3,8], inclusive)\n"
				<< "default of 1000,000 simulatios used" << endl;
	   }
	 }
	 G4cout << "using " << nsim << " Simulations per volume for CM MC" << endl;
   }

   MomentsEstimator::nsim = nsim;

   int iBatchFile = findArg(argc, argv, "m");
   G4String command ="/control/macroPath .";
   UImanager->ApplyCommand(command);
   if (iBatchFile > 0) {// batch mode
     G4String command = "/control/execute ";
     G4String fileName = getArg(argv, iBatchFile);
     UImanager->ApplyCommand(command+fileName);
	 if(iprintCM > 0) {
	   printGeometryMoments();
	 }
   }
   else           // interactive mode
   {
     G4UIExecutive* ui = new G4UIExecutive(argc, argv);
     UImanager->ApplyCommand("/control/execute initInter.mac");
	 if(iprintCM > 0) {
	   printGeometryMoments();
	 }
     ui->SessionStart();
     delete ui;
   }

   delete visManager;
   delete runManager;

   return 0;
}

