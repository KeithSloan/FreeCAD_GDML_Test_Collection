//#include <G4String.hh>
#include "G4GDMLParser.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include <iostream>
//#include <sstream>
#include "G4Types.hh"
#include "G4ios.hh"

//using namespace std;



void usage()
{
   //G4cout << G4endl;
   //G4cout << "Usage: load_gdml <intput_gdml_file:mandatory> \n"
   //G4cout << G4endl;
}


int main(int argc, char **argv)

{
   if (argc < 2)
   {
      G4cerr << "Error! Mandatory input file is not specified!" << G4endl;
      G4cerr << G4endl;
      usage();
      return -1;
   }

   if (argc > 3)
   {
      G4cerr << "Error! Too many arguments!" << G4endl;
      G4cerr << G4endl;
      usage();
      return -1;
   }

   //#G4String options [] = {"-o=", "-m=", "-help", "-printCM", "logNsim="};
   //#vector<G4String> invalidOptions = validateOptions(argc, argv, options, 5);
   //#if (invalidOptions.size() > 0)
   //# for(size_t i=0; i < invalidOptions.size(); i++) {
   //#    G4cerr << "Invalid option: " << invalidOptions[i] << endl;
   //#  }
   //#  usage();
   //#  return -1;
   //}

   //if(findArg(argc, argv, "help") > 0) {
   //  usage();
   //  return -1;
   }

    // Create the GDML parser
    G4GDMLParser parser;

    // File to parse
    //G4String gdmlFile = "geometry.gdml";

    //try {
    // Parse the GDML file
    parser.SetOverlapCheck(true);
    //parser.Read(gdmlFile);
    parser.Read(argv[1]);
    // If successful, get the world volume (or root volume)
    G4VPhysicalVolume* worldVolume = parser.GetWorldVolume();
        
    if (worldVolume) {
       std::cout << "GDML file parsed successfully. World volume is set." << std::endl;
    } else {
       std::cerr << "Failed to get world volume from the parsed GDML file." << std::endl;
        return 1;  // Return code to indicate failure
        }
    } catch (const std::exception& e) {
        // Handle exceptions (GDML parsing errors)
        std::cerr << "Error parsing GDML file: " << e.what() << std::endl;
        return 1;  // Return code to indicate failure
    }

    // Continue with the rest of your Geant4 setup if parsing succeeds
    // ...
    
    return 0;  // Return 0 to indicate success
}
