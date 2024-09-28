#include "G4GDMLParser.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include <iostream>

int main() {
    // Create the GDML parser
    G4GDMLParser parser;

    // File to parse
    G4String gdmlFile = "geometry.gdml";

    try {
        // Parse the GDML file
        parser.Read(gdmlFile);
        
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

