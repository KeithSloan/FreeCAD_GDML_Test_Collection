#include "myPhysicsList.hh"
#include <G4Gamma.hh>
#include <G4Geantino.hh>
#include <G4OpticalPhoton.hh>
#include <G4PhysicsListHelper.hh>
#include <G4VUserPhysicsList.hh>
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh"
// e-
// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

myPhysicsList::myPhysicsList() : G4VUserPhysicsList() {
  
}

myPhysicsList::~myPhysicsList()
{}


void myPhysicsList::ConstructProcess()
{
  AddTransportation();
  // ConstructEM();
}

void myPhysicsList::ConstructParticle()
{
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Geantino::GeantinoDefinition();
  G4Gamma::GammaDefinition();
  G4OpticalPhoton::OpticalPhotonDefinition();
 
}

void myPhysicsList::ConstructEM()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  //  Get pointer to gamma
  G4ParticleDefinition* particle = G4Gamma::GammaDefinition();

  // Construct and register processes for gamma
  ph->RegisterProcess(new G4PhotoElectricEffect(), particle);
  ph->RegisterProcess(new G4ComptonScattering(), particle);
  ph->RegisterProcess(new G4GammaConversion(), particle);
  ph->RegisterProcess(new G4RayleighScattering(), particle);
  
}

