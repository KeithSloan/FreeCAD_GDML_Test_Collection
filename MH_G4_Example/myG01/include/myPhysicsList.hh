#ifndef myPhysicsList_h
#define myPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class myPhysicsList: public G4VUserPhysicsList
{
public:
  myPhysicsList();
  ~myPhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void ConstructEM();

};


#endif
