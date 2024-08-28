#include "GeometryMoments.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4MultiUnion.hh"
#include "G4Transform3D.hh"
#include "G4ios.hh"

#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4AffineTransform.hh>
#include <G4RotationMatrix.hh>
#include <G4Sphere.hh>
#include <G4String.hh>
#include <G4Threading.hh>
#include <G4VMarker.hh>
#include <G4VSolid.hh>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <ios>
#include <utility>

#include "InertiaMatrix.hh"

#include <cmath>
#include <vector>

#include "G4LogicalVolumeStore.hh"
#include "Moments.hh"
#include "MomentsEstimator.hh"

using namespace std;
using namespace CLHEP;

void prettyPrint(G4String name, double volume,
				 G4ThreeVector cm0, InertiaMatrix M0,
				 G4String entityType) {
  G4cout << std::setprecision(2);

  G4String cmunit = "mm";
  if (cm0.mag() >= 1000.0*CLHEP::mm) {
	cm0 /= CLHEP::m;
	cmunit = " m";
  }
  else if( cm0.mag() >= 100.0*CLHEP::mm) {
	cm0 /= CLHEP::cm;
	cmunit = "cm";
  }
  
  G4String volunit = "mm^3";
  if (volume >= 1e9*CLHEP::mm3) {
	volume /= CLHEP::m3;
	volunit = " m^3";
  }
  else if (volume >= 1000*CLHEP::mm3) {
	volume /= CLHEP::cm3;
	volunit = "cm^3";
  }
  

  // inertia unit: either powers of cm^5, or powers of mm^5
  // If maxInertia > 0.1 cm^5, unit is cm^5
  // otherwise unit is mm^5
  // 0.1 cm^5 = 0.1*(10^5 mm^5) = 10^4 mm^5

  double maxInertia = M0.max();  // this is in base unit of freecad = mm^5
  G4String inertiaUnit;
  if (log10(maxInertia) >= 4){ 
	inertiaUnit = "cm^5";
	maxInertia /= 1e5;
  }
  else {
	inertiaUnit = "mm^5";
  }

  int k = 0;
  while (maxInertia > 10) {
	k += 1;
	maxInertia /= 10;
  }

  int power;
  if (inertiaUnit == "mm^5"){
	power = k;
  }
  else {
	power = 5 + k;
  }

  if (k > 0) {
	inertiaUnit = "x10^" + to_string(k) + " " + inertiaUnit;
  }
		  

  M0 /= pow(10, power);
    
  G4cout << left << setw(20) << name;
  G4cout << " Volume: " << right << setw(8) << volume << left << " " << volunit;
  G4cout << " --> center of \"mass\": ("
		 << right << setw(7) << cm0[0] << ","
		 << right << setw(7) << cm0[1] << ","
		 << right << setw(7) << cm0[2] << ")"
		 << left << " " << cmunit << ", ";

  char buffer[256];
  int n = sprintf(buffer, "moments of inertia: (%6.4g,%6.4g,%6.4g)", M0.xx(), M0.yy(), M0.zz());
  G4cout << buffer;
  
  /*
  G4cout << "moments of inertia: ( "
	"(" << right << setw(7) << M0.xx() << ","
		 << right << setw(7) << M0.xy() << ","
		 << right << setw(7) << M0.xz() << "), ";
  
  G4cout << 
	"(" << right << setw(7) << M0.yx() << ","
		 << right << setw(7) << M0.yy() << ","
		 << right << setw(7) << M0.yz()  << "), ";
  G4cout << 
	"(" << right << setw(7) << M0.zx() << ","
		 << right << setw(7) << M0.zy() << ","
		 << right << setw(7) << M0.zz() << ") )";
  
  */
  G4cout << " " << left << inertiaUnit; 

  if (entityType != "") {
	G4cout << " Entity type: " << entityType;
  }
  G4cout << endl;
}

void printLogicalVolume(G4LogicalVolume *lv, int level)
{
  for(int i=0; i < 4*level; i++) {
	G4cout << " ";
  }
  G4cout << lv->GetName();

  for(size_t i=0; i < lv->GetNoDaughters(); i++) {
	G4VPhysicalVolume *pv = lv->GetDaughter(i);
	G4LogicalVolume* daughterlv = pv->GetLogicalVolume();
	G4LogicalVolume* motherlv = pv->GetMotherLogical();
	if(motherlv != nullptr) {
	  G4cout << " mother: " << motherlv->GetName();
	}
	else {
	  G4cout << " No Mother";
	}
	
	if(daughterlv == nullptr) {
	  G4cout << " No daughter";
	  continue; 
	}
	else {
	  printLogicalVolume(daughterlv, level+1);
	}
  }
  G4cout << endl;
}

void printGeometryTree()
{
  /*
  G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume*>::iterator it;
  for(it = lvs->begin(); it != lvs->end(); it++) {
	
	printLogicalVolume((*it), 0);
  }
  */

  G4PhysicalVolumeStore *pvs = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume *>::iterator it;
  for(it = pvs->begin(); it != pvs->end(); it++) {
	G4String mother="none";
	if((*it)->GetMotherLogical() != nullptr) {
	  mother = (*it)->GetMotherLogical()->GetName();
	}
	G4String logName="none";
	if((*it)->GetLogicalVolume() != nullptr) {
	  logName = (*it)->GetLogicalVolume()->GetName();
	}

	G4cout << "physVol: "<< (*it)->GetName() << " logVol: " << logName << " mother: " << mother << endl;
  }

  G4cout << endl;
  
  G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume*>::iterator it1;
  for(it1 = lvs->begin(); it1 != lvs->end(); it1++) {
	G4cout << "logVol: "<< (*it1)->GetName() << " phys daughters: ";

	for(size_t i=0; i < (*it1)->GetNoDaughters(); i++) {
	  G4VPhysicalVolume *pv = (*it1)->GetDaughter(i);
	  G4cout << pv->GetName() << " ";
	}
	G4cout << endl;
  }

}

void getWorldMoments(G4VPhysicalVolume* pv, Moments& mom, G4AffineTransform parentTransform) {
  G4ThreeVector translation = pv->GetTranslation();
  G4RotationMatrix* rotation = pv->GetObjectRotation();

  G4AffineTransform thistransform(rotation, translation);
  auto transform = parentTransform*thistransform;
  
  G4ThreeVector cm0 = mom.cm;
  if (transform.IsRotated()) {
	// Rotate matrix of inertia
	G4RotationMatrix  R = transform.NetRotation();
	cm0 = R * mom.cm;
	mom.M = mom.M.rotated(R);	
  }

  cm0 += transform.NetTranslation();
  mom.cm = cm0;
  mom.M = mom.inertiaAboutOrigin();

}

std::vector<Moments> getMoments(G4VPhysicalVolume *pv, G4AffineTransform parentTransform) {

  std::vector<Moments> moments;
  
  G4LogicalVolume* lv0 = pv->GetLogicalVolume();
  G4LogicalVolume* lvmother = pv->GetMotherLogical();

  if (lvmother != nullptr) {
	G4VSolid* solid = lv0->GetSolid();
	// auto MCEstimator = MomentsEstimator(solid);
	auto MCEstimator = MomentsEstimator::getEstimator(solid);
	Moments mom0 = MCEstimator.estimate();
	getWorldMoments(pv, mom0, parentTransform);
  
	prettyPrint(pv->GetName(), mom0.Volume, mom0.cm, mom0.M, solid->GetEntityType()); 
	moments.push_back(mom0);
  }
  
  G4ThreeVector translation = pv->GetTranslation();
  G4RotationMatrix* rotation = pv->GetObjectRotation();

  G4AffineTransform thistransform(rotation, translation);

  for(size_t i=0; i < lv0->GetNoDaughters(); i++) {
	G4VPhysicalVolume *pvd = lv0->GetDaughter(i);
	auto daugherMoments = getMoments(pvd, thistransform);
	moments.insert(moments.end(), daugherMoments.begin(), daugherMoments.end());
  }
  
  return moments;
}


G4VPhysicalVolume* getWorldVolume()
{
  // Assumption is ther is only one volume that has no parent and that is the world volume
  G4PhysicalVolumeStore *pvs = G4PhysicalVolumeStore::GetInstance();
  map<G4String, std::vector<G4VPhysicalVolume*>> pvmap = pvs->GetMap();
  map<G4String, std::vector<G4VPhysicalVolume*>>::iterator it;
  for (it = pvmap.begin(); it != pvmap.end(); it++) {
	vector<G4VPhysicalVolume*> pv = it->second;
	vector<G4VPhysicalVolume*>::iterator itpv;

	for(itpv = pv.begin(); itpv != pv.end(); itpv++) {
	  G4LogicalVolume *lv = (*itpv)->GetLogicalVolume();
	  G4LogicalVolume* motherlv = (*itpv)->GetMotherLogical();
	  if(motherlv == nullptr) {
		G4cout << "World Logical Volume: " << lv->GetName() << endl;
		return (*itpv); 
	  }

	}
  }

  return nullptr;
}

void printGeometryMoments()
{

  // printGeometryTree();
  
  G4cout << std::fixed;
  G4cout << std::setprecision(2);
  G4cout << "***** geometry moments ****" << endl;
  
  G4PhysicalVolumeStore *pvs = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume *worldVol = getWorldVolume();
  G4AffineTransform worldTransform(G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0));

  std::vector<Moments> moms = getMoments(worldVol, worldTransform);
  
  double totVol = 0;
  G4ThreeVector cm = G4ThreeVector(0, 0, 0);  // center of gravity
  InertiaMatrix II = InertiaMatrix();
  
  for(size_t i=0; i < moms.size(); i++) {
	auto mom = moms[i];
	cm += mom.Volume * mom.cm;
	II += mom.M;
	totVol += mom.Volume;
  }
  
  G4cout << endl;
  prettyPrint("Total ", totVol, cm/totVol, II, "");
  G4cout << "***** geometry moments ****" << endl;
}



