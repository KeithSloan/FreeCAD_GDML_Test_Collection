#include "MomentsEstimator.hh"
#include "MultiUnionMomentsEstimator.hh"

#include "Moments.hh"
#include <CLHEP/Random/RandFlat.h>
#include <G4ios.hh>

using namespace CLHEP;

int MomentsEstimator::nsim = 1000000;

MomentsEstimator MomentsEstimator::getEstimator(G4VSolid *solid) {
  if (solid->GetEntityType() == "G4MultiUnion") {
	return MultiUnionMomentsEstimator(solid);
  }
  else {
	return MomentsEstimator(solid);
  }
}

bool MomentsEstimator::isBooleanSolid()
{
  G4String entityType = f_solid->GetEntityType();

  return entityType == "G4SubtractionSolid" ||
	entityType == "G4IntersectionSolid" ||
	entityType == "G4UnionSolid" ||
	entityType == "G4MultiUnion";
}


Moments MomentsEstimator::estimate()
{
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4ThreeVector pmin;
  G4ThreeVector pmax;
  f_solid->BoundingLimits(pmin, pmax);
  
  G4ThreeVector cm = G4ThreeVector(0, 0, 0);
  double Ixx = 0;
  double Ixy = 0;
  double Ixz = 0;
  double Iyy = 0;
  double Iyz = 0;
  double Izz = 0;
  
  int numIn = 0;
  for (int i=0; i < nsim; i++) {
	double x = RandFlat::shoot(pmin.x(), pmax.x());
	double y = RandFlat::shoot(pmin.y(), pmax.y());
	double z = RandFlat::shoot(pmin.z(), pmax.z());
	auto r = G4ThreeVector(x, y, z);
	EInside inside = f_solid->Inside(r);
	if (inside != kOutside) {
	  numIn++;
	  cm += r;

	  Ixx += z*z + y*y;
	  Ixy += -x*y;
	  Ixz += -x*z;
	  Iyy += x*x + z*z;
	  Iyz += -y*z;
	  Izz += x*x + y*y;
	}
  }

  // calculate moments about center of mass, using parallel axis
  // theorem.
  cm /= numIn;
  Ixx -= numIn*(cm.y()*cm.y() + cm.z()*cm.z());
  Iyy -= numIn*(cm.x()*cm.x() + cm.z()*cm.z());
  Izz -= numIn*(cm.x()*cm.x() + cm.y()*cm.y());
  Ixy += numIn*cm.x()*cm.y();
  Ixz += numIn*cm.x()*cm.z();
  Iyz += numIn*cm.y()*cm.z();


  moments.cm = cm;
  moments.M = InertiaMatrix(Ixx, Ixy, Ixz,
							Ixy, Iyy, Iyz,
							Ixz, Iyz, Izz);

  moments.M /= numIn;
  
  double MCVolume = (pmax.x() - pmin.x()) * (pmax.y() - pmin.y()) * (pmax.z() - pmin.z());
  MCVolume *= (double) numIn/nsim;
  moments.M *= MCVolume;
  if (isBooleanSolid()) {
	moments.Volume = MCVolume;
  }
  else {
	moments.Volume = f_solid->GetCubicVolume();
  }

  return moments;
}
