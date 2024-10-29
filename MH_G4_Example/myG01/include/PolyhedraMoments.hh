#ifndef POLYHEDRA_MOMENTS_CALCULATOR_HH
#define POLYHEDRA_MOMENTS_CALCULATOR_HH

#include "MomentsCalculator.hh"
#include "G4ThreeVector.hh"

Moments computeTriangleInertia(const G4ThreeVector &v1,
							   const G4ThreeVector &v2,
							   const G4ThreeVector &v3);

class PolyhedraMoments: public MomentsCalculator
{
public:
  inline PolyhedraMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~PolyhedraMoments() {};
  
  Moments calculate();

private:
  // return secton moment for the section between z(iz) and z(iz+1)
  Moments iMoment(int iz);

};

#endif
