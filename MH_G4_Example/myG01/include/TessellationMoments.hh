#ifndef TESSELLATION_MOMENTS_CALCULATOR_HH
#define TESSELLATION_MOMENTS_CALCULATOR_HH

#include "InertiaMatrix.hh"
#include "MomentsCalculator.hh"
#include "G4TessellatedSolid.hh"
#include "G4ThreeVector.hh"

Moments computeTriangleInertia(const G4ThreeVector &v1,
							   const G4ThreeVector &v2,
							   const G4ThreeVector &v3);

class TessellationMoments: public MomentsCalculator
{
public:
  inline TessellationMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~TessellationMoments() {};
  
  Moments calculate();
};

#endif
