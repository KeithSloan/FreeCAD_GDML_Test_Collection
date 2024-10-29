#ifndef EXTRUSION_MOMENTS_CALCULATOR_HH
#define EXTRUSION_MOMENTS_CALCULATOR_HH

#include "InertiaMatrix.hh"
#include "MomentsCalculator.hh"
#include "G4ExtrudedSolid.hh"
#include <G4TwoVector.hh>

class ExtrusionMoments: public MomentsCalculator
{
public:
  inline ExtrusionMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~ExtrusionMoments() {};
  
  Moments calculate();

};


#endif
