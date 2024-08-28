#ifndef MULTIUNIONMOMENTSESTIMATOR_H
#define MULTIUNIONMOMENTSESTIMATOR_H

#include "MomentsEstimator.hh"

class MultiUnionMomentsEstimator: public MomentsEstimator
{
public:
  inline MultiUnionMomentsEstimator(G4VSolid *solid): MomentsEstimator(solid) {};
  Moments estimate() override;
	
};

#endif
