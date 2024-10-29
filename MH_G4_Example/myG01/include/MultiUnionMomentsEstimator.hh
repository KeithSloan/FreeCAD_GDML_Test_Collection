#ifndef MULTIUNIONMOMENTSESTIMATOR_H
#define MULTIUNIONMOMENTSESTIMATOR_H

#include "MomentsEstimator.hh"
#include <G4MultiUnion.hh>

class MultiUnionMomentsEstimator: public MomentsEstimator
{
public:
  inline MultiUnionMomentsEstimator(G4VSolid *solid): MomentsEstimator(solid) {};
  ~MultiUnionMomentsEstimator() {};
							  
  Moments estimate() override;
	
};

#endif
