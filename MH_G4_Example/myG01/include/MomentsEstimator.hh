#ifndef MOMENTSESTIMATOR_H
#define MOMENTSESTIMATOR_H

#include "Moments.hh"
#include <G4VSolid.hh>

class MomentsEstimator
{
public:
  static int nsim;
  static MomentsEstimator* getEstimator(G4VSolid*);

  inline MomentsEstimator(G4VSolid *solid) {f_solid = solid; };
  virtual ~MomentsEstimator() {};
  
  virtual Moments estimate();

  bool isBooleanSolid();
  
  
protected:
  G4VSolid* f_solid;

	
};

#endif
