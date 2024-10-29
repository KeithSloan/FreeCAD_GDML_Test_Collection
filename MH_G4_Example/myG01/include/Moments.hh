#ifndef MOMENTS_H
#define MOMENTS_H

#include "InertiaMatrix.hh"
#include "G4ThreeVector.hh"

struct Moments {
  double Volume;
  G4ThreeVector cm;
  InertiaMatrix M;  // inertia matrix about the center of mass

  InertiaMatrix inertiaAboutOrigin();
  InertiaMatrix inertiaAboutCM();
};


#endif
