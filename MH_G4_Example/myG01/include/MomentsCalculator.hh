#ifndef MOMENTsCALCULATOR_H
#define MOMENTsCALCULATOR_H

#include "Moments.hh"
#include <G4Ellipsoid.hh>
#include <G4TwoVector.hh>
#include <G4VSolid.hh>

double polygonArea(const std::vector<G4TwoVector> vertices);
InertiaMatrix computeInertiaMatrix(std::vector<G4TwoVector> vertices);
G4TwoVector polygonCM(std::vector<G4TwoVector> vertices);

class MomentsCalculator
{
public:
  static MomentsCalculator* getCalculator(G4VSolid*);

  inline MomentsCalculator(G4VSolid *solid) {f_solid = solid; };
  virtual ~MomentsCalculator() {};

  virtual Moments calculate() = 0;

  bool isBooleanSolid();


protected:
  G4VSolid* f_solid;

};

class BoxMoments: public MomentsCalculator
{
public:
  inline BoxMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~BoxMoments() {};
  
  Moments calculate();
};

class EllipsoidMoments: public MomentsCalculator
{
public:
  inline EllipsoidMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~EllipsoidMoments() {};
					  
  Moments calculate();
};

class OrbMoments: public MomentsCalculator
{
public:
  inline OrbMoments(G4VSolid *solid) : MomentsCalculator(solid) {};

  ~OrbMoments() {};
  
  Moments calculate();

  
};

class SphereMoments: public MomentsCalculator
{
public:
  inline SphereMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~ SphereMoments() {};
  
  Moments calculate();

  
};

class ConeMoments: public MomentsCalculator
{
public:
  inline ConeMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~ConeMoments() {};
  
  Moments calculate();

  
};

class TubeMoments: public MomentsCalculator
{
public:
  inline TubeMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~TubeMoments() {};
  
  Moments calculate();

  
};

class TrdMoments: public MomentsCalculator
{
public:
  inline TrdMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~TrdMoments() {};
  
  Moments calculate();

  
};

class TorusMoments: public MomentsCalculator
{
public:
  inline TorusMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~TorusMoments() {};
  
  Moments calculate();

  
};

class ParaboloidMoments: public MomentsCalculator
{
public:
  inline ParaboloidMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~ParaboloidMoments() {};
  
  Moments calculate();

  
};

class HypeMoments: public MomentsCalculator
{
public:
  inline HypeMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~HypeMoments() {};
  
  Moments calculate();

};

class TetrahedronMoments: public MomentsCalculator
{
public:
  inline TetrahedronMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~TetrahedronMoments() {};
  
  Moments calculate();

};

class PolyconeMoments: public MomentsCalculator
{
public:
  inline PolyconeMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~PolyconeMoments() {};
  
  Moments calculate();

};

class EllipticalTubeMoments: public MomentsCalculator
{
public:
  inline EllipticalTubeMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~EllipticalTubeMoments() {};
  
  Moments calculate();

};

class EllipticalConeMoments: public MomentsCalculator
{
public:
  inline EllipticalConeMoments(G4VSolid *solid) : MomentsCalculator(solid) {};
  ~EllipticalConeMoments() {};
  
  Moments calculate();

};


#endif
