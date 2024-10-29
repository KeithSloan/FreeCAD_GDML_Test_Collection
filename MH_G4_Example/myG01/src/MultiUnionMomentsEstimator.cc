#include "MultiUnionMomentsEstimator.hh"

#include "G4MultiUnion.hh"
#include "Moments.hh"
#include "MomentsCalculator.hh"
#include "MultiUnionMomentsEstimator.hh"
#include <G4ThreeVector.hh>
#include <algorithm>

using namespace std;

typedef std::pair<G4ThreeVector, G4ThreeVector> bbox; // bounding box of a solid

bbox bounds(const G4MultiUnion* mu, int node)
{
  // Get transformed bounds of a node of a multiunion
	G4ThreeVector pmin, pmax;
	G4VSolid *solid = mu->GetSolid(node);
    G4Transform3D transform = mu->GetTransformation(node);
	solid->BoundingLimits(pmin, pmax);
    auto rot = transform.getRotation();
    auto translate = transform.getTranslation();
    pmin = rot * pmin;
    pmin = translate + pmin;
    pmax = rot * pmax;
    pmax = translate + pmax;

	double xmin = std::min(pmin[0], pmax[0]);
	double xmax = std::max(pmin[0], pmax[0]);
	double ymin = std::min(pmin[1], pmax[1]);
	double ymax = std::max(pmin[1], pmax[1]);
	double zmin = std::min(pmin[2], pmax[2]);
	double zmax = std::max(pmin[2], pmax[2]);

	pmin = G4ThreeVector(xmin, ymin , zmin);
	pmax = G4ThreeVector(xmax, ymax, zmax);
	
	return bbox(pmin, pmax);
}

bool boundsOverlap(bbox b0, bbox b1)
{
  auto pmin0 = b0.first;
  auto pmax0 = b0.second;
  auto pmin1 = b1.first;
  auto pmax1 = b1.second;

  for(int i=0; i < 3; i++) {
	if (pmin1[i] > pmax0[i] || pmin0[i] > pmax1[i])
	  return false;
  }

  return true;
}


bool multiUnionNodesOverlap(const G4MultiUnion* mu)
{
  // check if the nodes of the multiunion overlap
  int nodes = mu->GetNumberOfSolids();
  auto bbox0 = bounds(mu, 0);
  
  for(int i=1; i < nodes; i++) {
	auto bbox1 = bounds(mu, i);
	// we will assume there is overlap if the bounding boxes of consecutive
	// nodes overlap
	if (boundsOverlap(bbox0, bbox1))
	  return true;
	bbox0 = bbox1;
  }

  return false;
}


Moments MultiUnionMomentsEstimator::estimate()
{
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix() // momentsof Inertia
  };

  G4MultiUnion* mu = (G4MultiUnion *) f_solid;
  bool nodesOverlap = multiUnionNodesOverlap(mu);
  G4cout << "Nodes overlap: " << nodesOverlap << endl;
  if (nodesOverlap) {
	// if nodes overlap we treat the multiunion as a single solid
	// and use monte carlo to get its moments

	auto estimator = MomentsEstimator(mu);
	return estimator.estimate();
  }

  // otherwise, we do calculations on each individual node.
  // The reason is that by default Geant4 estimates volume for G4Multiuniomn
  // using only 10,000 events. This is far too low, especially for
  // a sparse G4Multiunion, which is the case for exports of arrays from
  // FreeCAD. There it would be more accurate to get volume of each
  // node and add them up.
  
  G4ThreeVector cm = G4ThreeVector(0, 0, 0);
  double totVol = 0;
  int nodes = mu->GetNumberOfSolids();
  G4cout << "nodes: " << nodes << endl;
  if(nodes > 0) {
	G4VSolid *solid = mu->GetSolid(0);
    double volume = solid->GetCubicVolume();
	G4cout << "node volume: " << volume << endl;
  }

  std::vector<Moments> nodeMoments;
  // In FreeCAD, the matrix of inertia for a multiunion is calculated about
  // the Center of Mass of the multiunion, so here we
  // have to do the same.
  
  for(int i=0; i < nodes; i++) {
	G4VSolid *solid = mu->GetSolid(i);
    G4Transform3D transform = mu->GetTransformation(i);
    G4ThreeVector translation = transform.getTranslation();
    double volume = solid->GetCubicVolume();
    G4ThreeVector center = G4ThreeVector(0, 0, 0);

	G4String entityType = solid->GetEntityType();

	auto calculator = MomentsCalculator::getCalculator(solid);
	if (calculator != nullptr) {
	  moments = calculator->calculate();
	  delete calculator;
	}
	else {
	  auto estimator = MomentsEstimator::getEstimator(solid);
	  moments = estimator->estimate();
	  delete estimator;
	}


    auto rot = transform.getRotation();
    auto translate = transform.getTranslation();
    center = rot * moments.cm;
    center = translate + center;
	// Make the center of the node, the translated center
	moments.cm = center;
	nodeMoments.push_back(moments);

    cm +=  moments.Volume * center;

    totVol += volume;
  }

  // recompute moments about the CM of the multiunion
  cm /= totVol;
  moments.cm = cm;
  moments.Volume = totVol;
  for(size_t i=0; i < nodeMoments.size(); i++) {
	Moments imoment = nodeMoments[i];
	double V = imoment.Volume;
	// the position of the node from the center of mass, i.e., relative to center of mass of the union
	G4ThreeVector cmnode = imoment.cm - cm;
	
	double Ixx = imoment.M.xx() + V*(cmnode.y()*cmnode.y() + cmnode.z()*cmnode.z());
	double Iyy = imoment.M.yy() + V*(cmnode.x()*cmnode.x() + cmnode.z()*cmnode.z());
	double Izz = imoment.M.zz() + V*(cmnode.x()*cmnode.x() + cmnode.y()*cmnode.y());
	
	double Ixy = imoment.M.xy() - V*cmnode.x()*cmnode.y();
	double Ixz = imoment.M.xz() - V*cmnode.x()*cmnode.z();
	double Iyz = imoment.M.yz() - V*cmnode.y()*cmnode.z();
	moments.M += InertiaMatrix(Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz); 
  }

  
  return moments;
}

