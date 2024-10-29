#include "TessellationMoments.hh"
#include "ExtrusionMoments.hh"
#include "InertiaMatrix.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "Moments.hh"
#include <CLHEP/Vector/AxisAngle.h>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4TwoVector.hh>
#include <G4ios.hh>
#include <vector>



G4RotationMatrix RotationOfTriangle(G4ThreeVector v1, G4ThreeVector v2, G4ThreeVector v3) {
  // return rotation matrix that takes the normal to the plane of the triangle formed
  // by v1, v2, v3 to the z-axis
  G4ThreeVector v12 = v2 - v1;
  G4ThreeVector v13 = v3 - v1;
  G4ThreeVector normal = v12.cross(v13);
  normal = normal.unit();
  G4ThreeVector axis = G4ThreeVector(normal.y(), -normal.x(), 0);  //normal cross k
  double thet = std::acos(normal.z());  // normal dot k
  
  return G4RotationMatrix(CLHEP::HepAxisAngle(axis, thet));
}

// Function to compute the inertia matrix of a triangle
Moments computeTriangleInertia0(const G4ThreeVector& v1, const G4ThreeVector& v2, const G4ThreeVector& v3) {
  double density = 1;
  // Mass of the triangle (area * density)
  G4ThreeVector v1v2 = v2 - v1;
  G4ThreeVector v1v3 = v3 - v1;
  G4ThreeVector cross = v1v2.cross(v1v3);
  double area = 0.5 * sqrt(cross.dot(cross));
  // Mass of the triangle (area * density)
  G4ThreeVector u = cross.unit();  
  G4ThreeVector trianglecm = (v1+v2+v3)/3;
  double h = u.dot(trianglecm);
  double volume = 1/3. * area * h;
  double mass = volume * density;

  // Vertices of the triangle (centered around the centroid)
  double x1 = v1.x(), y1 = v1.y(), z1 = v1.z();
  double x2 = v2.x(), y2 = v2.y(), z2 = v2.z();
  double x3 = v3.x(), y3 = v3.y(), z3 = v3.z();

  // Compute the inertia components
  // Moments of inertia relative to v1 (assumed to be at the origin)
  double Ixx = mass / 60.0 * (y1 * y1 + z1 * z1 + y2 * y2 + z2 * z2 + y3 * y3 + z3 * z3);
  double Iyy = mass / 60.0 * (x1 * x1 + z1 * z1 + x2 * x2 + z2 * z2 + x3 * x3 + z3 * z3);
  double Izz = mass / 60.0 * (x1 * x1 + y1 * y1 + x2 * x2 + y2 * y2 + x3 * x3 + y3 * y3);
  
  // Products of inertia
  double Ixy = mass / 120.0 * (x1 * y1 + x2 * y2 + x3 * y3);
  double Ixz = mass / 120.0 * (x1 * z1 + x2 * z2 + x3 * z3);
  double Iyz = mass / 120.0 * (y1 * z1 + y2 * z2 + y3 * z3);
  

  InertiaMatrix inertia(Ixx, Ixy, Ixz,
						Ixy, Iyy, Iyz,
						Ixz, Iyz, Izz);

  G4cout << "h: " << h << " area: " << area << " volume: " << volume << G4endl;
  

  Moments moments;
  moments.Volume = volume;
  moments.cm = (v1+v2+v3)/4;
  moments.M = inertia;
	
  return moments;
}

Moments computeTriangleInertia(const G4ThreeVector& v1, const G4ThreeVector& v2, const G4ThreeVector& v3) {
  double volume = v1.dot(v2.cross(v3))/6;
  G4ThreeVector cm = (v1+v2+v3)/4;

  
  double xx = (v3.x()*v3.x()+(v2.x()+v1.x())*v3.x()+v2.x()*v2.x()+v1.x()*v2.x()+v1.x()*v1.x())/60;
  double yy = (v3.y()*v3.y()+(v2.y()+v1.y())*v3.y()+v2.y()*v2.y()+v1.y()*v2.y()+v1.y()*v1.y())/60;
  double zz = (v3.z()*v3.z()+(v2.z()+v1.z())*v3.z()+v2.z()*v2.z()+v1.z()*v2.z()+v1.z()*v1.z())/60;

  double Ixx = yy + zz;
  double Iyy = xx + zz;
  double Izz = xx + yy;

  double Ixy = -((2*v3.x()+v2.x()+v1.x())*v3.y()+(v3.x()+2*v2.x()+v1.x())*v2.y()+(v3.x()+v2.x()+2*v1.x())*v1.y())/120;
  double Iyz = -((2*v3.y()+v2.y()+v1.y())*v3.z()+(v3.y()+2*v2.y()+v1.y())*v2.z()+(v3.y()+v2.y()+2*v1.y())*v1.z())/120;
  double Ixz = -((2*v3.x()+v2.x()+v1.x())*v3.z()+(v3.x()+2*v2.x()+v1.x())*v2.z()+(v3.x()+v2.x()+2*v1.x())*v1.z())/120;
  
  InertiaMatrix II = InertiaMatrix(Ixx, Ixy, Ixz,
								   Ixy, Iyy, Iyz,
								   Ixz, Iyz, Izz);

  II *= 6*volume;

  Moments moments;
  moments.Volume = volume;
  moments.cm = cm;
  moments.M = II;
	
  return moments;
  
}

/*
// Function to compute the inertia matrix of a triangle
Moments computeTriangleInertia0(const G4ThreeVector& v1, const G4ThreeVector& v2, const G4ThreeVector& v3) {

  // Rotation that takes the normal to the triangle to the z-axis
  G4RotationMatrix R = RotationOfTriangle(v1, v2, v3);
  G4ThreeVector v1R = R*v1;
  G4ThreeVector v2R = R*v2;
  G4ThreeVector v3R = R*v3;
  // All rotated vertecies should now have the same z. Check
  G4cout << "v1R: " << v1R << " v2R: " << v2R << " v3R: " << v3R << G4endl;
  
  std::vector<G4TwoVector> verts = {
	G4TwoVector(v1R.x(), v1R.y()),
	G4TwoVector(v2R.x(), v2R.y()),
	G4TwoVector(v3R.x(), v3R.y()),
  };	  
  G4TwoVector cm2D = polygonCM(verts);
  for (int i=0; i < 3; i++) {
	verts[i] -= cm2D;
  }
  
  double area = std::fabs(polygonArea(verts));
  // G4cout << "area: " << area << G4endl;

  InertiaMatrix dI_dz = computeInertiaMatrix(verts);
  double h = v1R.z();  // height of triangle, now parallel to x-y plane
  // dI_dz has units of length^4: dI_dz = inertia per unit length = (length^2) * volume/length
  // so dI_dz = area(top) * t^4 (t=0 at z=0, t = 1 at z=h)
  // So, ignoring the fact that the center of mass of the triangle is NOT necessarily
  // at the origin, the integrated moment of inertia of the tetrahedron is
  // Izz = integral ( dIzz_dz *dz); z = t *h ==> dz = h dt
  //       =  h * dIzz_dz (top) integral ( t^4 dt ) = h/5*dIzz_tz(top) -> triangle
  // But dII_dz(top) is computed about the 2D center of the triangle, so we have to adjust for
  // that using parallel axis theorem:
  // dIzz_z(z) = dII_z(top) * t^4 + area(t) * (cm2D_x(t)^2 + cm2D_y(t)^2
  // so I_zz = integral(dI_dz(top) * t^4 dz + area(top)*t^2 * (cm2D_x(top)^2 * t^2 + cm2D_y(top*t^2) dz
  // = h/5*dI_dz_zz(top) + area(top)*h/5 * (cm2D_x^2 + cm2D_y^2)
  double Izz = h/5*(dI_dz.zz() + area*(cm2D.x()*cm2D.x() + cm2D.y()*cm2D.y()));
  // The same scaling should work for Ixx and Iyy
  double Ixx = h/5*(dI_dz.xx() + area*(cm2D.y()*cm2D.y() + h*h));
  double Iyy = h/5*(dI_dz.yy() + area*(cm2D.x()*cm2D.x() + h*h));

  // Now the cross terms should be zer about the CM, so only contribution is
  // from parallel axis theorem:
  // Ixy = -integral(area(t) * cm2D_x(t)* cm2D_y(t) dz
  //     = -integral(area(top) * t^2 * cm2D_x(top) * t cm2D_y(top)*t * (dz=h*dt))
  double Ixy = -h/5.*area*cm2D.x()*cm2D.y();
  double Ixz = -h/5.*area*cm2D.x()*h;
  double Iyz = -h/5.*area*cm2D.y()*h;
  // Above are moments with the tetrahedron rotated so its normal is along the axis.
  // We should rotate back to moments for the oritginal tetrahedron orientation
  G4RotationMatrix Rinv = R.inverse();
  InertiaMatrix II = InertiaMatrix(Ixx, Ixy, Ixz,
								   Ixy, Iyy, Iyz,
								   Ixz, Iyz, Izz);
  					
  II = II.rotate(Rinv);
  
  
  // Mass of the triangle (area * density)
  G4ThreeVector v12 = v2 - v1;
  G4ThreeVector v13 = v3 - v1;
  G4ThreeVector normal = v12.cross(v13);
  G4ThreeVector u = normal.unit();
  
  G4ThreeVector cm = (v1+v2+v3)/4;
  
  double volume = 1/3. * area * h;
  
  
  Moments moments;
  moments.Volume = volume;
  moments.cm = cm;
  moments.M = II;
	
  return moments;
}
*/

Moments TessellationMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };
  
  G4TessellatedSolid *solid = (G4TessellatedSolid *) f_solid;
  InertiaMatrix II = InertiaMatrix();
  G4ThreeVector weightedCM = G4ThreeVector(0, 0, 0);

  double myVolume = 0;
  int nfaces = solid->GetNumberOfFacets();
  for (int i=0; i < nfaces; i++) {
	G4VFacet *f = solid->GetFacet(i);
	if (f->GetEntityType() == "G4TriangularFacet") {
	  G4TriangularFacet *tri = (G4TriangularFacet *) f;
	  if (tri->GetNumberOfVertices()!= 3) {
		G4cerr << "Hmm... A G4TriangularFacet with number of vertices " <<
		  tri->GetNumberOfVertices() << " not equal to 3" << G4endl;
		exit(0);
	  }
	  Moments moments_i = computeTriangleInertia(f->GetVertex(0), f->GetVertex(1), f->GetVertex(2)); 
	  II += moments_i.M;
	  weightedCM += moments_i.Volume * moments_i.cm;
	  myVolume += moments_i.Volume;
	}
	else { // must be quadrangular facet
	  G4QuadrangularFacet *quad = (G4QuadrangularFacet *) f;
	  if (quad->GetNumberOfVertices()!= 4) {
		G4cerr << "Hmm... A G4QuadrangularFacet with number of vertices " <<
		  quad->GetNumberOfVertices() << " not equal to 4" << G4endl;
		exit(0);
	  }
	  G4ThreeVector v1 = quad->GetVertex(0);
	  G4ThreeVector v2 = quad->GetVertex(1);
	  G4ThreeVector v3 = quad->GetVertex(2);
	  G4ThreeVector v4 = quad->GetVertex(3);
	  Moments moments_i = computeTriangleInertia(v1, v2, v3); 
	  II += moments_i.M;
	  weightedCM += moments_i.Volume * moments_i.cm;
	  myVolume += moments_i.Volume;
	  moments_i = computeTriangleInertia(v1, v3, v4); 
	  II += moments_i.M;
	  weightedCM += moments_i.Volume * moments_i.cm;
	  myVolume += moments_i.Volume;
	}
  }
  
  double volume = solid->GetCubicVolume();
  // Double check volumes
  G4cout << "geant4 vol: " << volume << " myVolume: " << myVolume << G4endl;
  
  moments.cm = 1/volume * weightedCM; 
  
  moments.M = II;
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;

}
