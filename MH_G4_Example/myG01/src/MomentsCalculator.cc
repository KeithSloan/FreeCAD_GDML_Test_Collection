#include "MomentsCalculator.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Torus.hh"
#include "G4Paraboloid.hh"
#include "G4Hype.hh"
#include "G4Tet.hh"
#include "G4ExtrudedSolid.hh"
#include "G4EllipticalTube.hh"
#include "ExtrusionMoments.hh"
#include "InertiaMatrix.hh"
#include "TessellationMoments.hh"
#include "PolyhedraMoments.hh"
#include "G4Polycone.hh"
#include "G4PolyconeHistorical.hh"
#include "G4TessellatedSolid.hh"
#include "G4EllipticalCone.hh"
#include "G4ThreeVector.hh"
#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"
#include "G4ios.hh"
#include <any>
#include <cmath>
#include <vector>

// ---------- some utility functions -----------------
double polygonArea(const std::vector<G4TwoVector> vertices) {
  // Note that this calculation of the area has a sign. The area is POSITIVE
  // if the points of the polygon are fiven in an anti-clockwise order
  // and NEGATIVE if the points are given in a clockwise direction.
  // It seems geant4 always orders the 2D vertexes of the xtrusion in a clockwise
  // direction, regardless of the order they are given in the gdml file.
  // I will maintain the sign here, because the CM uses the same area calculation
  // so it will only given the correct result if the sign of the area here keeps
  // its sign.
  double area = 0.0;
  int n = vertices.size();
  for (int i = 0; i < n; i++) {
	int j = (i + 1) % n;
	area += vertices[i].x() * vertices[j].y() - vertices[j].x() * vertices[i].y();
  }
  
  return area / 2.0;

}

G4TwoVector polygonCM(const std::vector<G4TwoVector> vertices) {
  double cx = 0.0, cy = 0.0;
  double area = polygonArea(vertices);
  double factor = 0.0;

  int n = vertices.size();
  for (int i = 0; i < n; i++) {
	int j = (i + 1) % n;
	factor = (vertices[i].x() * vertices[j].y() - vertices[j].x() * vertices[i].y());
	cx += (vertices[i].x() + vertices[j].x()) * factor;
	cy += (vertices[i].y() + vertices[j].y()) * factor;
  }

  cx /= (6.0 * area);
  cy /= (6.0 * area);

  return G4TwoVector(cx, cy);
}
										   
InertiaMatrix computeInertiaMatrix(std::vector<G4TwoVector> vertices) {
  // returns moments of inertia about the center of "mass" of the extrusions per unit height
  double Ixx2D = 0.0, Iyy2D = 0.0, Ixy2D = 0.0;

  G4TwoVector cm = polygonCM(vertices);
  
  int n = vertices.size();

  for (int i = 0; i < n; i++) {
	int j = (i + 1) % n;	
	G4TwoVector v1 = vertices[i] - cm;
	G4TwoVector v2 = vertices[j] - cm;

	double commonFactor = v1.x() * v2.y() - v2.x() * v1.y();

	Ixx2D += (v1.y() * v1.y() + v1.y() * v2.y() + v2.y() * v2.y()) * commonFactor;
	Iyy2D += (v1.x() * v1.x() + v1.x() * v2.x() + v2.x() * v2.x()) * commonFactor;
	Ixy2D += (v1.x() * (2 * v1.y() + v2.y()) + v2.x() * (2 * v2.y() + v1.y())) * commonFactor;
  }

  Ixx2D = std::fabs(Ixx2D) / 12.0;
  Iyy2D = std::fabs(Iyy2D) / 12.0;
  Ixy2D = std::fabs(Ixy2D) / 24.0;

  // The height adds to the z-axis moment of inertia
  double Izz = (Ixx2D + Iyy2D);

  // Moments in the x-y plane
  double Ixx = Ixx2D;
  double Iyy = Iyy2D;
  double Ixy = Ixy2D;

  // having calculate the moments about the center of mass
  // the extrusion is symmetric along the z-direction, Ixz and Iyz are 0
  double Ixz = 0;
  double Iyz = 0;

  return InertiaMatrix(Ixx, Ixy, Ixz,
					   Ixy, Iyy, Iyz,
					   Ixz, Iyz, Izz);

}

//------------------------------------------------------------------------------
MomentsCalculator* MomentsCalculator::getCalculator(G4VSolid *solid) {
  auto type = solid->GetEntityType();
  
  if (type == "G4Box") {
	return new BoxMoments(solid);
  }
  else if (type == "G4Ellipsoid") {
	return new EllipsoidMoments(solid);
  }
  else if (type == "G4Orb") {
	return new OrbMoments(solid);
  }
  else if (type == "G4Sphere") {
	return new SphereMoments(solid);
  }
  else if (type == "G4Cons") {
	return new ConeMoments(solid);
  }
  else if (type == "G4Tubs") {
	return new TubeMoments(solid);
  }
  else if (type == "G4Trd") {
	return new TrdMoments(solid);
  }
  else if (type == "G4Torus") {
	return new TorusMoments(solid);
  }
  else if (type == "G4Hype") {
	return new HypeMoments(solid);
  }
  else if (type == "G4Paraboloid") {
	return new ParaboloidMoments(solid);
  }
  else if (type == "G4ExtrudedSolid") {
	return new ExtrusionMoments(solid);
  }
  else if (type == "G4TessellatedSolid") {
	return new TessellationMoments(solid);
  }
  else if (type == "G4Tet") {
	return new TetrahedronMoments(solid);
  }
  else if (type == "G4Polycone") {
	return new PolyconeMoments(solid);
  }
  else if (type == "G4EllipticalTube") {
	return new EllipticalTubeMoments(solid);
  }
  else if (type == "G4EllipticalCone") {
	return new EllipticalConeMoments(solid);
  }
  else if (type == "G4Polyhedra") {
	return new PolyhedraMoments(solid);
  }

  return nullptr;
}

Moments BoxMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Box *solid = (G4Box *) f_solid;
  
  double x = 2 * solid->GetXHalfLength();
  double y = 2 * solid->GetYHalfLength();
  double z = 2 * solid->GetZHalfLength();

  double mass = solid->GetCubicVolume();
  double I_x = (1.0 / 12.0) * mass * (y * y + z * z);
  double I_y = (1.0 / 12.0) * mass * (x * x + z * z);
  double I_z = (1.0 / 12.0) * mass * (x * x + y * y);

  moments.cm = G4ThreeVector(0, 0, 0);
  moments.M = InertiaMatrix(I_x, 0, 0,
							0, I_y, 0,
							0, 0, I_z);
  moments.Volume = mass;

  return moments;
}

Moments EllipsoidMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Ellipsoid *solid = (G4Ellipsoid *) f_solid;
  
  double ax = solid->GetDx();
  double by = solid->GetDy();
  double cz = solid->GetDz();

  double mass = solid->GetCubicVolume();
  double I_x = (1.0 / 5.0) * mass * (by * by + cz * cz);
  double I_y = (1.0 / 5.0) * mass * (ax * ax + cz * cz);
  double I_z = (1.0 / 5.0) * mass * (ax * ax + by * by);

  moments.cm = G4ThreeVector(0, 0, 0);
  moments.M = InertiaMatrix(I_x, 0, 0,
							0, I_y, 0,
							0, 0, I_z);
  moments.Volume = mass;

  return moments;
}


Moments OrbMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Orb *solid = (G4Orb *) f_solid;

  double r = solid->GetRadius();

  double mass = solid->GetCubicVolume();
  double I_x = (2.0 / 5.0) * mass * r * r;
  double I_y = I_x;
  double I_z = I_x;

  moments.cm = G4ThreeVector(0, 0, 0);
  moments.M = InertiaMatrix(I_x, 0, 0,
							0, I_y, 0,
							0, 0, I_z);
  moments.Volume = mass;

  return moments;
}

Moments SphereMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Sphere *solid = (G4Sphere *) f_solid;

  double rmin = solid->GetInnerRadius();
  double rmax = solid->GetOuterRadius();
  double thet_s = solid->GetStartThetaAngle();
  double deltaThet = solid->GetDeltaThetaAngle();
  double thet_e =  deltaThet + thet_s;
  double phi_s = solid->GetStartPhiAngle();
  double deltaPhi = solid->GetDeltaPhiAngle();
  double phi_e = deltaPhi + phi_s;

  /*
	Diagonal moments of inertia
	x = r * sin(thet) * cos(phi)
	y = r * sin(thet) * sin(phi)
	z = r * cos(thet)

	dV = r^2 * dr * sin(thet) * dthet * dphi

	Izz = integral (x^2 + y^2) * dV
	    = integral r^2* ( sin^2(thet) * cos^2(phi) + sin^2(thet) * sin^2(phi)) * sin(thet)
		* r^2 dr * dthet * dphi

	radial_integral = integral r^4 dr = r^5/5.

	the angle term rdeuces to sin^2(thet) * (cos^2(phi) + sin^2(phi)) * sin(thet)*dthet*dphi
	= integral(sin^3(thet) dthet) * integral(dphi)
	theta_integral = integral (sin^2(thet) * sin(thet) dthet)
	               = integral (1 - cos^2(thet)) (-d(cos(thet))
				   = -cos(thet) + cos^3(thet) /3
	phi_integral = phi

	Will come to Ixx and Iyy in a second, but we introduce some notation here
	the indefine integral will be a lambda function:
	auto fr = [] (double r) return {r*r*r*r*r/5.;};
	and the integral will be:
	intfr = [fr] (double x1, double x2) { return fr(x2) - fr(x1); };
	and the pattern will repeat

	We first show Izz and then continute with Ixx, Iyy
   */

  auto fr = [] (double r) -> double  { return r*r*r*r*r/5.; };
  auto intfr = [fr] (double x1, double x2) -> double {return fr(x2) - fr(x1);};

  auto fthet_zz = [] (double thet) -> double { return std::cos(thet)*(-1 + std::cos(thet)*std::cos(thet)/3.); };
  auto intfthet_zz = [fthet_zz] (double x1, double x2) -> double {
	return fthet_zz(x2) - fthet_zz(x1); };

  auto fphi_zz = [] (double x) -> double {return x; };
  auto intfphi_zz = [fphi_zz] (double x1, double x2) -> double {return fphi_zz(x2) - fphi_zz(x1);};

  double I_zz = intfr(rmin, rmax) * intfthet_zz(thet_s, thet_e) * intfphi_zz(phi_s, phi_e);

  /*
	Now continue with Ixx
	Ixx = integral( (y^2 + z^2) * dV )
	    = radial_integral * integral ( (sin^2(thet) * sin^2(phi) + cos^2(thet) ) * sin(thet) * dthet * dphi)
		angular terms = (sin^2(thet) * ( 1 -cos^2(phi)) + cos^2(thet) ) * sin(thet)
		              = (1 - sin^2(thet) * cos^2(phi) * sin(thet)
					  = integral (sin(thet) dthet dphi) - integral(sin^3(thet) * cos^2(phi) dthet dphi
		 the only term we have not seen before is the cos^2(phi) dphi. Using the trig identify
		 cos^2(phi) = (1 + cos(2*phi))/2
		 the integral is phi/2 + sin(2*phi)/4;
		 so the integral becomes:
   */
  auto f1thet_xx = [] (double thet) -> double {return -std::cos(thet); };
  auto intf1thet_xx = [f1thet_xx] (double x1, double x2) -> double {return f1thet_xx(x2) - f1thet_xx(x1); };

  auto f1phi_xx = [] (double phi) -> double { return phi; }; // yes trivial lambda, but allows for the definite integral definition
  auto intf1phi_xx = [f1phi_xx] (double x1, double x2) -> double {return f1phi_xx(x2) - f1phi_xx(x1); };

  auto f2thet_xx = [fthet_zz] (double thet) -> double { return fthet_zz(thet); };
  auto intf2thet_xx = [f2thet_xx] (double x1, double x2) -> double {return f2thet_xx(x2) - f2thet_xx(x1); };

  auto f2phi_xx = [] (double phi) -> double {return phi/2 + std::sin(2*phi)/4.; };
  auto intf2phi_xx = [f2phi_xx] (double x1, double x2) -> double {return f2phi_xx(x2) - f2phi_xx(x1); };

  double I_xx = intfr(rmin, rmax) * (intf1thet_xx(thet_s, thet_e) * intf1phi_xx(phi_s, phi_e)
									 - intf2thet_xx(thet_s, thet_e) * intf2phi_xx(phi_s, phi_e));

  /*
	Now Iyy
	Iyy = integral( (x^2 + z^2) * dV )
	    = radial_integral * integral ( (sin^2(thet) * cos^2(phi) + cos^2(thet) ) * sin(thet) * dthet * dphi)
		angular terms = (sin^2(thet) * ( 1 -sin^2(phi)) + cos^2(thet) ) * sin(thet)
		              = (1 - sin^2(thet) * sin^2(phi) * sin(thet)
					  = integral (sin(thet) dthet dphi) - integral(sin^3(thet) * sin^2(phi) dthet dphi
		 the only term we have not seen before is the sin^2(phi) dphi. Using the trig identify
		 sin^2(phi) = (1 - cos(2*phi))/2
		 the integral is phi/2 - sin(2*phi)/4;
		 so the integral becomes:
   */
  
  auto f1thet_yy = [] (double thet) -> double {return -std::cos(thet); };
  auto intf1thet_yy = [f1thet_yy] (double x1, double x2) -> double {return f1thet_yy(x2) - f1thet_yy(x1); };

  auto f1phi_yy = [] (double phi) -> double { return phi; }; // yes trivial lambda, but allows for the definite integral definition
  auto intf1phi_yy = [f1phi_yy] (double x1, double x2) -> double {return f1phi_yy(x2) - f1phi_yy(x1); };

  auto f2thet_yy = [fthet_zz] (double thet) -> double { return fthet_zz(thet); };
  auto intf2thet_yy = [f2thet_yy] (double x1, double x2) -> double {return f2thet_yy(x2) - f2thet_yy(x1); };

  auto f2phi_yy = [] (double phi) -> double {return phi/2 - std::sin(2*phi)/4.; };
  auto intf2phi_yy = [f2phi_yy] (double x1, double x2) -> double {return f2phi_yy(x2) - f2phi_yy(x1); };

  double I_yy = intfr(rmin, rmax) * (intf1thet_yy(thet_s, thet_e) * intf1phi_yy(phi_s, phi_e)
									 - intf2thet_yy(thet_s, thet_e) * intf2phi_yy(phi_s, phi_e));
 

  /*
	The off diagonal terms ought to be slightly easier:
	Ixy = integral ( -x*y * dV)
	   = radial_integral * integral (sin(thet)*cos(phi) * sin(thet)*sin(phi) * sin(thet) *dthet* dphi
	   the angular integral are a product of the theta integral:
	   theta integral ( sin^3(thet) dthet )  // we've done this before (fthet_zz)
	   The phi integral is integral( cos(phi)*sin(phi) dphi = sin^2(phi)/2
   */
  auto fthet_xy = [fthet_zz] (double thet) -> double {return fthet_zz(thet); };
  auto intfthet_xy = [fthet_xy]  (double x1, double x2) -> double {return fthet_xy(x2) - fthet_xy(x1); };

  auto fphi_xy = [] (double phi) -> double {return std::sin(phi)*std::sin(phi)/2.; };
  auto intfphi_xy = [fphi_xy]  (double x1, double x2) -> double {return fphi_xy(x2) - fphi_xy(x1); };

  double I_xy = -intfr(rmin, rmax) * intfthet_xy(thet_s, thet_e) * intfphi_xy(phi_s, phi_e);
  /*
	Ixz = integral ( -x*z * dV)
	   = radial_integral * integral (sin(thet)*cos(phi) * cos(thet) * sin(thet) *dthet* dphi
	   the angular integral are a product of the theta integral:
	   theta integral ( sin^2(thet) * cos(thet) * dthet )
	   = sin^3(thet)/3
	   The phi integral is integral( cos(phi)*dphi = sin(phi)
   */
  auto fthet_xz = [] (double thet) -> double {return std::sin(thet)*std::sin(thet)*std::sin(thet)/3.; };
  auto intfthet_xz = [fthet_xz]  (double x1, double x2) -> double {return fthet_xz(x2) - fthet_xz(x1); };

  auto fphi_xz = [] (double phi) -> double {return std::sin(phi); };
  auto intfphi_xz = [fphi_xz]  (double x1, double x2) -> double {return fphi_xz(x2) - fphi_xz(x1); };

  double I_xz = -intfr(rmin, rmax) * intfthet_xz(thet_s, thet_e) * intfphi_xz(phi_s, phi_e);

  /*
	Iyz = integral ( -y*z * dV)
	= radial_integral * integral (sin(thet)*sin(phi) * cos(thet) * sin(thet) *dthet* dphi
	the angular integral are a product of the theta integral:
	theta integral ( sin^2(thet) * cos(thet) * dthet )
	= sin^3(thet)/3
	The phi integral is integral( sin(phi)*dphi = -cos(phi)
  */
  auto fthet_yz = [] (double thet) -> double {return std::sin(thet)*std::sin(thet)*std::sin(thet)/3.; };
  auto intfthet_yz = [fthet_yz]  (double x1, double x2) -> double {return fthet_yz(x2) - fthet_yz(x1); };

  auto fphi_yz = [] (double phi) -> double {return -std::cos(phi); };
  auto intfphi_yz = [fphi_yz]  (double x1, double x2) -> double {return fphi_yz(x2) - fphi_yz(x1); };

  double I_yz = -intfr(rmin, rmax) * intfthet_yz(thet_s, thet_e) * intfphi_yz(phi_s, phi_e);

  // center of mass
  double volume = solid->GetCubicVolume();

  G4double R4 = rmax*rmax*rmax*rmax;
  G4double r4 = rmin*rmin*rmin*rmin;
  G4double radial_part = 1/4.*(R4 - r4) / volume;
  
  G4double xcm = radial_part * ( deltaThet/2. - 1/4. * (std::sin(2*thet_e) - std::sin(2*thet_s)) ) *
	(std::sin(phi_e) - std::sin(phi_s));
  
  G4double ycm = radial_part * ( deltaThet/2. - 1/4. * (std::sin(2*thet_e) - std::sin(2*thet_s)) ) *
	(std::cos(phi_s) - std::cos(phi_e));
  
  G4double zcm = radial_part * 1/2. * (std::cos(thet_s)*std::cos(thet_s) - std::cos(thet_e)*std::cos(thet_e)) * deltaPhi;
  
  

  moments.cm = G4ThreeVector(xcm, ycm, zcm);
  moments.M = InertiaMatrix(I_xx, I_xy, I_xz,
							I_xy, I_yy, I_yz,
							I_xz, I_yz, I_zz);
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;
}

Moments ConeMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Cons *solid = (G4Cons *) f_solid;

  double rmin1 = solid->GetInnerRadiusMinusZ();
  double rmax1 = solid->GetOuterRadiusMinusZ();
  double rmin2 = solid->GetInnerRadiusPlusZ();
  double rmax2 = solid->GetOuterRadiusPlusZ();
  double h = solid->GetZHalfLength();
  h *= 2;

  double deltaPhi = solid->GetDeltaPhiAngle();
  double phi_s = solid->GetStartPhiAngle();
  double phi_e = solid->GetDeltaPhiAngle() + phi_s;
 

  /*
	For cone use cylindrical coordinates: r, phi, z
	x = r cos(phi)
	y = r sin(phi)

	In what follows we need to parametrize the inner and out radius of the cone as
	a function of height along the cone z:

	rmin(t) = rmin1 + t*(rmin2 - rmin1), t=0..1. This give rmin=rmin1 at t=0, rmin=rmin2 at t = 1
	rmax(t) = rmax1 + t*(rmax2 - rmax1), t=0..1. This give rmax=rmax1 at t=0, rmax=rmax2 at t = 1

	t = (z+h/2)/h  gives t = 0 at z = -h/2 and t = 1 at z = h/2
	above implies dt = dz/h and z = h*t - h/2

	dV = r dr dphi dz

	We will start with Izz:

	Izz = integral( x^2 + y2  ) dV
	    = integral( r^2*cos^2(phi) + r^2*sin^2(phi) * dV
		= integral( r^3 dr dPhi dz
		= deltaPhi * 1/4*integral( (rmax(z)^4 - rmin(z)^4) dz
		Now change variables from z to t:
		= deltaPhi * h/4 * integral( rmax(t)^4 - rmin(t)^4) dt, t=0..1)

	According to wxMaxima
	Ir = integral (r1 + t*(r2-r1))^4 dt, t=0..1 = 1/5*(r2^4 + r1*r2^3 + r1^2*r2^2 + r1^3*r2 + r1^4) 

	so finally
	Izz = deltaPhi * h/20 * (Ir(r1=r1max, r2=r2max) - Ir(r1=r1min, r2=r2min)
  */

  auto Ir = [] (double r1, double r2) -> double {return r2*r2*r2*r2 + r2*r2*r2*r1 + r2*r2*r1*r1 + r2*r1*r1*r1 + r1*r1*r1*r1;};
  double I_zz = deltaPhi * h/20. * (Ir(rmax1, rmax2) - Ir(rmin1, rmin2));

  /*
	For Ixx we need integral(y^2 + z^2) * dV
	and for Iyy we need integral(x^2 + z^2) * dV

	We will calculate the integrals over x^2, y^2, z^2 individually and later add the appropriate terms.

	integral (x^2) dV = integral (r^3 dr dz) * cos^2(phi) dphi
	so the radial and z integral are the same as obove, but the phi integral becomes integral(cos^2(phi) dphi) = (phi/2 + sin(2*phi)/4) 
   */

  auto fphi_xx = [] (double phi) -> double {return phi/2 + std::sin(2*phi)/4;};
  auto intfphi_xx = [fphi_xx] (double a1, double a2) -> double {return fphi_xx(a2) - fphi_xx(a1); };
  double intxx = h/20. * (Ir(rmax1, rmax2) - Ir(rmin1, rmin2)) * intfphi_xx(phi_s, phi_e);

  /*
	for y^2 the only difference is is the phi term: integral(sin^2(phi) dphi) which gives (phi/2 - sin(2*phi)/4) (Note the minus sign)
   */

  auto fphi_yy = [] (double phi) -> double {return phi/2 - std::sin(2*phi)/4;};
  auto intfphi_yy = [fphi_yy] (double a1, double a2) -> double {return fphi_yy(a2) - fphi_yy(a1); };
  double intyy = h/20. * (Ir(rmax1, rmax2) - Ir(rmin1, rmin2)) * intfphi_yy(phi_s, phi_e);

  /*
	Now the z^2 term
	integral(z^2 dV) = integral( r dr z^2 dz) integral(dPhi)

	integral (dPhi) = phi

	integral( r(z)^2/2 z^2 dz )
	changing variables from z to t:
	1/2* r(t)^2 * (h*t - h/2)^2 h dt = h/2 *1/60 * h^2 * (2* r2^2 + r1*r2 + 2*r1^2) (from wxMaxima)
   */

  auto Ir_zz = [] (double r1, double r2) -> double {return 2*r2*r2 + r1*r2 + 2*r1*r1; };
  double intzz = deltaPhi * h*h*h/120. * (Ir_zz(rmax1, rmax2) - Ir_zz(rmin1, rmin2));

  double I_xx = intyy + intzz;
  double I_yy = intxx + intzz;

  /*
	Now the off-diagonal elements:
	Ixy = -integral(x*y * dV)
	    = -integral( r^3 dr dz) * integral(cos(phi) * sin(phi * dphi)
		The radial term is as above (Ir) and the integral over phi is sin^2(phi)/2
  */
  auto fphi_xy = [] (double a) -> double {return std::sin(a) * std::sin(a) / 2.; };
  auto intfphi_xy = [fphi_xy] (double a1, double a2) -> double {return fphi_xy(a2) - fphi_xy(a1);}; 
  double I_xy = -h/20. * (Ir(rmax1, rmax2) - Ir(rmin1, rmin2)) * intfphi_xy(phi_s, phi_e);

  /*
	Ixz = -integral( x*z dV)
	    = -integral( r^2 dr dz) * integral( (cos_phi) dphi)
		= - h/3 * integral (r(t)^3 * (t*h - h/2) dt) * sin(phi)
		= -h^2/3 * 1/40*(3*r2^3+ r1*r2^2 - r1^2*r2 - 3*r1^3) * sin(phi)
   */
  auto Ir_xz = [] (double r1, double r2) -> double  {return 1/40.*(3*r2*r2*r2 + r1*r2*r2 - r1*r1*r2 - 3*r1*r1*r1); };
  auto fphi_xz = [] (double phi) -> double {return std::sin(phi);};
  auto intfphi_xz = [fphi_xz] (double a1, double a2) -> double {return fphi_xz(a2) - fphi_xz(a1); };
  
  double I_xz = -h*h/3. * (Ir_xz(rmax1, rmax2) - Ir_xz(rmin1, rmin2)) * intfphi_xz(phi_s, phi_e);

  /*
	Iyz = -integral( y*z dV)
	= -integral( r^2 dr dz) * integral( (sin_phi) dphi)
		= - h/3 * integral (r(t)^3 * (t*h - h/2) dt) * sin(phi)
		= -h^2/3 * 1/40*(3*r2^3+ r1*r2^2 - r1^2*r2 - 3*r1^3) * (-cos(phi))
	 SO same as I_xz, except for the phi integration
  */
  auto Ir_yz = Ir_xz;
  auto fphi_yz = [] (double phi) -> double {return -std::cos(phi);};
  auto intfphi_yz = [fphi_yz] (double a1, double a2) -> double {return fphi_yz(a2) - fphi_yz(a1); };
  
  double I_yz = -h*h/3. * (Ir_yz(rmax1, rmax2) - Ir_yz(rmin1, rmin2)) * intfphi_yz(phi_s, phi_e);

  double volume = solid->GetCubicVolume();

  // center of mass
  G4ThreeVector cm = G4ThreeVector(0, 0, 0);
  if (volume != 0) {
	double radial_xy = 1./12.*(-rmin2*rmin2*rmin2 - rmin1*rmin2*rmin2 - rmin2*rmin1*rmin1
                             - rmin1*rmin1*rmin1 + rmax2*rmax2*rmax2 
                             + rmax1*rmax2*rmax2 + rmax2*rmax1*rmax1 + rmax1*rmax1*rmax1) * h;
						 
	double xcm = radial_xy * (std::sin(phi_e) - std::sin(phi_s));
	double ycm = radial_xy * (std::cos(phi_s) - std::cos(phi_e));
	double zcm = 1/24. * (-rmin2*rmin2 + rmin1*rmin1 + rmax2*rmax2 - rmax1*rmax1) * h*h * deltaPhi;
	
	cm = G4ThreeVector(xcm, ycm, zcm) / volume;
  }

  moments.cm = cm;
  moments.M = InertiaMatrix(I_xx, I_xy, I_xz,
							I_xy, I_yy, I_yz,
							I_xz, I_yz, I_zz);
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;

}

Moments TubeMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Tubs *solid = (G4Tubs *) f_solid;

  double rmin = solid->GetInnerRadius();
  double rmax = solid->GetOuterRadius();
  double h = solid->GetZHalfLength();
  h *= 2;

  double deltaPhi = solid->GetDeltaPhiAngle();
  double phi_s = solid->GetStartPhiAngle();
  double phi_e = solid->GetDeltaPhiAngle() + phi_s;

  /*
	Cylindrical coordinates:
	x = r cos(phi)
	z = r sin(phi)
	z = z

	dV = r dr dz dphi

	as in Cone section, we calculate integrals of binomials first

	xx = integral( x^2 dV ) = integral( r^2 cos^2(phi) r dr dz dPhi

	   = (r^4/4) * (z) * (phi/2 + sin(2*phi)/4)

	yy = integral( y^2 dV) = integral( r^2 sin^2(phi) r dr dz dPhi
	   = (r^4/4) * (z) * (phi/2 - sin(2*phi)/4)

	zz = integral( z^2 dV) = integral( r dr z^2 dz dPhi
	   = (r^2/2) * (z^3/3) * deltaPhi

	xy = integral( r^2 r dr dz sin(phi) cos(phi) dphi
	   = (r^4/4) * (z) * (sin^2(phi)/2)

	xz = integral( r^2 dr z dz cos(phi) dphi
	   = (r^3/3) * (z^2/2) * (sin(phi))

	yz = integral( r^2 dr z dz sin(phi) dphi
	   = (r^3/3) * (z^2/2) * (-cos(phi))
   */

  // for I_zz = xx + yy, the phi term results in deltaPhi
  auto r4 = [] (double r) -> double {return r*r*r*r/4; };
  auto intr4 = [r4] (double r1, double r2) -> double {return r4(r2) - r4(r1); };

  auto r3 = [] (double r) -> double {return r*r*r/3; };
  auto intr3 = [r3] (double r1, double r2) -> double {return r3(r2) - r3(r1); };
  
  auto r2 = [] (double r) -> double {return r*r/2; };
  auto intr2 = [r2] (double rs, double re) -> double {return r2(re) - r2(rs); };

  auto z2 = [] (double z) -> double {return z*z/2; };
  auto intz2 = [z2] (double zs, double ze) -> double {return z2(ze) - z2(zs); };

  auto z3 = [] (double z) -> double {return z*z*z/3; };
  auto intz3 = [z3] (double zs, double ze) -> double {return z3(ze) - z3(zs); };

  auto fphi_xx = [] (double phi) -> double {return phi/2 + std::sin(2*phi)/4;};
  auto intfphi_xx = [fphi_xx] (double a1, double a2) -> double {return fphi_xx(a2) - fphi_xx(a1); };
  
  auto fphi_yy = [] (double phi) -> double {return phi/2 - std::sin(2*phi)/4;};
  auto intfphi_yy = [fphi_yy] (double a1, double a2) -> double {return fphi_yy(a2) - fphi_yy(a1); };

  auto fphi_xy = [] (double phi) -> double {return std::sin(phi)*std::sin(phi)/2;};
  auto intfphi_xy = [fphi_xy] (double a1, double a2) -> double {return fphi_xy(a2) - fphi_xy(a1); };

  double I_zz = deltaPhi * h * intr4(rmin, rmax);

  double xx = intr4(rmin, rmax) * h * intfphi_xx(phi_s, phi_e);
  double yy = intr4(rmin, rmax) * h * intfphi_yy(phi_s, phi_e);
  double zz = intz3(-h/2, h/2) * intr2(rmin, rmax) * deltaPhi;
  double xy = intr4(rmin, rmax) * h * intfphi_xy(phi_s, phi_e);
  double xz = intr3(rmin, rmax) * intz2(-h/2, h/2) * (std::sin(phi_e) - std::sin(phi_s));
  double yz = intr3(rmin, rmax) * intz2(-h/2, h/2) * (-std::cos(phi_e) + std::cos(phi_s));

  double I_xx = yy + zz;
  double I_yy = xx + zz;
  double I_xy = -xy;
  double I_xz = -xz;
  double I_yz = -yz;
  
  /*
	I found out after testing that G4Tubs gives 0 for the center of mass. Clearly wrong!
	So I add here the calculation of the center of mass.

	V*xcm = integral(x dV) = integral( r^2 dr dz cos(phi) dphi)
	      = (r^3/3) * h * sin(phi)

	V*ycm = integral(y dV) = integral( r^2 dr dz sin(phi) dphi)
	      = (r^3/3) * h * (-cos(phi))

  
	V*zcm = 0, since the Tube is symmetric about the x-y plane (goes from z=-h/2 to z=h/2
		  
   */
  
  double volume = solid->GetCubicVolume();

  double xcm = intr3(rmin, rmax) * h * (std::sin(phi_e) - std::sin(phi_s));
  double ycm = intr3(rmin, rmax) * h * (-std::cos(phi_e) + std::cos(phi_s));

  // moments.cm = solid->GetCenterOfMass();
  moments.cm = 1./volume*G4ThreeVector(xcm , ycm ,0);
  
  moments.M = InertiaMatrix(I_xx, I_xy, I_xz,
							I_xy, I_yy, I_yz,
							I_xz, I_yz, I_zz);
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;

}

Moments TrdMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Trd *solid = (G4Trd *) f_solid;

  double x1 = solid->GetXHalfLength1();
  double x2 = solid->GetXHalfLength2();
  double y1 = solid->GetYHalfLength1();
  double y2 = solid->GetYHalfLength2();
  double h = solid->GetZHalfLength();
  h *= 2;


  /* ChatGPT gives the moments as
	     // Function to calculate the principal moments of inertia
    void calculateMomentsOfInertia(double& I_x, double& I_y, double& I_z) {
        // Moment of inertia about z-axis
        I_z = (1.0 / 12.0) * mass * (x1 * x1 + x2 * x2 + y1 * y1 + y2 * y2);
        
        // Moment of inertia about x-axis
        I_x = (1.0 / 12.0) * mass * (height * height + ((y1 + y2) * (y1 + y2)) / 2.0);
        
        // Moment of inertia about y-axis
        I_y = (1.0 / 12.0) * mass * (height * height + ((x1 + x2) * (x1 + x2)) / 2.0);
    }

	They look right, but I think I'll calculate these myself.
	dV = dx dy dz

	xx = integral(x^2 dx dy dz)

	Since the x, y edges vary along z we need to parametrize:
	x(t) = x1 + t(x2 - x1), t=0..1
	y(t) = y1 + t(y2 - y1), t=0..1
	z = -h/2 + t*(h) ==> t = (z+h/2)/h ==> dz = h dt

	so xx = integral(dz (x^3/3) * y) =
	      = h * integral(dt x(t)^3/3 - (-x)^3/3) * (y(t) - (-y(t))
		  = h * integral(dt 2/3*dx^3(t) * 2 *y(t))
		  = 4/3 * h * integral(dt (x1 + t(x2 - x1))^3 * (y1 + t(y2 - y1)))
		  = 4/3/20 * h * ((4 x2^3 + 3 x1 x2^2 + 2 x1^2 x2 + x1^3) * y2
		                  + (x2^3 + 2 x1 x2^2 + 3 x1^2 x2 + 4 x1^3) * y1 ) (from wxMaxima)

	yy = same as above, but replace x with y:

	   = 1/15 * h * ((4 y2^3 + 3 y1 y2^2 + 2 y1^2 y2 + y1^3) * x2
		                  + (y2^3 + 2 y1 y2^2 + 3 y1^2 y2 + 4 y1^3) * x1 ) (from wxMaxima)

	zz = integral( z^2 dx dy dz )
	   = 4 * h * integral( x(t) * y(t) * z^2(t) dt
	   = 4 * h * h^2/120 * ((4 x2 + x1) * y1 + (x2 + 4 x1) * y1) 

	The cross terms (xy, xz, yz) should all be zero because of reflection symmetry about the axes.
   */

  double xx = 1/15. * h * ((4*x2*x2*x2 + 3 *x1*x2*x2 + 2*x1*x1*x2 +x1*x1*x1)*y2
						   + (x2*x2*x2 + 2 *x1*x2*x2 + 3*x1*x1*x2 +4*x1*x1*x1)*y1);

  double yy = 1/15. * h * ((4*y2*y2*y2 + 3 *y1*y2*y2 + 2*y1*y1*y2 +y1*y1*y1)*x2
						   + (y2*y2*y2 + 2 *y1*y2*y2 + 3*y1*y1*y2 +4*y1*y1*y1)*x1);

  double zz = 1/30. * h*h*h * ((4*x1 + x1)*y2 + (x2 + 4*x1)*y1);

  double I_xx = yy + zz;
  double I_yy = xx + zz;
  double I_zz = xx + yy;
  
  double volume = solid->GetCubicVolume();

  double zcm = 4*h*h/12. * (x2*y2 - x1 * y1);
  moments.cm = 1/volume * G4ThreeVector(0, 0, zcm); 
  
  moments.M = InertiaMatrix(I_xx, 0, 0,
							0, I_yy, 0,
							0, 0, I_zz);
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;

}

Moments TorusMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Torus *solid = (G4Torus *) f_solid;

  double rmin = solid->GetRmin();
  double rmax = solid->GetRmax();
  double rtor = solid->GetRtor();
  double phi_s = solid->GetSPhi();
  double deltaPhi = solid->GetDPhi();
  double phi_e = phi_s + deltaPhi;


  double I_zz = deltaPhi * M_PI /4 * rtor*((4*(rmax*rmax - rmin*rmin)*rtor*rtor) +
										   3*(std::pow(rmax, 4) - std::pow(rmin, 4)));

  auto fphi_xx = [] (double phi) -> double {return phi/2 + std::sin(2*phi)/4;};
  auto intfphi_xx = [fphi_xx] (double a1, double a2) -> double {return fphi_xx(a2) - fphi_xx(a1);};

  auto fphi_yy = [] (double phi) -> double {return phi/2 - std::sin(2*phi)/4;};
  auto intfphi_yy = [fphi_yy] (double a1, double a2) -> double {return fphi_yy(a2) - fphi_yy(a1);};

  auto Iradial = [rtor] (double r) -> double {return M_PI/4. * rtor*r*r*(4*rtor*rtor + 3*r*r); };
  auto intradial = [Iradial] (double r1, double r2) -> double {return Iradial(r2) - Iradial(r1);};

  double xx = intradial(rmin, rmax) * intfphi_xx(phi_s, phi_e);
  double yy = intradial(rmin, rmax) * intfphi_yy(phi_s, phi_e);
  double zz = deltaPhi * M_PI/4. * rtor * (rmax*rmax*rmax*rmax - rmin*rmin*rmin*rmin);

  double I_xx = yy + zz;
  double I_yy = xx + zz;
  
  double I_xy = -1/4. * M_PI * rtor * (4*(rmin*rmin - rmax*rmax)*rtor*rtor + 3*(std::pow(rmin,4)-std::pow(rmax, 4)))
	* (std::cos(phi_e)*std::cos(phi_e) - std::cos(phi_s)*std::cos(phi_s))/2.;
	
  double volume = solid->GetCubicVolume();
  
  // center of mass
  double zcm = 0;
  double xcm = M_PI/4. * (4*(rmax*rmax - rmin*rmin)*rtor*rtor + (std::pow(rmax, 4) - std::pow(rmin, 4)))
	* (std::sin(phi_e) - std::sin(phi_s));

  double ycm = M_PI/4. * (4*(rmax*rmax - rmin*rmin)*rtor*rtor + (std::pow(rmax, 4) - std::pow(rmin, 4)))
	* (std::cos(phi_s) - std::cos(phi_e));
  
  G4ThreeVector CM = G4ThreeVector(xcm, ycm, zcm);
  if ( volume != 0) {
	CM /= volume;
  }
  moments.cm = CM;

  
  moments.M = InertiaMatrix(I_xx, I_xy, 0,
							I_xy, I_yy, 0,
							0, 0, I_zz);
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;

}

Moments ParaboloidMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Paraboloid *solid = (G4Paraboloid *) f_solid;

  double r1 = solid->GetRadiusMinusZ();
  double r2 = solid->GetRadiusPlusZ();
  double h = 2*solid->GetZHalfLength();

  /*
	wxMaxima had some trouble with this, but it is relatively easy to do by hand,
	so i will do so here.
	From geant4 solids web page:
	https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
	definition of parabolid is given by the parameteric equations:
	rho**2 = k1* z + k2
	k1 and k2 obtained from:
	r1**2 = k1(-h/2) + k2
	r2**2 = k1(h/2) + k2
	From which k2 = 1/2(r1^2 + r2^2)
	and k1 = (r2^2 - r1^2)/h

	so rho^2(z) = (r2^2 - r1^2) * z/h + (1/2 (r1^2 + r2^2)

	dV = rho drho dphi dz
	x = rho * cos(phi)
	y = rho * sin(phi)

	xx = integral( x^2 dV) = integral( rho^3 drho cos^2(phi) dphi dz)
	   = ( phi/2 + sin(2Phi)/2)[0, 2pi] integral( 1/4*rho(z)^4 dz )   // note square brackets indicate limits of integration     = pi integral ( (k1*z + k2)^2 dz = pi * integral (k1^2*z^2 + 2k1 *k2 *z + k2^2) dz
	   = pi * (k1^2 * z^3[-h/2, h/2] + k2^2 * z[-h/2, h/2]
	   
	= pi * h * ( 1/8 * (r2^2 - r1^2)^2  + 1/4 * (r1^2 + r2^2)^2)

	yy = xx

	so Izz = xx + yy = pi * h * (1/4 * ( r2^2 - r1^2)^2 + 1/2 (r1^2 + r2^2)^2)

	zz = integral( z^2 dV) = integral( rho drho z^2 dphi dz)
	  = 2*pi * integral( rho^2/2 z^2 dz)
	  = pi * integral ( k1*z + k2) * z^2 dz
	  = pi * k2 * z^3/3[-h/2, h/2]
	  = pi * (1/2*(r1^2 + r2^2) * h^3/24)

	we also need the zcm:
	zcm = integral(z dV)
	    = 2*pi* integral(rho drho * z dz)
		= 2 * pi integral( rho^2/2 * z dz)
		= 2 * pi * integral (k1*z + k2)/2 * z dz
		= pi * k1*z^3/3[-h/2, h/2]
		= 2*pi * (r2^2 - r1^2)/h * h^3/8
		= 1/4*pi*(r2^2 - r1^2) * h^2
   */


  double k1 = (r2*r2 - r1*r1)/h;
  double k2 = (r2*r2 + r1*r1)/2;
  double Izz = M_PI/24*h*(12*k2*k2 + k1*k1*h*h);
  double Ixx = M_PI*h/48*(12*k2*k2 + 4*h*h*k2 + k1*k1*h);
  
  double volume = solid->GetCubicVolume();
  double zcm = M_PI/12*h*h*h*k1/volume;
  
  // above are calculaetd about origin
  // moment of inertia about z-axis is same
  Ixx -= volume*zcm*zcm;
  double Iyy = Ixx;
  
  moments.Volume = volume;

  moments.cm = G4ThreeVector(0, 0, zcm);
  
  moments.M = InertiaMatrix(Ixx, 0, 0,
							0, Iyy, 0,
							0, 0, Izz);

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;

}

Moments HypeMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  G4Hype *solid = (G4Hype *) f_solid;

  double Rin = solid->GetInnerRadius();
  double Rout = solid->GetOuterRadius();
  double h = 2 * solid->GetZHalfLength();

  double inst = solid->GetInnerStereo();
  double outst = solid->GetOuterStereo();

  double ain = std::tan(inst);
  double aout = std::tan(outst);
  
  /*
	From my GDMl Workbench python implementation it seems the geometry is as follows:
	rin = inner radius at z=0
	rout = outer radius at z = 0
	The radius at any aother z is:
	rout(z) = sqrt(aout*aout*z*z + Rout*Rout) and
	rin(z) = sqrt(ain*ain*z*z + Rin*Rin)

	x = r cos(phi)
	y = r sin(phi)
	z = [-h/2, h/2]
	dV = r dr dphi dz

	xxout = integral( x*x dV) = integral( r3 dr cos(phi)^2 dphi dz)
	   = pi integral( r(z)^4/4 dz ) = pi*integral( (aout*aout*z*z + Rout*Rout)^2 dz )
	   = pi/4 integral(aout^4 z^4 + 2*aout*aout*z*z*Rout*Rout + Rout^4) dx
	   = pi/4 * (aout^4* z^5/5[-h/2,h/2] + 2*aout^2*z^3/3[-h/2, h/2] + Rout^4 z[-h/2,h/2])
	   = pi/4 * ( 2/(5*2^5) * aout^4 * h^5 + 4/(3*2^3)*aout^2*Rout^2 * h^3 +  Rout^4 * h)

	yyout = xxout

	zz = integral (z^2 * r dr * dphi * dz)
	   = 2*pi integral( r^2/2 * z^2 dz)
	   = pi * integral ( (aout*aout*z*z + Rout*Rout) * z^2 dz
	   = pi * (aout*aout* z^5/5)[-h/2, h/2] + Rout*Rout * z^3/3[-h/2, h/2])
	   = pi * ((2/(5*2^5) * aout*aout*h^5 + 2/(3*2^3) Rout*Rout*h^3)
   */
  double xxout = M_PI/960. * h*(3.*aout*aout*aout*aout*h*h*h*h +40.*Rout*Rout*aout*aout*h*h+240.*Rout*Rout*Rout*Rout);  
  double yyout = xxout;

  double zzout = M_PI * ( 2./(5*32) * aout*aout*std::pow(h, 5) + 2./(3*8)*Rout*Rout*h*h*h);

  double xxin = M_PI/960. * h*(3.*ain*ain*ain*ain*h*h*h*h +40.*Rin*Rin*ain*ain*h*h+240.*Rin*Rin*Rin*Rin);
  double yyin = xxin;

  double zzin = M_PI * ( 2./(5*32) * ain*ain*std::pow(h, 5) + 2./(3*8)*Rin*Rin*h*h*h);

  double I_zz = (xxout+yyout) - (xxin+yyin);

  double I_xx = (yyout + zzout) - (yyin + zzin);

  double I_yy = (xxout + zzout) - (xxin + zzin);

  
  double volume = solid->GetCubicVolume();
  
  moments.Volume = volume;

  moments.cm = G4ThreeVector(0, 0, 0);

  
  moments.M = InertiaMatrix(I_xx, 0, 0,
							0, I_yy, 0,
							0, 0, I_zz);

  return moments;

}

Moments TetrahedronMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };
  
  G4Tet *solid = (G4Tet *) f_solid;

  G4ThreeVector anchor, v1, v2, v3;
  auto verts = solid->GetVertices();

  std::vector<G4TriangularFacet> facets;
  G4TriangularFacet *f1 = new G4TriangularFacet(verts[0], verts[2], verts[1], ABSOLUTE);
  G4TriangularFacet *f2 = new G4TriangularFacet(verts[1], verts[2], verts[3], ABSOLUTE);
  G4TriangularFacet *f3 = new G4TriangularFacet(verts[0], verts[1], verts[3], ABSOLUTE);
  G4TriangularFacet *f4 = new G4TriangularFacet(verts[0], verts[3], verts[2], ABSOLUTE);
  G4TessellatedSolid* tess = new G4TessellatedSolid("tetrahedron");
  tess->AddFacet((G4VFacet *) f1);
  tess->AddFacet((G4VFacet *) f2);
  tess->AddFacet((G4VFacet *) f3);
  tess->AddFacet((G4VFacet *) f4);

  auto calculator = MomentsCalculator::getCalculator(tess);
  moments =  calculator->calculate();
  if (moments.Volume < 0) {
	moments.M *= -1;
	moments.Volume *= -1;
  }

  delete f1;
  delete f2;
  delete f3;
  delete f4;
  // delete tess;

  return moments;
}


Moments PolyconeMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };
  
  G4Polycone *solid = (G4Polycone *) f_solid;

  G4PolyconeHistorical *hp = solid->GetOriginalParameters();

  double totVol = 0;
  G4ThreeVector CM(0, 0, 0);
  
  for(int i=0; i < hp->Num_z_planes-1; i++) {
	G4Cons *cons = new G4Cons("cons", hp->Rmin[i], hp->Rmax[i],
							  hp->Rmin[i+1], hp->Rmax[i+1],
							  std::fabs(hp->Z_values[i+1] - hp->Z_values[i])/2.0,
							  hp->Start_angle, hp->Opening_angle);
	
	auto calculator = MomentsCalculator::getCalculator(cons);
	auto section_moments =  calculator->calculate();
	section_moments.cm += G4ThreeVector(0, 0, std::fabs(hp->Z_values[i+1]+hp->Z_values[i])/2.0);
	section_moments.inertiaAboutOrigin();
	moments.M += section_moments.Volume * section_moments.M;
	totVol += section_moments.Volume;
	CM += section_moments.Volume * section_moments.cm;

	delete cons;
  }

  // double check volume calculation
  // G4cout << "myPolycone volume: " << totVol << " geant4's: " << solid->GetCubicVolume() << G4endl;

  moments.cm = CM/totVol;
  moments.M = 1./totVol * moments.M;
  moments.Volume = totVol;

  if (CM != G4ThreeVector(0, 0, 0)) {
	moments.inertiaAboutCM();
  }

  return moments;
}

Moments EllipticalTubeMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };
  
  G4EllipticalTube *solid = (G4EllipticalTube *) f_solid;
  double volume = solid->GetCubicVolume();
  
  
  double a = solid->GetDx();
  double b = solid->GetDy();
  double h = 2*solid->GetDz();

  double Ixx = 1/4. * volume * (b*b + h*h/3.0);
  double Iyy = 1/4. * volume * (a*a + h*h/3.0);
  double Izz = 1/4. * volume * (a*a + b*b);

  InertiaMatrix II(Ixx, 0, 0,
				   0, Iyy, 0,
				   0, 0, Izz);
  
  
  moments.cm = G4ThreeVector(0, 0, 0);
  moments.M = II;
  moments.Volume = volume;

  return moments;
}

Moments EllipticalConeMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };
  
  G4EllipticalCone *solid = (G4EllipticalCone *) f_solid;
  double volume = solid->GetCubicVolume();
  
  
  double dx = solid->GetSemiAxisX();
  double dy = solid->GetSemiAxisY();
  double zmax = solid->GetZMax();
  double zcut = solid->GetZTopCut();
  double h = 2*zcut;

  // The semi-major at the bottom a0 is
  // a0 = dx*(zmax + zcut)
  // The semi major axis at the top is
  // a1 = dx*(zmax - zcut)
  // Similarly, bor the semiminro axes b0, b1:
  // b0 = dy*(zmax + zcut)
  // b1 = dy*(zmax - zcut)
  //
  // The semi axes at "height" t are given by
  // a(t) = a0 + t*(a1 - a0)  t = [0,1]
  // b(t) = b0 + t*(b1 - b0)
  // z(t) = -zcut +2*t*zcut
  // From the Elliptical tube calculation, the moments of inertial of a thin elliptical
  // slab of thickness dz are
  // dIzz = dV/4 * (a^2 + b^2),
  // dV = area*dz = (pi*a*b) * dz
  // dz = 2*zcut = h
  // ==> dIzz = h*pi/4 * a(t)*b(t)*(a(t)^2 + b(t)^2)*dt
  // The above can be integated using wxMaxima
  // dIxx = dV/4*b(t)^2 + dV*z(t)*2;   // the second piece is from the parallel axis theorem
  //      = h*pi * (1/4*a(t)*b(t)*b(t)^2 + a(t)*b(t)*z(t)^2) * dt
  // Similarly
  // dIyy = = h*pi * (1/4*a(t)*b(t)*a(t)^2 + a(t)*b(t)*z(t)^2) * dt
  //
  // Now the center of mass:
  // xcm = ycm =0
  // but zcm = 1/volume*integral (z(t) * dV) = integral (z(t) * pi * a(t) * b(t) * h * dt)
  //           = 1/volume*pi * h * integral (z(t) * a(t) * b(t) * dt)
  //           = 1/volume*pi*(a1*b1 - a0 *b0) * h * zcut)/6

  double a0 = dx*(zmax + zcut);
  double a1 = dx*(zmax - zcut);
  double b0 = dy*(zmax + zcut);
  double b1 = dy*(zmax - zcut);

  // integrating in wxMaxima
  double Izz=(M_PI*((4*a1+a0)*b1*b1*b1+(3*a1+2*a0)*b0*b1*b1
					+((2*a1+3*a0)*b0*b0+4*a1*a1*a1
					  +3*a0*a1*a1+2*a0*a0*a1+a0*a0*a0)*b1
					+(a1+4*a0)*b0*b0*b0+(a1*a1*a1+2*a0*a1*a1+3*a0*a0*a1+4*a0*a0*a0)*b0)*h)/80;

  double Ixx = (M_PI*h*(((32*a1+8*a0)*b1+(8*a1+32*a0)*b0)*zcut*zcut+(12*a1+3*a0)*b1*b1*b1
						+(9*a1+6*a0)*b0*b1*b1+(6*a1+9*a0)*b0*b0*b1+(3*a1+12*a0)*b0*b0*b0))/240;

  double Iyy = (M_PI*h*(((32*a1+8*a0)*b1+(8*a1+32*a0)*b0)*zcut*zcut
						+(12*a1*a1*a1+9*a0*a1*a1+6*a0*a0*a1+3*a0*a0*a0)*b1
						+(3*a1*a1*a1+6*a0*a1*a1+9*a0*a0*a1+12*a0*a0*a0)*b0))/240;

  double zcm = 1/volume*(M_PI*(a1*b1-a0*b0)*h*zcut)/6;

  // above is inertia about origin; we need to calculate aobut the
  // center of mass. Izz stays same, Iyy and Ixx change
  Ixx -= volume*zcm*zcm;
  Iyy -= volume*zcm*zcm;

  InertiaMatrix II(Ixx, 0, 0,
				   0, Iyy, 0,
				   0, 0, Izz);
  
  moments.cm = G4ThreeVector(0, 0, zcm);
  moments.M = II;
  moments.Volume = volume;

  if (zcm != 0) {
	moments.inertiaAboutCM();
  }

  return moments;
}
