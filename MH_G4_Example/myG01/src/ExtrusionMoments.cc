#include "ExtrusionMoments.hh"
#include "InertiaMatrix.hh"
#include <G4ExtrudedSolid.hh>
#include <G4TwoVector.hh>
#include <G4ios.hh>
#include <vector>

Moments ExtrusionMoments::calculate() {
  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };
  
  G4ExtrudedSolid *solid = (G4ExtrudedSolid *) f_solid;

  std::vector<G4TwoVector> polygon = solid->GetPolygon();

  int nsections = solid->GetNofZSections();
  InertiaMatrix dI_dz = computeInertiaMatrix(polygon);

  double area = std::fabs(polygonArea(polygon));
  
  InertiaMatrix II = InertiaMatrix();
  G4TwoVector cm2D = polygonCM(polygon);

  // G4cout << "cm2D: " << cm2D << G4endl;
  // G4cout << "area: " << area << G4endl;
  // G4cout << "I_per_z: " << dI_dz << G4endl;
  
  double Izz = 0;
  double Ixx = 0;
  double Iyy = 0;
  double Ixy = 0;
  double Ixz = 0;
  double Iyz = 0;

  G4ThreeVector cm = G4ThreeVector(0, 0, 0);
  double totalVolume = 0;
  
  // step 1, just take scaling into account
  for (int isection = 0; isection < nsections - 1; isection++) {
	G4ExtrudedSolid::ZSection zsection0 = solid->GetZSection(isection);
	G4ExtrudedSolid::ZSection zsection1 = solid->GetZSection(isection+1);
	double f0 = zsection0.fScale;
	double f1 = zsection1.fScale;
	double z0 = zsection0.fZ;
	double z1 = zsection1.fZ;
	G4TwoVector offset0 = zsection0.fOffset;
	G4TwoVector offset1 = zsection1.fOffset;
	
	auto meancm = [f0, f1] (double x0, double x1) -> double {
	  return 1./12. * ((3*f1*f1 + 2*f0*f1 + f0*f0) * x1 + (f1*f1 + 2*f0*f1 + 3*f0*f0)*x0);};

	/*
	 Take height of zsection from z0 to z1 into account, plus scaling factor
	 scalingFactor(t) = scalingFactor(z=z0) + t*(scalingfactor(z=z1) - scalingFactor(z=z0)
	 f(t) = f0 + t*(f1 - f0)
	 z(t) = z0 + t*(z1-z0)
	 Offsets in x, y also vary linearly:
	 ofx = ofx0 + t*(ofx1 - ofx0)
	 ofy = ofy0 + t*(ofy1 - ofy0)
	 since the area scales as the square of the scaling factor,
	 I_per_z = computeInertiaMatrix() computes the inertia matrix per unit length along z
	 For Izz the calculation is:
	 dI_dz(z) = dI_dz(0) * f^2(z) + area0*f^2(z) * ((xcm + ofx(z))^2 + (ycm + ofy(z))^2)
	 Now z(t) = z0 + t*(z1-z0) ==> dz = (z1-z0)*dt
	 So for a section between z=z0 and z=z1 we have
	 Izz = integral(dI_dz(z) * dz)[z0, z1]
	     = dI_dz * (z1-z0) integral( f^2(t) dt)[0,1]
		        + (z1-z0)* area0 * integral( f^2(t) * ((xcm + ofx(t) )^2 + (ycm + ofy(t)^2) dt )) [0,1]
	 From wxMaxima:
	 integral(f^4(t) dt)[0,1] = (f1^4+f0*f1^3+f0^2*f1^2+f0^3*f1+f0^4)/5
	 Each term in the parallel axis theorem integration has the form:
	*/
	auto intc2 = [f0, f1] (double xcm, double ofx0, double ofx1) -> double { return
		1./30*((10*f1*f1+10*f0*f1+10*f0*f0)*xcm*xcm
			   +((15*f1*f1+10*f0*f1+5*f0*f0)*ofx1 + (5*f1*f1+10*f0*f1+15*f0*f0)*ofx0)*xcm
			   +(6*f1*f1+3*f0*f1+f0*f0)*ofx1*ofx1
			   +(3*f1*f1+4*f0*f1+3*f0*f0)*ofx0*ofx1
			   +(f1*f1+3*f0*f1+6*f0*f0)*ofx0*ofx0); };
	  
	// double fScale = integral( f(t)^4 dt)
	double fScale = (std::pow(f1, 4) + f0*std::pow(f1, 3) + f0*f0*f1*f1 + f1*std::pow(f0, 3) + std::pow(f0, 4))/5.;
	Izz += dI_dz.zz() * (z1-z0) * fScale +
	  (z1-z0) * area * (intc2(cm2D.x(), offset0.x(), offset1.x())
						+ intc2(cm2D.y(), offset0.y(), offset1.y()));
	
	// G4cout << " fscale: " << fScale << G4endl;
	/*
	  Now the XX and yy moments:
	  dI_dz calculated by computeInertiaMatrix is for an infinitely thin slab lying in
	  the x-y plane. To take into account the xtrusion into the z-direction and the offset
	  in the x-y direction, we use the parallel axis theorem:
	  Ixx(z) = integral(dI_dz(0)*f^4(z) dz  + area(z)*dz* ((ycm+ofy(z))^2 + z^2))[z=z0, z1]
	   = (z1-z0) * (dI_dz.xx integral(f^2(t) dt) + area0 integral(f^2(t) * ((ycm+ofy(t))^2 + z(t)^2)[0,1]

	  The same holds for Iyy
	 */

	double intz2sf2 = ((6*f1*f1+3*f0*f1+f0*f0)*z1*z1+(3*f1*f1+4*f0*f1+3*f0*f0)*z0*z1+(f1*f1+3*f0*f1+6*f0*f0)*z0*z0)/30.;
	
	// double intz2sf2 = 1./30.*((6*f1*f1+3*f0*f1)*z1*z1 +
	//						  (3*f1*f1 + 4*f0*f1+3*f0*f0)*z0*z1 +
	//						  (f1*f1+3*f0*f1+6*f0*f0)*z0*z0);

	Ixx += (z1-z0)*dI_dz.xx()*fScale +
	  (z1-z0)*area*(intc2(cm2D.y(), offset0.y(), offset1.y()) + intz2sf2);

	Iyy += (z1-z0)*dI_dz.yy()*fScale +
	  (z1-z0)*area*(intc2(cm2D.x(), offset0.x(), offset1.x()) + intz2sf2);
	
	double volume = (z1-z0)*area*fScale;
	// zcm = 1/volume * (z1-z0) * integal(z dV)
	//     = 1/volume * (z1-z0) * integral( z(t) * sf(t)^2 * area0 dt)
	//       = area0/volume*(z1-z0) integral( (z0 + t*(z1-z0) * (sf0 + t(sf1-sf0))^2 dt)
	totalVolume += volume;

	// G4cout << "z0: " << z0 << " z1: " << z1 << G4endl;
	// double zcm_i = meancm(z0, z1);
	double Dz = (z1-z0);
	double Df = (f1-f0);
	
	double zcm_i = area/volume*Dz * (z0* (f0*f0 + f0*Df + Df*Df/3) + f0*f0*Dz/2 + 2./3.*f0*Df*Dz + Df*Df*Dz/4);
	double xcm_i = cm2D.x() + meancm(offset0.x(), offset1.x());
	double ycm_i = cm2D.y() + meancm(offset0.y(), offset1.y());
	// G4cout << "Layer cm: " << G4ThreeVector(xcm_i, ycm_i, zcm_i) << G4endl;
	cm += volume*G4ThreeVector(xcm_i, ycm_i, zcm_i);


	double dx0 = cm2D.x() + offset0.x();
	double dx = offset1.x() - offset0.x();
	// dx(t) = dx0 + t*dx; // x-center of layer at z = z0 + t*(z1-z0) from the origin
	double dy0 = cm2D.y() + offset0.y();
	double dy = offset1.y() - offset0.y();
	// dy(t) = dy0 + t*dy; // y-center of layer at z = z0 + t*(z1-z0) from the origin

	// integral(dx(t)*dy(t)*f(t) * f^2(t) dt)  // The integral of f^2(t) is to account for the change of area as the scale factor changes along z
	
	double int_xyf2 = (((20*dx0+15*dx)*dy0+(15*dx0+12*dx)*dy)*f1*f1
					   +((20*dx0+10*dx)*dy0+(10*dx0+6*dx)*dy)*f0*f1
					   +((20*dx0+5*dx)*dy0+(5*dx0+2*dx)*dy)*f0*f0)/60;
	  
	Ixy -= Dz*area*int_xyf2;  // parallel axis theorem applied to off diagonal elements	
	// integral (dx(t), z(t), f(t)) [0,1]
	double int_xzf2 = (((15*dx0+12*dx)*f1*f1+(10*dx0+6*dx)*f0*f1
						+(5*dx0+2*dx)*f0*f0)*z1
					   +((5*dx0+3*dx)*f1*f1
						 +(10*dx0+4*dx)*f0*f1+(15*dx0+3*dx)*f0*f0)*z0)/60;
	
	Ixz -= Dz*area*int_xzf2;  // parallel axis theorem applied to off diagonal elements

	double int_yzf2 = (((15*dy0+12*dy)*f1*f1
						+(10*dy0+6*dy)*f0*f1+(5*dy0+2*dy)*f0*f0)*z1
					   +((5*dy0+3*dy)*f1*f1+(10*dy0+4*dy)*f0*f1
						 +(15*dy0+3*dy)*f0*f0)*z0)/60;
	
	Iyz -= Dz*area*int_yzf2;  // parallel axis theorem applied to off diagonal elements	
  }

  

  double volume = solid->GetCubicVolume();
  // Double check:
  // G4cout << "geant volume: " << volume << " my volume: " << totalVolume << G4endl;

  moments.cm = 1/volume * cm; 
  
  moments.M = InertiaMatrix(Ixx, Ixy, Ixz,
							Ixy, Iyy, Iyz,
							Ixz, Iyz, Izz);
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;


}
