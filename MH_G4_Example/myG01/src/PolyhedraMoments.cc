#include "PolyhedraMoments.hh"
#include "G4Polyhedra.hh"
#include "InertiaMatrix.hh"
#include "Moments.hh"
#include <G4PolyhedraHistorical.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <cmath>

Moments PolyhedraMoments::iMoment(int iz)
{
  G4Polyhedra *solid = (G4Polyhedra *) f_solid;
  G4PolyhedraHistorical * phistortical = solid->GetOriginalParameters();
  double phistart = solid->GetStartPhi();
  
  double deltaPhi = phistortical->Opening_angle/phistortical->numSide;
  double tanDphi = tan(deltaPhi/2);
  // It turns out that the G4PolyhedraHistorical stores to corner distances, not the
  // tangent distances that are described in the manual and eneterd into the gdml file
  // rather than rederive all the equations in terms of the corner radii, I will
  // simply compute the tangent radii by multiplying by cos(deltaPhi/2)
  double cosDphi = cos(deltaPhi/2);
  double rmin0 = phistortical->Rmin[iz] * cosDphi;
  double rmin1 = phistortical->Rmin[iz+1] * cosDphi;
  double rmax0 = phistortical->Rmax[iz] * cosDphi;
  double rmax1 = phistortical->Rmax[iz+1] * cosDphi;
  double z0 = phistortical->Z_values[iz];
  double z1 = phistortical->Z_values[iz+1];
  double h = z1-z0;

  double rmin0_2 = rmin0*rmin0;
  double rmin0_3 = rmin0_2*rmin0;
  double rmin0_4 = rmin0_3*rmin0;

  double rmax0_2 = rmax0*rmax0;
  double rmax0_3 = rmax0_2*rmax0;
  double rmax0_4 = rmax0_3*rmax0;

  double rmin1_2 = rmin1*rmin1;
  double rmin1_3 = rmin1_2*rmin1;
  double rmin1_4 = rmin1_3*rmin1;

  double rmax1_2 = rmax1*rmax1;
  double rmax1_3 = rmax1_2*rmax1;
  double rmax1_4 = rmax1_3*rmax1;

  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  // The following expressions come
  // from considering a trapezoid in the x-y plane. The trapezoid corners have coordinates
  // x1 = rmin, y1 = -rmin/tan(Dphi/2)
  // x2 = rmin, y2 = rmin/tan(Dphi/2)
  // x3 = rmax, y3 = -rmax/tan(Dphi/2)
  // x4 = rmax, y4 = rmax/tan(Dphi/2)
  // we make basic use of the fact that the moment of inertia of a thin rod of length l about
  // its center is 1/2 m*l^2
  // The lenght l at x is l = 2*x*tan(Dphi/2)
  // dm = density * dV = 1 * l dx dz
  // For a section between z0=z[i] and z1=z[i+1] the coordinates can be parameterized
  // as follows
  // rmin(t) = rmin0 + t*(rmin1 - rmin0)
  // rmax(t) = rmax0 + t*(rmax1 - rmax0)
  // z(t) = z0 + t*(z1-z0) = z0 + h*t
  // dz = h * dt
  // for Izz, d^2Izz = 1/12 * l^2 * dV + dV * x^2  (second term from parallel axis theorem)
  //    = 2/3 inetegral ((x^3 * tan^3(Dphi/2) dx ) + integral(x^3 * tan(Dphi/2) * dx)) * dz)
  //    = (2/3*tan^3(Dphi/2) + 2*tan(Dphi/2)) * 1/4 integral(rmax^4(t) - rmin^4(t), dt)
  // Ixx = integral( 1/2 dm * l^2 + dm * z^2) // second piece from parallel axis theorem
  //     = h/6*tan^3(Dphi/2) * integral(rmax^4(t) - rmin^4(t)) dt +
  //              h*tan(Dphi/2) integral(rmax^2(t) - rmin^2(t) dt
  // Iyy = integral dm * (x^2 + y^2)
  //     = h * tan(Dphi/2) * integral1/2*(rmax^4(t) - rmin^4(t) dt)
  //             + h * integral(rmax^2(t) - rmin^2(t)) *z^2(t) dt
  //
  // Volume:
  // dV = Area * dz
  // trapezoid area = (lmax + lmin)/2 * (rmax - rmin) = tan(Dphi/2) * (rmax + rmin)*(rmax - rmin)
  //                = tan(Dphi/2) * (rmax^2(t) - rmin^2(t)
  // so volume = tan(Dphi/2) * h * integral((rmax^2(t) - rmin^2(t) dt)
  //
  // center of mass:
  // zcm = 1/volume * integral (z dV) = h * integral (z(t) * Area(t) ) dt
  //                                  = h * tan(Dphi/2) * integral(z(t) * (rmax^2(t) - rmin^2(t)) dt
  //xcm = integral (x dm) / volume
  //      = 1/volume integral( x l dx dz) = integral (2 x^2 * tan(Dphi/2) dx dz)
  //      = 2*h/3*tan(Dphi/2) * integar( rmax^3(t) - rmin^3(t)) dt
  // ycm = 0
  //
  // Ixy = 0
  // Ixz = 0.
  // Iyx = 0
  // wxMaxim was used to do the symbolic integrals
  // The above is for a segment symmetric about the y axis
  // TO get the inertia of the rotated section we use:
  // I(rotated) = R * I R^T,
  // where R is a rotation matrix to azimuthal angle phi
  //
  double II_integral = -((h*(rmin1_4+rmin0*rmin1_3+rmin0_2*rmin1_2+rmin0_3*rmin1+rmin0_4-rmax1_4-rmax0*rmax1_3-rmax0_2*rmax1_2-rmax0_3*rmax1-rmax0_4))/20);

  double Izz =  (2./3.*tanDphi*tanDphi*tanDphi + 2 * tanDphi) * II_integral;

  double Ixx1 = -((h*(rmin1_4+rmin0*rmin1_3+rmin0_2*rmin1_2+rmin0_3*rmin1+rmin0_4-rmax1_4-rmax0*rmax1_3-rmax0_2*rmax1_2-rmax0_3*rmax1-rmax0_4)*tanDphi*tanDphi*tanDphi)/30);

  double h_2 = h * h;
  double Ixx2 = -((h*tanDphi*((10*rmin1_2+10*rmin0*rmin1+10*rmin0_2-10*rmax1_2-10*rmax0*rmax1-10*rmax0_2)*z0*z0+(15*h*rmin1_2+10*h*rmin0*rmin1+5*h*rmin0_2-15*h*rmax1_2-10*h*rmax0*rmax1-5*h*rmax0_2)*z0+6*h_2*rmin1_2+3*h_2*rmin0*rmin1+h_2*rmin0_2-6*h_2*rmax1_2-3*h_2*rmax0*rmax1-h_2*rmax0_2))/30);

  double Ixx = Ixx1 + Ixx2;

  double Iyy1 = -((h*(rmin1_4+rmin0*rmin1_3+rmin0_2*rmin1_2+rmin0_3*rmin1+rmin0_4-rmax1_4-rmax0*rmax1_3-rmax0_2*rmax1_2-rmax0_3*rmax1-rmax0_4)*tanDphi)/10);

  double  Iyy2 = -((h*tanDphi*((10*rmin1_2+10*rmin0*rmin1+10*rmin0_2-10*rmax1_2-10*rmax0*rmax1-10*rmax0_2)*z0*z0+(15*h*rmin1_2+10*h*rmin0*rmin1+5*h*rmin0_2-15*h*rmax1_2-10*h*rmax0*rmax1-5*h*rmax0_2)*z0+6*h_2*rmin1_2+3*h_2*rmin0*rmin1+h_2*rmin0_2-6*h_2*rmax1_2-3*h_2*rmax0*rmax1-h_2*rmax0_2))/30);

  double Iyy = Iyy1 + Iyy2;

  double Ixz = (h*tanDphi*((5*rmin1_3+5*rmin0*rmin1_2+5*rmin0_2*rmin1+5*rmin0_3-5*rmax1_3-5*rmax0*rmax1_2-5*rmax0_2*rmax1-5*rmax0_3)*z0+4*h*rmin1_3+3*h*rmin0*rmin1_2+2*h*rmin0_2*rmin1+h*rmin0_3-4*h*rmax1_3-3*h*rmax0*rmax1_2-2*h*rmax0_2*rmax1-h*rmax0_3))/30;
  
  double volume_integral = -((h*(rmin1_2+rmin0*rmin1+rmin0_2-rmax1_2-rmax0*rmax1-rmax0_2))/3);
  moments.Volume = phistortical->numSide * tanDphi * volume_integral;

  double zcm_integral = -((h*((4*rmin1_2+4*rmin0*rmin1+4*rmin0_2-4*rmax1_2-4*rmax0*rmax1-4*rmax0_2)*z0+3*h*rmin1_2+2*h*rmin0*rmin1+h*rmin0_2-3*h*rmax1_2-2*h*rmax0*rmax1-h*rmax0_2))/12);
  
  double zcm = tanDphi * zcm_integral;

  // note this is NOT xcm, but integral( x dV);
  double xcm = -((h*(rmin1_3+rmin0*rmin1_2+rmin0_2*rmin1+rmin0_3-rmax1_3-rmax0*rmax1_2-rmax0_2*rmax1-rmax0_3)*tanDphi)/6);
  
  auto II = InertiaMatrix(Ixx, 0, Ixz,
						  0, Iyy, 0,
						  Ixz, 0, Izz);

  G4ThreeVector axis = G4ThreeVector(0, 0, 1);

  // The above calculations are for a section centerd at phi = 0; To get the mements
  // for the numSide sections we need to rotate to phi = phiStart + isiDe * deltaPhi;
  for(int j = 0; j < phistortical->numSide; j++) {
	double phi = phistart + deltaPhi/2 + j*deltaPhi;
	G4RotationMatrix R = G4RotationMatrix(axis, phi);
	InertiaMatrix IIRotated = II.rotated(R);
	moments.M += IIRotated;
	moments.cm += R*G4ThreeVector(xcm, 0, zcm);
  }


  moments.cm /= moments.Volume;

  return moments;
  
}

Moments PolyhedraMoments::calculate()
{
  G4Polyhedra *solid = (G4Polyhedra *) f_solid;

  Moments moments = {
	0,  //vol
	G4ThreeVector(0, 0, 0), // cm
	InertiaMatrix()
  };

  double myVolume = 0;

  G4ThreeVector weighted_cm = G4ThreeVector(0, 0, 0);
  
  G4PolyhedraHistorical *phistortical = solid->GetOriginalParameters();
  for(int i=0; i < phistortical->Num_z_planes-1; i++) {
	Moments section_moments = iMoment(i);
	moments.M += section_moments.M;
	weighted_cm += section_moments.Volume * section_moments.cm;
	myVolume += section_moments.Volume;
  }


  double volume = f_solid->GetCubicVolume();
  
  moments.cm = 1/volume * weighted_cm;
  
  moments.Volume = volume;

  if (moments.cm != G4ThreeVector(0, 0, 0)) {
	moments.M = moments.inertiaAboutCM();
  }

  return moments;


}
