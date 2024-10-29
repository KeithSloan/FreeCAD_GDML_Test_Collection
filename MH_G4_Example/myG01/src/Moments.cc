#include "Moments.hh"
#include <G4ThreeVector.hh>

InertiaMatrix Moments::inertiaAboutOrigin()
{
  if (cm != G4ThreeVector(0, 0, 0)) {
	double Ixx = M.xx() + Volume * (cm.z()*cm.z() + cm.y()*cm.y());
	double Ixy = M.xy() - Volume * cm.x() * cm.y();
	double Ixz = M.xz() - Volume * cm.x() * cm.z();
	double Iyy = M.yy() + Volume * (cm.x()*cm.x() + cm.z()*cm.z());
	double Iyz = M.yz() - Volume * cm.y() * cm.z();
	double Izz = M.zz() + Volume * (cm.x()*cm.x() + cm.y()*cm.y());

    InertiaMatrix mat = InertiaMatrix(Ixx, Ixy, Ixz,
									Ixy, Iyy, Iyz,
									Ixz, Iyz, Izz);

	return mat;
  }
  else {
	return M;
  }
}

InertiaMatrix Moments::inertiaAboutCM()
{
  if (cm != G4ThreeVector(0, 0, 0)) {
	double Ixx = M.xx() - Volume * (cm.z()*cm.z() + cm.y()*cm.y());
	double Ixy = M.xy() + Volume * cm.x() * cm.y();
	double Ixz = M.xz() + Volume * cm.x() * cm.z();
	double Iyy = M.yy() - Volume * (cm.x()*cm.x() + cm.z()*cm.z());
	double Iyz = M.yz() + Volume * cm.y() * cm.z();
	double Izz = M.zz() - Volume * (cm.x()*cm.x() + cm.y()*cm.y());

    InertiaMatrix mat = InertiaMatrix(Ixx, Ixy, Ixz,
									  Ixy, Iyy, Iyz,
									  Ixz, Iyz, Izz);

	return mat;
  }
  else {
	return M;
  }
}
