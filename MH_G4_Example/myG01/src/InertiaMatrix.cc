#include "InertiaMatrix.hh"
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

InertiaMatrix::InertiaMatrix() :
  Ixx(0), Ixy(0), Ixz(0),
  Iyx(0), Iyy(0), Iyz(0),
    Izx(0), Izy(0), Izz(0)
{
}

InertiaMatrix::InertiaMatrix(double a, double b, double c,
							 double d, double e, double f,
							 double g, double h, double i) :
  Ixx(a), Ixy(b), Ixz(c),
  Iyx(d), Iyy(e), Iyz(f),
    Izx(g), Izy(h), Izz(i)
{
}

InertiaMatrix::InertiaMatrix(const InertiaMatrix& M):
  Ixx(M.xx()), Ixy(M.xy()), Ixz(M.xz()),
  Iyx(M.yx()), Iyy(M.yy()), Iyz(M.yz()),
  Izx(M.zx()), Izy(M.zy()), Izz(M.zz())
{
}


InertiaMatrix::InertiaMatrix(G4ThreeVector Ix, G4ThreeVector Iy, G4ThreeVector Iz):
  Ixx(Ix.x()), Ixy(Ix.y()), Ixz(Ix.z()),
  Iyx(Iy.x()), Iyy(Iy.y()), Iyz(Iy.z()),
  Izx(Iz.x()), Izy(Iy.y()), Izz(Iz.z())
{
  /// Symmetrize (arbitrarily);
  Iyx = Ixy;
  Izx = Ixz;
  Izy = Iyz;
}

InertiaMatrix::~InertiaMatrix() {}

InertiaMatrix& InertiaMatrix::operator=(const InertiaMatrix& other) {
  if (this != &other) {
	Ixx = other.Ixx;
	Ixy = other.Ixy;
	Ixz = other.Ixz;

	Iyx = other.Iyx;
	Iyy = other.Iyy;
	Iyz = other.Iyz;

	Izx = other.Izx;
	Izy = other.Izy;
	Izz = other.Izz;
  }
  return *this;
}

InertiaMatrix InertiaMatrix::rotated(G4RotationMatrix &R)
{
  InertiaMatrix m;
  m.Ixx =
	R.xx() * R.xx() * Ixx +
	R.xy() * R.xx() * Iyx +
	R.xz() * R.xx() * Izx +

	R.xx() * R.xy() * Ixy +
	R.xy() * R.xy() * Iyy +
	R.xz() * R.xy() * Izy +
	
	R.xx() * R.xz() * Ixz +
	R.xy() * R.xz() * Iyz +
	R.xz() * R.xz() * Izz;

  m.Ixy =
	R.xx() * R.yx() * Ixx +
	R.xy() * R.yx() * Iyx +
	R.xz() * R.yx() * Izx +

	R.xx() * R.yy() * Ixy +
	R.xy() * R.yy() * Iyy +
	R.xz() * R.yy() * Izy +
	
	R.xx() * R.yz() * Ixz +
	R.xy() * R.yz() * Iyz +
	R.xz() * R.yz() * Izz;

  m.Ixz =
	R.xx() * R.zx() * Ixx +
	R.xy() * R.zx() * Iyx +
	R.xz() * R.zx() * Izx +

	R.xx() * R.zy() * Ixy +
	R.xy() * R.zy() * Iyy +
	R.xz() * R.zy() * Izy +
	
	R.xx() * R.zz() * Ixz +
	R.xy() * R.zz() * Iyz +
	R.xz() * R.zz() * Izz;

  m.Iyx =
	R.yx() * R.xx() * Ixx +
	R.yy() * R.xx() * Iyx +
	R.yz() * R.xx() * Izx +

	R.yx() * R.xy() * Ixy +
	R.yy() * R.xy() * Iyy +
	R.yz() * R.xy() * Izy +
	
	R.yx() * R.xz() * Ixz +
	R.yy() * R.xz() * Iyz +
	R.yz() * R.xz() * Izz;

  m.Iyy =
	R.yx() * R.yx() * Ixx +
	R.yy() * R.yx() * Iyx +
	R.yz() * R.yx() * Izx +

	R.yx() * R.yy() * Ixy +
	R.yy() * R.yy() * Iyy +
	R.yz() * R.yy() * Izy +
	
	R.yx() * R.yz() * Ixz +
	R.yy() * R.yz() * Iyz +
	R.yz() * R.yz() * Izz;

  m.Iyz =
	R.yx() * R.zx() * Ixx +
	R.yy() * R.zx() * Iyx +
	R.yz() * R.zx() * Izx +

	R.yx() * R.zy() * Ixy +
	R.yy() * R.zy() * Iyy +
	R.yz() * R.zy() * Izy +
	
	R.yx() * R.zz() * Ixz +
	R.yy() * R.zz() * Iyz +
	R.yz() * R.zz() * Izz;
  
  m.Izx =
	R.zx() * R.xx() * Ixx +
	R.zy() * R.xx() * Iyx +
	R.zz() * R.xx() * Izx +

	R.zx() * R.xy() * Ixy +
	R.zy() * R.xy() * Iyy +
	R.zz() * R.xy() * Izy +
	
	R.zx() * R.xz() * Ixz +
	R.zy() * R.xz() * Iyz +
	R.zz() * R.xz() * Izz;

  m.Izy =
	R.zx() * R.yx() * Ixx +
	R.zy() * R.yx() * Iyx +
	R.zz() * R.yx() * Izx +

	R.zx() * R.yy() * Ixy +
	R.zy() * R.yy() * Iyy +
	R.zz() * R.yy() * Izy +
	
	R.zx() * R.yz() * Ixz +
	R.zy() * R.yz() * Iyz +
	R.zz() * R.yz() * Izz;

  m.Izz =
	R.zx() * R.zx() * Ixx +
	R.zy() * R.zx() * Iyx +
	R.zz() * R.zx() * Izx +

	R.zx() * R.zy() * Ixy +
	R.zy() * R.zy() * Iyy +
	R.zz() * R.zy() * Izy +
	
	R.zx() * R.zz() * Ixz +
	R.zy() * R.zz() * Iyz +
	R.zz() * R.zz() * Izz;
  
  return m;
}

InertiaMatrix& InertiaMatrix::rotate(G4RotationMatrix &R)
{
  auto m = (*this).rotated(R);
  (*this) = m;

  return *this;
}

InertiaMatrix InertiaMatrix::aboutCM(G4ThreeVector cm, double Volume)
{
  if (cm != G4ThreeVector(0, 0, 0)) {
	double Ixxc = Ixx - Volume * (cm.z()*cm.z() + cm.y()*cm.y());
	double Ixyc = Ixy + Volume * cm.x() * cm.y();
	double Ixzc = Ixz + Volume * cm.x() * cm.z();
	double Iyyc = Iyy - Volume * (cm.x()*cm.x() + cm.z()*cm.z());
	double Iyzc = Iyz + Volume * cm.y() * cm.z();
	double Izzc = Izz - Volume * (cm.x()*cm.x() + cm.y()*cm.y());

    InertiaMatrix mat = InertiaMatrix(Ixxc, Ixyc, Ixzc,
									  Ixyc, Iyyc, Iyzc,
									  Ixzc, Iyzc, Izzc);
	return mat;
  }
  else {
	return InertiaMatrix(*this);
  }
}

InertiaMatrix InertiaMatrix::aboutOrigin(G4ThreeVector cm, double Volume)
{
  if (cm != G4ThreeVector(0, 0, 0)) {
	double Ixxc = Ixx + Volume * (cm.z()*cm.z() + cm.y()*cm.y());
	double Ixyc = Ixy - Volume * cm.x() * cm.y();
	double Ixzc = Ixz - Volume * cm.x() * cm.z();
	double Iyyc = Iyy + Volume * (cm.x()*cm.x() + cm.z()*cm.z());
	double Iyzc = Iyz - Volume * cm.y() * cm.z();
	double Izzc = Izz + Volume * (cm.x()*cm.x() + cm.y()*cm.y());

    InertiaMatrix mat = InertiaMatrix(Ixxc, Ixyc, Ixzc,
									  Ixyc, Iyyc, Iyzc,
									  Ixzc, Iyzc, Izzc);
	return mat;
  }
  else {
	return InertiaMatrix(*this);
  }
}

InertiaMatrix & InertiaMatrix::operator*=(double a)
{
  (*this).Ixx *= a;
  (*this).Ixy *= a;
  (*this).Ixz *= a;

  (*this).Iyx *= a;
  (*this).Iyy *= a;
  (*this).Iyz *= a;

  (*this).Izx *= a;
  (*this).Izy *= a;
  (*this).Izz *= a;

  return *this;
}

InertiaMatrix& InertiaMatrix::operator+=(const InertiaMatrix &I)
{
  (*this).Ixx += I.Ixx;
  (*this).Ixy += I.Ixy;
  (*this).Ixz += I.Ixz;

  (*this).Iyx += I.Iyx;
  (*this).Iyy += I.Iyy;
  (*this).Iyz += I.Iyz;

  (*this).Izx += I.Izx;
  (*this).Izy += I.Izy;
  (*this).Izz += I.Izz;

  return *this;
}

InertiaMatrix& InertiaMatrix::operator/=(double a)
{
  (*this).Ixx /= a;
  (*this).Ixy /= a;
  (*this).Ixz /= a;

  (*this).Iyx /= a;
  (*this).Iyy /= a;
  (*this).Iyz /= a;

  (*this).Izx /= a;
  (*this).Izy /= a;
  (*this).Izz /= a;

  return *this;
}

InertiaMatrix operator* (const InertiaMatrix& I, double a)
{
  return InertiaMatrix(I.xx()*a, I.xy()*a, I.xz()*a,
					   I.yx()*a, I.yy()*a, I.yz()*a,
					   I.zx()*a, I.zy()*a, I.zz()*a					   
					   );
}

InertiaMatrix operator* (double a, const InertiaMatrix& I)
{
  return InertiaMatrix(I.xx()*a, I.xy()*a, I.xz()*a,
					   I.yx()*a, I.yy()*a, I.yz()*a,
					   I.zx()*a, I.zy()*a, I.zz()*a					   
					   );
}

InertiaMatrix operator/ (const InertiaMatrix& I, double a)
{
  return InertiaMatrix(I.xx()/a, I.xy()/a, I.xz()/a,
					   I.yx()/a, I.yy()/a, I.yz()/a,
					   I.zx()/a, I.zy()/a, I.zz()/a					   
					   );
}

std::ostream& operator<<(std::ostream& os, const InertiaMatrix& obj)
{
  os << "( ("<< obj.Ixx << ", " << obj.Ixy << ", " << obj.Ixz << ")" << std::endl;
  os << "  ("<< obj.Iyx << ", " << obj.Iyy << ", " << obj.Iyz << ")" << std::endl;
  os << "  ("<< obj.Izx << ", " << obj.Izy << ", " << obj.Izz << ") )" << std::endl;

  return os;
}
