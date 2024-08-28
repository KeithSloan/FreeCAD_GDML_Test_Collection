#ifndef INERTIA_MATRIX_H
#define INERTIA_MATRIX_H

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

class InertiaMatrix {

  
public:
  InertiaMatrix();
  InertiaMatrix(G4ThreeVector Ix, G4ThreeVector Iy, G4ThreeVector Iz);
  InertiaMatrix(double Ixx, double Ixy, double Ixz,
				double Iyx, double Iyy, double Iyz,
				double Izx, double Izy, double Izz);

  InertiaMatrix(const InertiaMatrix& M);
  // InertiaMatrix(InertiaMatrix&& M) = default;

  ~InertiaMatrix();

  InertiaMatrix rotated(G4RotationMatrix &rot);
  // return a rotated copy
  
  InertiaMatrix& rotate(G4RotationMatrix &rot);
  // rotate the inertia matrix
    
  InertiaMatrix & operator *= (double);
  InertiaMatrix & operator /= (double);
  InertiaMatrix & operator += (const InertiaMatrix &);

  inline double xx() const;
  inline double xy() const;
  inline double xz() const;
  
  inline double yx() const;
  inline double yy() const;
  inline double yz() const;

  inline double zx() const;
  inline double zy() const;
  inline double zz() const;

  inline double max() const;

  friend std::ostream& operator<<(std::ostream& os, const InertiaMatrix& obj);
  
private:
  double Ixx, Ixy, Ixz,
	Iyx, Iyy, Iyz,
	Izx, Izy, Izz;
  
  
};

// Global methods
InertiaMatrix operator * (const InertiaMatrix &, double a);
InertiaMatrix operator * (double a, const InertiaMatrix &);
InertiaMatrix operator / (const InertiaMatrix &, double a);


#include "InertiaMatrix.icc"
#endif
