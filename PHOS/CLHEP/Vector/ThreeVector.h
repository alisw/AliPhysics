// -*- C++ -*-
// CLASSDOC OFF
// $Id$
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// Hep3Vector is a general 3-vector class defining vectors in three
// dimension using HepDouble components. Rotations of these vectors are
// performed by multiplying with an object of the HepRotation class.
//
// .SS See Also
// LorentzVector.h, Rotation.h, LorentzRotation.h
//
// .SS Authors
// Leif Lonnblad and Anders Nilsson.
//

#ifndef HEP_THREEVECTOR_H
#define HEP_THREEVECTOR_H

#ifdef GNUPRAGMA
#pragma interface
#endif

#include "CLHEP/config/CLHEP.h"

#ifdef HEP_NO_INLINE_IN_DECLARATION
#define inline
#endif

#include <iostream.h>
class HepRotation;

class Hep3Vector {

public:

  inline Hep3Vector(HepDouble x = 0.0, HepDouble y = 0.0, HepDouble z = 0.0);
  // The constructor.

  inline Hep3Vector(const Hep3Vector &);
  // The copy constructor.

  inline ~Hep3Vector();
  // The destructor.

  inline HepDouble x() const;
  inline HepDouble y() const;
  inline HepDouble z() const;
  // The components in cartesian coordinate system.

  HepDouble operator () (int) const;
  // Get components by index (Geant4).

  inline void setX(HepDouble);
  inline void setY(HepDouble);
  inline void setZ(HepDouble);
  // Set the components in cartesian coordinate system.

  inline HepDouble phi() const;
  // The azimuth angle.

  inline HepDouble theta() const;
  // The polar angle.

  inline HepDouble cosTheta() const;
  // Cosine of the polar angle.

  inline HepDouble mag2() const;
  // The magnitude squared (rho^2 in spherical coordinate system).

  inline HepDouble mag() const;
  // The magnitude (rho in spherical coordinate system).

  inline void setPhi(HepDouble);
  // Set phi keeping mag and theta constant (BaBar).

  inline void setTheta(HepDouble);
  // Set theta keeping mag and phi constant (BaBar).

  inline void setMag(HepDouble);
  // Set magnitude keeping theta and phi constant (BaBar).

  inline HepDouble perp2() const;
  // The transverse component squared (R^2 in cylindrical coordinate system).

  inline HepDouble perp() const;
  // The transverse component (R in cylindrical coordinate system).

  inline HepDouble perp2(const Hep3Vector &) const;
  // The transverse component w.r.t. given axis squared.

  inline HepDouble perp(const Hep3Vector &) const;
  // The transverse component w.r.t. given axis.

  inline Hep3Vector & operator = (const Hep3Vector &);
  // Assignment.

  inline HepBoolean operator == (const Hep3Vector &) const;
  inline HepBoolean operator != (const Hep3Vector &) const;
  // Comparisons (Geant4). 

  inline Hep3Vector & operator += (const Hep3Vector &);
  // Addition.

  inline Hep3Vector & operator -= (const Hep3Vector &);
  // Subtraction.

  inline Hep3Vector operator - () const;
  // Unary minus.

  inline Hep3Vector & operator *= (HepDouble);
  // Scaling with real numbers.

  inline Hep3Vector unit() const;
  // Unit vector parallel to this.

  inline Hep3Vector orthogonal() const;
  // Vector orthogonal to this (Geant4).

  inline HepDouble dot(const Hep3Vector &) const;
  // Scalar product.

  inline Hep3Vector cross(const Hep3Vector &) const;
  // Cross product.

  inline HepDouble angle(const Hep3Vector &) const;
  // The angle w.r.t. another 3-vector.

  void rotateX(HepDouble);
  // Rotates the Hep3Vector around the x-axis.

  void rotateY(HepDouble);
  // Rotates the Hep3Vector around the y-axis.

  void rotateZ(HepDouble);
  // Rotates the Hep3Vector around the z-axis.

  void rotateUz(Hep3Vector&);
  // Rotates reference frame from Uz to newUz (unit vector) (Geant4).

  void rotate(HepDouble, const Hep3Vector &);
  // Rotates around the axis specified by another Hep3Vector.

  Hep3Vector & operator *= (const HepRotation &);
  Hep3Vector & transform(const HepRotation &);
  // Transformation with a Rotation matrix.

private:

  HepDouble dx, dy, dz;
  // The components.
};

#ifdef HEP_NO_INLINE_IN_DECLARATION
#undef inline
#endif

ostream & operator << (ostream &, const Hep3Vector &);
// Output to a stream.

extern const Hep3Vector HepXHat, HepYHat, HepZHat;

#ifdef HEP_SHORT_NAMES
typedef Hep3Vector Vector3;
typedef Hep3Vector DVector3;
typedef Hep3Vector FVector3;
static const Hep3Vector & xhat = HepXHat;
static const Hep3Vector & yhat = HepYHat;
static const Hep3Vector & zhat = HepZHat;
#endif
typedef Hep3Vector HepThreeVectorD;
typedef Hep3Vector HepThreeVectorF;

#ifdef HEP_DEBUG_INLINE

Hep3Vector operator + (const Hep3Vector &, const Hep3Vector &);
// Addition of 3-vectors.

Hep3Vector operator - (const Hep3Vector &, const Hep3Vector &);
// Subtraction of 3-vectors.

HepDouble operator * (const Hep3Vector &, const Hep3Vector &);
// Scalar product of 3-vectors.

Hep3Vector operator * (const Hep3Vector &, HepDouble a);
Hep3Vector operator * (HepDouble a, const Hep3Vector &);
// Scaling of 3-vectors with a real number

#else
#include "CLHEP/Vector/ThreeVector.icc"
#endif

#endif /* HEP_THREEVECTOR_H */
