// -*- C++ -*-
// CLASSDOC OFF
// $Id$
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// HepLorentzVector is a Lorentz vector consisting of Hep3Vector and
// HepDouble components. Lorentz transformations (rotations and boosts)
// of these vectors are perfomed by multiplying with objects of
// the HepLorenzRotation class.
//
// .SS See Also
// ThreeVector.h, Rotation.h, LorentzRotation.h
//
// .SS Authors
// Leif Lonnblad and Anders Nilsson. Modified by Evgueni Tcherniaev.
//

#ifndef HEP_LORENTZVECTOR_H
#define HEP_LORENTZVECTOR_H

#ifdef GNUPRAGMA
#pragma interface
#endif

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#ifdef HEP_NO_INLINE_IN_DECLARATION
#define inline
#endif

#include <iostream.h>
class HepLorentzRotation;

class HepLorentzVector {

public:

  inline HepLorentzVector(HepDouble x = 0.0, HepDouble y = 0.0,
			  HepDouble z = 0.0, HepDouble t = 0.0);
  // Constructor giving the components x, y, z, t.

  inline HepLorentzVector(const Hep3Vector &, HepDouble);
  // Constructor giving a 3-Vector and a time component.

  inline HepLorentzVector(const HepLorentzVector &);
  // Copy constructor.

  inline ~HepLorentzVector();
  // The destructor.

  inline operator Hep3Vector () const;
  inline operator Hep3Vector & ();
  // Conversion (cast) to Hep3Vector.

  inline HepDouble x() const;
  inline HepDouble y() const;
  inline HepDouble z() const;
  inline HepDouble t() const;
  // Get position and time.

  inline void setX(HepDouble);
  inline void setY(HepDouble);
  inline void setZ(HepDouble);
  inline void setT(HepDouble);
  // Set position and time.

  inline HepDouble px() const;
  inline HepDouble py() const;
  inline HepDouble pz() const;
  inline HepDouble e() const;
  // Get momentum and energy.

  inline void setPx(HepDouble);
  inline void setPy(HepDouble);
  inline void setPz(HepDouble);
  inline void setE(HepDouble);
  // Set momentum and energy.

  inline Hep3Vector vect() const;
  // Get spatial component. 

  inline void setVect(const Hep3Vector &);
  // Set spatial component. 

  inline HepDouble theta() const;
  inline HepDouble cosTheta() const;
  inline HepDouble phi() const;
  inline HepDouble rho() const;
  // Get spatial vector components in spherical coordinate system.

  inline void setTheta(HepDouble);
  inline void setPhi(HepDouble);
  inline void setRho(HepDouble);
  // Set spatial vector components in spherical coordinate system.

  inline HepDouble operator () (int) const;
  // Get components by index.

  inline HepLorentzVector & operator = (const HepLorentzVector &);
  // Assignment. 

  inline HepLorentzVector   operator +  (const HepLorentzVector &) const;
  inline HepLorentzVector & operator += (const HepLorentzVector &);
  // Additions.

  inline HepLorentzVector   operator -  (const HepLorentzVector &) const;
  inline HepLorentzVector & operator -= (const HepLorentzVector &);
  // Subtractions.

  inline HepLorentzVector operator - () const;
  // Unary minus.

  inline HepBoolean operator == (const HepLorentzVector &) const;
  inline HepBoolean operator != (const HepLorentzVector &) const;
  // Comparisons.

  inline HepDouble perp2() const;
  // Transverse component of the spatial vector squared.

  inline HepDouble perp() const;
  // Transverse component of the spatial vector (R in cylindrical system).

  inline HepDouble perp2(const Hep3Vector &) const;
  // Transverse component of the spatial vector w.r.t. given axis squared.

  inline HepDouble perp(const Hep3Vector &) const;
  // Transverse component of the spatial vector w.r.t. given axis.

  inline HepDouble angle(const Hep3Vector &) const;
  // Angle wrt. another vector.

  inline HepDouble mag2() const;
  inline HepDouble m2() const;
  // Invariant mass squared.

  inline HepDouble mag() const;
  inline HepDouble m() const;
  // Invariant mass. If mag2() is negative then -sqrt(-mag2()) is returned.

  inline HepDouble dot(const HepLorentzVector &) const;
  inline HepDouble operator * (const HepLorentzVector &) const;
  // Scalar product.

  inline HepDouble plus() const;
  inline HepDouble minus() const;
  // Returns the positive/negative light-cone component t +/- z.


  inline Hep3Vector boostVector() const ;
  // Returns the spatial components divided by the time component.

  void boost(HepDouble, HepDouble, HepDouble);
  inline void boost(const Hep3Vector &);
  // Lorentz boost.

  inline void rotateX(HepDouble);
  // Rotate the spatial component around the x-axis.

  inline void rotateY(HepDouble);
  // Rotate the spatial component around the y-axis.

  inline void rotateZ(HepDouble);
  // Rotate the spatial component around the z-axis.

  inline void rotateUz(Hep3Vector &);
  // Rotates the reference frame from Uz to newUz (unit vector).

  inline void rotate(HepDouble, const Hep3Vector &);
  // Rotate the spatial component around specified axis.

  inline HepLorentzVector & operator *= (const HepRotation &);
  inline HepLorentzVector & transform(const HepRotation &);
  // Transformation with HepRotation.

  HepLorentzVector & operator *= (const HepLorentzRotation &);
  HepLorentzVector & transform(const HepLorentzRotation &);
  // Transformation with HepLorenzRotation.

private:

  Hep3Vector pp;
  HepDouble  ee;

};

ostream & operator << (ostream &, const HepLorentzVector &);
// Output to a stream.

#ifdef HEP_NO_INLINE_IN_DECLARATION
#undef inline
#endif

#ifdef HEP_SHORT_NAMES
typedef HepLorentzVector VectorL;
typedef HepLorentzVector Vector4;
typedef HepLorentzVector DVectorL;
typedef HepLorentzVector DVector4;
typedef HepLorentzVector FVectorL;
typedef HepLorentzVector FVector4;
#endif

typedef HepLorentzVector HepLorentzVectorD;
typedef HepLorentzVector HepLorentzVectorF;

#ifndef HEP_DEBUG_INLINE
#include "CLHEP/Vector/LorentzVector.icc"
#endif

#endif /* HEP_LORENTZVECTOR_H */
