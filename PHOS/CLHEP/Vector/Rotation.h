// -*- C++ -*-
// CLASSDOC OFF
// $Id$
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepRotation class for performing rotations
// on objects of the Hep3Vector class.
//
// .SS See Also
// ThreeVector.h, LorentzVector.h, LorentzRotation.h
//
// .SS Author
// Leif Lonnblad

#ifndef HEP_ROTATION_H
#define HEP_ROTATION_H

#ifdef GNUPRAGMA
#pragma interface
#endif

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"

#ifdef HEP_NO_INLINE_IN_DECLARATION
#define inline
#endif

class HepRotation {

public:

  inline HepRotation();
  // Default constructor. Gives a unit matrix.

  inline HepRotation(const HepRotation &);
  // Copy constructor.

  inline HepDouble xx() const;
  inline HepDouble xy() const;
  inline HepDouble xz() const;
  inline HepDouble yx() const;
  inline HepDouble yy() const;
  inline HepDouble yz() const;
  inline HepDouble zx() const;
  inline HepDouble zy() const;
  inline HepDouble zz() const;
  // Elements of the rotation matrix (Geant4).

  HepDouble operator () (int, int) const;
  // Returns (i,j) element of the rotation matrix.

  inline HepRotation & operator = (const HepRotation &);
  // Assignment.

  inline HepBoolean operator == (const HepRotation &) const;
  inline HepBoolean operator != (const HepRotation &) const;
  // Comparisons (Geant4).

  inline HepBoolean isIdentity() const;
  // Returns true if the identity matrix (Geant4).

  inline Hep3Vector operator * (const Hep3Vector &) const;
  // Multiplication with a Hep3Vector.

  HepRotation operator * (const HepRotation &) const;
  inline HepRotation & operator *= (const HepRotation &);
  inline HepRotation & transform(const HepRotation &);
  // Matrix multiplication.
  // Note a *= b; <=> a = a * b; while a.transform(b); <=> a = b * a;

  inline HepRotation inverse() const;
  // Returns the inverse.

  inline HepRotation & invert();
  // Inverts the Rotation matrix.

  HepRotation & rotateX(HepDouble);
  // Rotation around the x-axis.

  HepRotation & rotateY(HepDouble);
  // Rotation around the y-axis.

  HepRotation & rotateZ(HepDouble);
  // Rotation around the z-axis.

  HepRotation & rotate(HepDouble, const Hep3Vector &);
  inline HepRotation & rotate(HepDouble, const Hep3Vector *);
  // Rotation around a specified vector.

  HepRotation & rotateAxes(const Hep3Vector & newX,
                           const Hep3Vector & newY,
                           const Hep3Vector & newZ);
  // Rotation of local axes (Geant4).

  HepDouble phiX() const;
  HepDouble phiY() const;
  HepDouble phiZ() const;
  HepDouble thetaX() const;
  HepDouble thetaY() const;
  HepDouble thetaZ() const;
  // Return angles (RADS) made by rotated axes against original axes (Geant4).

  void getAngleAxis(HepDouble &, Hep3Vector &) const;
  // Returns the rotation angle and rotation axis (Geant4).

protected:

  inline HepRotation(HepDouble, HepDouble, HepDouble, HepDouble, HepDouble,
		     HepDouble, HepDouble, HepDouble, HepDouble);
  // Protected constructor.

  HepDouble rxx, rxy, rxz, ryx, ryy, ryz, rzx, rzy, rzz;
  // The matrix elements.
};

#ifdef HEP_NO_INLINE_IN_DECLARATION
#undef inline
#endif

#ifdef HEP_SHORT_NAMES
typedef HepRotation Rotation;
#endif

#ifndef HEP_DEBUG_INLINE
#include "CLHEP/Vector/Rotation.icc"
#endif

#endif /* HEP_ROTATION_H */

