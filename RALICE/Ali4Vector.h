#ifndef ALI4VECTOR_H
#define ALI4VECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class Ali4Vector
// Handling of Lorentz 4-vectors in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have Lorentz 4-vector characteristics.
//
// All 4-vectors are treated in the contravariant form and the convention
// for the metric and the 4-vector components is according to the one
// used in the book "Classical Electrodynamics" by J.D. Jackson.
//
// The dotproduct is defined such that p.Dot(p) yields the Lorentz invariant
// scalar of the 4-vector p (i.e. m**2 in case p is a 4-momentum).   
//
// Note :
// ------
// Vectors (v) and reference frames (f) are specified via
// SetVector(Float_t* v,TString f) under the following conventions :
//
// f="car" ==> 3-vector part of v in Cartesian coordinates   (x,y,z)
// f="sph" ==> 3-vector part of v in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> 3-vector part of v in Cylindrical coordinates (rho,phi,z)
//
// All angles are in radians.
//
// Example :
// ---------
//
// Ali4Vector a;
//
// Float_t v[4]={25,-1,3,7};
// a.SetVector(v,"car");
//
// Float_t vec[4];
// a.GetVector(vec,"sph");
//
// Ali4Vector b;
// Float_t v2[4]={33,6,-18,2};
// b.SetVector(v2,"car");
//
// Float_t dotpro=a.Dot(b);
//
// Float_t x0=16;
// Ali3Vector x;
// Float_t vec2[3]={1,2,3};
// x.SetVector(vec2,"car");
//
// Ali4Vector c;
// c.SetVector(x0,x);
// c.GetVector(vec,"car");
// c.Info("cyl");
// c=a+b;
// c=a-b;
// c=a*5;
//
//--- NvE 01-apr-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
 
#include "Ali3Vector.h"
 
class Ali4Vector
{
 public:
  Ali4Vector();                                     // Default constructor for contravariant vector
  virtual ~Ali4Vector();                            // Destructor
  virtual void SetVector(Double_t v0,Ali3Vector v); // Store contravariant vector
  virtual void SetVector(Double_t* v,TString f);    // Store contravariant vector v^i in frame f
  virtual void GetVector(Double_t* v,TString f);    // Provide contravariant vector v^i in frame f
  virtual void SetVector(Float_t*  v,TString f);    // Store contravariant vector v^i in frame f
  virtual void GetVector(Float_t*  v,TString f);    // Provide contravariant vector v^i in frame f
  Double_t GetScalar();                             // Provide the scalar part of v
  Ali3Vector Get3Vector();                          // Provide the 3-vector part of v
  virtual void Info(TString f="car");               // Print contravariant components in frame f
  Double_t Dot(Ali4Vector& q);                      // Provide dot product v^i*q_i
  Ali4Vector operator+(Ali4Vector& q);              // Add contravariant vector q
  Ali4Vector operator-(Ali4Vector& q);              // Subtract contravariant vector q
  Ali4Vector operator*(Double_t s);                 // Multiply contravariant vector with scalar s
  Ali4Vector operator/(Double_t s);                 // Divide contravariant vector by scalar s
  Ali4Vector& operator+=(Ali4Vector& q);            // Add contravariant vector q
  Ali4Vector& operator-=(Ali4Vector& q);            // Subtract contravariant vector q
  Ali4Vector& operator*=(Double_t s);               // Multiply with scalar s
  Ali4Vector& operator/=(Double_t s);               // Divide by scalar s

 protected:
  Double_t fV0;  // The scalar part
  Ali3Vector fV; // The 3-vector part

 ClassDef(Ali4Vector,1) // Class definition to enable ROOT I/O
};
#endif
