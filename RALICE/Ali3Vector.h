#ifndef ALI3VECTOR_H
#define ALI3VECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class Ali3Vector
// Handling of 3-vectors in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have 3-dimensional vector characteristics.
//
// Note :
// ------
// Vectors (v) and reference frames (f) are specified via
// SetVector(Float_t* v,TString f) under the following conventions :
//
// f="car" ==> v in Cartesian coordinates   (x,y,z)
// f="sph" ==> v in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> v in Cylindrical coordinates (rho,phi,z)
//
// All angles are in radians.
//
// Example :
// ---------
//
// Ali3Vector a;
// Float_t v[3]={-1,25,7};
// a.SetVector(v,"car");
// a.Info();
//
// Float_t vec[3];
// a.GetVector(vec,"sph");
//
// Ali3Vector b;
// Float_t v2[3]={6,-18,33};
// b.SetVector(v2,"car");
//
// Float_t dotpro=a.Dot(b);
//
// Ali3Vector c=a.Cross(b);
// c.Info("sph");
// c.GetVector(vec,"cyl");
// Float_t norm=c.GetNorm();
// c=a+b;
// c=a-b;
// c=a*5;
//
//--- NvE 30-mar-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TString.h"
 
class Ali3Vector
{
 public:
  Ali3Vector();                                  // Default constructor
  virtual ~Ali3Vector();                         // Destructor
  virtual void SetVector(Double_t* v,TString f); // Store vector v in frame f
  virtual void GetVector(Double_t* v,TString f); // Provide vector v in frame f
  virtual void SetVector(Float_t*  v,TString f); // Store vector v in frame f
  virtual void GetVector(Float_t*  v,TString f); // Provide vector v in frame f
  virtual void Info(TString f="car");            // Print vector components in frame f
  Double_t GetNorm();                            // Provide norm of the vector
  Double_t Dot(Ali3Vector& q);                   // Provide dot product with q
  Ali3Vector Cross(Ali3Vector& q);               // Provide cross product with q
  Ali3Vector operator+(Ali3Vector& q);           // Add vector q
  Ali3Vector operator-(Ali3Vector& q);           // Subtract vector q
  Ali3Vector operator*(Double_t s);              // Multiply vector with scalar s
  Ali3Vector operator/(Double_t s);              // Divide vector by scalar s
  Ali3Vector& operator+=(Ali3Vector& q);         // Add vector q
  Ali3Vector& operator-=(Ali3Vector& q);         // Subtract vector q
  Ali3Vector& operator*=(Double_t s);            // Multiply with scalar s
  Ali3Vector& operator/=(Double_t s);            // Divide by scalar s

 protected:
  Double_t fV,fTheta,fPhi; // Vector in spherical coordinates

 ClassDef(Ali3Vector,1) // Class definition to enable ROOT I/O
};
#endif
