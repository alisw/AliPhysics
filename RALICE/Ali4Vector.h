#ifndef ALI4VECTOR_H
#define ALI4VECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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
  virtual void SetScalar(Double_t v0,Double_t dv0=0); // Set the scalar part (with error) of v
  virtual void SetScalarError(Double_t dv0);        // Set error on the scalar part of v
  Double_t GetScalar();                             // Provide the scalar part of v
  virtual void Set3Vector(Ali3Vector v);            // Set the 3-vector part of v
  virtual void Set3Vector(Double_t* v,TString f);   // Set the 3-vector part of v in frame f
  virtual void Set3Vector(Float_t*  v,TString f);   // Set the 3-vector part of v in frame f
  Ali3Vector Get3Vector();                          // Provide the 3-vector part of v
  virtual void SetInvariant(Double_t v2,Double_t dv2=0); // Set the Lorentz invariant (with error)
  virtual void SetInvariantError(Double_t dv2);     // Set error on the Lorentz invariant
  Double_t GetInvariant();                          // Provide the Lorentz invariant
  virtual void SetErrors(Double_t* v,TString f);    // Store errors of vector v^i in frame f
  virtual void GetErrors(Double_t* v,TString f);    // Provide errors of vector v^i in frame f
  virtual void SetErrors(Float_t*  v,TString f);    // Store errors of vector v^i in frame f
  virtual void GetErrors(Float_t*  v,TString f);    // Provide errors of vector v^i in frame f
  virtual void Info(TString f="car");               // Print contravariant components in frame f
  Double_t Dot(Ali4Vector& q);                      // Provide dot product v^i*q_i
  Double_t GetResultError();                        // Provide error on scalar result (e.g. Dot)
  Ali4Vector operator+(Ali4Vector& q);              // Add contravariant vector q
  Ali4Vector operator-(Ali4Vector& q);              // Subtract contravariant vector q
  Ali4Vector operator*(Double_t s);                 // Multiply contravariant vector with scalar s
  Ali4Vector operator/(Double_t s);                 // Divide contravariant vector by scalar s
  Ali4Vector& operator+=(Ali4Vector& q);            // Add contravariant vector q
  Ali4Vector& operator-=(Ali4Vector& q);            // Subtract contravariant vector q
  Ali4Vector& operator*=(Double_t s);               // Multiply with scalar s
  Ali4Vector& operator/=(Double_t s);               // Divide by scalar s
  Int_t GetScalarFlag();                            // Provide the fScalar flag value

 protected:
  Double_t fV2;      // The Lorentz invariant (v^i*v_i)
  Double_t fV0;      // The scalar part
  Ali3Vector fV;     // The 3-vector part
  Double_t fDv2;     // The error on the Lorentz invariant
  Double_t fDv0;     // The error on the scalar part
  Double_t fDresult; // The error on the scalar result of an operation (e.g. dotproduct) 
  Int_t fScalar;     // Flag denoting scalar mode

 ClassDef(Ali4Vector,1) // Handling of Lorentz 4-vectors in various reference frames.
};
#endif
