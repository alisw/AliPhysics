#ifndef ALI4VECTOR_H
#define ALI4VECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "Ali3Vector.h"
 
class Ali4Vector
{
 public:
  Ali4Vector();                                     // Default constructor for contravariant vector
  virtual ~Ali4Vector();                            // Destructor
  Ali4Vector(const Ali4Vector& v);                  // Copy constructor
  virtual void Load(Ali4Vector& q);                 // Load all attributes of input Ali4Vector
  virtual void SetZero();                           // (Re)set all attributes to zero
  void SetVector(Double_t v0,Ali3Vector& v);        // Store contravariant vector
  void SetVector(Double_t* v,TString f,TString u="rad"); // Store contravariant vector v^i in frame f with ang units u
  void GetVector(Double_t* v,TString f,TString u="rad"); // Provide contravariant vector v^i in frame f in ang units u
  void SetVector(Float_t*  v,TString f,TString u="rad"); // Store contravariant vector v^i in frame f with ang units u
  void GetVector(Float_t*  v,TString f,TString u="rad"); // Provide contravariant vector v^i in frame f in ang units u
  void SetScalar(Double_t v0,Double_t dv0=0);       // Set the scalar part (with error) of v
  void SetScalarError(Double_t dv0);                // Set error on the scalar part of v
  Double_t GetScalar();                             // Provide the scalar part of v
  void Set3Vector(Ali3Vector& v);                   // Set the 3-vector part of v
  void Set3Vector(Double_t* v,TString f,TString u="rad"); // Set the 3-vector part of v in frame f with ang units u
  void Set3Vector(Float_t*  v,TString f,TString u="rad"); // Set the 3-vector part of v in frame f with ang units u
  Ali3Vector Get3Vector() const;                    // Provide the 3-vector part of v
  void SetInvariant(Double_t v2,Double_t dv2=0);    // Set the Lorentz invariant (with error)
  void SetInvariantError(Double_t dv2);             // Set error on the Lorentz invariant
  Double_t GetInvariant();                          // Provide the Lorentz invariant
  void SetErrors(Double_t* v,TString f,TString u="rad"); // Store errors of vector v^i in frame f with ang units u
  void GetErrors(Double_t* v,TString f,TString u="rad"); // Provide errors of vector v^i in frame f in ang units u
  void SetErrors(Float_t*  v,TString f,TString u="rad"); // Store errors of vector v^i in frame f with ang units u
  void GetErrors(Float_t*  v,TString f,TString u="rad"); // Provide errors of vector v^i in frame f in ang units u
  virtual void Data(TString f="car",TString u="rad");    // Print contravariant components in frame f in ang units u
  Double_t Dot(Ali4Vector& q);                      // Provide dot product v^i*q_i
  Double_t GetResultError() const;                  // Provide error on scalar result (e.g. Dot)
  Ali4Vector operator+(Ali4Vector& q);              // Add contravariant vector q
  Ali4Vector operator-(Ali4Vector& q);              // Subtract contravariant vector q
  Ali4Vector operator*(Double_t s);                 // Multiply contravariant vector with scalar s
  Ali4Vector operator/(Double_t s);                 // Divide contravariant vector by scalar s
  Ali4Vector& operator+=(Ali4Vector& q);            // Add contravariant vector q
  Ali4Vector& operator-=(Ali4Vector& q);            // Subtract contravariant vector q
  Ali4Vector& operator*=(Double_t s);               // Multiply with scalar s
  Ali4Vector& operator/=(Double_t s);               // Divide by scalar s
  Int_t GetScalarFlag() const;                      // Provide the fScalar flag value
  Ali3Vector GetVecTrans() const;                   // Provide transverse vector part w.r.t. z-axis
  Ali3Vector GetVecLong() const;                    // Provide longitudinal vector part w.r.t. z-axis
  Double_t GetPseudoRapidity();                     // Provide pseudorapidity of vector part w.r.t z-axis
  Ali3Vector GetBetaVector() const;                 // Provide the beta 3-vector
  Double_t GetBeta();                               // Provide the norm of the beta 3-vector, i.e. v/c
  Double_t GetGamma();                              // Provide the Lorentz gamma factor
  Double_t GetX(Int_t i,TString f,TString u="rad"); // Provide i-th vector component in frame f in units u
  virtual Double_t GetOpeningAngle(Ali4Vector& q,TString u="rad"); // Opening angle between 3-vector parts in units u
  virtual Double_t GetOpeningAngle(Ali3Vector& q,TString u="rad"); // Opening angle with 3-vector q in units u

 protected:
  Double32_t fV2;      // The Lorentz invariant (v^i*v_i)
  Double32_t fV0;      // The scalar part
  Ali3Vector fV;       // The 3-vector part
  Double32_t fDv2;     // The error on the Lorentz invariant
  Double32_t fDv0;     // The error on the scalar part
  Double32_t fDresult; //! The error on the scalar result of an operation (e.g. dotproduct) 
  Int_t fScalar;       // Flag denoting scalar mode
  Double_t GetScaTrans(); // Provide "transverse value" of scalar part w.r.t. z-axis
  Double_t GetScaLong();  // Provide "longitudinal value" of scalar part w.r.t. z-axis

 ClassDef(Ali4Vector,11) // Handling of Lorentz 4-vectors in various reference frames.
};
#endif
