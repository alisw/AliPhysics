#ifndef ALI3VECTOR_H
#define ALI3VECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "TObject.h"
#include "TString.h"
#include "TRotMatrix.h"
 
class Ali3Vector
{
 public:
  Ali3Vector();                                  // Default constructor
  virtual ~Ali3Vector();                         // Destructor
  Ali3Vector(const Ali3Vector& v);               // Copy constructor
  virtual void Load(Ali3Vector& q);              // Load all attributes of input Ali3Vector
  virtual void SetZero();                        // (Re)set all attributes to zero.
  void SetVector(Double_t* v,TString f);         // Store vector v in frame f
  void GetVector(Double_t* v,TString f) const;   // Provide vector v in frame f
  void SetVector(Float_t*  v,TString f);         // Store vector v in frame f
  void GetVector(Float_t*  v,TString f) const;   // Provide vector v in frame f
  void SetErrors(Double_t* e,TString f);         // Store errors of vector in frame f
  void GetErrors(Double_t* e,TString f) const;   // Provide errors of vector in frame f
  void SetErrors(Float_t*  e,TString f);         // Store errors of vector in frame f
  void GetErrors(Float_t*  e,TString f) const;   // Provide errors of vector in frame f
  virtual void Data(TString f="car") const;      // Print vector components in frame f
  Double_t GetNorm();                            // Provide norm of the vector
  Double_t Dot(Ali3Vector& q);                   // Provide dot product with q
  Double_t GetPseudoRapidity();                  // Provide the pseudorapidity w.r.t z-axis
  Double_t GetResultError() const;               // Provide error on scalar result (e.g. norm)
  Ali3Vector Cross(Ali3Vector& q) const;         // Provide cross product with q
  Ali3Vector operator+(Ali3Vector& q) const;     // Add vector q
  Ali3Vector operator-(Ali3Vector& q) const;     // Subtract vector q
  Ali3Vector operator*(Double_t s) const;        // Multiply vector with scalar s
  Ali3Vector operator/(Double_t s) const;        // Divide vector by scalar s
  Ali3Vector& operator+=(Ali3Vector& q);         // Add vector q
  Ali3Vector& operator-=(Ali3Vector& q);         // Subtract vector q
  Ali3Vector& operator*=(Double_t s);            // Multiply with scalar s
  Ali3Vector& operator/=(Double_t s);            // Divide by scalar s
  Ali3Vector GetVecTrans() const;                // Provide transverse vector w.r.t. z-axis
  Ali3Vector GetVecLong() const;                 // Provide longitudinal vector w.r.t. z-axis
  Ali3Vector GetPrimed(TRotMatrix* m) const;     // Provide vector components in a rotated frame
  Ali3Vector GetUnprimed(TRotMatrix* m) const;   // Provide original vector components from a rotated one

 protected:
  Double_t fV,fTheta,fPhi;    // Vector in spherical coordinates
  Double_t fDx,fDy,fDz;       // Errors on Cartesian coordinates
  Double_t fDresult;          // Error on scalar result (e.g. norm or dotproduct)

 ClassDef(Ali3Vector,8) // Handling of 3-vectors in various reference frames.
};
#endif
