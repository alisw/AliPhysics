#ifndef ALI3VECTOR_H
#define ALI3VECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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
  virtual void SetErrors(Double_t* e,TString f); // Store errors of vector in frame f
  virtual void GetErrors(Double_t* e,TString f); // Provide errors of vector in frame f
  virtual void SetErrors(Float_t*  e,TString f); // Store errors of vector in frame f
  virtual void GetErrors(Float_t*  e,TString f); // Provide errors of vector in frame f
  virtual void Info(TString f="car");            // Print vector components in frame f
  Double_t GetNorm();                            // Provide norm of the vector
  Double_t Dot(Ali3Vector& q);                   // Provide dot product with q
  Double_t GetPseudoRapidity();                  // Provide the pseudorapidity w.r.t z-axis
  Double_t GetResultError();                     // Provide error on scalar result (e.g. norm)
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
  Double_t fV,fTheta,fPhi;    // Vector in spherical coordinates
  Double_t fDx,fDy,fDz;       // Errors on Cartesian coordinates
  Double_t fDresult;          // Error on scalar result (e.g. norm or dotproduct)

 ClassDef(Ali3Vector,1) // Handling of 3-vectors in various reference frames.
};
#endif
