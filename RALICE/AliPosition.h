#ifndef ALIPOSITION_H
#define ALIPOSITION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "TObject.h"
#include "TString.h"

#include "Ali3Vector.h"
 
class AliPosition : public Ali3Vector
{
 public:
  AliPosition();                                         // Default constructor
  virtual ~AliPosition();                                // Destructor
  AliPosition(const AliPosition& p);                     // Copy constructor
  virtual void SetPosition(Double_t* r,TString f);       // Store position r in frame f
  virtual void GetPosition(Double_t* r,TString f);       // Provide position r in frame f
  virtual void SetPosition(Float_t*  r,TString f);       // Store position r in frame f
  virtual void GetPosition(Float_t*  r,TString f);       // Provide position r in frame f
  AliPosition& GetPosition();                            // Provide position
  virtual void SetPosition(Ali3Vector& r);               // Store position r
  Double_t GetDistance(AliPosition& p);                  // Provide distance to position p
  Double_t GetDistance(AliPosition* p) { return GetDistance(*p); }
  virtual void SetPositionErrors(Double_t* r,TString f); // Store position r in frame f
  virtual void GetPositionErrors(Double_t* r,TString f); // Provide position r in frame f
  virtual void SetPositionErrors(Float_t*  r,TString f); // Store position r in frame f
  virtual void GetPositionErrors(Float_t*  r,TString f); // Provide position r in frame f

 ClassDef(AliPosition,2) // Handling of positions in various reference frames.
};
#endif
