#ifndef ALIPOSITION_H
#define ALIPOSITION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TString.h"

#include "Ali3Vector.h"
 
class AliPosition : public Ali3Vector
{
 public:
  AliPosition();                                         // Default constructor
  virtual ~AliPosition();                                // Destructor
  virtual void SetPosition(Double_t* r,TString f);       // Store position r in frame f
  virtual void GetPosition(Double_t* r,TString f);       // Provide position r in frame f
  virtual void SetPosition(Float_t*  r,TString f);       // Store position r in frame f
  virtual void GetPosition(Float_t*  r,TString f);       // Provide position r in frame f
  AliPosition& GetPosition();                            // Provide position
  virtual void SetPosition(Ali3Vector& r);               // Store position r

  virtual void SetPositionErrors(Double_t* r,TString f); // Store position r in frame f
  virtual void GetPositionErrors(Double_t* r,TString f); // Provide position r in frame f
  virtual void SetPositionErrors(Float_t*  r,TString f); // Store position r in frame f
  virtual void GetPositionErrors(Float_t*  r,TString f); // Provide position r in frame f

 ClassDef(AliPosition,1) // Handling of positions in various reference frames.
};
#endif
