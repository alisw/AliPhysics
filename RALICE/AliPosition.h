#ifndef ALIPOSITION_H
#define ALIPOSITION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class AliPosition
// Handling of positions in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have a unique position in 3-dimensional space.
//
// Note :
// ------
// Positions (r) and reference frames (f) are specified via
// SetPosition(Float_t* r,TString f) under the following conventions :
//
// f="car" ==> r in Cartesian coordinates   (x,y,z)
// f="sph" ==> r in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> r in Cylindrical coordinates (rho,phi,z)
//
// All angles are in radians.
//
// Example :
// ---------
//
// AliPosition q;
// Float_t pos[3]={-1,25,7};
// q.SetPosition(pos,"car");
// Float_t loc[3];
// q.GetPosition(loc,"sph");
//
//--- NvE 06-feb-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TString.h"

#include "Ali3Vector.h"
 
class AliPosition : public Ali3Vector
{
 public:
  AliPosition();                                   // Default constructor
  virtual ~AliPosition();                          // Destructor
  virtual void SetPosition(Double_t* r,TString f); // Store position r in frame f
  virtual void GetPosition(Double_t* r,TString f); // Provide position r in frame f
  virtual void SetPosition(Float_t*  r,TString f); // Store position r in frame f
  virtual void GetPosition(Float_t*  r,TString f); // Provide position r in frame f
  AliPosition& GetPosition();                      // Provide position
  virtual void SetPosition(Ali3Vector& r);         // Store position r

 ClassDef(AliPosition,1) // Class definition to enable ROOT I/O
};
#endif
