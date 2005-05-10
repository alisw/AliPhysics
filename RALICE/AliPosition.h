#ifndef ALIPOSITION_H
#define ALIPOSITION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "TObject.h"
#include "TString.h"

#include "Ali3Vector.h"
#include "AliTimestamp.h"
 
class AliPosition : public Ali3Vector
{
 public:
  AliPosition();                                         // Default constructor
  virtual ~AliPosition();                                // Destructor
  AliPosition(const AliPosition& p);                     // Copy constructor
  void SetPosition(Double_t* r,TString f);               // Store position r in frame f
  void GetPosition(Double_t* r,TString f) const;         // Provide position r in frame f
  void SetPosition(Float_t*  r,TString f);               // Store position r in frame f
  void GetPosition(Float_t*  r,TString f) const;         // Provide position r in frame f
  AliPosition& GetPosition();                            // Provide position
  void SetPosition(Ali3Vector& r);                       // Store position r
  Double_t GetDistance(AliPosition& p);                  // Provide distance to position p
  Double_t GetDistance(AliPosition* p) { return GetDistance(*p); }
  void SetPositionErrors(Double_t* r,TString f);         // Store position r in frame f
  void GetPositionErrors(Double_t* r,TString f) const;   // Provide position r in frame f
  void SetPositionErrors(Float_t*  r,TString f);         // Store position r in frame f
  void GetPositionErrors(Float_t*  r,TString f) const;   // Provide position r in frame f
  void ResetPosition();                                  // Reset position and errors to 0
  void SetUnitScale(Float_t s);                          // Set unit scale for the position coordinates
  Float_t GetUnitScale() const;                          // Provide unit scale for the position coordinates
  void SetTimestamp(AliTimestamp& t);                    // Set the timestamp for this position
  AliTimestamp* GetTimestamp();                          // Provide the timestamp for this position
  void RemoveTimestamp();                                // Remove the timestamp from this position
  virtual void Data(TString f="car") const;              // Print all position/time info for coord. frame f

 protected:
  Float_t fScale;        // The unit scale used for the position coordinates
  AliTimestamp* fTstamp; // The timestamp for this position

 ClassDef(AliPosition,7) // Handling of positions (with timestamps) in various reference frames.
};
#endif
