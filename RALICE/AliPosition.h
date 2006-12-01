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
  void SetPosition(Double_t* r,TString f,TString u="rad");       // Store position r in frame f with ang units u
  void GetPosition(Double_t* r,TString f,TString u="rad",Float_t s=-1) const; // Provide position r in frame f in ang units u
  void SetPosition(Float_t*  r,TString f,TString u="rad");       // Store position r in frame f with ang units u
  void GetPosition(Float_t*  r,TString f,TString u="rad",Float_t s=-1) const; // Provide position r in frame f in ang units u
  AliPosition& GetPosition();                            // Provide position
  void SetPosition(Ali3Vector& r);                       // Store position r
  Double_t GetDistance(AliPosition& p,Float_t scale=-1); // Provide distance to position p
  Double_t GetDistance(AliPosition* p,Float_t scale=-1) { return GetDistance(*p,scale); }
  void SetPositionErrors(Double_t* r,TString f,TString u="rad");       // Store position r in frame f with ang units u
  void GetPositionErrors(Double_t* r,TString f,TString u="rad",Float_t s=-1) const; // Provide position r in frame f in ang units u
  void SetPositionErrors(Float_t*  r,TString f,TString u="rad");       // Store position r in frame f with ang units u
  void GetPositionErrors(Float_t*  r,TString f,TString u="rad",Float_t s=-1) const; // Provide position r in frame f in ang units u
  void ResetPosition();                                  // Reset position and errors to 0
  void SetUnitScale(Float_t s);                          // Set metric unit scale for the position coordinates
  Float_t GetUnitScale() const;                          // Provide metric unit scale for the position coordinates
  void SetTimestamp(AliTimestamp& t);                    // Set the timestamp for this position
  AliTimestamp* GetTimestamp();                          // Provide the timestamp for this position
  void RemoveTimestamp();                                // Remove the timestamp from this position
  virtual void Data(TString f="car",TString u="rad") const; // Print position/time info for frame f and ang units u

 protected:
  Float_t fScale;        // The unit scale used for the position coordinates
  AliTimestamp* fTstamp; // The timestamp for this position

 ClassDef(AliPosition,9) // Handling of positions (with timestamps) in various reference frames.
};
#endif
