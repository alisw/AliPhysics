#ifndef AliTOFFormatDCS_H
#define AliTOFFormatDCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

// AliTOFFormatDCS class
// describing the format of the calibration data
// coming from DCS

class AliTOFFormatDCS : public TObject 
{
 public:
  
  AliTOFFormatDCS();
  AliTOFFormatDCS(const AliTOFFormatDCS & format);
  AliTOFFormatDCS& operator=(const AliTOFFormatDCS & format);
  virtual ~AliTOFFormatDCS();
  
  Float_t GetFloat(Int_t i) const {return fFloats[i];}
  Float_t GetTimeStampFloat(Int_t i) const {return fTimeStampsFloat[i];}
  Float_t GetDelta(Int_t i) const {return fDeltas[i];}
  Float_t GetTimeStampDelta(Int_t i) const {return fTimeStampsDelta[i];}
  Short_t GetShort() const {return fShort;}
  
  void SetFloat(Int_t i, Float_t valfloat) {fFloats[i]=valfloat;}
  void SetTimeStampFloat(Int_t i, Float_t timestampfloat) {fTimeStampsFloat[i]=timestampfloat;}
  void SetDelta(Int_t i, Float_t valdelta) {fDeltas[i]=valdelta;}
  void SetTimeStampDelta(Int_t i, Float_t timestampdelta) {fTimeStampsDelta[i]=timestampdelta;}
  void SetShort(Short_t valshort) {fShort=valshort;}
  
 private:
  
  Float_t fFloats[3];           // Floats for values at determined 
                                // time intervals
  Float_t fTimeStampsFloat[3];  // Time stamps of the Floats 
                                // for values at determined time intervals
  Float_t fDeltas[2];           // Significant Increments
  Float_t fTimeStampsDelta[2];  // Time stamps for the significant increments
  Short_t fShort;               // Short Integer
  
  ClassDef(AliTOFFormatDCS, 1);
};

#endif
