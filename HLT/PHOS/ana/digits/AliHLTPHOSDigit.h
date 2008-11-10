//-*- Mode: C++ -*-
// $Id$

//insert copyright

#ifndef ALIHLTPHOSDIGIT_H
#define ALIHLTPHOSDIGIT_H

#include "TObject.h"
//#include "AliHLTPHOSAltroConfig.h"
#include "AliHLTPHOSBase.h"

//class AliHLTPHOSDigit : public TObject, public AliHLTPHOSAltroConfig
class AliHLTPHOSDigit : public TObject, public AliHLTPHOSBase
{
   
public: 
  AliHLTPHOSDigit();
  virtual ~AliHLTPHOSDigit();

  void SetX(Int_t x) { fX = x; }
  void SetZ(Int_t z) { fZ = z; }
  void SetAmplitude(Float_t amp) { fAmplitude = amp; }
  void SetTime(Float_t time) { fTime = time; }
  void SetEnergy(Float_t energy) { fEnergy = energy; }
  void SetGain(Int_t gain) { fGain = gain; }

  void SetRawData(Int_t*);

  void SetCrazyness(Int_t crazyness) { fCrazyness = crazyness; }
  void SetBaseline(Float_t baseline) { fBaseline = baseline; }
  
  void SetSamples(Int_t samples) { fSamples = samples; }
  void SetPreSamples(Int_t presamples) { fPreSamples = presamples; }
 
  void SetDebugVar(Int_t val) { fDebugVar = val; }
  
  Int_t GetX() { return fX; }
  Int_t GetZ() { return fZ; }
  Float_t GetAmplitude() { return fAmplitude; }
  Float_t GetTime() { return fTime; }
  Float_t GetEnergy() { return fEnergy; }
  Int_t GetGain() { return fGain; }

  Int_t* GetRawData() { return fData; }
 
  Int_t GetCrazyness() {return fCrazyness; }
  Float_t GetBaseline() { return fBaseline; }
  
  Int_t GetSamples() { return fSamples; }
  Int_t GetPreSamples() { return fPreSamples; }
  
  Int_t GetDebugVar() { return fDebugVar; }
  

private:
  
  Int_t fX;
  Int_t fZ;
  Float_t fAmplitude;
  Float_t fTime;
  Float_t fEnergy;
  Int_t fGain;
  Int_t fSamples;
  Int_t fPreSamples;
  Int_t fTotalSamples;
  
  Int_t fDebugVar; //can be anything, not intended for use in analysis
  
  Int_t *fData;   //[fTotalSamples]

  Int_t fCrazyness; 
  Float_t fBaseline;

  ClassDef(AliHLTPHOSDigit, 1);
  
};

#endif
