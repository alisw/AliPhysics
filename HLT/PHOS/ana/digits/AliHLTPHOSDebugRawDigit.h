//-*- Mode: C++ -*-
// $Id$

//insert copyright

#ifndef ALIHLTPHOSDEBUGRAWDIGIT_H
#define ALIHLTPHOSDEBUGRAWDIGIT_H

#include "TObject.h"

class AliHLTPHOSDebugRawDigit : public TObject
{
  
public: 
  AliHLTPHOSDebugRawDigit();
  virtual ~AliHLTPHOSDebugRawDigit();

  void SetX(Int_t x) { fX = x; }
  void SetZ(Int_t z) { fZ = z; }
  void SetAmplitude(Float_t amp) { fAmplitude = amp; }
  void SetTime(Float_t time) { fTime = time; }
  void SetEnergy(Float_t energy) { fEnergy = energy; }
  void SetGain(Int_t gain) { fGain = gain; }

  void SetRawData(UInt_t*);

  void SetCrazyness(Int_t crazyness) { fCrazyness = crazyness; }
  void SetBaseline(Int_t baseline) { fBaseline = baseline; }
  
  Int_t GetX() { return fX; }
  Int_t GetZ() { return fZ; }
  Float_t GetAmplitude() { return fAmplitude; }
  Float_t GetTime() { return fTime; }
  Float_t GetEnergy() { return fEnergy; }
  Int_t GetGain() { return fGain; }

  UInt_t* GetRawData() { return fData; }
 
  Int_t GetCrazyness() {return fCrazyness; }
  Int_t GetBaseline() { return fBaseline; }
  

private:
  
  Int_t fX;
  Int_t fZ;
  Float_t fAmplitude;
  Float_t fTime;
  Float_t fEnergy;
  Int_t fGain;
  
  UInt_t fData[70];
  Int_t fCrazyness; 
  Int_t fBaseline;

  ClassDef(AliHLTPHOSDebugRawDigit, 1);
  
};

#endif
