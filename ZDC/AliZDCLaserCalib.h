#ifndef ALIZDCLASERCALIB_H
#define ALIZDCLASERCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////
//  Class for ZDC signal stability monitor      //
//  takes into account PTMs ageing effects      //
//  used for fine adjustments in reconstruction //
//////////////////////////////////////////////////

#include "TNamed.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCLaserCalib: public TNamed {

 public:
  AliZDCLaserCalib();
  AliZDCLaserCalib(const char* name);
  AliZDCLaserCalib(const AliZDCLaserCalib &calibda);
  AliZDCLaserCalib& operator= (const AliZDCLaserCalib &calibda);
  virtual ~AliZDCLaserCalib();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Float_t GetReferenceValue(Int_t i) const {return fReferenceValue[i];}
  Float_t GetSector(Int_t i) const {return fSector[i];}
  Float_t GetGain(Int_t i)   const {return fGain[i];}
  Float_t GetPMRefValue(Int_t i) const {return fPMRefValue[i];}
  Float_t GetPMRefWidth(Int_t i) const {return fPMRefWidth[i];}
  
  void  SetReferenceValue(Int_t i, Float_t ival) {fReferenceValue[i] = ival;}
  void  SetSector(Int_t i, Float_t ival) {fSector[i] = ival;}
  void  SetGain(Int_t i, Float_t ival)   {fGain[i] = ival;}
  void  SetfPMRefValue(Int_t i, Float_t ival){fPMRefValue[i] = ival;}
  void  SetfPMRefWidth(Int_t i, Float_t ival){fPMRefWidth[i] = ival;}
  
 protected:
  // 2 reference ch. x 2 gain chain
  Float_t fReferenceValue[2]; // Reference (inital) value for PM Ref
  Float_t fSector[2];     // sector fSector=1(side C), 4(sideA)
  Float_t fGain[2];	  // fGain=0 (high gain chain), 1 (low gain chain) 
  Float_t fPMRefValue[2]; // ADC meanvalue
  Float_t fPMRefWidth[2]; // ADC spectrum width
  //
  ClassDef(AliZDCLaserCalib,1)    // ZDC LASER calibration data
};

#endif
