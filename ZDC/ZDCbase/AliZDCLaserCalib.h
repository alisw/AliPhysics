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
  Int_t GetDetector(Int_t i) const {return fDetector[i];}
  Int_t GetSector(Int_t i)   const {return fSector[i];}
  Float_t GetPMValue(Int_t i)  const {return fPMValue[i];}
  Float_t GetPMWidth(Int_t i)  const {return fPMWidth[i];}
  
  void  SetDetector(Int_t i, Int_t ival) {fDetector[i] = ival;}
  void  SetSector(Int_t i, Int_t ival)   {fSector[i] = ival;}
  void  SetfPMValue(Int_t i, Float_t ival) {fPMValue[i] = ival;}
  void  SetfPMWidth(Int_t i, Float_t ival) {fPMWidth[i] = ival;}
  
 protected:
  Int_t fDetector[24];// detector code
  Int_t fSector[24];  // sector in detector (=5 for reference PMs)
  Float_t fPMValue[24]; // ADC spectrum mean value
  Float_t fPMWidth[24]; // ADC spectrum width
  //
  ClassDef(AliZDCLaserCalib,4)    // ZDC LASER calibration data
};

#endif
