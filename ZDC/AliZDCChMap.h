#ifndef ALIZDCCHMAP_H
#define ALIZDCCHMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Class for ZDC calibration containing      //
//    the map ADC ch. <-> physics signal      //
//       needed for reconstruction            //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCChMap: public TNamed {

 public:
  AliZDCChMap();
  AliZDCChMap(const char* name);
  AliZDCChMap(const AliZDCChMap &calibda);
  AliZDCChMap& operator= (const AliZDCChMap &calibda);
  virtual ~AliZDCChMap();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Int_t GetADCModule(Int_t i)  const {return fADCModule[i];}
  Int_t GetADCChannel(Int_t i) const {return fADCChannel[i];}
  Int_t GetDetector(Int_t i)   const {return fDetector[i];}
  Int_t GetSector(Int_t i)     const {return fSector[i];}

  void  SetADCModule(Int_t i, Int_t mod)  {fADCModule[i] = mod;}
  void  SetADCChannel(Int_t i, Int_t ich) {fADCChannel[i] = ich;}
  void  SetDetector(Int_t i, Int_t ival)  {fDetector[i] = ival;}
  void  SetSector(Int_t i, Int_t ival)    {fSector[i] = ival;}
  
 protected:
  // 22 signal ch. + 2 reference ch.
  // in-time + out-of-time signals
  // -> 48 channels to be mapped
  Int_t  fADCModule[48];  // ADC module
  Int_t  fADCChannel[48]; // ADC channel
  Int_t  fDetector[48];   // detector
  Int_t  fSector[48];     // sector
  //
  ClassDef(AliZDCChMap,1)    // ZDC pedestal calibration data
};

#endif
