#ifndef ALIZDCCHMAP_H
#define ALIZDCCHMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Class for ZDC calibration containing      //
//    the map ADC ch. <-> physics signal      //
//    the scaler map <-> counted signal       //
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
  Int_t* GetModuleMap() const {return (int*)fModuleMap;}
  Int_t GetModuleMap(Int_t iModType, Int_t iMapEntry) const {return fModuleMap[iModType][iMapEntry];}
  //
  Int_t GetADCModule(Int_t i)  const {return fADCModule[i];}
  Int_t GetADCChannel(Int_t i) const {return fADCChannel[i];}
  Int_t GetDetector(Int_t i)   const {return fDetector[i];}
  Int_t GetSector(Int_t i)     const {return fSector[i];}
  Int_t GetADCSignalCode(Int_t i) const {return fADCSignalCode[i];}
  //
  Int_t GetScChannel(Int_t i)    const {return fScalerChannel[i];}
  Int_t GetScDetector(Int_t i)   const {return fScDetector[i];}
  Int_t GetScSector(Int_t i)     const {return fScSector[i];}
  Int_t GetScSignalCode(Int_t i) const {return fScSignalCode[i];}
  //
  Int_t GetTDCChannel(Int_t i)    const {return fTDCChannel[i];}
  Int_t GetTDCSignalCode(Int_t i) const {return fTDCSignalCode[i];}
  
  void  SetModuleMap(Int_t iEntry, Int_t iGeoAdd, Int_t iModType, Int_t iNCh)
        {fModuleMap[iEntry][0] = iGeoAdd;
	 fModuleMap[iEntry][1] = iModType;
	 fModuleMap[iEntry][2] = iNCh;}
  
  void  SetADCModule(Int_t i, Int_t mod)  {fADCModule[i] = mod;}
  void  SetADCChannel(Int_t i, Int_t ich) {fADCChannel[i] = ich;}
  void  SetDetector(Int_t i, Int_t ival)  {fDetector[i] = ival;}
  void  SetSector(Int_t i, Int_t ival)    {fSector[i] = ival;}
  void  SetADCSignalCode(Int_t i, Int_t ival) {fADCSignalCode[i] = ival;}
  //
  void  SetScChannel(Int_t i, Int_t ich)     {fScalerChannel[i] = ich;}
  void  SetScDetector(Int_t i, Int_t ival)   {fScDetector[i] = ival;}
  void  SetScSector(Int_t i, Int_t ival)     {fScSector[i] = ival;}
  void  SetScSignalCode(Int_t i, Int_t ival) {fScSignalCode[i] = ival;}
  //
  void  SetTDCChannel(Int_t i, Int_t ich)     {fTDCChannel[i] = ich;}
  void  SetTDCSignalCode(Int_t i, Int_t ival) {fTDCSignalCode[i] = ival;}
  
 protected:
  Int_t  fModuleMap[10][3]; // 10 module maps: GEO, mod. type, no. ch.
  // ************ ADC ************
  // 22 signal ch. + 2 reference ch.
  // in-time + out-of-time signals
  // -> 48 channels to be mapped
  Int_t  fADCModule[48];     // ADC module
  Int_t  fADCChannel[48];    // ADC channel
  Int_t  fDetector[48];      // detector
  Int_t  fSector[48];        // sector
  Int_t  fADCSignalCode[48]; // ADC signal code
  //
  // ************ VME scaler ************
  Int_t  fScalerChannel[32]; // Scaler channel
  Int_t  fScDetector[32];    // detector
  Int_t  fScSector[32];	     // sector
  Int_t  fScSignalCode[32];  // scaler signal code
  //
  // ************ ZDC TDC ************
  Int_t  fTDCChannel[32];    // TDC channel
  Int_t  fTDCSignalCode[32]; // TDC signal code
  
  ClassDef(AliZDCChMap,4)    // ZDC pedestal calibration data
};

#endif
