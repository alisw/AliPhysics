#ifndef ALITPCCALIBDB_H
#define ALITPCCALIBDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
class AliTPCCalPad;
class AliCDBEntry;
class AliTPCParam;
//class AliCDBStorage;

class AliTPCcalibDB : public TObject
{
 public: 
  static AliTPCcalibDB* Instance();
  AliTPCcalibDB();
  virtual ~AliTPCcalibDB();
  static void Terminate();
  void   SetRun(Long64_t run);   
  //
  AliTPCCalPad* GetPadGainFactor();
  AliTPCCalPad* GetPadTime0();
  AliTPCCalPad* GetPadPRFWidth();
  AliTPCCalPad* GetPadNoise();
  AliTPCCalPad* GetPedestals();
  AliTPCParam*  GetParameters();
  //
protected:
  void         Update();  //update entries
  AliCDBEntry* GetCDBEntry(const char* cdbPath);   
  Long64_t        fRun;         // current run number            
//  AliCDBStorage* fLocator;      // Storage locator retrieved from AliCDBManager
  //
  // calibration parameters per pad
  //
  AliTPCCalPad* fPadGainFactor;
  AliTPCCalPad* fPadTime0;
  AliTPCCalPad* fPadPRFWidth;
  AliTPCCalPad* fPadNoise;
  AliTPCCalPad* fPedestals;
  //
  //
  AliTPCParam * fParam;
  //
  static AliTPCcalibDB* fgInstance;
  static Bool_t       fgTerminated;
 private:
  AliTPCcalibDB(const AliTPCcalibDB &param); // copy constructor
  AliTPCcalibDB &operator = (const AliTPCcalibDB & param);
  ClassDef(AliTPCcalibDB, 0)
};


#endif
