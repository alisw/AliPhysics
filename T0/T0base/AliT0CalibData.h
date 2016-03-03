#ifndef AliT0CalibData_H
#define AliT0CalibData_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TMap.h"

class AliT0CalibData: public TNamed {

 public:
  AliT0CalibData();
  AliT0CalibData(const char* name);
  AliT0CalibData(const AliT0CalibData &calibda);
  AliT0CalibData& operator= (const AliT0CalibData &calibda);
  virtual ~AliT0CalibData();
 
  void     ReadAsciiLookup(const Char_t *filename);
  Int_t    GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel);
  void     PrintLookup(Option_t* option= "") const;
  TMap    *GetMapLookup(void) {return &fLookup;}
  Int_t    GetNumberOfTRMs() const {return fNumberOfTRMs;}
  void     SetNumberOfTRMs(Int_t ntrms=2) {fNumberOfTRMs = ntrms;}


 protected:

  TMap fLookup;           //lookup table
  Int_t fNumberOfTRMs;    // number of TRMs in setup

  //
  ClassDef(AliT0CalibData,8)    // T0 Sensor Calibration data
};

typedef AliT0CalibData AliSTARTCalibData; // for backward compatibility

#endif

