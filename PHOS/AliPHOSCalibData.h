#ifndef ALIPHOSCALIBDATA_H
#define ALIPHOSCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  class for PHOS calibration                //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TString.h"
#include "AliPHOSEmcCalibData.h"
#include "AliPHOSCpvCalibData.h"
#include "AliCDBMetaData.h"

class AliPHOSCalibData: public TNamed {

 public:
  AliPHOSCalibData();
  AliPHOSCalibData(Int_t runNumber);
  virtual ~AliPHOSCalibData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 
  
  void CreateNew();
  void RandomEmc();
  void RandomCpv();

  Float_t GetADCchannelEmc(Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const;
  
  void SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value);
  void SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value);

  Float_t GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const;
  
  void SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value);
  void SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value);

  void SetDB(const char* db) {fDB=db;}
  void SetEmcDataPath(const char* emcPath) {fEmcDataPath=emcPath;}
  void SetCpvDataPath(const char* cpvPath) {fCpvDataPath=cpvPath;}

  Bool_t WriteEmc(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md);
  Bool_t WriteCpv(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md);

 private:

  AliPHOSEmcCalibData* fCalibDataEmc; // EMC calibration data
  AliPHOSCpvCalibData* fCalibDataCpv; // CPV calibration data
  
  TString fDB;
  TString fEmcDataPath; // path to EMC calibration data
  TString fCpvDataPath; // path to CPV calibration data

  //
  ClassDef(AliPHOSCalibData,1)    // PHOS Calibration data
};

#endif
