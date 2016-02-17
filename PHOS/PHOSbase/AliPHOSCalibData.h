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

class AliPHOSEmcCalibData;
class AliPHOSCpvCalibData;
class AliPHOSEmcBadChannelsMap;
class AliPHOSCpvBadChannelsMap;
class AliCDBMetaData;

class AliPHOSCalibData: public TNamed {

 public:
  AliPHOSCalibData();
  AliPHOSCalibData(Int_t runNumber);
  AliPHOSCalibData(AliPHOSCalibData & phosCDB);
  virtual ~AliPHOSCalibData();

  AliPHOSCalibData & operator = (const AliPHOSCalibData & rhs);

  void Reset();
  virtual void Print(Option_t *option = "") const; 

  AliPHOSEmcCalibData *GetCalibDataEmc() const {return fCalibDataEmc;}
  AliPHOSCpvCalibData *GetCalibDataCpv() const {return fCalibDataCpv;}
  
  void CreateNew();
  void RandomEmc(Float_t ccMin=0.5   , Float_t ccMax=1.5);
  void RandomCpv(Float_t ccMin=0.5, Float_t ccMax=2.);

  //----First EMC parameters---------
  Float_t GetADCchannelEmc(Int_t module, Int_t column, Int_t row) const;
  void    SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value);

  Float_t GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const;
  void    SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value);

  Float_t GetHighLowRatioEmc(Int_t module, Int_t column, Int_t row) const ;
  void    SetHighLowRatioEmc(Int_t module, Int_t column, Int_t row, Float_t value) ;
  
  Float_t GetTimeShiftEmc(Int_t module, Int_t column, Int_t row) const;
  void    SetTimeShiftEmc(Int_t module, Int_t column, Int_t row, Float_t value) ;

  Float_t GetLGTimeShiftEmc(Int_t module, Int_t column, Int_t row) const;
  void    SetLGTimeShiftEmc(Int_t module, Int_t column, Int_t row, Float_t value) ;

  Int_t  GetAltroOffsetEmc(Int_t module, Int_t column, Int_t row) const;
  void   SetAltroOffsetEmc(Int_t module, Int_t column, Int_t row, Int_t value) ;

  Float_t GetSampleTimeStep() const ;
  void    SetSampleTimeStep(Float_t step) ;

  //----Now CPV parameters-----------
  Float_t GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const;
  void    SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value);

  Float_t GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const;
  void    SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value);

  //----Bad channels map-------------
  Int_t  GetNumOfEmcBadChannels() const;
  Bool_t IsBadChannelEmc(Int_t module, Int_t col, Int_t row) const; 
  void   EmcBadChannelIds(Int_t *badIds=0); 

  Int_t  GetNumOfCpvBadChannels() const;
  Bool_t IsBadChannelCpv(Int_t module, Int_t col, Int_t row) const; 
  void   CpvBadChannelIds(Int_t *badIds=0); 


  void SetEmcDataPath(const char* emcPath) {fEmcDataPath=emcPath;}
  void SetCpvDataPath(const char* cpvPath) {fCpvDataPath=cpvPath;}

  Bool_t WriteEmc(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md);
  Bool_t WriteCpv(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md);
  Bool_t WriteEmcBadChannelsMap(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md);
  Bool_t WriteCpvBadChannelsMap(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md);


  //----Decalibration factors for simulation-------------
  Float_t GetADCchannelEmcDecalib(Int_t module, Int_t column, Int_t row) const;
  void    SetADCchannelEmcDecalib(Int_t module, Int_t column, Int_t row, Float_t value);  
  
 private:

  AliPHOSEmcCalibData* fCalibDataEmc; // EMC calibration data
  AliPHOSCpvCalibData* fCalibDataCpv; // CPV calibration data
  AliPHOSEmcBadChannelsMap* fEmcBadChannelsMap; // EMC bad channels map
  AliPHOSCpvBadChannelsMap* fCpvBadChannelsMap; // CPV bad channels map

  TString fEmcDataPath; // path to EMC calibration data
  TString fCpvDataPath; // path to CPV calibration data
  TString fEmcBadChannelsMapPath; // path to bad channels map
  TString fCpvBadChannelsMapPath; // path to bad channels map

  ClassDef(AliPHOSCalibData,7)    // PHOS Calibration data
};

#endif
