#ifndef ALITRDCALIBDB_H
#define ALITRDCALIBDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ALITRDPIDUTIL_H
#include "AliTRDpidUtil.h"
#endif

#ifndef ALITRDPIDRESPONSE_H
#include "AliTRDPIDResponse.h"
#endif

#include "AliTRDCalTrapConfig.h"
#include "AliTRDtrapConfig.h"

class TString;

class AliCDBEntry;

class AliTRDrecoParam;
class AliTRDCalTrkAttach;
class AliTRDCalPID;
class AliTRDCalMonitoring;
class AliTRDCalROC;
class AliTRDCalDet;
class AliTRDCalSingleChamberStatus;
class AliTRDCalOnlineGainTableROC;

class AliTRDcalibDB : public TObject {

 public:

  enum { kNlayer  =   6
       , kNstack  =   5
       , kNsector =  18
       , kNdet    = 540 };
  
  enum { kFltrSet = 1
       , kReadout
       , kTimebin
       , kTrkMode
       , kTrigSet
       , kAddOpti };

  static AliTRDcalibDB               *Instance();
  static void                         Terminate();
  void                                SetRun(Long64_t run);
  Long64_t                            GetRun() const { return fRun; }

  Float_t                             GetNoise(Int_t det, Int_t col, Int_t row);
        AliTRDCalROC                 *GetNoiseROC(Int_t det);
  const AliTRDCalDet                 *GetNoiseDet();

  Float_t                             GetVdrift(Int_t det, Int_t col, Int_t row);
  Float_t                             GetVdriftAverage(Int_t det);
        AliTRDCalROC                 *GetVdriftROC(Int_t det);
  const AliTRDCalDet                 *GetVdriftDet();
  TObjArray                          *GetPHQ();
  const AliTRDCalDet                 *GetExBDet();

  Float_t                             GetT0(Int_t det, Int_t col, Int_t row);
  Float_t                             GetT0Average(Int_t det);
        AliTRDCalROC                 *GetT0ROC(Int_t det);
  const AliTRDCalDet                 *GetT0Det();

  Float_t                             GetGainFactor(Int_t det, Int_t col, Int_t row);
  Float_t                             GetGainFactorAverage(Int_t det);
  AliTRDCalROC                       *GetGainFactorROC(Int_t det);
  const AliTRDCalDet                 *GetGainFactorDet();

  Float_t                             GetOnlineGainFactor(Int_t det, Int_t col, Int_t row);
  AliTRDCalOnlineGainTableROC        *GetOnlineGainTableROC(Int_t det);

  AliTRDCalROC                       *GetPRFROC(Int_t det);
  Float_t                             GetPRFWidth(Int_t det, Int_t col, Int_t row);

  Float_t*                            GetSampledPRF() const { return fPRFsmp; };
  Int_t                               GetPRFbin() const     { return fPRFbin; };
  Float_t                             GetPRFlo() const      { return fPRFlo;  };
  Float_t                             GetPRFhi() const      { return fPRFhi;  };

  static Int_t                        ExtractTimeBinsFromString(TString tbstr);
  Int_t                               GetNumberOfTimeBinsDCS();
  void                                GetFilterType(TString &filterType);
  void                                GetGlobalConfiguration(TString &config);
  void                                GetGlobalConfigurationByChamber(TString &config,Int_t par, Int_t opt=0);
  void                                GetGlobalConfigurationVersion(TString &version);
  static Int_t                        GetNumberOfParsDCS(TString cname, Char_t delimiter='_');
  static Int_t                        GetNumberOfOptsDCS(TString cname, Int_t cfgType);
  static void                         GetDCSConfigParOption(TString cname, Int_t cfgType, Int_t option, TString &cfgo);

  Int_t                               GetOnlineGainTableID();

  Bool_t                              HasOnlineFilterPedestal();
  Bool_t                              HasOnlineFilterGain();
  Bool_t                              HasOnlineTailCancellation();

  Char_t                              GetPadStatus(Int_t det, Int_t col, Int_t row);
  AliTRDCalSingleChamberStatus       *GetPadStatusROC(Int_t det);
  AliTRDrecoParam*                    GetRecoParam(Int_t *eventtype);
  AliTRDPIDResponse                  *GetPIDResponse(AliTRDPIDResponse::ETRDPIDMethod m);

  Char_t                              GetChamberStatus(Int_t det);

  Bool_t                              IsPadMasked(Int_t det, Int_t col, Int_t row);
  Bool_t                              IsPadBridgedLeft(Int_t det, Int_t col, Int_t row);
  Bool_t                              IsPadBridgedRight(Int_t det, Int_t col, Int_t row);
  Bool_t                              IsPadNotConnected(Int_t det, Int_t col, Int_t row);
  
  Bool_t                              IsChamberGood(Int_t det);
  Bool_t                              IsChamberNoData(Int_t det);
  Bool_t                              IsHalfChamberNoData(Int_t det, Int_t side);
  Bool_t                              IsChamberBadCalibrated(Int_t det);
  Bool_t                              IsChamberNotCalibrated(Int_t det);

  const AliTRDCalMonitoring          *GetMonitoringObject();
  const AliTRDCalPID                 *GetPIDObject(AliTRDpidUtil::ETRDPIDMethod m);
  const AliTRDCalTrkAttach           *GetAttachObject();

  // Related functions, these depend on calibration data
  Int_t                               PadResponse(Double_t signal, Double_t dist
                                                , Int_t layer, Double_t *pad) const;

  AliTRDtrapConfig*                   GetTrapConfig();
  void                                GetTrapConfig(TString &name, TString &version) { name = fTrapConfigName; version = fTrapConfigVersion; }
  void                                SetTrapConfig(const TString name, const TString version) { fTrapConfigName = name; fTrapConfigVersion = version; }
  void                                SetTrapConfig(AliTRDtrapConfig *trapcfg) { fTrapConfig = trapcfg; }
  void                                Invalidate();

 protected:

  AliTRDtrapConfig*                   LoadTrapConfig(const TString &name = "", const TString &version = "");
  Int_t                               GetNumberOfTimeBinsDCSBoard(); // Old method as fallback for patched OCDB 
  virtual                            ~AliTRDcalibDB();
  // For caching see also implentation of GetCachedCDBObject in the .cxx file
  enum { kIDVdriftPad = 0
       , kIDVdriftChamber
       , kIDExBChamber
       , kIDT0Pad
       , kIDT0Chamber
       , kIDGainFactorPad
       , kIDGainFactorChamber
       , kIDOnlineGainFactor
       , kIDNoiseChamber
       , kIDNoisePad
       , kIDPRFWidth
       , kIDFEE
       , kIDTrapConfig
       , kIDChamberPos
       , kIDStackPos
       , kIDSuperModulePos
       , kIDPIDNN
       , kIDPIDLQ
       , kIDPIDLQ1D
       , kIDRecoParam
       , kIDMonitoringData
       , kIDChamberStatus
       , kIDPadStatus
       , kIDDCS
       , kIDAttach
       , kIDPHQ
       , kCDBCacheSize };         // IDs of cached objects

  const TObject *GetCachedCDBObject(Int_t id);
  
  void           SamplePRF();
  
  AliCDBEntry   *GetCDBEntry(const Char_t *cdbPath);
  const TObject *CacheCDBEntry(Int_t id, const Char_t *cdbPath);

  static AliTRDcalibDB *fgInstance;                 //  Instance of this class (singleton implementation)
  static Bool_t         fgTerminated;               //  Defines if this class has already been terminated

  AliCDBEntry          *fCDBEntries[kCDBCacheSize]; //  Cache for CDB entries
  TObject              *fCDBCache[kCDBCacheSize];   //  Cache for calibration objects.

  Long64_t              fRun;                       //  Run Number

  Float_t              *fPRFsmp;                    //! Sampled pad response
  Int_t                 fPRFbin;                    //  Number of bins for the PRF
  Float_t               fPRFlo;                     //  Lower boundary of the PRF
  Float_t               fPRFhi;                     //  Higher boundary of the PRF
  Float_t               fPRFwid;                    //  Bin width of the sampled PRF
  Int_t                 fPRFpad;                    //  Distance to next pad in PRF

  AliTRDPIDResponse    *fPIDResponse;               //  TRD PID Response function

  Int_t                 fOnlineGainTableID;         //  ID for online gain table 

  AliTRDtrapConfig*     fTrapConfig;                //  TRAP configuration
  TString               fTrapConfigName;            //  name of the TRAPconfig
  TString               fTrapConfigVersion;         //  version of the TRAPconfig
  
 private:

  AliTRDcalibDB();                                  //  This is a singleton, constructor is private!  
  AliTRDcalibDB(const AliTRDcalibDB &c);   
  AliTRDcalibDB &operator=(const AliTRDcalibDB &c); 

  ClassDef(AliTRDcalibDB, 8)                        //  Provides central access to the CDB

};

#endif

