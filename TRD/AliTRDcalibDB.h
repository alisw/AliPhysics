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

class AliCDBEntry;

class AliTRDrecoParam;
class AliTRDCalTrkAttach;
class AliTRDCalPID;
class AliTRDCalMonitoring;
class AliTRDCalROC;
class AliTRDCalDet;
class AliTRDCalSingleChamberStatus;
class TString;
class AliTRDcalibDB : public TObject {

 public:

  enum { kNlayer  =   6
       , kNstack  =   5
       , kNsector =  18
       , kNdet    = 540 };
  
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

  Float_t                             GetT0(Int_t det, Int_t col, Int_t row);
  Float_t                             GetT0Average(Int_t det);
        AliTRDCalROC                 *GetT0ROC(Int_t det);
  const AliTRDCalDet                 *GetT0Det();

  Float_t                             GetGainFactor(Int_t det, Int_t col, Int_t row);
  Float_t                             GetGainFactorAverage(Int_t det);
  AliTRDCalROC                       *GetGainFactorROC(Int_t det);
  const AliTRDCalDet                 *GetGainFactorDet();

  AliTRDCalROC                       *GetPRFROC(Int_t det);
  Float_t                             GetPRFWidth(Int_t det, Int_t col, Int_t row);

  Float_t*                            GetSampledPRF() const { return fPRFsmp; };
  Int_t                               GetPRFbin() const     { return fPRFbin; };
  Float_t                             GetPRFlo() const      { return fPRFlo;  };
  Float_t                             GetPRFhi() const      { return fPRFhi;  };

  Int_t                               GetNumberOfTimeBinsDCS();
  void                                GetFilterType(TString &filterType);
  void                                GetGlobalConfiguration(TString &config);
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
  
  Bool_t                              IsChamberInstalled(Int_t det);
  Bool_t                              IsChamberMasked(Int_t det);
  Bool_t                              IsHalfChamberMasked(Int_t det, Int_t side);

  const AliTRDCalMonitoring          *GetMonitoringObject();
  const AliTRDCalPID                 *GetPIDObject(AliTRDpidUtil::ETRDPIDMethod m);
  const AliTRDCalTrkAttach           *GetAttachObject();

  // Related functions, these depend on calibration data
         Int_t                        PadResponse(Double_t signal, Double_t dist
                                                , Int_t layer, Double_t *pad) const;

 protected:

  // For caching see also implentation of GetCachedCDBObject in the .cxx file
  enum { kIDVdriftPad = 0
       , kIDVdriftChamber
       , kIDT0Pad
       , kIDT0Chamber
       , kIDGainFactorPad
       , kIDGainFactorChamber
       , kIDNoiseChamber
       , kIDNoisePad
       , kIDPRFWidth
       , kIDFEE
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
       , kCDBCacheSize };         // IDs of cached objects

  const TObject *GetCachedCDBObject(Int_t id);
  
  void           Invalidate();
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
  
 private:

  AliTRDcalibDB();                                  //  This is a singleton, constructor is private!  
  AliTRDcalibDB(const AliTRDcalibDB &c);   
  AliTRDcalibDB &operator=(const AliTRDcalibDB &c); 
  virtual ~AliTRDcalibDB();

  ClassDef(AliTRDcalibDB, 5)                         //  Provides central access to the CDB

};

#endif

