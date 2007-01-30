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

#include "TObject.h"

class AliCDBEntry;

class AliTRDCalPIDLQ;
class AliTRDCalMonitoring;

class AliTRDcalibDB : public TObject {

 public:

  enum { kNplan =   6
       , kNcham =   5
       , kNsect =  18
       , kNdet  = 540 };
  
  static AliTRDcalibDB *Instance();
  static void           Terminate();

  void     SetRun(Long64_t run);
  Long64_t GetRun() { return fRun; }

  Float_t  GetVdrift(Int_t det, Int_t col, Int_t row);
  Float_t  GetVdriftAverage(Int_t det);

  Float_t  GetT0(Int_t det, Int_t col, Int_t row);
  Float_t  GetT0Average(Int_t det);

  Float_t  GetGainFactor(Int_t det, Int_t col, Int_t row);
  Float_t  GetGainFactorAverage(Int_t det);

  Float_t  GetPRFWidth(Int_t det, Int_t col, Int_t row);

  Int_t    GetNumberOfTimeBins();

  Char_t   GetPadStatus(Int_t det, Int_t col, Int_t row);
  Char_t   GetChamberStatus(Int_t det);

  Bool_t   IsPadMasked(Int_t det, Int_t col, Int_t row);
  Bool_t   IsPadBridgedLeft(Int_t det, Int_t col, Int_t row);
  Bool_t   IsPadBridgedRight(Int_t det, Int_t col, Int_t row);
  
  Bool_t   IsChamberInstalled(Int_t det);
  Bool_t   IsChamberMasked(Int_t det);

  const AliTRDCalMonitoring *GetMonitoringObject();
  const AliTRDCalPIDLQ      *GetPIDLQObject();

  // Related functions, these depend on calibration data
  static Float_t GetOmegaTau(Float_t vdrift, Float_t bz);
         Int_t   PadResponse(Double_t signal, Double_t dist, Int_t plane, Double_t *pad) const;
  
 protected:

  // For caching see also implentation of GetCachedCDBObject in the .cxx file
  enum { kCDBCacheSize = 15 };   // Number of cached objects
  enum { kIDVdriftPad = 0
       , kIDVdriftChamber
       , kIDT0Pad
       , kIDT0Chamber
       , kIDGainFactorPad
       , kIDGainFactorChamber
       , kIDPRFWidth
       , kIDGlobals
       , kIDChamberPos
       , kIDStackPos
       , kIDSuperModulePos
       , kIDPIDLQ
       , kIDMonitoringData
       , kIDChamberStatus
       , kIDPadStatus };         // IDs of cached objects

  const TObject* GetCachedCDBObject(Int_t id);
  
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
  
 private:

  AliTRDcalibDB();                                  //  This is a singleton, constructor is private!  
  AliTRDcalibDB(const AliTRDcalibDB &c);   
  AliTRDcalibDB &operator=(const AliTRDcalibDB &c); 
  virtual ~AliTRDcalibDB();

  ClassDef(AliTRDcalibDB,3)                         //  Provides central access to the CDB

};

#endif
