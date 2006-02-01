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

#include <iostream>
#include "TObject.h"

#include "AliLog.h"
#include "AliTRDgeometry.h"

#include <AliCDBStorage.h>
#include <AliCDBEntry.h>

class AliTRDCalPIDLQ;

class AliTRDcalibDB : public TObject
{
public:
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
  
  static AliTRDcalibDB* Instance();
  static void Terminate();

  void SetRun(Long64_t run);
  
  Bool_t GetChamberPos(Int_t det, Float_t* xyz);
  Bool_t GetChamberPos(Int_t plane, Int_t chamber, Int_t sector, Float_t* xyz) { return GetChamberPos(AliTRDgeometry::GetDetector(plane, chamber, sector), xyz); }  
  
  Bool_t GetChamberRot(Int_t det, Float_t* xyz);
  Bool_t GetChamberRot(Int_t plane, Int_t chamber, Int_t sector, Float_t* xyz) { return GetChamberRot(AliTRDgeometry::GetDetector(plane, chamber, sector), xyz); }  
  
  Bool_t GetStackPos(Int_t chamber, Int_t sector, Float_t* xyz);
  Bool_t GetStackRot(Int_t chamber, Int_t sector, Float_t* xyz);
   
  Float_t GetVdrift(Int_t det, Int_t col, Int_t row);
  Float_t GetVdrift(Int_t plane, Int_t chamber, Int_t sector, Int_t col, Int_t row) { return GetVdrift(AliTRDgeometry::GetDetector(plane, chamber, sector), col, row); }
    
  Float_t GetT0(Int_t det, Int_t col, Int_t row);
  Float_t GetT0(Int_t plane, Int_t chamber, Int_t sector, Int_t col, Int_t row) { return GetT0(AliTRDgeometry::GetDetector(plane, chamber, sector), col, row); }
  
  Float_t GetGainFactor(Int_t det, Int_t col, Int_t row);
  Float_t GetGainFactor(Int_t plane, Int_t chamber, Int_t sector, Int_t col, Int_t row) { return GetGainFactor(AliTRDgeometry::GetDetector(plane, chamber, sector), col, row); }

  Float_t GetPRFWidth(Int_t det, Int_t col, Int_t row);
  Float_t GetPRFWidth(Int_t plane, Int_t chamber, Int_t sector, Int_t col, Int_t row) { return GetPRFWidth(AliTRDgeometry::GetDetector(plane, chamber, sector), col, row); }
  
  Float_t GetSamplingFrequency(); 
  Int_t GetNumberOfTimeBins();
  
  const AliTRDCalPIDLQ* GetPIDLQObject();
  
  //Related functions, these depend on calibration data
  static Float_t GetOmegaTau(Float_t vdrift);
  Int_t PadResponse(Double_t signal, Double_t dist, Int_t plane, Double_t *pad) const;
  
protected:
  // for caching see also implentation of GetCachedCDBObject in the .cxx file
  enum { kCDBCacheSize = 8 };   // Number of cached objects
  enum { kIDVdrift = 0, kIDT0 = 1, kIDGainFactor = 2, kIDPRFWidth = 3, kIDGlobals = 4, 
         kIDChamber = 5, kIDStack = 6, kIDPIDLQ = 7 };    // IDs of cached objects
  
  const TObject* GetCachedCDBObject(Int_t id);
  
  void Invalidate();
  void SamplePRF();
  
  AliCDBEntry* GetCDBEntry(const char* cdbPath);
  const TObject* CacheCDBEntry(Int_t id, const char* cdbPath);
  const TObject* CacheMergeCDBEntry(Int_t id, const char* cdbPadPath, const char* cdbChamberPath);
  
  static AliTRDcalibDB* fgInstance;     // Instance of this class (singleton implementation)
  static Bool_t fgTerminated;               // Defines if this class has already been terminated and therefore does not return instances in GetInstance anymore

  AliCDBStorage* fLocator;                  // Storage locator retrieved from AliCDBManager
  
  AliCDBEntry* fCDBEntries[kCDBCacheSize];    // Cache for CDB entries 
  TObject* fCDBCache[kCDBCacheSize];          // Cache for calibration objects.
      
  Long64_t fRun;
  
  struct 
  {
    Float_t             *fPRFsmp;                             //! Sampled pad response
    Int_t                fPRFbin;                             //  Number of bins for the PRF
    Float_t              fPRFlo;                              //  Lower boundary of the PRF
    Float_t              fPRFhi;                              //  Higher boundary of the PRF
    Float_t              fPRFwid;                             //  Bin width of the sampled PRF
    Int_t                fPRFpad;                             //  Distance to next pad in PRF
  } fPadResponse;
  
private:
  // this is a singleton, constructor is private!  
  AliTRDcalibDB();
  virtual ~AliTRDcalibDB();

  ClassDef(AliTRDcalibDB, 0)
};

#endif
