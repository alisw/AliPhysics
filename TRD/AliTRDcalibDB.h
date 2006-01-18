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

//includes neccessary here for compiliation of dynamic_cast
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"

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
  
  AliTRDCalPIDLQ* GetPIDLQObject();
  
  //Related functions, these depend on calibration data
  static Float_t GetOmegaTau(Float_t vdrift);
  Int_t PadResponse(Double_t signal, Double_t dist, Int_t plane, Double_t *pad) const;
  
protected:
  enum { kCDBCacheSize = 8 };   // Number of cached objects
  enum { kIDVdrift = 0, kIDT0 = 1, kIDGainFactor = 2, kIDPRFWidth = 3, kIDGlobals = 4, 
         kIDChamber = 5, kIDStack = 6, kIDPIDLQ = 7 };    // IDs of cached objects
  
  TObject* GetCachedCDBObject(Int_t id)
  {
    //
    // Retrieves a cdb object with the given id. The objects are cached as long as the run number is not changed.
    //
    // Put together the available objects here by using the lines
    //   a) For usual calibration objects:
    //      ase kID<Name> : return CacheCDBEntry(kID<Name>, "TRD/Calib/<Path>"); break;
    //      See function CacheCDBEntry for details.
    //   and
    //   b) For calibration data which depends on two objects: One containing a value per detector and one the local fluctuations per pad:
    //      case kID<Name> : return CacheMergeCDBEntry(kID<Name>, "TRD/Calib/<padPath>", "TRD/Calib/<chamberPath>"); break;
    //      See function CacheMergeCDBEntry for details.
    //
    
    switch (id)
    {
      // parameters defined per pad and chamber
      case kIDVdrift : return CacheMergeCDBEntry(kIDVdrift, "TRD/Calib/LocalVdrift", "TRD/Calib/ChamberVdrift"); break;
      case kIDT0 : return CacheMergeCDBEntry(kIDT0, "TRD/Calib/LocalT0", "TRD/Calib/ChamberT0"); break;
      
      // parameters defined per pad
      case kIDGainFactor : return CacheCDBEntry(kIDGainFactor, "TRD/Calib/GainFactor"); break;
      case kIDPRFWidth : return CacheCDBEntry(kIDPRFWidth, "TRD/Calib/PRFWidth"); break;
    
      // global parameters
      case kIDGlobals : return CacheCDBEntry(kIDGlobals, "TRD/Calib/Globals"); break;
      case kIDChamber : return CacheCDBEntry(kIDChamber, "TRD/Calib/Chamber"); break;
      case kIDStack : return CacheCDBEntry(kIDStack, "TRD/Calib/Stack"); break;
      case kIDPIDLQ : return CacheCDBEntry(kIDPIDLQ, "TRD/Calib/PIDLQ"); break;
    }
    return 0;
  }
  
  void Invalidate();
  void SamplePRF();
  
  inline AliCDBEntry* GetCDBEntry(const char* cdbPath);
  inline TObject* CacheCDBEntry(Int_t id, const char* cdbPath);
  inline TObject* CacheMergeCDBEntry(Int_t id, const char* cdbPadPath, const char* cdbChamberPath);
  
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

AliCDBEntry* AliTRDcalibDB::GetCDBEntry(const char* cdbPath)
{
  // 
  // Retrieves an entry with path <cdbPath> from the CDB.
  //
    
  if (fRun < 0)
  {
    AliFatal("AliTRDcalibDB: Run number not set! Use AliTRDcalibDB::SetRun.");
    //std::cerr << "AliTRDcalibDB: Run number not set! Use AliTRDcalibDB::SetRun." << std::endl;
    return 0;
  }
  if (!fLocator) 
  { 
    std::cerr << "AliTRDcalibDB: Storage Locator not available." << std::endl; 
    return 0; 
  } 
  AliCDBEntry* entry = fLocator->Get(cdbPath, fRun); 
  if (!entry) 
  { 
    std::cerr << "AliTRDcalibDB: Failed to get entry: " << cdbPath << std::endl; 
    return 0; 
  }
  
  std::cout << "AliTRDcalibDB: Retrieved object: " << cdbPath << std::endl;
  return entry;
}

TObject* AliTRDcalibDB::CacheCDBEntry(Int_t id, const char* cdbPath)
{
  //
  // Caches the entry <id> with cdb path <cdbPath>
  //
  
  if (!fCDBCache[id])
  {
    fCDBEntries[id] = GetCDBEntry(cdbPath);
    if (fCDBEntries[id])
      fCDBCache[id] = fCDBEntries[id]->GetObject();
  }
  return fCDBCache[id];
}

TObject* AliTRDcalibDB::CacheMergeCDBEntry(Int_t id, const char* cdbPadPath, const char* cdbChamberPath)
{
  //
  // Retrieves and caches an object (id <id>) from the CDB. This function is specialized for parameters which are stored
  // as local variation at pad level of a global variable defined per detector chamber. It uses the classes AliTRDCalPad and AliTRDCalDet.
  // Before storing the object it retrieves the local variations (cdbPadPath) and the global variable (cdbChamberPath) and merges them using
  // the AliTRDCalPad::ScaleROCs.
  //
    
  if (!fCDBCache[id]) 
  {
    AliTRDCalPad* padObject = 0;
    AliTRDCalDet* detObject = 0;
   
    fCDBEntries[id] = GetCDBEntry(cdbPadPath);
    if (fCDBEntries[id])
      padObject = dynamic_cast<AliTRDCalPad*>(fCDBEntries[id]->GetObject());
   
    AliCDBEntry* mergeEntry = GetCDBEntry(cdbChamberPath);
    if (mergeEntry)
      detObject = dynamic_cast<AliTRDCalDet*>(mergeEntry->GetObject());
    
    if (!padObject || !detObject) 
    {
      if (fCDBEntries[id]) {
        if (fCDBEntries[id]->IsOwner() == kFALSE && padObject)
          delete padObject;
        delete fCDBEntries[id];
        fCDBEntries[id] = 0;
      }
      if (mergeEntry) 
      {
        if (mergeEntry->IsOwner() == kFALSE && detObject)
          delete detObject;
        delete mergeEntry;
      }
      return 0;
    }
    
    padObject->ScaleROCs(detObject);
    if (mergeEntry->IsOwner() == kFALSE)
      delete detObject;
    delete mergeEntry;
    
    fCDBCache[id] = padObject;
  }
  
  return fCDBCache[id];
}

#endif
