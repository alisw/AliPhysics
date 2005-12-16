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

class AliTRDCalChamber;
class AliTRDCalStack;
class AliTRDCalGlobals;

  // defines call to function by providing plane, chamber, sector instead of detector
#define HEADER_PAD(funcname) \
  Float_t funcname(Int_t plane, Int_t chamber, Int_t sector, Int_t col, Int_t row) \
  { return funcname(AliTRDgeometry::GetDetector(plane, chamber, sector), col, row); }

  // defines call to function by providing plane, chamber, sector instead of detector
#define HEADER_CHAMBER(funcname) \
  Bool_t funcname(Int_t plane, Int_t chamber, Int_t sector, Float_t* xyz) \
  { return funcname(AliTRDgeometry::GetDetector(plane, chamber, sector), xyz); }

  // This macro creates a piece of code in the function GetCDBObject which retrieves the object with path <cdbPath> from
  // the CDB when GetCDBObject(cdbEntryId) is called.
  /* EXAMPLE: When called e.g. CACHE_CDB_OBJECT(kIDVdrift, "TRD/Calib/LocalVdrift") it produces the following piece of code:
      case kIDVdrift : { \
        if (!fCDBCache[kIDVdrift]) \
        { \
          fCDBEntries[kIDVdrift] = GetCDBEntry(TRD/Calib/LocalVdrift");  \
          if (fCDBEntries[kIDVdrift]) \
            fCDBCache[kIDVdrift] = fCDBEntries[kIDVdrift]->GetObject(); \
        } \
        return fCDBCache[kIDVdrift]; \
      } \
      break;
   */
#define CACHE_CDB_OBJECT(cdbEntryId, cdbPath) \
      case cdbEntryId : { \
        if (!fCDBCache[cdbEntryId]) \
        { \
          fCDBEntries[cdbEntryId] = GetCDBEntry(cdbPath);  \
          if (fCDBEntries[cdbEntryId]) \
            fCDBCache[cdbEntryId] = fCDBEntries[cdbEntryId]->GetObject(); \
        } \
        return fCDBCache[cdbEntryId]; \
      } \
      break;  

// The following macro retrieves and caches an object from the CDB. The macro is specialized for parameters which are stored
// as local variation at pad level of a global variable defined per detector chamber. It uses the classes AliTRDCalPad and AliTRDCalDet.
// Before storing the object it retrieves the local variations (cdbPath) and the global variable (cdbMergePath) and merges them using
// the AliTRDCalPad::ScaleROCs.
#define CACHE_MERGE_CDB_OBJECT(cdbEntryId, cdbPath, cdbMergePath) \
      case cdbEntryId : { \
        if (!fCDBCache[cdbEntryId]) { \
          AliTRDCalDet* detObject = 0; \
          AliTRDCalPad* padObject = 0; \
          \
          fCDBEntries[cdbEntryId] = GetCDBEntry(cdbPath); \
          if (fCDBEntries[cdbEntryId]) \
            padObject = dynamic_cast<AliTRDCalPad*>(fCDBEntries[cdbEntryId]->GetObject()); \
          \
          AliCDBEntry* mergeEntry = GetCDBEntry(cdbMergePath); \
          if (mergeEntry) \
            detObject = dynamic_cast<AliTRDCalDet*>(mergeEntry->GetObject()); \
          if (!padObject || !detObject) { \
            if (fCDBEntries[cdbEntryId]) { \
              delete fCDBEntries[cdbEntryId]; \
              fCDBEntries[cdbEntryId] = 0; \
            } \
            if (mergeEntry) \
              delete mergeEntry; \
            return 0; \
          } \
          padObject->ScaleROCs(detObject); \
          delete mergeEntry; \
          fCDBCache[cdbEntryId] = padObject; \
        } \
        return fCDBCache[cdbEntryId]; \
      } \
      break;  
      
class AliTRDcalibDB : public TObject
{
public:
  static AliTRDcalibDB* Instance();
  static void Terminate();

  void SetRun(Long64_t run);
  
  Bool_t GetChamberPos(Int_t det, Float_t* xyz);
  HEADER_CHAMBER(GetChamberPos);
  
  Bool_t GetChamberRot(Int_t det, Float_t* xyz);
  HEADER_CHAMBER(GetChamberRot);
  
  Bool_t GetStackPos(Int_t chamber, Int_t sector, Float_t* xyz);
  Bool_t GetStackRot(Int_t chamber, Int_t sector, Float_t* xyz);
   
  Float_t GetVdrift(Int_t det, Int_t col, Int_t row);
  HEADER_PAD(GetVdrift);
    
  Float_t GetT0(Int_t det, Int_t col, Int_t row);
  HEADER_PAD(GetT0);
  
  Float_t GetGainFactor(Int_t det, Int_t col, Int_t row);
  HEADER_PAD(GetGainFactor);

  Float_t GetSamplingFrequency(); 
  Int_t GetNumberOfTimeBins();
  
  //Related functions, these depend on calibration data
  virtual Float_t GetOmegaTau(Float_t vdrift);
  
protected:
  void Invalidate();
  
  static AliTRDcalibDB* fgInstance;     // Instance of this class (singleton implementation)
  static Bool_t fgTerminated;               // Defines if this class has already been terminated and therefore does not return instances in GetInstance anymore

  AliCDBStorage* fLocator;                  // Storage locator retrieved from AliCDBManager
  
  enum { kCDBCacheSize = 7 };   // Number of cached objects
  enum { kIDVdrift = 0, kIDT0 = 1, kIDGainFactor = 2, kIDPRFWidth = 3, kIDGlobals = 4, kIDChamber = 5, kIDStack = 6 };    // IDs of cached objects
  
  AliCDBEntry* fCDBEntries[kCDBCacheSize];    // Cache for CDB entries 
  TObject* fCDBCache[kCDBCacheSize];          // Cache for calibration objects.
      
  AliCDBEntry* GetCDBEntry(const char* cdbPath)
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
  
  TObject* GetCachedCDBObject(Int_t id)
  {
    //
    // Retrieves a cdb object with the given id. The objects are cached as long as the run number is not changed.
    //
    // Put together the available objects by using the macros CACHE_CDB_OBJECT and CACHE_MERGE_CDB_OBJECT.
    // The maximal number of cached objects is given above in kCDBCacheSize, the ids of the objects in the enum below kCDBCacheSize.
    // CACHE_CDB_OBJECT can be used for usual calibration objects.
    // CACHE_MERGE_CDB_OBJECT for calibration data which depends on two objects: One containing a value per detector and one the local
    //   fluctuations per pad.
    //
    
    switch (id)
    {
      // parameters defined per pad
      CACHE_MERGE_CDB_OBJECT(kIDVdrift, "TRD/Calib/LocalVdrift", "TRD/Calib/ChamberVdrift");
      CACHE_MERGE_CDB_OBJECT(kIDT0, "TRD/Calib/LocalT0", "TRD/Calib/ChamberT0");
      CACHE_CDB_OBJECT(kIDGainFactor, "TRD/Calib/GainFactor");
      CACHE_CDB_OBJECT(kIDPRFWidth, "TRD/Calib/PRFWidth");
    
      // global parameters
      CACHE_CDB_OBJECT(kIDGlobals, "TRD/Calib/Globals");
      CACHE_CDB_OBJECT(kIDChamber, "TRD/Calib/Chamber");
      CACHE_CDB_OBJECT(kIDStack, "TRD/Calib/Stack");
    }
    return 0;
  }
  
  Long64_t fRun;
  
private:
  // this is a singleton, constructor is private!  
  AliTRDcalibDB();
  virtual ~AliTRDcalibDB();

  ClassDef(AliTRDcalibDB, 0)
};

#undef HEADER_PAD
#undef HEADER_CHAMBER
#undef CACHE_CDB_OBJECT
#undef CACHE_MERGE_CDB_OBJECT

#endif
