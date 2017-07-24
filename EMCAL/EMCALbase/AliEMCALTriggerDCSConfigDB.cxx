/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

#include "AliEMCALTriggerDCSConfigDB.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerDCSConfigDB) ;
/// \endcond

AliEMCALTriggerDCSConfigDB* AliEMCALTriggerDCSConfigDB::fgInstance   = 0;
Bool_t                      AliEMCALTriggerDCSConfigDB::fgTerminated = kFALSE;

///
/// Singleton implementation
/// \return an instance of this class, it is created if neccessary
//_____________________________________________________________________________
AliEMCALTriggerDCSConfigDB* AliEMCALTriggerDCSConfigDB::Instance()
{  
  if (fgTerminated != kFALSE) 
  {
    return 0;
  }
  
  if (fgInstance == 0) 
  {
    fgInstance = new AliEMCALTriggerDCSConfigDB();
  }
  
  return fgInstance;
}

///
/// Singleton implementation
/// Deletes the instance of this class and sets the terminated flag,
/// instances cannot be requested anymore
/// This function can be called several times.
//_____________________________________________________________________________
void AliEMCALTriggerDCSConfigDB::Terminate()
{
  fgTerminated = kTRUE;
  
  if (fgInstance != 0) 
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

///
/// Default constructor
///
/// TODO Default runnumber is set to 0, this should be changed later
///      to an invalid value (e.g. -1) to prevent
/// TODO invalid calibration data to be used.
//_____________________________________________________________________________
AliEMCALTriggerDCSConfigDB::AliEMCALTriggerDCSConfigDB() : TObject()
,fRun(-1)
{
  for (Int_t i = 0; i < kCDBCacheSize; ++i) 
  {
    fCDBCache[i]   = 0;
    fCDBEntries[i] = 0;
  }
}

///
/// Copy constructor (not that it make any sense for a singleton...)
//_____________________________________________________________________________
AliEMCALTriggerDCSConfigDB::AliEMCALTriggerDCSConfigDB(const AliEMCALTriggerDCSConfigDB &c) : TObject(c)
,fRun(-1)
{
  for (Int_t i = 0; i < kCDBCacheSize; ++i) 
  {
    fCDBCache[i]   = 0;
    fCDBEntries[i] = 0;
  }
}

//
// Assignment operator (not that it make any sense for a singleton...)
//_____________________________________________________________________________
AliEMCALTriggerDCSConfigDB &AliEMCALTriggerDCSConfigDB::operator=(const AliEMCALTriggerDCSConfigDB &c) 
{
  if (this != &c) 
  {
    AliFatal("No assignment operator defined");
  }
  
  return *this;
}

///
/// Destructor
//_____________________________________________________________________________
AliEMCALTriggerDCSConfigDB::~AliEMCALTriggerDCSConfigDB() 
{
	Invalidate();
}

///
/// Retrieves a cdb object with the given id. The objects are cached as
/// long as the run number is not changed.
//_____________________________________________________________________________
const TObject *AliEMCALTriggerDCSConfigDB::GetCachedCDBObject(Int_t id)
{
  switch (id) 
  {
      // Parameters defined per pad and chamber
    case kIDTriggerConfig : 
      return CacheCDBEntry(kIDTriggerConfig, "EMCAL/Calib/Trigger"); 
      break;
    default:			
      AliError("Object not found!");
      break;
  }
  
  return 0x0;
}

/// 
/// Retrieves an entry with path "cdbPath" from the CDB.
//_____________________________________________________________________________
AliCDBEntry* AliEMCALTriggerDCSConfigDB::GetCDBEntry(const char *cdbPath)
{
  AliCDBEntry *entry = AliCDBManager::Instance()->Get(cdbPath,fRun);
  
  if (!entry) 
  { 
    AliError(Form("Failed to get entry: %s",cdbPath));
    return 0; 
  }
  
  return entry;
}

///
/// Caches the entry "id" with cdb path "cdbPath"
//_____________________________________________________________________________
const TObject *AliEMCALTriggerDCSConfigDB::CacheCDBEntry(Int_t id, const char *cdbPath)
{
  if (!fCDBCache[id]) 
  {
    fCDBEntries[id] = GetCDBEntry(cdbPath);
    
    if (fCDBEntries[id]) fCDBCache[id] = fCDBEntries[id]->GetObject();
  }
  
  return fCDBCache[id];
}

///
/// Sets current run number. Calibration data is read from the corresponding file.
/// When the run number changes the caching is invalidated.
//_____________________________________________________________________________
void AliEMCALTriggerDCSConfigDB::SetRun(Long64_t run)
{
  if (fRun == run) return;
  
  fRun = run;
  
  Invalidate();
}

///
/// Invalidates cache (when run number is changed).
//_____________________________________________________________________________
void AliEMCALTriggerDCSConfigDB::Invalidate()
{
  for (Int_t i = 0; i < kCDBCacheSize; ++i) 
  {
    if (fCDBEntries[i]) 
    {
      if (AliCDBManager::Instance()->GetCacheFlag() == kFALSE) 
      {
        if ((fCDBEntries[i]->IsOwner() == kFALSE) && (fCDBCache[i])) delete fCDBCache[i];
        
        delete fCDBEntries[i];
      }
      
      fCDBEntries[i] = 0;
      fCDBCache[i]   = 0;
    }
  }
}

///
/// Get DCS config
//_____________________________________________________________________________
const AliEMCALTriggerDCSConfig* AliEMCALTriggerDCSConfigDB::GetTriggerDCSConfig()
{
  const AliEMCALTriggerDCSConfig* dcsConf = dynamic_cast<const AliEMCALTriggerDCSConfig*>(GetCachedCDBObject(kIDTriggerConfig));
  
  if (!dcsConf) 
  {
    AliError("Trigger DCS configuration not found!");
    return 0x0;
  }
  else
    return dcsConf;
}
