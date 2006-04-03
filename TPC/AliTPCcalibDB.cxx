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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
// Request an instance with AliTPCcalibDB::Instance()                        //
// If a new event is processed set the event number with SetRun              //
// Then request the calibration data                                         // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliLog.h>

#include "AliTPCcalibDB.h"

#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalDet.h"

ClassImp(AliTPCcalibDB)

AliTPCcalibDB* AliTPCcalibDB::fgInstance = 0;
Bool_t AliTPCcalibDB::fgTerminated = kFALSE;


//_ singleton implementation __________________________________________________
AliTPCcalibDB* AliTPCcalibDB::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  
  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTPCcalibDB();
  
  return fgInstance;
}

void AliTPCcalibDB::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0)
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

//_____________________________________________________________________________
AliTPCcalibDB::AliTPCcalibDB()
{
  //
  // constructor
  //  
   // TODO Default runnumber is set to 0, this should be changed later to an invalid value (e.g. -1) to prevent
  // TODO invalid calibration data to be used.
  fRun = 0;
      
  AliCDBManager* manager = AliCDBManager::Instance();
  if (!manager)
  {
    AliFatal("AliTRDcalibDB: CRITICAL: Failed to get instance of AliCDBManager.");
    fLocator = 0;
  }
  else
    fLocator = manager->GetStorage("local://$ALICE_ROOT");    

  //
  //
  //
  fPadGainFactor = 0;
  fPadTime0      = 0;
  fPadPRFWidth   = 0;
  fPadNoise      = 0;
  Update();    // temporary
}

//_____________________________________________________________________________
AliTPCcalibDB::~AliTPCcalibDB() 
{
  //
  // destructor
  //
  if (fPadGainFactor) delete fPadGainFactor;
  if (fPadTime0) delete fPadTime0;
  if (fPadPRFWidth) delete fPadPRFWidth;
  if (fPadNoise) delete fPadNoise;
}


//_____________________________________________________________________________
AliCDBEntry* AliTPCcalibDB::GetCDBEntry(const char* cdbPath)
{
  // 
  // Retrieves an entry with path <cdbPath> from the CDB.
  //
  char chinfo[1000];
    
  if (fRun < 0)
  {
    AliFatal("AliTPCcalibDB: Run number not set! Use AliTPCcalibDB::SetRun.");
    return 0;
  }
  if (!fLocator) 
  { 
    AliError("AliTPCcalibDB: Storage Locator not available."); 
    return 0; 
  } 
  AliCDBEntry* entry = fLocator->Get(cdbPath, fRun); 
  if (!entry) 
  { 
    sprintf(chinfo,"AliTPCcalibDB: Failed to get entry:\t%s ", cdbPath);
    AliError(chinfo); 
    return 0; 
  }
  return entry;
}


//_____________________________________________________________________________
void AliTPCcalibDB::SetRun(Long64_t run)
{
  //
  // Sets current run number. Calibration data is read from the corresponding file. 
  //  
  if (fRun == run)
    return;  
  fRun = run;
  Update();
}
  


void AliTPCcalibDB::Update(){
  //
  AliCDBEntry * entry=0;
  //
  entry          = GetCDBEntry("TPC/Calib/PadGainFactor");
  if (entry){
    if (fPadGainFactor) delete fPadGainFactor;
    entry->SetOwner(kTRUE);
    fPadGainFactor = (AliTPCCalPad*)entry->GetObject();
  }
  //
  entry          = GetCDBEntry("TPC/Calib/PadTime0");
  if (entry){
    if (fPadTime0) delete fPadTime0;
    entry->SetOwner(kTRUE);
    fPadTime0 = (AliTPCCalPad*)entry->GetObject();
  }
  //
  entry          = GetCDBEntry("TPC/Calib/PadPRF");
  if (entry){
    if (fPadPRFWidth) delete fPadPRFWidth;
    entry->SetOwner(kTRUE);
    fPadPRFWidth = (AliTPCCalPad*)entry->GetObject();
  }
  //
  entry          = GetCDBEntry("TPC/Calib/PadNoise");
  if (entry){
    if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fPadNoise = (AliTPCCalPad*)entry->GetObject();
  }
  //

}
