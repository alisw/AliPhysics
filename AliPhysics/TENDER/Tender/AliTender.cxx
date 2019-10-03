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

/* $Id$ */

#include <TChain.h>
#include <TFile.h>
 
#include "AliTender.h"
#include "AliTenderSupply.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliLog.h"


ClassImp(AliTender)

//______________________________________________________________________________
AliTender::AliTender():
           AliAnalysisTaskSE(),
           fRun(0),
           fRunChanged(kFALSE),
           fHandleCDB(kFALSE),
           fCDBkey(0),
           fDefaultStorage(),
           fCDB(NULL),
           fESDhandler(NULL),
           fESD(NULL),
           fSupplies(NULL),
           fCDBSettings(NULL)
{
// Dummy constructor
}

//______________________________________________________________________________
AliTender::AliTender(const char* name):
           AliAnalysisTaskSE(name),
           fRun(0),
           fRunChanged(kFALSE),
           fHandleCDB(kFALSE),
           fCDBkey(0),
           fDefaultStorage(),
           fCDB(NULL),
           fESDhandler(NULL),
           fESD(NULL),
           fSupplies(NULL),
           fCDBSettings(NULL)
{
// Default constructor
  DefineOutput(1,  AliESDEvent::Class());
}

//______________________________________________________________________________
AliTender::~AliTender()
{
// Destructor
  if (fSupplies) {
    fSupplies->Delete();
    delete fSupplies;
  }
}

//______________________________________________________________________________
void AliTender::AddSupply(AliTenderSupply *supply)
{
// Addition of supplies.
  if (!fSupplies) fSupplies = new TObjArray();
  if (fSupplies->FindObject(supply)) {
     Error("AddSupply", "Tender supply %s already connected.", supply->GetName());
     return;
  }   
  fSupplies->Add(supply);
  supply->SetTender(this);
}
   
//______________________________________________________________________________
void AliTender::ConnectInputData(Option_t* option)
{
// Connect the input data, create CDB manager.
  if (fDebug > 1) Printf("AliTender::ConnectInputData()\n");
  
  if (!fESDhandler) { 
    AliAnalysisTaskSE::ConnectInputData(option);
    fESDhandler = dynamic_cast<AliESDInputHandler *>(fInputHandler);
  }
  
  if (fESDhandler) {
    fESD = (AliESDEvent*)fESDhandler->GetEvent();
  } else {
     AliFatal("No ESD input event handler connected") ; 
  }
  
  // Obtain run number
  Int_t run = AliAnalysisManager::GetAnalysisManager()->GetRunFromPath();
  if (!run) {
     AliWarning("AliTaskCDBconnect: Could not set run from path");
  } else {
     fRun = run;
     fRunChanged = kTRUE;
     printf("AliTender: #### Setting run to: %d\n", fRun);
  }   

  fCDB = AliCDBManager::Instance();
  // Initialize OCDB (only done when explicitly requested)
  if(fHandleCDB){
    // Create CDB manager
    if (!fDefaultStorage.Length()) AliFatal("Default CDB storage not set.");
    // SetDefault storage. Specific storages must be set by AliTenderSupply::Init()
    fCDB->SetDefaultStorage(fDefaultStorage);
    // Unlock CDB
    fCDBkey = fCDB->SetLock(kFALSE, fCDBkey);
    if(run){ fCDB->SetRun(fRun); }
    // Lock CDB
    fCDBkey = fCDB->SetLock(kTRUE, fCDBkey);
  }
  TIter next(fSupplies);
  AliTenderSupply *supply;
  while ((supply=(AliTenderSupply*)next())) supply->Init();
}

//______________________________________________________________________________
void AliTender::UserCreateOutputObjects()
{
// Nothing for the moment, but we may need ESD event replication here.
  if (fDebug > 1) Printf("AliTender::CreateOutputObjects()\n");
  fESDhandler = dynamic_cast<AliESDInputHandler *>(fInputHandler);
  if (fESDhandler && TObject::TestBit(kCheckEventSelection)) {
     fESDhandler->SetUserCallSelectionMask(kTRUE);
     Info("UserCreateOutputObjects","The TENDER will check the event selection. Make sure you add the tender as FIRST wagon!");
  }   
}

//______________________________________________________________________________
void AliTender::UserExec(Option_t* option)
{
//
// Execute all supplied analysis of one event. Notify run change via RunChanged().
  if (fDebug > 1) {
    Long64_t entry = fESDhandler->GetReadEntry();
    Printf("AliTender::Exec() %s ==> processing event %lld\n", fESDhandler->GetTree()->GetCurrentFile()->GetName(),entry);
  }  
  fESD = (AliESDEvent*)fESDhandler->GetEvent();

// Call the user analysis

  // Intercept when the run number changed
  if (fRun != fESD->GetRunNumber()) {
    fRunChanged = kTRUE;
    fRun = fESD->GetRunNumber();
    fCDB = AliCDBManager::Instance();
    if(fHandleCDB){
      // Unlock CDB
      fCDBkey = fCDB->SetLock(kFALSE, fCDBkey);
      fCDB->SetRun(fRun);
      // Lock CDB
      fCDBkey = fCDB->SetLock(kTRUE, fCDBkey);
    } 
  }
  TIter next(fSupplies);
  AliTenderSupply *supply;
  while ((supply=(AliTenderSupply*)next())) supply->ProcessEvent();
  fRunChanged = kFALSE;

  if (TObject::TestBit(kCheckEventSelection)) fESDhandler->CheckSelectionMask();

  TString opt = option;
  if (!opt.Contains("NoPost")) PostData(1, fESD);
}

//______________________________________________________________________________
void AliTender::SetDefaultCDBStorage(const char *dbString)
{
// Set default CDB storage
   fDefaultStorage = dbString;
}
