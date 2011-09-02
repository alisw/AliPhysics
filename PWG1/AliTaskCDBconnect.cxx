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

#include "AliTaskCDBconnect.h"

#include <TChain.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include "TGeoManager.h"
 
#include "AliAnalysisManager.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliLog.h"

ClassImp(AliTaskCDBconnect)

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect():
           AliAnalysisTask(),
           fRun(0),
           fRunChanged(kFALSE),
           fESDhandler(NULL),
           fESD(NULL),
           fGRPManager(NULL)
{
// Dummy constructor
}

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect(const char* name):
           AliAnalysisTask(name, "ESD analysis tender car"),
           fRun(0),
           fRunChanged(kFALSE),
           fESDhandler(NULL),
           fESD(NULL),
           fGRPManager(NULL)
{
// Default constructor
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  DefineInput (0, TChain::Class());
//  DefineOutput(0,  AliESDEvent::Class());
}

//______________________________________________________________________________
AliTaskCDBconnect::~AliTaskCDBconnect()
{
// Destructor
  if (fGRPManager) delete fGRPManager;
//  if (gGeoManager) delete gGeoManager;
}  

//______________________________________________________________________________
void AliTaskCDBconnect::LocalInit()
{
// Init CDB locally if run number is defined.
}
  
//______________________________________________________________________________
void AliTaskCDBconnect::ConnectInputData(Option_t* /*option*/)
{
// Connect the input data, create CDB manager.
}

//______________________________________________________________________________
void AliTaskCDBconnect::InitGRP()
{
// Initialize geometry and mag. field
  if (!fGRPManager) {
  // magnetic field
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      printf("AliCDBconnect: #### Loading field map...\n");
      fGRPManager = new AliGRPManager();
      if(!fGRPManager->ReadGRPEntry()) { 
        AliError("Cannot get GRP entry"); 
      }
      if( !fGRPManager->SetMagField() ) { 
        AliError("Problem with magnetic field setup"); 
      }
    }

    // geometry
    if (!gGeoManager) {
      printf("AliCDBconnect: #### Loading geometry...\n");
      AliGeomManager::LoadGeometry();
      if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD") ) {
        AliError("Problem with align objects"); 
      }
    }  
  }  
}

//______________________________________________________________________________
void AliTaskCDBconnect::CreateOutputObjects()
{
// Init CDB locally if run number is defined.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  fESDhandler = dynamic_cast<AliESDInputHandler *>(mgr->GetInputEventHandler());
    
  if (!fESDhandler) {
     AliFatal("No ESD input event handler connected");
     return;
  }   
  // Try to get event number before the first event is read (this has precedence
  // over existing fRun)
  Int_t run = mgr->GetRunFromPath();
  if (!run && !fRun) {
     AliError("AliTaskCDBconnect: Run not set - no CDB connection");
     return;
  }
  // Create CDB manager
  AliCDBManager *cdb = AliCDBManager::Instance();
  // If CDB is already locked, return
  if (cdb->GetLock()) return;
  // SetDefault storage. Specific storages must be set by TaskCDBconnectSupply::Init()
  //  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if (!cdb->GetRaw()) {
     cdb->SetDefaultStorage("raw://");
  }   
  if (run && (run != fRun)) {
     fRunChanged = kTRUE;
     fRun = run;
  } else {
     fRunChanged = kFALSE;
  }
  // Set run
  if (fRunChanged || !fGRPManager) {
     printf("AliCDBconnect: #### Setting run to: %d\n", fRun);
     cdb->SetRun(fRun);
     // Initialize GRP manager only once
     if (fRun) InitGRP();
  }   
}

//______________________________________________________________________________
Bool_t AliTaskCDBconnect::Notify()
{
// Init CDB locally if run number is defined.
  CreateOutputObjects();
  return kTRUE;
}

//______________________________________________________________________________
void AliTaskCDBconnect::Exec(Option_t* /*option*/)
{
//
// Execute all supplied analysis of one event. Notify run change via RunChanged().
  fESD = fESDhandler->GetEvent();
  // Intercept when the run number changed
  if (fRun != fESD->GetRunNumber()) {
    fRunChanged = kTRUE;
    fRun = fESD->GetRunNumber();
    CreateOutputObjects();
  }
}

//______________________________________________________________________________
void AliTaskCDBconnect::Terminate(Option_t *)
{
// Initialize CDB also in Terminate
//   CreateOutputObjects();
}
