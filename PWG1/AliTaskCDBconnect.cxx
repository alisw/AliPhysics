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
  if (!fRun) {
    AliWarning("AliTaskCDBconnect::LocalInit Run number not defined yet");
    return;
  }      
  // Create CDB manager
  AliCDBManager *cdb = AliCDBManager::Instance();
  // SetDefault storage. Specific storages must be set by TaskCDBconnectSupply::Init()
//  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetDefaultStorage("raw://");
  // Set run
  cdb->SetRun(fRun);
  // Initialize GRP manager only once
  InitGRP();
}
  
//______________________________________________________________________________
void AliTaskCDBconnect::ConnectInputData(Option_t* /*option*/)
{
// Connect the input data, create CDB manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  fESDhandler = dynamic_cast<AliESDInputHandler *>(mgr->GetInputEventHandler());
    
  if (fESDhandler) {
     fESD = fESDhandler->GetEvent();
  } else {
     AliFatal("No ESD input event handler connected") ; 
  }
  // Create CDB manager
  AliCDBManager *cdb = AliCDBManager::Instance();
  // SetDefault storage. Specific storages must be set by TaskCDBconnectSupply::Init()
//  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetDefaultStorage("raw://");
  Int_t run = fESD->GetRunNumber();
  if (run != fRun) fRunChanged = kTRUE;
  else fRunChanged = kFALSE;
  // Set run
  if (fRunChanged) cdb->SetRun(fRun);
  // Initialize GRP manager only once
  InitGRP();
}

//______________________________________________________________________________
void AliTaskCDBconnect::InitGRP()
{
// Initialize geometry and mag. field
  if (!fGRPManager) {
  // magnetic field
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      printf("Loading field map...\n");
      fGRPManager = new AliGRPManager();
      if(!fGRPManager->ReadGRPEntry()) { 
        printf("Cannot get GRP entry\n"); 
      }
      if( !fGRPManager->SetMagField() ) { 
        printf("Problem with magnetic field setup\n"); 
      }
    }

    // geometry
    printf("Loading geometry...\n");
    AliGeomManager::LoadGeometry();
    if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC") ) {
      printf("Problem with align objects\n"); 
    }
  }  
}

//______________________________________________________________________________
void AliTaskCDBconnect::CreateOutputObjects()
{
// Nothing for the moment, but we may need ESD event replication here.
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
    AliCDBManager::Instance()->SetRun(fRun);
  }
//  PostData(0, fESD);
}
