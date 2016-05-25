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
#include "AliVEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

ClassImp(AliTaskCDBconnect)

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect():
           AliAnalysisTask(),
           fRun(0),
           fGRPManager(NULL)
{
// Dummy constructor
}

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect(const char* name, const char *storage, Int_t run)
          :AliAnalysisTask(name, "ESD analysis tender car"),
           fRun(run),
           fGRPManager(NULL)
{
// Default constructor
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(storage);
  if (gSystem->AccessPathName("OCDB.root",kFileExists)==0) cdb->SetSnapshotMode("OCDB.root"); 
  DefineInput (0, TChain::Class());
  if (run>0) InitGRP();
}

//______________________________________________________________________________
AliTaskCDBconnect::~AliTaskCDBconnect()
{
  // Destructor
  delete fGRPManager;
}  

//______________________________________________________________________________
void AliTaskCDBconnect::InitGRP()
{
  // Initialize geometry and mag. field
  AliCDBManager *cdb = AliCDBManager::Instance();
  if (!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("raw://");
  cdb->SetRun(fRun);
  if (!fGRPManager) fGRPManager = new AliGRPManager();
  AliInfo("AliCDBconnect: #### Loading GRP to init B-field...");
  if(!fGRPManager->ReadGRPEntry()) AliFatal("Cannot get GRP entry"); 
  if(!fGRPManager->SetMagField())  AliFatal("Problem with magnetic field setup"); 
  //
  // geometry
  if (!gGeoManager) {
    AliInfo("AliCDBconnect: #### Loading geometry...");
    AliGeomManager::LoadGeometry("geometry.root");
    if(!AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD MUON")) AliWarning("Problem with align objects");
  }  
}

//______________________________________________________________________________
void AliTaskCDBconnect::CreateOutputObjects()
{
  // Init CDB locally if run number is defined.
  //
  //  try to init before the analysis set
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  if ( fRun>0 && !fGRPManager) {
    // in the proof or plugin mode the initialization done in the constructor is not available
    InitGRP();
  }
  else {
    AliInfo("Run number is not available at this stage, InitGRP will be called in the execution loop");
  }
  //
}

//______________________________________________________________________________
void AliTaskCDBconnect::Exec(Option_t* /*option*/)
{
//
// Execute all supplied analysis of one event. Notify run change via RunChanged().
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  AliInputEventHandler* inp = (AliInputEventHandler*)mgr->GetInputEventHandler();
  if (!inp) AliFatal("No input event handler connected");
  //
  AliVEvent* ev = inp->GetEvent();
  if (!ev) AliFatal("No event returned");
  int run = ev->GetRunNumber();
  // Intercept when the run number changed
  if (fRun != run) {
    fRun = run;
    InitGRP();
  }
}

//______________________________________________________________________________
void AliTaskCDBconnect::SetSpecificStorage(const char* calibType, const char* dbString, Int_t version, Int_t subVersion)
{
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetSpecificStorage(calibType,dbString,version,subVersion);
 }
