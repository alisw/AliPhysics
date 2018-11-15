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
#include <TRegexp.h>
#include <TGeoGlobalMagField.h>
#include "TGeoManager.h"
#include <TInterpreter.h>
 
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
  fFallBackToRaw(kFALSE),
  fRun(0),
  fLock(0),
  fStorage(),
  fSpecCDBUri(),
  fGRPManager(NULL)
{
  // Dummy constructor
  fSpecCDBUri.SetOwner();
}

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect(const char* name, const char *storage, Int_t run, Bool_t fallback)
  :AliAnalysisTask(name, "ESD analysis tender car"),
   fFallBackToRaw(fallback),
   fRun(run),
   fLock(0),
   fStorage(storage),
   fSpecCDBUri(),
   fGRPManager(NULL)
{
  // Default constructor
  fSpecCDBUri.SetOwner();
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
  
  if (!cdb->IsDefaultStorageSet()) {
    if (gSystem->AccessPathName("OCDB.root",kFileExists)==0) {
      cdb->SetDefaultStorage("local://");
      cdb->SetRaw(kFALSE);
      cdb->SetSnapshotMode("OCDB.root");
    }
    else {
      //  
      // automatic setting of year
      Int_t year = -1;
      if      (fRun<139674) year = 2010;
      else if (fRun<170718) year = 2011;
      else if (fRun<194479) year = 2012;
      else if (fRun<199999) year = 2013;
      else if (fRun<208504) year = 2014;
      else if (fRun<247170) year = 2015;
      else if (fRun<267255) year = 2016;
      else if (fRun<282901) year = 2017;
      else {
	year = 2018;
	// no more updates needed
	// TDatime today;
	// if (today.GetYear()!=year) AliErrorF("Adjust CDB connect, we are now in %d!",today.GetYear());
      }
      //
      Bool_t useCVMFS = kFALSE;
      TString inpStor = fStorage.Strip(TString::kTrailing,'/');
      if (inpStor == "cvmfs:") {
	fStorage = Form("local:///cvmfs/alice.cern.ch/calibration/data/%4d/OCDB",year);
	useCVMFS = kTRUE;
      }
      else if (inpStor.BeginsWith("local:///cvmfs") && inpStor.EndsWith("/OCDB")) {
	TString tmp = inpStor;
	tmp.ReplaceAll("local://","");
	TString strYold = tmp(TRegexp("/[0-9][0-9][0-9][0-9]/OCDB"));
	TString strYnew = Form("/%4d/OCDB",year);
	if (strYold.IsNull()) tmp += strYnew;
	else tmp.ReplaceAll(strYold,strYnew);
	useCVMFS = kTRUE;
	fStorage = Form("local://%s", tmp.Data());
      }
      // check if cvfms is linked    
      if (useCVMFS) {
	TString cvmfspath = fStorage;
	cvmfspath = cvmfspath.ReplaceAll("local://", "");
	if(gSystem->AccessPathName(cvmfspath.Data(),kFileExists)) {
	  if (fFallBackToRaw) {
	    AliErrorF("could not access %s, switching to raw://",fStorage.Data());
	    fStorage = "raw://";
	  }
	  else AliFatalF("could not access %s, fallback to raw:// disabled",fStorage.Data());
	}
      }
      AliInfoF("Setting default storage to %s",fStorage.Data());
      cdb->SetDefaultStorage(fStorage);
    }
    //
    // set specific storages
    for (Int_t i = 0; i < fSpecCDBUri.GetEntriesFast(); i++) {
      TNamed* obj = (TNamed*)fSpecCDBUri[i];
      if (!obj) continue;
      UInt_t vsv = obj->GetUniqueID();
      Int_t ver    = int(vsv>>16)-1;
      Int_t subver = int(vsv&0xffff)-1;
      cdb->SetSpecificStorage(obj->GetName(), obj->GetTitle(), ver, subver);
    }
  }
  if (cdb->GetRun()!=fRun) {    
    fLock = cdb->SetLock(kFALSE,fLock);
    cdb->SetRun(fRun);
    if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {
      // If we are using a specific storage override, this is the place to load it
      gInterpreter->ProcessLine("localOCDBaccessConfig()");
    }
    fLock = cdb->SetLock(kTRUE,fLock);
  }
  //
  if (!fGRPManager) fGRPManager = new AliGRPManager();
  AliInfo("AliCDBconnect: #### Loading GRP to init B-field...");
  if(!fGRPManager->ReadGRPEntry()) AliFatal("Cannot get GRP entry"); 
  if(!fGRPManager->SetMagField())  AliFatal("Problem with magnetic field setup"); 
  //
  // geometry
  if (!gGeoManager) {
    AliInfo("AliCDBconnect: #### Loading geometry...");
    AliGeomManager::LoadGeometry("geometry.root");
    if(!AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD TOF EMCAL PHOS MUON")) AliWarning("Problem with align objects");
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
  if (fRun>0 && !fGRPManager) {
    // in the proof or plugin mode the initialization done in the constructor is not available
    InitGRP();
  }
  else {
    AliInfo("Run number is not available at this stage, InitGRP will be called in the execution loop");
  }
}

//______________________________________________________________________________
void AliTaskCDBconnect::ConnectInputData(Option_t* option)
{
  // Connect the input data, create CDB manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");  
  AliAnalysisTask::ConnectInputData(option);
  Int_t run = AliAnalysisManager::GetAnalysisManager()->GetRunFromPath();
  if (run<=0) {
    AliWarning("AliTaskCDBconnect: Could not set run from path");
    return;
  }
  if (fRun != run) {
    fRun = run;
    InitGRP();
  }
}


//______________________________________________________________________________
void AliTaskCDBconnect::Exec(Option_t* /*option*/)
{
  // Execute all supplied analysis of one event. Notify run change via RunChanged().
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  AliInputEventHandler* inp = (AliInputEventHandler*)mgr->GetInputEventHandler();
  if (!inp) AliFatal("No input event handler connected");
  AliVEvent* ev = inp->GetEvent();
  if (!ev) AliFatal("No event returned");
  int run = ev->GetRunNumber();
  // Intercept when the run number changed
  if (fRun != run && run>0) {
    fRun = run;
    InitGRP();
  }
}

//______________________________________________________________________________
void AliTaskCDBconnect::SetSpecificStorage(const char* calibType, const char* dbString, Int_t version, Int_t subVersion)
{
  // Set a specific storage
  TNamed *nmpath = new TNamed(calibType,dbString);
  if (version<0) version = -1;
  if (subVersion<0) subVersion = -1;
  nmpath->SetUniqueID((UInt_t(version+1)<<16)+UInt_t(subVersion+1));
  fSpecCDBUri.AddLast(nmpath);
}

//______________________________________________________________________________
AliTaskCDBconnect* AliTaskCDBconnect::AddTaskCDBconnect(const char* path/*="raw://"*/, Int_t run) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCDBconnect", "No analysis manager to connect to.");
    return NULL;
  }
  AliTaskCDBconnect *task= new AliTaskCDBconnect("CDBconnect", path, run);
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,  0, cinput1);
  return task;
}

//______________________________________________________________________________
AliTaskCDBconnect* AliTaskCDBconnect::AddTaskCDBconnect() {
  return AliTaskCDBconnect::AddTaskCDBconnect("cvmfs://");
}
