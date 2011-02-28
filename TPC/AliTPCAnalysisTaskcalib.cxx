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
// ANALYSIS task to perrorm TPC calibration                                  //

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPCAnalysisTaskcalib.h"
#include "TChain.h"
#include "AliTPCcalibBase.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTimeStamp.h"

ClassImp(AliTPCAnalysisTaskcalib)


AliTPCAnalysisTaskcalib::AliTPCAnalysisTaskcalib()
  :AliAnalysisTask(),
   fCalibJobs(0),
   fESD(0),
   fESDfriend(0),
   fDebugOutputPath("")
{
  //
  // default constructor
  // 
  
}


AliTPCAnalysisTaskcalib::AliTPCAnalysisTaskcalib(const char *name) 
  :AliAnalysisTask(name,""),
   fCalibJobs(0),
   fESD(0),
   fESDfriend(0),
   fDebugOutputPath("")
{
  //
  // Constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TObjArray::Class());
  fCalibJobs = new TObjArray(0);
  fCalibJobs->SetOwner(kTRUE);
}

AliTPCAnalysisTaskcalib::~AliTPCAnalysisTaskcalib() {
  //
  // destructor
  //
  printf("AliTPCAnalysisTaskcalib::~AliTPCAnalysisTaskcalib");
  fCalibJobs->Delete();
}

void AliTPCAnalysisTaskcalib::Exec(Option_t *) {
  //
  // Exec function
  // Loop over tracks and call  Process function
  if (!fESD) {
    //Printf("ERROR: fESD not available");
    return;
  }
  fESDfriend=static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) {
    //Printf("ERROR: fESDfriend not available");
    return;
  }
  if (fESDfriend->TestSkipBit()) return;
  Int_t n=fESD->GetNumberOfTracks();
  Process(fESD);
  Int_t run = fESD->GetRunNumber();
  for (Int_t i=0;i<n;++i) {
    AliESDfriendTrack *friendTrack=fESDfriend->GetTrack(i);
    AliESDtrack *track=fESD->GetTrack(i);
    TObject *calibObject=0;
    AliTPCseed *seed=0;
    if (!friendTrack) continue;
    for (Int_t j=0;(calibObject=friendTrack->GetCalibObject(j));++j)
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject)))
	break;
    if (track) Process(track, run);
    if (seed)
      Process(seed);
  }
  PostData(0,fCalibJobs);
}

void AliTPCAnalysisTaskcalib::ConnectInputData(Option_t *) {
  //
  //
  //
  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    //Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      //Printf("ERROR: Could not get ESDInputHandler");
    } 
    else {
      fESD = esdH->GetEvent();
      //Printf("*** CONNECTED NEW EVENT ****");
    }
  }
}

void AliTPCAnalysisTaskcalib::CreateOutputObjects() {
  //
  //
  //
  //OpenFile(0, "RECREATE");

  PostData(0,fCalibJobs);
}

void AliTPCAnalysisTaskcalib::Terminate(Option_t */*option*/) {
  //
  // Terminate
  //
  AliTPCcalibBase *job=0;
  Int_t njobs = fCalibJobs->GetEntriesFast();
  for (Int_t i=0;i<njobs;i++){
    job = (AliTPCcalibBase*)fCalibJobs->UncheckedAt(i);
    if (job) job->Terminate();
  }
  
}

void AliTPCAnalysisTaskcalib::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call 
  // on the slaves before sending data
  //
  Terminate("slave");
  if(!fDebugOutputPath.IsNull()) { 
    RegisterDebugOutput();
  }
  
}


void AliTPCAnalysisTaskcalib::Process(AliESDEvent *event) {
  //
  // Process ESD event
  //
  AliTPCcalibBase *job=0;
  Int_t njobs = fCalibJobs->GetEntriesFast();
  for (Int_t i=0;i<njobs;i++){
    job = (AliTPCcalibBase*)fCalibJobs->UncheckedAt(i);
    if (job) {
      job->UpdateEventInfo(event);
      if (job->AcceptTrigger())
	job->Process(event);
    }
  }
}

void AliTPCAnalysisTaskcalib::Process(AliTPCseed *track) {
  //
  // Process TPC track
  //
  AliTPCcalibBase *job=0;
  Int_t njobs = fCalibJobs->GetEntriesFast();
  for (Int_t i=0;i<njobs;i++){
    job = (AliTPCcalibBase*)fCalibJobs->UncheckedAt(i);
    if (job)  
      if (job->AcceptTrigger())
	job->Process(track);
  }
}

void AliTPCAnalysisTaskcalib::Process(AliESDtrack *track, Int_t run) {
  //
  // Process ESD track
  //
  AliTPCcalibBase *job=0;
  Int_t njobs = fCalibJobs->GetEntriesFast();
  for (Int_t i=0;i<njobs;i++){
    job = (AliTPCcalibBase*)fCalibJobs->UncheckedAt(i);
    if (job) 
      if (job->AcceptTrigger())
	job->Process(track,run);
  }
}

Long64_t AliTPCAnalysisTaskcalib::Merge(TCollection *li) {
  TIterator *i=fCalibJobs->MakeIterator();
  AliTPCcalibBase *job;
  Long64_t n=0;
  while ((job=dynamic_cast<AliTPCcalibBase*>(i->Next())))
    n+=job->Merge(li);
  return n;
}

void AliTPCAnalysisTaskcalib::Analyze() {
  //
  // Analyze the content of the task
  //
  AliTPCcalibBase *job=0;
  Int_t njobs = fCalibJobs->GetEntriesFast();
  for (Int_t i=0;i<njobs;i++){
    job = (AliTPCcalibBase*)fCalibJobs->UncheckedAt(i);
    if (job) job->Analyze();
  }
}


void AliTPCAnalysisTaskcalib::RegisterDebugOutput(){
  //
  //
  //
  AliTPCcalibBase *job=0;
  Int_t njobs = fCalibJobs->GetEntriesFast();
  for (Int_t i=0;i<njobs;i++){
    job = (AliTPCcalibBase*)fCalibJobs->UncheckedAt(i);
    if (job) job->RegisterDebugOutput(fDebugOutputPath.Data());
  }
  TString dsName=GetName();
  dsName+=".root";
  TFile fff(dsName.Data(),"recreate");
  fCalibJobs->Write("TPCCalib",TObject::kSingleKey);
  fff.Close();
  //
  // store  - copy debug output to the destination position
  // currently ONLY for local copy
  TString dsName2=fDebugOutputPath.Data();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=gSystem->HostName();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+="/";
  TTimeStamp s;
  dsName2+=Int_t(s.GetNanoSec());
  dsName2+="/";
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=dsName;
  AliInfo(Form("copy %s\t%s\n",dsName.Data(),dsName2.Data()));
  printf("copy %s\t%s\n",dsName.Data(),dsName2.Data());
  TFile::Cp(dsName.Data(),dsName2.Data());

}
