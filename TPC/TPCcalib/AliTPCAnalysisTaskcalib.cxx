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
#include "AliVEvent.h"
#include "AliVfriendEvent.h"
#include "AliVTrack.h"
#include "AliVfriendTrack.h"
#include "AliTPCseed.h"
#include "AliVEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTimeStamp.h"

ClassImp(AliTPCAnalysisTaskcalib)


AliTPCAnalysisTaskcalib::AliTPCAnalysisTaskcalib()
  :AliAnalysisTask(),
   fCalibJobs(0),
   fV(0),
   fVfriend(0),
   fDebugOutputPath("")
{
  //
  // default constructor
  // 
  
}


AliTPCAnalysisTaskcalib::AliTPCAnalysisTaskcalib(const char *name) 
  :AliAnalysisTask(name,""),
   fCalibJobs(0),
   fV(0),
   fVfriend(0),
   fDebugOutputPath("")
{
  //
  // Constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, AliTPCcalibBase::Class());
  DefineOutput(1, AliTPCcalibBase::Class());
  DefineOutput(2, AliTPCcalibBase::Class());
  DefineOutput(3, AliTPCcalibBase::Class());
  DefineOutput(4, AliTPCcalibBase::Class());
  DefineOutput(5, AliTPCcalibBase::Class());
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
    //Printf(" **************** AliTPCAnalysisTaskcalib::Exec() **************** ");
  if (!fV) {
    //Printf("ERROR: fV not available");
    return;
  }
  //fVfriend=fV->FindFriend();
  Int_t n=fV->GetNumberOfTracks();
  Process(fV);
  if (!fVfriend) {
    Printf("ERROR: fVfriend not available");
    return;
  }
  if (fVfriend->TestSkipBit()) {
      //Printf("Skipping Event...");
      return;}
  //else if (!fVfriend->TestSkipBit()){
      //Printf("continue with event...");
  //}
  //
  Int_t run = fV->GetRunNumber();
  for (Int_t i=0;i<n;++i) {
    AliVfriendTrack *friendTrack=const_cast<AliVfriendTrack*>(fVfriend->GetTrack(i));
    AliVTrack *track=fV->GetVTrack(i);
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
}

void AliTPCAnalysisTaskcalib::ConnectInputData(Option_t *) {
  //
  //
  //
  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliVEventHandler *vH = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    TString classInputHandler = vH->ClassName();
      if (!vH) {
      Printf("ERROR: Could not get VEventHandler");
    }
    else {
      fV = vH->GetEvent();
      //if (fV) {Printf("*** CONNECTED NEW EVENT ****");}
      if (classInputHandler.Contains("HLT")) { // we are running in HLT
        fVfriend = vH->GetVfriendEvent();
        //if (fVfriend) Printf("Connected friend Event from V manager!");
      }
      else { /// we are running offline
        if (vH && vH->GetTree()) {
          //Printf("...We got the tree...");
          if (vH->GetTree()->GetBranch("ESDfriend.")){
            //Printf("friend branch found, use AliESDInputHandler");
            fVfriend = ((AliESDInputHandler*)vH)->GetESDfriend();
          }
        }
      }
    }
  }
}

void AliTPCAnalysisTaskcalib::CreateOutputObjects() {
  //
  //
  //
  //OpenFile(0, "RECREATE");

  for (Int_t i=0; i<fCalibJobs->GetEntries(); i++)
  {
    if (fCalibJobs->At(i))
      PostData(i,(AliTPCcalibBase*)fCalibJobs->At(i));
  }
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


void AliTPCAnalysisTaskcalib::Process(AliVEvent *event) {
  //
  // Process V event
  //

    //Printf("AliTPCAnalysisTaskcalib::Process event");
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

void AliTPCAnalysisTaskcalib::Process(AliVTrack *track, Int_t run) {
  //
  // Process V track
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
