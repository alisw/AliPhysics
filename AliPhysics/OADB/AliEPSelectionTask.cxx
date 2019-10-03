/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*****************************************************
//   Class AliEPSelectionTask
//   Class to determine event plane
//   author: Alberica Toia, Johanna Gramling
//*****************************************************

#include "AliEPSelectionTask.h"

#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TFile.h>
#include <TObjString.h>
#include <TString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <iostream>
#include <TRandom2.h>
#include <TArrayF.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliBackgroundSelection.h"
#include "AliESDUtils.h"
#include "AliOADBContainer.h"
#include "AliAODMCHeader.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliMultSelection.h"
#include "AliEventplane.h"

using std::cout;
using std::endl;
ClassImp(AliEPSelectionTask)

//________________________________________________________________________
AliEPSelectionTask::AliEPSelectionTask():
AliAnalysisTaskSE(),
  fAnalysisInput("ESD"),
  fTrackType("TPC"),
  fPeriod(""),
  fCentrality(-1.),
  fUseMCRP(kFALSE),
  fUsePhiWeight(kFALSE),
  fUsePtWeight(kFALSE),
  fUseRecentering(kFALSE),
  fSaveTrackContribution(kFALSE),
  fUserphidist(kFALSE),
  fUsercuts(kFALSE),
  fRunNumber(-15),
  fAODfilterbit(1),
  fEtaGap(0.),
  fSplitMethod(0),
  fESDtrackCuts(0),
  fEPContainer(0),
  fQxContainer(0),
  fQyContainer(0),
  fSparseDist(0),
  fHruns(0),
  fQVector(0),
  fQContributionX(0),
  fQContributionY(0),
  fEventplaneQ(0),
  fQsub1(0),
  fQsub2(0),
  fQsubRes(0),
  fOutputList(0),
  fHOutEventplaneQ(0),
  fHOutPhi(0),
  fHOutPhiCorr(0),
  fHOutsub1sub2(0),
  fHOutNTEPRes(0),
  fHOutPTPsi(0),
  fHOutDiff(0),
  fHOutleadPTPsi(0)
{
  // Default constructor
  AliInfo("Event Plane Selection enabled.");
  for(Int_t i = 0; i < 4; ++i) {
     fPhiDist[i] = 0;
  }
  for(Int_t i = 0; i < 2; ++i) {
     fQDist[i] = 0;
  }
}

//________________________________________________________________________
AliEPSelectionTask::AliEPSelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisInput("ESD"),
  fTrackType("TPC"),
  fPeriod(""),
  fCentrality(-1.),
  fUseMCRP(kFALSE),
  fUsePhiWeight(kFALSE),
  fUsePtWeight(kFALSE),
  fUseRecentering(kFALSE),
  fSaveTrackContribution(kFALSE),
  fUserphidist(kFALSE),
  fUsercuts(kFALSE),
  fRunNumber(-15),
  fAODfilterbit(1),
  fEtaGap(0.),
  fSplitMethod(0),
  fESDtrackCuts(0),
  fEPContainer(0),
  fQxContainer(0),
  fQyContainer(0),
  fSparseDist(0),
  fHruns(0),
  fQVector(0),
  fQContributionX(0),
  fQContributionY(0),
  fEventplaneQ(0),
  fQsub1(0),
  fQsub2(0),
  fQsubRes(0),
  fOutputList(0),
  fHOutEventplaneQ(0),
  fHOutPhi(0),
  fHOutPhiCorr(0),
  fHOutsub1sub2(0),
  fHOutNTEPRes(0),
  fHOutPTPsi(0),
  fHOutDiff(0),
  fHOutleadPTPsi(0)
{
  // Default constructor
  AliInfo("Event Plane Selection enabled.");
  DefineOutput(1, TList::Class());
  for(Int_t i = 0; i < 4; i++) {
     fPhiDist[i] = 0;
  }
  for(Int_t i = 0; i < 2; ++i) {
     fQDist[i] = 0;
  }
}

//________________________________________________________________________
AliEPSelectionTask::~AliEPSelectionTask()
{
  // Destructor
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
      delete fOutputList;
      fOutputList = 0;
  }
  if (fESDtrackCuts){
      delete fESDtrackCuts;
      fESDtrackCuts = 0;
  }
  if (fUserphidist) {
    if (fPhiDist[0]) {
      delete fPhiDist[0];
      fPhiDist[0] = 0;
    }
  }
  if (fEPContainer){
      delete fEPContainer;
      fEPContainer = 0;
  }
  if (fPeriod.CompareTo("LHC11h")==0){
      for(Int_t i = 0; i < 4; i++) {
        if(fPhiDist[i]){
          delete fPhiDist[i];
          fPhiDist[i] = 0;
        }
      }
      if(fHruns) delete fHruns;
  }
  if(fQDist[0] && fQDist[1]) {
    for(Int_t i = 0; i < 2; i++) {
      if(fQDist[i]){
	delete fQDist[i];
	fQDist[i] = 0;
      }
    }
  }
}

//________________________________________________________________________
void AliEPSelectionTask::UserCreateOutputObjects()
{
  // Create the output containers
  if (fDebug>1) printf("AliEPSelectionTask::UserCreateOutputObjects() \n");
  AliLog::SetClassDebugLevel("AliEPSelectionTask", AliLog::kInfo);

  fOutputList = new TList();
  fOutputList->SetOwner();
  fHOutEventplaneQ = new TH1F("fHOutEventplaneQ","fHOutEventplaneQ; Event Plane Q",100,0,TMath::Pi());
  fHOutPhi = new TH1F("fHOutPhi","fHOutPhi; Phi Distribution",100,0,TMath::TwoPi());
  fHOutPhiCorr = new TH1F("fHOutPhiCorr","fHOutPhiCorr; Corrected Phi Distribution",100,0,TMath::TwoPi());
  fHOutsub1sub2 = new TH2F("fHOutsub1sub2","fHOutsub1sub2; EP1; EP2",100,0,TMath::Pi(),100,0,TMath::Pi());
  fHOutNTEPRes = new TH2F("fHOutNTEPRes","fHOutNTEPRes; Number of Tracks; Event Plane Resolution",100,0,5000,100,-TMath::Pi(),TMath::Pi());
  fHOutPTPsi = new TH2F("fHOutPTPsi","fHOutPTPsi; PT; Phi-EP",100,0,20,100,0,TMath::Pi());
  fHOutDiff = new TH2F("fHOutDiff","fHOutDiff; EP; MCEP",100,0,TMath::Pi(),100,0,TMath::Pi());
  fHOutleadPTPsi = new TH2F("fHOutleadPTPsi","fHOutleadPTPsi; leadPT; EP",100,0,TMath::Pi(),100,0,TMath::Pi());

  fOutputList->Add(fHOutEventplaneQ);
  fOutputList->Add(fHOutPhi);
  fOutputList->Add(fHOutPhiCorr);
  fOutputList->Add(fHOutsub1sub2);
  fOutputList->Add(fHOutNTEPRes);
  fOutputList->Add(fHOutPTPsi);
  fOutputList->Add(fHOutleadPTPsi);
  fOutputList->Add(fHOutDiff);

  PostData(1, fOutputList);


}

//________________________________________________________________________
void AliEPSelectionTask::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  if (fDebug>1) printf(" **** AliEPSelectionTask::UserExec() \n");

//   fRunNumber = -15;

  AliEventplane *esdEP;
  TVector2 qq1;
  TVector2 qq2;
  Double_t fRP = 0.; // monte carlo reaction plane angle

  if (fAnalysisInput.CompareTo("ESD")==0){

    AliVEvent* event = InputEvent();
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (esd){
      if (!(fRunNumber == esd->GetRunNumber())) {
	  fRunNumber = esd->GetRunNumber();
          AliInfo(Form("Changing Phi-distribution to run %d",fRunNumber));
          SetOADBandPeriod();
	  SetPhiDist();
	  SetQvectorDist(); // recentring objects
      }


      if (fUseMCRP) {
	  AliMCEvent* mcEvent  = MCEvent();
	  if (mcEvent && mcEvent->GenEventHeader()) {
	      AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
	      if (headerH) fRP = headerH->ReactionPlaneAngle();
	  }
      }

      esdEP = esd->GetEventplane();
      if (fSaveTrackContribution) {
	esdEP->GetQContributionXArray()->Set(esd->GetNumberOfTracks());
	esdEP->GetQContributionYArray()->Set(esd->GetNumberOfTracks());
        esdEP->GetQContributionXArraysub1()->Set(esd->GetNumberOfTracks());
	esdEP->GetQContributionYArraysub1()->Set(esd->GetNumberOfTracks());
        esdEP->GetQContributionXArraysub2()->Set(esd->GetNumberOfTracks());
	esdEP->GetQContributionYArraysub2()->Set(esd->GetNumberOfTracks());
      }

      TObjArray* tracklist = new TObjArray;
      if (fTrackType.CompareTo("GLOBAL")==0) tracklist = fESDtrackCuts->GetAcceptedTracks(esd,kFALSE);
      if (fTrackType.CompareTo("TPC")==0 && fPeriod.CompareTo("LHC10h")==0) tracklist = fESDtrackCuts->GetAcceptedTracks(esd,kTRUE);
      else if (fTrackType.CompareTo("TPC")==0 && fPeriod.CompareTo("LHC11h")==0) tracklist = GetTracksForLHC11h(esd);
      const int nt = tracklist->GetEntries();

      if (nt>4){

	// qvector full event
	fQVector = new TVector2(GetQ(esdEP,tracklist));
	fEventplaneQ = fQVector->Phi()/2;
	// q vector subevents
	GetQsub(qq1, qq2, tracklist, esdEP);
	fQsub1 = new TVector2(qq1);
	fQsub2 = new TVector2(qq2);
	fQsubRes = (fQsub1->Phi()/2 - fQsub2->Phi()/2);

	esdEP->SetQVector(fQVector);
	esdEP->SetEventplaneQ(fEventplaneQ);
	esdEP->SetQsub(fQsub1,fQsub2);
	esdEP->SetQsubRes(fQsubRes);

	fHOutEventplaneQ->Fill(fEventplaneQ);
	fHOutsub1sub2->Fill(fQsub1->Phi()/2,fQsub2->Phi()/2);
	fHOutNTEPRes->Fill(nt,fQsubRes);

	if (fUseMCRP) fHOutDiff->Fill(fEventplaneQ, fRP);

	for (int iter = 0; iter<nt;iter++){
	  AliESDtrack* track = dynamic_cast<AliESDtrack*> (tracklist->At(iter));
	  if (track) {
	    float delta = track->Phi()-fEventplaneQ;
	    while (delta < 0) delta += TMath::Pi();
	    while (delta > TMath::Pi()) delta -= TMath::Pi();
	    fHOutPTPsi->Fill(track->Pt(),delta);
	    fHOutPhi->Fill(track->Phi());
	    fHOutPhiCorr->Fill(track->Phi(),GetPhiWeight(track));
          }
	}

	AliESDtrack* trmax = esd->GetTrack(0);
	for (int iter = 1; iter<nt;iter++){
	  AliESDtrack* track = dynamic_cast<AliESDtrack*> (tracklist->At(iter));
	  if (track && (track->Pt() > trmax->Pt())) trmax = track;
	}
	fHOutleadPTPsi->Fill(trmax->Phi(),fEventplaneQ);
      }
      tracklist->Clear();
      delete tracklist;
      tracklist = 0;
    }
  }

    else if (fAnalysisInput.CompareTo("AOD")==0){
    AliVEvent* event = InputEvent();
    AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);

    if (aod){
      // get centrality of the event
      AliAODHeader *header=dynamic_cast<AliAODHeader*>(aod->GetHeader());
      if(!header) AliFatal("Not a standard AOD");
      if(AliMultSelection *multSelection = (AliMultSelection*)event->FindListObject("MultSelection")){
        fCentrality = multSelection->GetMultiplicityPercentile("V0M",kFALSE);
      }
      else{
        AliCentrality *centrality=header->GetCentralityP();
        if(!centrality)  AliError(Form("No AliCentrality attached to AOD header"));
        fCentrality = centrality->GetCentralityPercentile("V0M");
      }
      if (!(fRunNumber == aod->GetRunNumber())) {
        fRunNumber = aod->GetRunNumber();
        AliInfo(Form("Changing Phi-distribution to run %d",fRunNumber));
        SetOADBandPeriod();
        SetPhiDist();
        SetQvectorDist();
      }

      if (fUseMCRP) {
	AliAODMCHeader *headerH = dynamic_cast<AliAODMCHeader*>(aod->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
	if (headerH) fRP = headerH->GetReactionPlaneAngle();
      }

      esdEP = header->GetEventplaneP();
      if(!esdEP) return; // protection against missing EP branch (nanoAODs)
      esdEP->Reset();

      Int_t maxID = 0;
      TObjArray* tracklist = GetAODTracksAndMaxID(aod,maxID);

      if (fSaveTrackContribution) {
	esdEP->GetQContributionXArray()->Set(maxID+1);
	esdEP->GetQContributionYArray()->Set(maxID+1);
	esdEP->GetQContributionXArraysub1()->Set(maxID+1);
	esdEP->GetQContributionYArraysub1()->Set(maxID+1);
	esdEP->GetQContributionXArraysub2()->Set(maxID+1);
	esdEP->GetQContributionYArraysub2()->Set(maxID+1);
      }

      const int NT = tracklist->GetEntries();

      if (NT>4){

	// qvector full event
	fQVector = new TVector2(GetQ(esdEP,tracklist));
	fEventplaneQ = fQVector->Phi()/2;
	// q vector subevents
	GetQsub(qq1, qq2, tracklist, esdEP);
	fQsub1 = new TVector2(qq1);
	fQsub2 = new TVector2(qq2);
	fQsubRes = (fQsub1->Phi()/2 - fQsub2->Phi()/2);

	esdEP->SetQVector(fQVector);
	esdEP->SetEventplaneQ(fEventplaneQ);
	esdEP->SetQsub(fQsub1,fQsub2);
	esdEP->SetQsubRes(fQsubRes);

	fHOutEventplaneQ->Fill(fEventplaneQ);
	fHOutsub1sub2->Fill(fQsub1->Phi()/2.,fQsub2->Phi()/2.);
	fHOutNTEPRes->Fill(NT,fQsubRes);

	if (fUseMCRP) fHOutDiff->Fill(fEventplaneQ, fRP);

	for (int iter = 0; iter<NT;iter++){
	  AliAODTrack* track = dynamic_cast<AliAODTrack*> (tracklist->At(iter));
	  if (track) {
	    float delta = track->Phi()-fEventplaneQ;
	    while (delta < 0) delta += TMath::Pi();
	    while (delta > TMath::Pi()) delta -= TMath::Pi();
	    fHOutPTPsi->Fill(track->Pt(),delta);
	    fHOutPhi->Fill(track->Phi());
	    fHOutPhiCorr->Fill(track->Phi(),GetPhiWeight(track));
	  }
	}

	AliAODTrack* trmax = dynamic_cast<AliAODTrack*>(aod->GetTrack(0));
	if(!trmax) AliFatal("Not a standard AOD");
	for (int iter = 1; iter<NT;iter++){
	  AliAODTrack* track = dynamic_cast<AliAODTrack*> (tracklist->At(iter));
	  if (track && (track->Pt() > trmax->Pt())) trmax = track;
	}
	fHOutleadPTPsi->Fill(trmax->Phi(),fEventplaneQ);
      }
      delete tracklist;
      tracklist = 0;
    }


  }


  else {
    printf(" Analysis Input not known!\n\n ");
    return;
  }
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliEPSelectionTask::Terminate(Option_t */*option*/)
{
  // Terminate analysis
}

//__________________________________________________________________________
TVector2 AliEPSelectionTask::GetQ(AliEventplane* EP, TObjArray* tracklist)
{
  // Get the Q vector
  TVector2 mQ;
  float mQx=0, mQy=0;
  AliVTrack* track;
  Double_t weight;
  Int_t idtemp = -1;
  // get recentering values
  Double_t mean[2], rms[2];
  Recenter(0, mean);
  Recenter(1, rms);

  int nt = tracklist->GetEntries();

  for (int i=0; i<nt; i++){
    weight = 1;
    track = dynamic_cast<AliVTrack*> (tracklist->At(i));
    if (track) {
      weight = GetWeight(track);
    if (fSaveTrackContribution){
      idtemp = track->GetID();
      if ((fAnalysisInput.CompareTo("AOD")==0) && (fAODfilterbit == 128)) idtemp = idtemp*(-1) - 1;
      EP->GetQContributionXArray()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
      EP->GetQContributionYArray()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
     }
     mQx += (weight*cos(2*track->Phi())/rms[0]);
     mQy += (weight*sin(2*track->Phi())/rms[1]);
    }
  }
  mQ.Set(mQx-(mean[0]/rms[0]), mQy-(mean[1]/rms[1]));
  return mQ;
}

  //________________________________________________________________________
void AliEPSelectionTask::GetQsub(TVector2 &Q1, TVector2 &Q2, TObjArray* tracklist,AliEventplane* EP)
{
  // Get Qsub
  TVector2 mQ[2];
  float mQx1=0, mQy1=0, mQx2=0, mQy2=0;
  Double_t weight;
  // get recentering values
  Double_t mean[2], rms[2];
  Recenter(0, mean);
  Recenter(1, rms);

  AliVTrack* track;
  TRandom2 rn = 0;

  int nt = tracklist->GetEntries();
  int trackcounter1=0, trackcounter2=0;
  int idtemp = 0;

  if (fSplitMethod == AliEPSelectionTask::kRandom){

    for (Int_t i = 0; i < nt; i++) {
      weight = 1;
      track = dynamic_cast<AliVTrack*> (tracklist->At(i));
      if (!track) continue;
      weight = GetWeight(track);
      idtemp = track->GetID();
      if ((fAnalysisInput.CompareTo("AOD")==0) && (fAODfilterbit == 128)) idtemp = idtemp*(-1) - 1;

      // This loop splits the track set into 2 random subsets
      if( trackcounter1 < int(nt/2.) && trackcounter2 < int(nt/2.)){
        float random = rn.Rndm();
        if(random < .5){
          mQx1 += (weight*cos(2*track->Phi())/rms[0]);
          mQy1 += (weight*sin(2*track->Phi())/rms[1]);
          if (fSaveTrackContribution){
            EP->GetQContributionXArraysub1()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
            EP->GetQContributionYArraysub1()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
          }
          trackcounter1++;
        }
        else {
          mQx2 += (weight*cos(2*track->Phi())/rms[0]);
          mQy2 += (weight*sin(2*track->Phi())/rms[1]);
          if (fSaveTrackContribution){
            EP->GetQContributionXArraysub2()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
            EP->GetQContributionYArraysub2()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
          }
          trackcounter2++;
        }
      }
      else if( trackcounter1 >= int(nt/2.)){
        mQx2 += (weight*cos(2*track->Phi())/rms[0]);
        mQy2 += (weight*sin(2*track->Phi())/rms[1]);
        if (fSaveTrackContribution){
          EP->GetQContributionXArraysub2()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
          EP->GetQContributionYArraysub2()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
        }
        trackcounter2++;
      }
      else {
        mQx1 += (weight*cos(2*track->Phi())/rms[0]);
        mQy1 += (weight*sin(2*track->Phi())/rms[1]);
        if (fSaveTrackContribution){
          EP->GetQContributionXArraysub1()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
          EP->GetQContributionYArraysub1()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
        }
        trackcounter1++;
      }
    }
  } else if (fSplitMethod == AliEPSelectionTask::kEta) {

    for (Int_t i = 0; i < nt; i++) {
      weight = 1;
      track = dynamic_cast<AliVTrack*> (tracklist->At(i));
      if (!track) continue;
      weight = GetWeight(track);
      Double_t eta = track->Eta();
      idtemp = track->GetID();
      if ((fAnalysisInput.CompareTo("AOD")==0) && (fAODfilterbit == 128)) idtemp = idtemp*(-1) - 1;

      if (eta > fEtaGap/2.) {
        mQx1 += (weight*cos(2*track->Phi())/rms[0]);
        mQy1 += (weight*sin(2*track->Phi())/rms[1]);
        if (fSaveTrackContribution){
          EP->GetQContributionXArraysub1()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
          EP->GetQContributionYArraysub1()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
        }
      } else if (eta < -1.*fEtaGap/2.) {
        mQx2 += (weight*cos(2*track->Phi())/rms[0]);
        mQy2 += (weight*sin(2*track->Phi())/rms[1]);
        if (fSaveTrackContribution){
          EP->GetQContributionXArraysub2()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
          EP->GetQContributionYArraysub2()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
        }
      }
    }
  } else if (fSplitMethod == AliEPSelectionTask::kCharge) {

    for (Int_t i = 0; i < nt; i++) {
      weight = 1;
      track = dynamic_cast<AliVTrack*> (tracklist->At(i));
      if (!track) continue;
      weight = GetWeight(track);
      Short_t cha = track->Charge();
      idtemp = track->GetID();
      if ((fAnalysisInput.CompareTo("AOD")==0) && (fAODfilterbit == 128)) idtemp = idtemp*(-1) - 1;

      if (cha > 0) {
        mQx1 += (weight*cos(2*track->Phi())/rms[0]);
        mQy1 += (weight*sin(2*track->Phi())/rms[1]);
        if (fSaveTrackContribution){
          EP->GetQContributionXArraysub1()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
          EP->GetQContributionYArraysub1()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
        }
      } else if (cha < 0) {
        mQx2 += (weight*cos(2*track->Phi())/rms[0]);
        mQy2 += (weight*sin(2*track->Phi())/rms[1]);
        if (fSaveTrackContribution){
          EP->GetQContributionXArraysub2()->AddAt(weight*cos(2*track->Phi())/rms[0],idtemp);
          EP->GetQContributionYArraysub2()->AddAt(weight*sin(2*track->Phi())/rms[1],idtemp);
        }
      }
    }
  } else {
    printf("plane resolution determination method not available!\n\n ");
    return;
  }
  // apply recenetering
  mQ[0].Set(mQx1-(mean[0]/rms[0]), mQy1-(mean[1]/rms[1]));
  mQ[1].Set(mQx2-(mean[0]/rms[0]), mQy2-(mean[1]/rms[1]));
  Q1 = mQ[0];
  Q2 = mQ[1];
}

//________________________________________________________________________
void AliEPSelectionTask::SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts){

  if(fESDtrackCuts){
    delete fESDtrackCuts;
    fESDtrackCuts = 0;
  }
  if (fAnalysisInput.CompareTo("AOD")==0){
    AliInfo("ESD track cuts not possible for AOD analysis; please use SetPersonalAODtrackCuts(); using TPC only track cuts");
    fUsercuts = kFALSE;
    SetTrackType("TPC");
    return;
  }
  fUsercuts = kTRUE;
  fESDtrackCuts = trackcuts;
}

//________________________________________________________________________
void AliEPSelectionTask::SetPersonalAODtrackCuts(UInt_t filterbit, Float_t etalow, Float_t etaup, Float_t ptlow, Float_t ptup, Int_t ntpc){

  if(fESDtrackCuts){
    delete fESDtrackCuts;
    fESDtrackCuts = 0;
  }
  if (fAnalysisInput.CompareTo("ESD")==0){
    AliInfo("AOD track cuts not possible for ESD analysis; please use SetPersonalESDtrackCuts(); using TPC only track cuts");
    fUsercuts = kFALSE;
    SetTrackType("TPC");
    return;
  }
  fUsercuts = kTRUE;
  fESDtrackCuts = new AliESDtrackCuts();
  fESDtrackCuts->SetPtRange(ptlow,ptup);
  fESDtrackCuts->SetMinNClustersTPC(ntpc);
  fESDtrackCuts->SetEtaRange(etalow,etaup);
  fAODfilterbit = filterbit;
}

//_____________________________________________________________________________

void AliEPSelectionTask::SetTrackType(TString tracktype){
  fTrackType = tracktype;
  if (!fUsercuts) {
  if (fTrackType.CompareTo("GLOBAL")==0){
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    fAODfilterbit = 32;
  }
  if (fTrackType.CompareTo("TPC")==0){
    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fAODfilterbit = 128;
  }
  fESDtrackCuts->SetPtRange(0.15,20.);
  fESDtrackCuts->SetEtaRange(-0.8,0.8);
  }
}

//________________________________________________________________________
Double_t AliEPSelectionTask::GetWeight(TObject* track1)
{
  Double_t ptweight=1;
  AliVTrack* track = dynamic_cast<AliVTrack*>(track1);
  if (fUsePtWeight && track) {
    if (track->Pt()<2) ptweight=track->Pt();
    else ptweight=2;
  }
  return ptweight*GetPhiWeight(track);
}

//________________________________________________________________________
Double_t AliEPSelectionTask::GetPhiWeight(TObject* track1)
{
  Double_t phiweight=1;
  AliVTrack* track = dynamic_cast<AliVTrack*>(track1);

  TH1F *phiDist = 0x0;
  if(track) phiDist = SelectPhiDist(track);

  if (fUsePhiWeight && phiDist && track) {
    Double_t nParticles = phiDist->Integral();
    Double_t nPhibins = phiDist->GetNbinsX();

    Double_t Phi = track->Phi();

    while (Phi<0) Phi += TMath::TwoPi();
    while (Phi>TMath::TwoPi()) Phi -= TMath::TwoPi();

    Double_t PhiDistValue = phiDist->GetBinContent(1+TMath::FloorNint((track->Phi())*nPhibins/TMath::TwoPi()));

    if (PhiDistValue > 0) phiweight = nParticles/nPhibins/PhiDistValue;
  }
  return phiweight;
}

//________________________________________________________________________
void AliEPSelectionTask::Recenter(Int_t var, Double_t * values)
{

  if (fUseRecentering && fQDist[0] && fQDist[1] && fCentrality!=-1.) {
    Int_t centbin = fQDist[0]->FindBin(fCentrality);

    if(var==0) { // fill mean
      values[0] = fQDist[0]->GetBinContent(centbin);
      values[1] = fQDist[1]->GetBinContent(centbin);
    }
    else if(var==1) { // fill rms
      values[0] = fQDist[0]->GetBinError(centbin);
      values[1] = fQDist[1]->GetBinError(centbin);
      // protection against division by zero
      if(values[0]==0.0) values[0]=1.0;
      if(values[1]==0.0) values[1]=1.0;
    }
  }
  else { //default (no recentering)
    values[0] = var;
    values[1] = var;
  }
  return;
}

//__________________________________________________________________________
void AliEPSelectionTask::SetPhiDist()
{
  if(!fUserphidist && (fPeriod.CompareTo("LHC10h") == 0 || fPeriod.CompareTo("LHC11h") == 0)) { // if it's already set and custom class is required, we use the one provided by the user

    if (fPeriod.CompareTo("LHC10h")==0)
       {
        fPhiDist[0] = (TH1F*) fEPContainer->GetObject(fRunNumber, "Default");}
        else if(fPeriod.CompareTo("LHC11h")==0){
            Int_t runbin=fHruns->FindBin(fRunNumber);
            if (fHruns->GetBinContent(runbin) > 1){
              fSparseDist->GetAxis(0)->SetRange(runbin,runbin);
              }
            else if(fHruns->GetBinContent(runbin) < 2){
               fSparseDist->GetAxis(0)->SetRange(1,2901); // not calibrated run, use integrated phi-weights
               AliInfo("Using integrated Phi-weights for this run");
               }
            for (Int_t i = 0; i<4 ;i++)
            {
	     if(fPhiDist[i]){
               delete fPhiDist[i];
               fPhiDist[i] = 0x0;
               }
             if(i == 0){
               fSparseDist->GetAxis(1)->SetRange(1,1);  // neg charge
               fSparseDist->GetAxis(2)->SetRange(1,1);} // neg eta
             if(i == 1){
               fSparseDist->GetAxis(1)->SetRange(2,2);  // pos charge
               fSparseDist->GetAxis(2)->SetRange(1,1);} // neg eta
             if(i == 2){
               fSparseDist->GetAxis(1)->SetRange(1,1);  // neg charge
               fSparseDist->GetAxis(2)->SetRange(2,2);} // pos eta
             if(i == 3){
               fSparseDist->GetAxis(1)->SetRange(2,2);  // pos charge
               fSparseDist->GetAxis(2)->SetRange(2,2);} // pos eta
             fPhiDist[i] = (TH1F*)fSparseDist->Projection(3); // Projection on Phi
             fPhiDist[i]->SetName(Form("phidist%d%d",i,fRunNumber));
             fSparseDist->GetAxis(1)->SetRange(1,2); // reset axes
             fSparseDist->GetAxis(2)->SetRange(1,2);
             }
             fSparseDist->GetAxis(0)->SetRange(1,2901);// reset run axis
       }

    if (!fPhiDist[0]) AliFatal(Form("Cannot find OADB phi distribution for run %d", fRunNumber));

  }


  if (fPeriod.CompareTo("LHC10h")==0 || fUserphidist){
     Bool_t emptybins;

     int iter = 0;
     while (iter<3){
         emptybins = kFALSE;

         for (int i=1; i<fPhiDist[0]->GetNbinsX(); i++){
	   if (!((fPhiDist[0]->GetBinContent(i))>0)) {
	     emptybins = kTRUE;
	   }
         }
         if (emptybins) {
	   cout << "empty bins - rebinning!" << endl;
	   fPhiDist[0]->Rebin();
	   iter++;
         }
         else iter = 3;
     }

     if (emptybins) {
       AliError("After Maximum of rebinning still empty Phi-bins!!!");
     }
  }
  if (fPeriod.CompareTo("LHC10h") != 0 && fPeriod.CompareTo("LHC11h") != 0 && !fUserphidist){
  AliInfo("No Phi-weights available. All Phi weights set to 1");
  SetUsePhiWeight(kFALSE);
  }
}

//__________________________________________________________________________
void AliEPSelectionTask::SetQvectorDist()
{
  if(!fUseRecentering) return;
  AliInfo(Form("Setting q vector distributions"));
  fQDist[0] = (TProfile*) fQxContainer->GetObject(fRunNumber, "Default");
  fQDist[1] = (TProfile*) fQyContainer->GetObject(fRunNumber, "Default");

  if (!fQDist[0] || !fQDist[1]) {
    AliError(Form("Cannot find OADB q-vector distributions for run %d. Using default values (mean=0,rms=1).", fRunNumber));
    return;
  }

  Bool_t emptybins;

  int iter = 0;
  while (iter<3){
    emptybins = kFALSE;

    for (int i=1; i<=fQDist[0]->GetNbinsX()-2; i++){
      if (!((fQDist[0]->GetBinEntries(i))>0)) {
	emptybins = kTRUE;
      }
    }
    if (emptybins) {
      cout << "empty bins - rebinning!" << endl;
      fQDist[0]->Rebin();
      fQDist[1]->Rebin();
      iter++;
    }
    else iter = 3;
  }

  if (emptybins) {
    AliError("After Maximum of rebinning still empty Qxy-bins!!!");
  }
}

//__________________________________________________________________________
void AliEPSelectionTask::SetPersonalPhiDistribution(const char* infilename, char* listname)
{

  fUserphidist = kTRUE;

  TFile f(infilename);
  TObject* list = f.Get(listname);
  fPhiDist[0] = (TH1F*)list->FindObject("fHOutPhi");
  if (!fPhiDist[0]) AliFatal("Phi Distribution not found!!!");

  f.Close();
}

//_________________________________________________________________________
TObjArray* AliEPSelectionTask::GetAODTracksAndMaxID(AliAODEvent* aod, Int_t& maxid)
{
  TObjArray *acctracks = new TObjArray();

  AliAODTrack *tr = 0;
  Int_t maxid1 = 0;
  Int_t maxidtemp = -1;
  Float_t ptlow = 0;
  Float_t ptup = 0;
  Float_t etalow = 0;
  Float_t etaup = 0;
  fESDtrackCuts->GetPtRange(ptlow,ptup);
  fESDtrackCuts->GetEtaRange(etalow,etaup);
  Int_t ntpc = fESDtrackCuts->GetMinNClusterTPC();

  for (Int_t i = 0; i < aod->GetNumberOfTracks() ; i++){
     tr = dynamic_cast<AliAODTrack*>(aod->GetTrack(i));
     if(!tr) AliFatal("Not a standard AOD");
     maxidtemp = tr->GetID();
     if(maxidtemp < 0 && fAODfilterbit != 128) continue;
     if(maxidtemp > -1 && fAODfilterbit == 128) continue;
     if (fAODfilterbit == 128) maxidtemp = maxidtemp*(-1) - 1;
     if (maxidtemp > maxid1) maxid1 = maxidtemp;
     if(tr->TestFilterBit(fAODfilterbit) && tr->Pt() < ptup && tr->Pt() > ptlow && tr->Eta() < etaup && tr->Eta() > etalow && tr->GetTPCNcls() > ntpc){
     acctracks->Add(tr);
     }
  }

  maxid = maxid1;

  return acctracks;

}

//_________________________________________________________________________
void AliEPSelectionTask::SetOADBandPeriod()
{
  TString oadbfilename;

  if (fRunNumber >= 136851 && fRunNumber <= 139515) // LHC10h
    {fPeriod = "LHC10h";
     if (!fUserphidist) { // if it's already set and custom class is required, we use the one provided by the user

        if (fAnalysisInput.CompareTo("AOD")==0){
           oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.aod.root", AliAnalysisManager::GetOADBPath()));
           } else if (fAnalysisInput.CompareTo("ESD")==0){
           oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.root", AliAnalysisManager::GetOADBPath()));
           }

       TFile foadb(oadbfilename);
       if(!foadb.IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

       AliInfo("Using Standard OADB");
       fEPContainer = (AliOADBContainer*) foadb.Get("epphidist");
       if (!fEPContainer) AliFatal("Cannot fetch OADB container for EP selection");
       foadb.Close();
       }
     }

  if (fRunNumber >= 166529 && fRunNumber <= 170593) // LHC11h
     {fPeriod = "LHC11h";
      if (!fUserphidist) {
      // if it's already set and custom class is required, we use the one provided by the user

      oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist2011.root", AliAnalysisManager::GetOADBPath()));
      TFile *foadb = TFile::Open(oadbfilename);
      if(!foadb->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

      AliInfo("Using Standard OADB");
      fSparseDist = (THnSparse*) foadb->Get("Default");
      if (!fSparseDist) AliFatal("Cannot fetch OADB container for EP selection");
      foadb->Close();
      if(!fHruns){
           fHruns = (TH1F*)fSparseDist->Projection(0); //projection on run axis;
           fHruns->SetName("runsHisto");
      }
      }

      if(fUseRecentering) {
	oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/eprecentering.root", AliAnalysisManager::GetOADBPath()));
	TFile *foadb = TFile::Open(oadbfilename);
	if(!foadb->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

	AliInfo("Using Standard OADB");
	fQxContainer = (AliOADBContainer*) foadb->Get("eprecentering.Qx");
	fQyContainer = (AliOADBContainer*) foadb->Get("eprecentering.Qy");
	if (!fQxContainer || !fQyContainer) AliFatal("Cannot fetch OADB container for EP recentering");
	foadb->Close();
      }

     }
}

//_________________________________________________________________________
TH1F* AliEPSelectionTask::SelectPhiDist(AliVTrack *track)
{
  if (fPeriod.CompareTo("LHC10h")==0  || fUserphidist) return fPhiDist[0];
  else if(fPeriod.CompareTo("LHC11h")==0)
    {
     if (track->Charge() < 0)
       {
        if(track->Eta() < 0.)       return fPhiDist[0];
        else if (track->Eta() > 0.) return fPhiDist[2];
       }
      else if (track->Charge() > 0)
       {
        if(track->Eta() < 0.)       return fPhiDist[1];
        else if (track->Eta() > 0.) return fPhiDist[3];
       }

    }
  return 0;
}

TObjArray* AliEPSelectionTask::GetTracksForLHC11h(AliESDEvent* esd)
{
  // Need to do this hack beacuse only this type of TPC only tracks in AOD is available and one type of Phi-weights is used
  TObjArray *acctracks = new TObjArray();
  acctracks->SetOwner(kTRUE);

  const AliESDVertex *vtxSPD = esd->GetPrimaryVertexSPD();

  Float_t ptlow = 0;
  Float_t ptup = 0;
  Float_t etalow = 0;
  Float_t etaup = 0;
  fESDtrackCuts->GetPtRange(ptlow,ptup);
  fESDtrackCuts->GetEtaRange(etalow,etaup);

  Double_t pDCA[3] = { 0. }; // momentum at DCA
  Double_t rDCA[3] = { 0. }; // position at DCA
  Float_t  dDCA[2] = {0.};    // DCA to the vertex d and z
  Float_t  cDCA[3] = {0.};    // covariance of impact parameters



  for (Int_t i = 0; i < esd->GetNumberOfTracks() ; i++){

    AliESDtrack* esdTrack = esd->GetTrack(i); //carefull do not modify it othwise  need to work with a copy
    //
    // Track selection
    if (!fESDtrackCuts->AcceptTrack(esdTrack)) continue;

    // create a tpc only tracl
    AliESDtrack *track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esd),esdTrack->GetID());
    if(!track) continue;

    if(track->Pt()>0.)
    {
      // only constrain tracks above threshold
      AliExternalTrackParam exParam;
      // take the B-field from the ESD, no 3D fieldMap available at this point
      Bool_t relate = false;
      relate = track->RelateToVertexTPC(vtxSPD,esd->GetMagneticField(),kVeryBig,&exParam);
      if(!relate){
        delete track;
        continue;
      }
      // fetch the track parameters at the DCA (unconstraint)
      if(track->GetTPCInnerParam()){
	track->GetTPCInnerParam()->GetPxPyPz(pDCA);
	track->GetTPCInnerParam()->GetXYZ(rDCA);
      }
      // get the DCA to the vertex:
      track->GetImpactParametersTPC(dDCA,cDCA);
      // set the constrained parameters to the track
      track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
    }


    Float_t eta = track->Eta();
    Float_t pT = track->Pt();

    if(pT < ptlow || pT > ptup || eta < etalow || eta > etaup){
      delete track;
      continue;
    }

    acctracks->Add(track);
   }


  return acctracks;

}
