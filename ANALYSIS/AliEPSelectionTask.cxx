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
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliESDUtils.h"
#include "AliOADBContainer.h"

#include "AliEventplane.h"

ClassImp(AliEPSelectionTask)

//________________________________________________________________________
AliEPSelectionTask::AliEPSelectionTask():
AliAnalysisTaskSE(),
  fAnalysisInput("ESD"),
  fTrackType("TPC"),
  fUseMCRP(kFALSE),
  fUsePhiWeight(kFALSE),
  fUsePtWeight(kFALSE),
  fSaveTrackContribution(kFALSE),
  fuserphidist(kFALSE),
  fusercuts(kFALSE),
  frunNumber(-15),
  fESDtrackCuts(0),
  ftracklist(0),
  fEPContainer(0),
  fPhiDist(0),
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
}   

//________________________________________________________________________
AliEPSelectionTask::AliEPSelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisInput("ESD"),
  fTrackType("TPC"),
  fUseMCRP(kFALSE),
  fUsePhiWeight(kFALSE),
  fUsePtWeight(kFALSE),  
  fSaveTrackContribution(kFALSE),
  fuserphidist(kFALSE),
  fusercuts(kFALSE),
  frunNumber(-15),
  fESDtrackCuts(0),
  ftracklist(0),
  fEPContainer(0),
  fPhiDist(0),
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
  if (fQVector){
      delete fQVector;
      fQVector = 0;
  }
  if (fQsub1){
      delete fQsub1;
      fQsub1 = 0;
  }
  if (fQsub2){
      delete fQsub2;
      fQsub2 = 0;
  }
  if (fPhiDist){
      delete fPhiDist;
      fPhiDist = 0;
  }
  if (ftracklist){
      delete ftracklist;
      ftracklist = 0;
  }
  if (fEPContainer){
      delete fEPContainer;
      fEPContainer = 0;
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
  
    if(!fuserphidist) { // if it's already set and custom class is required, we use the one provided by the user

    TString oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.root", AliAnalysisManager::GetOADBPath()));

    TFile foadb(oadbfilename); 
    if(!foadb.IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

    AliInfo("Using Standard OADB");
    fEPContainer = (AliOADBContainer*) foadb.Get("epphidist");    
    if (!fEPContainer) AliFatal("Cannot fetch OADB container for EP selection");
    foadb.Close();
    }
}

//________________________________________________________________________
void AliEPSelectionTask::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if (fDebug>1) printf(" **** AliEPSelectionTask::UserExec() \n");
  
//   frunNumber = -15;
 
  AliEventplane* esdEP = 0;
  TVector2 qq1;
  TVector2 qq2;
  Double_t fRP = 0.; // the monte carlo reaction plane angle
    
  if (fAnalysisInput.CompareTo("ESD")==0){

    AliVEvent* event = InputEvent();
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (esd){    
      if (!(frunNumber == esd->GetRunNumber())) {
	  frunNumber = esd->GetRunNumber();
	    SetPhiDist();      
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
      }
      
      if (fTrackType.CompareTo("GLOBAL")==0) ftracklist = fESDtrackCuts->GetAcceptedTracks(esd,kFALSE);
      if (fTrackType.CompareTo("TPC")==0) ftracklist = fESDtrackCuts->GetAcceptedTracks(esd,kTRUE);
      const int nt = ftracklist->GetEntries();
      
      if (nt>4){
	fQVector = new TVector2(GetQ(esdEP,ftracklist));
	fEventplaneQ = fQVector->Phi()/2; 
	GetQsub(qq1, qq2, ftracklist);
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
	  AliESDtrack* track = dynamic_cast<AliESDtrack*> (ftracklist->At(iter));
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
	  AliESDtrack* track = dynamic_cast<AliESDtrack*> (ftracklist->At(iter));
	  if (track && (track->Pt() > trmax->Pt())) trmax = track;
	}
	fHOutleadPTPsi->Fill(trmax->Phi(),fEventplaneQ);      
      }
      delete ftracklist;
    }
  }
  
  else if (fAnalysisInput.CompareTo("AOD")==0){
    //AliAODEvent *aod =  dynamic_cast<AliAODEvent*> (InputEvent());
    // to be implemented
    printf("  AOD analysis not yet implemented!!!\n\n");
    return;
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
  AliESDtrack* track;
  Double_t weight;
  
  int nt = tracklist->GetEntries();

  for (int i=0; i<nt; i++){
    weight = 1;
    track = dynamic_cast<AliESDtrack*> (tracklist->At(i));
    if (track) {
      weight = GetWeight(track);
      if (fSaveTrackContribution){
	EP->GetQContributionXArray()->AddAt(weight*cos(2*track->Phi()),track->GetID());
	EP->GetQContributionYArray()->AddAt(weight*sin(2*track->Phi()),track->GetID());
      }
      mQx += (weight*cos(2*track->Phi()));
      mQy += (weight*sin(2*track->Phi()));
    }
  }
  mQ.Set(mQx,mQy);
  return mQ;
}
  
  //________________________________________________________________________
void AliEPSelectionTask::GetQsub(TVector2 &Q1, TVector2 &Q2, TObjArray* tracklist)
{
// Get Qsub
  TVector2 mQ[2];
  float mQx1=0, mQy1=0, mQx2=0, mQy2=0;
  Double_t weight;

  AliESDtrack* track;
  TRandom2 rn = 0;
  
  int nt = tracklist->GetEntries();
  int trackcounter1=0, trackcounter2=0;
  
  for (Int_t i = 0; i < nt; i++) {
    weight = 1;
    track = dynamic_cast<AliESDtrack*> (tracklist->At(i));
    if (!track) continue;
    weight = GetWeight(track);
    
    // This loop splits the track set into 2 random subsets
    if( trackcounter1 < int(nt/2.) && trackcounter2 < int(nt/2.)){
      float random = rn.Rndm();
      if(random < .5){
        mQx1 += (weight*cos(2*track->Phi()));
        mQy1 += (weight*sin(2*track->Phi()));
        trackcounter1++;
      }
      else {
        mQx2 += (weight*cos(2*track->Phi()));
        mQy2 += (weight*sin(2*track->Phi()));
        trackcounter2++;
      }
    }
    else if( trackcounter1 >= int(nt/2.)){
      mQx2 += (weight*cos(2*track->Phi()));
      mQy2 += (weight*sin(2*track->Phi()));
      trackcounter2++;
    }
    else {
      mQx1 += (weight*cos(2*track->Phi()));
      mQy1 += (weight*sin(2*track->Phi()));
      trackcounter1++;
    }
  }
  mQ[0].Set(mQx1,mQy1);
  mQ[1].Set(mQx2,mQy2);
  Q1 = mQ[0];
  Q2 = mQ[1];
}

//________________________________________________________________________
void AliEPSelectionTask::SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts){
  
  fusercuts = kTRUE;
  fESDtrackCuts = trackcuts;
}

//_____________________________________________________________________________
void AliEPSelectionTask::SetTrackType(TString tracktype){
// Set the track type
  fTrackType = tracktype;
  if (!fusercuts) {
  if (fTrackType.CompareTo("GLOBAL")==0) fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  if (fTrackType.CompareTo("TPC")==0)    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fESDtrackCuts->SetPtRange(0.15,20);
  fESDtrackCuts->SetEtaRange(-0.8,0.8);
  }
}

//________________________________________________________________________
Double_t AliEPSelectionTask::GetWeight(AliESDtrack* track)
{
// Get weight for track
  Double_t ptweight=1;

  if (fUsePtWeight) {      
    if (track->Pt()<2) ptweight=track->Pt();
    else ptweight=2;
  }
  return ptweight*GetPhiWeight(track);
}

//________________________________________________________________________
Double_t AliEPSelectionTask::GetPhiWeight(AliESDtrack* track)
{
// Get phi weight for track
  Double_t phiweight=1;
  
  if (fUsePhiWeight && fPhiDist && track) {
    Double_t nParticles = fPhiDist->Integral();
    Double_t nPhibins = fPhiDist->GetNbinsX();
  
    Double_t phi = track->Phi();
    
    while (phi<0) phi += TMath::TwoPi();
    while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();
      
    Double_t phiDistValue = fPhiDist->GetBinContent(1+TMath::FloorNint((track->Phi())*nPhibins/TMath::TwoPi()));
    
    if (phiDistValue > 0) phiweight = nParticles/nPhibins/phiDistValue;
  }
  return phiweight;
}

//__________________________________________________________________________
void AliEPSelectionTask::SetPhiDist() 
{
// Set the phi distribution
  if(!fuserphidist) { // if it's already set and custom class is required, we use the one provided by the user

    fPhiDist = (TH1F*) fEPContainer->GetObject(frunNumber, "Default");
    if (!fPhiDist) AliFatal(Form("Cannot find OADB phi distribution for run %d", frunNumber));

  } 
  else {
    AliInfo("Using Custom Phi Distribution");
  }
    
  Bool_t emptybins;

  int iter = 0;  
  while (iter<3){
      emptybins = kFALSE;
   
      for (int i=1; i<fPhiDist->GetNbinsX(); i++){
	if (!((fPhiDist->GetBinContent(i))>0)) {
	  emptybins = kTRUE;
	}
      }  
      if (emptybins) {
	cout << "empty bins - rebinning!" << endl;
	fPhiDist->Rebin();
	iter++;
      }      
      else iter = 3;
  }
  
  if (emptybins) {
    AliError("After Maximum of rebinning still empty Phi-bins!!!");
  }
}

//__________________________________________________________________________
void AliEPSelectionTask::SetPersonalPhiDistribution(const char* infilename, char* listname)
{
    // Set a personal phi distribution
  fuserphidist = kTRUE;
  
  TFile f(infilename);
  TObject* list = f.Get(listname);
  fPhiDist = (TH1F*)list->FindObject("fHOutPhi");
  if (!fPhiDist) AliFatal("Phi Distribution not found!!!");

  f.Close();
} 
