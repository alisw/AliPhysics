/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//                  Base class for DStar Analysis
//
//
//  The D* spectra study is done in pt bins:
//  [0.7,1] [1,2] [2,3] [3,5] [5,8] [8,12],[12,18]
//
//  Optimized cuts used and TPC PID is on request (flag in che .C) 
//  Cuts option of analysis: 0 Heidelberg ; 1 Utrecht
//  Side Band and like sign background are implemented in the macro
//
//-----------------------------------------------------------------------
//
//                         Author A.Grelli 
//              ERC-QGP Utrecht University - a.grelli@uu.nl,
//                         Author Y.Wang
//        University of Heidelberg - yifei@physi.uni-heidelberg.de
//                         Author C.Ivan 
//             ERC-QGP Utrecht University - c.ivan@uu.nl,
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TH1I.h>
#include "TROOT.h"


#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"
//#include "AliAODPidHF.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
//#include "AliAODJet.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
//#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDStarSpectra.h"

ClassImp(AliAnalysisTaskSEDStarSpectra)


// I like pink

//__________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra():
  AliAnalysisTaskSE(),
  fEvents(0),
  fAnalysis(0),
  fD0Window(0),
  fPeakWindow(0),
  fUseMCInfo(kFALSE), 
  fOutput(0),
  fOutputSpectrum(0),
  fOutputAll(0),
  fOutputPID3(0),
  fOutputPID2(0),
  fOutputPID1(0),
  fNSigma(3),
  fPID(kTRUE),
  fCuts(0),
  fCEvents(0),     
  fTrueDiff2(0)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra(const Char_t* name, AliRDHFCutsDStartoKpipi* cuts) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fAnalysis(0),
  fD0Window(0),
  fPeakWindow(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputSpectrum(0),
  fOutputAll(0),
  fOutputPID3(0),
  fOutputPID2(0),
  fOutputPID1(0),
  fNSigma(3),
  fPID(kTRUE),
  fCuts(0),
  fCEvents(0),     
  fTrueDiff2(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarSpectra","Calling Constructor");

  fCuts=cuts;

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());  //Spectrum output
  DefineOutput(3,TList::Class());  //All Entries output
  DefineOutput(4,TList::Class());  //3sigma PID output
  DefineOutput(5,TList::Class());  //2sigma PID output
  DefineOutput(6,TList::Class());  //1sigma PID output
  DefineOutput(7,AliRDHFCutsDStartoKpipi::Class());  //My private output

}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::~AliAnalysisTaskSEDStarSpectra() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEDStarSpectra","Calling Destructor");
  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fOutputSpectrum) {
    delete fOutputSpectrum;
    fOutputSpectrum = 0;
  }
  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }
  if (fOutputPID3) {
    delete fOutputPID3;
    fOutputPID3 = 0;
  }
  if (fOutputPID2) {
    delete fOutputPID2;
    fOutputPID2 = 0;
  }
  if (fOutputPID1) {
    delete fOutputPID1;
    fOutputPID1 = 0;
  }
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  if(fCEvents){
    delete fCEvents;
    fCEvents =0;
  }
}
//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::Init(){
  //
  // Initialization
  //

  if(fDebug > 1) printf("AnalysisTaskSEDStarSpectra::Init() \n");
   AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
  // Post the data
  PostData(7,copyfCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }
  
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayDStartoD0pi=0;
 
  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
    }
  } else {
    arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
  }


  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
 

  fCEvents->Fill(1);
  // Load the event
  fEvents++;
  AliInfo(Form("Event %d",fEvents));
  if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));

  // counters for efficiencies
  Int_t icountReco = 0;
  
  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if(!vtx1) return;
  if(vtx1->GetNContributors()<1) return;

  if (!arrayDStartoD0pi){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast())); 

  // loop over the tracks to search for candidates soft pion
  
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {

    // D* candidates and D0 from D*
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    if(!dstarD0pi->GetSecondaryVtx()) continue;
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    if (!theD0particle) continue;
    
    Int_t isDStar = 0;
    
    TClonesArray *mcArray = 0; // fix coverity

    // mc analysis 
    if(fUseMCInfo){
    //MC array need for maching
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
	AliError("Could not find Monte-Carlo in AOD");
	return;
      }
      // find associated MC particle for D* ->D0toKpi
      Int_t mcLabel = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      if(mcLabel>=0) isDStar = 1;
    }

    // soft pion candidate
    AliAODTrack *track2 = (AliAODTrack*)dstarD0pi->GetBachelor(); 
    
    Double_t pt = dstarD0pi->Pt();
    Int_t isTkSelected = fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kTracks); // quality cuts on tracks
    if(!isTkSelected) continue;

    // cut in acceptance for the soft pion and for the D0 daughters  - TO BE REMOVED ONCE WILL BE IN THE CUTTING CLASS        
    Bool_t okTracks = SingleTrackSelections(theD0particle, track2);
    if (!okTracks) continue;
    
    Int_t ptbin=fCuts->PtBin(dstarD0pi->Pt());
    
    // set the D0 search window bin by bin
    if (ptbin==0){
      if(fAnalysis==1){
	fD0Window=0.015;
	fPeakWindow=0.0018;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (ptbin==1){
      if(fAnalysis==1){
	fD0Window=0.015;
	fPeakWindow=0.0018;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (ptbin==2){
      if(fAnalysis==1){
	fD0Window=0.018;
	fPeakWindow=0.0018;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (ptbin==3){
      if(fAnalysis==1){
	fD0Window=0.036;
	fPeakWindow=0.0018;
      }else{
	fD0Window=0.022;
	fPeakWindow=0.0016;
      }
    }
    if (ptbin==4){ 
      if(fAnalysis==1){
	fD0Window=0.036;
	fPeakWindow=0.0016;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.0014;
      }
    }
    if (ptbin>=5){
      if(fAnalysis==1){
	fD0Window=0.062;
	fPeakWindow=0.0014;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.0014;
      }
    } 

    FillSpectrum(dstarD0pi,isDStar,1,3,fCuts,fOutputPID3);
    FillSpectrum(dstarD0pi,isDStar,1,2,fCuts,fOutputPID2);
    FillSpectrum(dstarD0pi,isDStar,1,1,fCuts,fOutputPID1);
    FillSpectrum(dstarD0pi,isDStar,fPID,fNSigma,fCuts,fOutputSpectrum);
    //FillSpectrum(dstarD0pi,isDStar,0,0,fCuts,fOutputAll);

    SideBandBackground(dstarD0pi,1,3,fCuts,fOutputPID3);
    SideBandBackground(dstarD0pi,1,2,fCuts,fOutputPID2);
    SideBandBackground(dstarD0pi,1,1,fCuts,fOutputPID1);
    SideBandBackground(dstarD0pi,fPID,fNSigma,fCuts,fOutputSpectrum);
    //SideBandBackground(dstarD0pi,0,0,fCuts,fOutputAll);

    WrongSignForDStar(dstarD0pi,1,3,fCuts,fOutputPID3);
    WrongSignForDStar(dstarD0pi,1,2,fCuts,fOutputPID2);
    WrongSignForDStar(dstarD0pi,1,1,fCuts,fOutputPID1);
    WrongSignForDStar(dstarD0pi,fPID,fNSigma,fCuts,fOutputSpectrum);
    //WrongSignForDStar(dstarD0pi,0,0,fCuts,fOutputAll);
   
    if(isDStar == 1) { 
      fTrueDiff2->Fill(pt,dstarD0pi->DeltaInvMass());
    }
    
  }
  
  AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));
  
  PostData(1,fOutput);
  PostData(2,fOutputSpectrum);
  PostData(3,fOutputAll);
  PostData(4,fOutputPID3);
  PostData(5,fOutputPID2);
  PostData(6,fOutputPID1);


}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarSpectra::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  fTrueDiff2      = dynamic_cast<TH2F*>(fOutput->FindObject("fTrueDiff2"));

  fOutputSpectrum = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputSpectrum) {
    printf("ERROR: fOutputSpectrum not available\n");
    return;
  }
  fOutputAll = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputAll) {
    printf("ERROR: fOutputAll not available\n");
    return;
  }
  fOutputPID3 = dynamic_cast<TList*> (GetOutputData(4));
  if (!fOutputPID3) {
    printf("ERROR: fOutputPID3 not available\n");
    return;
  }
  fOutputPID2 = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputPID2) {
    printf("ERROR: fOutputPID2 not available\n");
    return;
  }
  fOutputPID1 = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputPID1) {
    printf("ERROR: fOutputPID1 not available\n");
    return;
  }

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserCreateOutputObjects() { 
 // output
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");

 
  fOutputSpectrum = new TList();
  fOutputSpectrum->SetOwner();
  fOutputSpectrum->SetName("listSpectrum");

  fOutputPID3 = new TList();
  fOutputPID3->SetOwner();
  fOutputPID3->SetName("listPID3");

  fOutputPID2 = new TList();
  fOutputPID2->SetOwner();
  fOutputPID2->SetName("listPID2");
    
  fOutputPID1 = new TList();
  fOutputPID1->SetOwner();
  fOutputPID1->SetName("listPID1");
    
  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("listAll");
 
  // define histograms
  DefineHistograms();

  PostData(1,fOutput);
  PostData(2,fOutputSpectrum);
  PostData(3,fOutputAll);
  PostData(4,fOutputPID3);
  PostData(5,fOutputPID2);
  PostData(6,fOutputPID1);

  return;
}
//___________________________________ hiostograms _______________________________________
void  AliAnalysisTaskSEDStarSpectra::DefineHistograms(){

  fCEvents = new TH1F("fCEvents","conter",10,0,10);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");
  fOutput->Add(fCEvents);

  fTrueDiff2 = new TH2F("DiffDstar_pt","True Reco diff vs pt",200,0,15,900,0,0.3);
  fOutput->Add(fTrueDiff2);

  const Int_t nhist=9;
  TString nameMass=" ", nameSgn=" ", nameBkg=" ";

  for(Int_t i=-2;i<nhist;i++){
    nameMass="histDeltaMass_";
    nameMass+=i+1;
    nameSgn="histDeltaSgn_";
    nameSgn+=i+1;
    nameBkg="histDeltaBkg_";
    nameBkg+=i+1; 
    
    if (i==-2) {
      nameMass="histDeltaMass";
      nameSgn="histDeltaSgn";
      nameBkg="histDeltaBkg";
    }
    
    TH1F* spectrumMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} invariant mass; #DeltaM [GeV/c^{2}]; Entries",200,0.1,0.2);
    TH1F* spectrumSgn = new TH1F(nameSgn.Data(), "D^{*}-D^{0} Signal invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries",200,0.1,0.2);
    TH1F* spectrumBkg = new TH1F(nameBkg.Data(), "D^{*}-D^{0} Background invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries",200,0.1,0.2);
    
    nameMass="histD0Mass_";
    nameMass+=i+1;
    nameSgn="histD0Sgn_";
    nameSgn+=i+1;
    nameBkg="histD0Bkg_";
    nameBkg+=i+1;
    
    if (i==-2) {
      nameMass="histD0Mass";
      nameSgn="histD0Sgn";
      nameBkg="histD0Bkg";
    }

    TH1F* spectrumD0Mass = new TH1F(nameMass.Data(),"D^{0} invariant mass; M(D^{0}) [GeV/c^{2}]; Entries",200,1.75,1.95);
    TH1F* spectrumD0Sgn = new TH1F(nameSgn.Data(), "D^{0} Signal invariant mass - MC; M(D^{0}) [GeV/c^{2}]; Entries",200,1.75,1.95);
    TH1F* spectrumD0Bkg = new TH1F(nameBkg.Data(), "D^{0} Background invariant mass - MC; M(D^{0}) [GeV/c^{2}]; Entries",200,1.75,1.95);

    nameMass="histDstarMass_";
    nameMass+=i+1;
    nameSgn="histDstarSgn_";
    nameSgn+=i+1;
    nameBkg="histDstarBkg_";
    nameBkg+=i+1;

    if (i==-2) {
      nameMass="histDstarMass";
      nameSgn="histDstarSgn";
      nameBkg="histDstarBkg";
    }

    TH1F* spectrumDstarMass = new TH1F(nameMass.Data(),"D^{*} invariant mass; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);
    TH1F* spectrumDstarSgn = new TH1F(nameSgn.Data(), "D^{*} Signal invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);
    TH1F* spectrumDstarBkg = new TH1F(nameBkg.Data(), "D^{*} Background invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);

    nameMass="histSideBandMass_";
    nameMass+=i+1;
    if (i==-2) { 
      nameMass="histSideBandMass";
    }
    
    TH1F* spectrumSideBandMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} sideband mass; M(D^{*}) [GeV/c^{2}]; Entries",200,0.1,0.2);

    nameMass="histWrongSignMass_";
    nameMass+=i+1;
    if (i==-2) { 
      nameMass="histWrongSignMass";
    }
    
    TH1F* spectrumWrongSignMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} wrongsign mass; M(D^{*}) [GeV/c^{2}]; Entries",200,0.1,0.2);


    spectrumMass->Sumw2();
    spectrumSgn->Sumw2();
    spectrumBkg->Sumw2();
    
    spectrumMass->SetLineColor(6);
    spectrumSgn->SetLineColor(2);
    spectrumBkg->SetLineColor(4);
    
    spectrumMass->SetMarkerStyle(20);
    spectrumSgn->SetMarkerStyle(20);
    spectrumBkg->SetMarkerStyle(20);
    spectrumMass->SetMarkerSize(0.6);
    spectrumSgn->SetMarkerSize(0.6);
    spectrumBkg->SetMarkerSize(0.6);
    spectrumMass->SetMarkerColor(6);
    spectrumSgn->SetMarkerColor(2);
    spectrumBkg->SetMarkerColor(4);

    spectrumD0Mass->Sumw2();
    spectrumD0Sgn->Sumw2();
    spectrumD0Bkg->Sumw2();

    spectrumD0Mass->SetLineColor(6);
    spectrumD0Sgn->SetLineColor(2);
    spectrumD0Bkg->SetLineColor(4);

    spectrumD0Mass->SetMarkerStyle(20);
    spectrumD0Sgn->SetMarkerStyle(20);
    spectrumD0Bkg->SetMarkerStyle(20);
    spectrumD0Mass->SetMarkerSize(0.6);
    spectrumD0Sgn->SetMarkerSize(0.6);
    spectrumD0Bkg->SetMarkerSize(0.6);
    spectrumD0Mass->SetMarkerColor(6);
    spectrumD0Sgn->SetMarkerColor(2);
    spectrumD0Bkg->SetMarkerColor(4);

    spectrumDstarMass->Sumw2();
    spectrumDstarSgn->Sumw2();
    spectrumDstarBkg->Sumw2();

    spectrumDstarMass->SetLineColor(6);
    spectrumDstarSgn->SetLineColor(2);
    spectrumDstarBkg->SetLineColor(4);

    spectrumDstarMass->SetMarkerStyle(20);
    spectrumDstarSgn->SetMarkerStyle(20);
    spectrumDstarBkg->SetMarkerStyle(20);
    spectrumDstarMass->SetMarkerSize(0.6);
    spectrumDstarSgn->SetMarkerSize(0.6);
    spectrumDstarBkg->SetMarkerSize(0.6);
    spectrumDstarMass->SetMarkerColor(6);
    spectrumDstarSgn->SetMarkerColor(2);
    spectrumDstarBkg->SetMarkerColor(4);

    spectrumSideBandMass->Sumw2();
    spectrumSideBandMass->SetLineColor(4);
    spectrumSideBandMass->SetMarkerStyle(20);
    spectrumSideBandMass->SetMarkerSize(0.6);
    spectrumSideBandMass->SetMarkerColor(4);

    spectrumWrongSignMass->Sumw2();
    spectrumWrongSignMass->SetLineColor(4);
    spectrumWrongSignMass->SetMarkerStyle(20);
    spectrumWrongSignMass->SetMarkerSize(0.6);
    spectrumWrongSignMass->SetMarkerColor(4);

    TH1F* allMass = (TH1F*)spectrumMass->Clone();
    TH1F* allSgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* allBkg  = (TH1F*)spectrumBkg->Clone();

    TH1F* pid3Mass = (TH1F*)spectrumMass->Clone();
    TH1F* pid3Sgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* pid3Bkg  = (TH1F*)spectrumBkg->Clone();

    TH1F* pid2Mass = (TH1F*)spectrumMass->Clone();
    TH1F* pid2Sgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* pid2Bkg  = (TH1F*)spectrumBkg->Clone();

    TH1F* pid1Mass = (TH1F*)spectrumMass->Clone();
    TH1F* pid1Sgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* pid1Bkg  = (TH1F*)spectrumBkg->Clone();

    fOutputSpectrum->Add(spectrumMass);
    fOutputSpectrum->Add(spectrumSgn);
    fOutputSpectrum->Add(spectrumBkg);

    fOutputAll->Add(allMass);
    fOutputAll->Add(allSgn);
    fOutputAll->Add(allBkg);

    fOutputPID3->Add(pid3Mass);
    fOutputPID3->Add(pid3Sgn);
    fOutputPID3->Add(pid3Bkg);

    fOutputPID2->Add(pid2Mass);
    fOutputPID2->Add(pid2Sgn);
    fOutputPID2->Add(pid2Bkg);

    fOutputPID1->Add(pid1Mass);
    fOutputPID1->Add(pid1Sgn);
    fOutputPID1->Add(pid1Bkg);

    TH1F* allD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* allD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* allD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    TH1F* pid3D0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* pid3D0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* pid3D0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    TH1F* pid2D0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* pid2D0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* pid2D0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    TH1F* pid1D0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* pid1D0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* pid1D0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    fOutputSpectrum->Add(spectrumD0Mass);
    fOutputSpectrum->Add(spectrumD0Sgn);
    fOutputSpectrum->Add(spectrumD0Bkg);

    fOutputAll->Add(allD0Mass);
    fOutputAll->Add(allD0Sgn);
    fOutputAll->Add(allD0Bkg);

    fOutputPID3->Add(pid3D0Mass);
    fOutputPID3->Add(pid3D0Sgn);
    fOutputPID3->Add(pid3D0Bkg);

    fOutputPID2->Add(pid2D0Mass);
    fOutputPID2->Add(pid2D0Sgn);
    fOutputPID2->Add(pid2D0Bkg);

    fOutputPID1->Add(pid1D0Mass);
    fOutputPID1->Add(pid1D0Sgn);
    fOutputPID1->Add(pid1D0Bkg);
  
    TH1F* allDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* allDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* allDstarBkg = (TH1F*)spectrumDstarBkg->Clone();
    
    TH1F* pid3DstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* pid3DstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* pid3DstarBkg = (TH1F*)spectrumDstarBkg->Clone();
    
    TH1F* pid2DstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* pid2DstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* pid2DstarBkg = (TH1F*)spectrumDstarBkg->Clone();

    TH1F* pid1DstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* pid1DstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* pid1DstarBkg = (TH1F*)spectrumDstarBkg->Clone();

    fOutputSpectrum->Add(spectrumDstarMass);
    fOutputSpectrum->Add(spectrumDstarSgn);
    fOutputSpectrum->Add(spectrumDstarBkg);

    fOutputAll->Add(allDstarMass);
    fOutputAll->Add(allDstarSgn);
    fOutputAll->Add(allDstarBkg);

    fOutputPID3->Add(pid3DstarMass);
    fOutputPID3->Add(pid3DstarSgn);
    fOutputPID3->Add(pid3DstarBkg);

    fOutputPID2->Add(pid2DstarMass);
    fOutputPID2->Add(pid2DstarSgn);
    fOutputPID2->Add(pid2DstarBkg);

    fOutputPID1->Add(pid1DstarMass);
    fOutputPID1->Add(pid1DstarSgn);
    fOutputPID1->Add(pid1DstarBkg);

    TH1F* allSideBandMass = (TH1F*)spectrumSideBandMass->Clone();
    TH1F* pid3SideBandMass = (TH1F*)spectrumSideBandMass->Clone();
    TH1F* pid2SideBandMass = (TH1F*)spectrumSideBandMass->Clone();
    TH1F* pid1SideBandMass = (TH1F*)spectrumSideBandMass->Clone();

    fOutputSpectrum->Add(spectrumSideBandMass);
    fOutputAll->Add(allSideBandMass);
    fOutputPID3->Add(pid3SideBandMass);
    fOutputPID2->Add(pid2SideBandMass);
    fOutputPID1->Add(pid1SideBandMass);

    TH1F* allWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();
    TH1F* pid3WrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();
    TH1F* pid2WrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();
    TH1F* pid1WrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();

    fOutputSpectrum->Add(spectrumWrongSignMass);
    fOutputAll->Add(allWrongSignMass);
    fOutputPID3->Add(pid3WrongSignMass);
    fOutputPID2->Add(pid2WrongSignMass);
    fOutputPID1->Add(pid1WrongSignMass);
    
  }
  
  // pt spectra
  nameMass="ptMass";
  nameSgn="ptSgn";
  nameBkg="ptBkg";
  
  TH1F* ptspectrumMass = new TH1F(nameMass.Data(),"D^{*} p_{T}; p_{T} [GeV]; Entries",200,0,10);
  TH1F* ptspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal p_{T} - MC; p_{T} [GeV]; Entries",200,0,10);
  TH1F* ptspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background p_{T} - MC; p_{T} [GeV]; Entries",200,0,10);
  
  ptspectrumMass->Sumw2();
  ptspectrumSgn->Sumw2();
  ptspectrumBkg->Sumw2();
  
  ptspectrumMass->SetLineColor(6);
  ptspectrumSgn->SetLineColor(2);
  ptspectrumBkg->SetLineColor(4);
  
  ptspectrumMass->SetMarkerStyle(20);
  ptspectrumSgn->SetMarkerStyle(20);
  ptspectrumBkg->SetMarkerStyle(20);
  ptspectrumMass->SetMarkerSize(0.6);
  ptspectrumSgn->SetMarkerSize(0.6);
  ptspectrumBkg->SetMarkerSize(0.6);
  ptspectrumMass->SetMarkerColor(6);
  ptspectrumSgn->SetMarkerColor(2);
  ptspectrumBkg->SetMarkerColor(4);
  
  TH1F* ptallMass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptallSgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptallBkg = (TH1F*)ptspectrumBkg->Clone();
  
  TH1F* ptpid3Mass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptpid3Sgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptpid3Bkg = (TH1F*)ptspectrumBkg->Clone();
  
  TH1F* ptpid2Mass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptpid2Sgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptpid2Bkg = (TH1F*)ptspectrumBkg->Clone();
  
  TH1F* ptpid1Mass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptpid1Sgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptpid1Bkg = (TH1F*)ptspectrumBkg->Clone();
  
  fOutputSpectrum->Add(ptspectrumMass);
  fOutputSpectrum->Add(ptspectrumSgn);
  fOutputSpectrum->Add(ptspectrumBkg);
  
  fOutputAll->Add(ptallMass);
  fOutputAll->Add(ptallSgn);
  fOutputAll->Add(ptallBkg);
  
  fOutputPID3->Add(ptpid3Mass);
  fOutputPID3->Add(ptpid3Sgn);
  fOutputPID3->Add(ptpid3Bkg);
  
  fOutputPID2->Add(ptpid2Mass);
  fOutputPID2->Add(ptpid2Sgn);
  fOutputPID2->Add(ptpid2Bkg);
  
  fOutputPID1->Add(ptpid1Mass);
  fOutputPID1->Add(ptpid1Sgn);
  fOutputPID1->Add(ptpid1Bkg);
  
  // eta spectra
  nameMass="etaMass";
  nameSgn="etaSgn";
  nameBkg="etaBkg";
  
  TH1F* etaspectrumMass = new TH1F(nameMass.Data(),"D^{*} #eta; #eta; Entries",200,-1,1);
  TH1F* etaspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal #eta - MC; #eta; Entries",200,-1,1);
  TH1F* etaspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background #eta - MC; #eta; Entries",200,-1,1);
  
  etaspectrumMass->Sumw2();
  etaspectrumSgn->Sumw2();
  etaspectrumBkg->Sumw2();
  
  etaspectrumMass->SetLineColor(6);
  etaspectrumSgn->SetLineColor(2);
  etaspectrumBkg->SetLineColor(4);
  
  etaspectrumMass->SetMarkerStyle(20);
  etaspectrumSgn->SetMarkerStyle(20);
  etaspectrumBkg->SetMarkerStyle(20);
  etaspectrumMass->SetMarkerSize(0.6);
  etaspectrumSgn->SetMarkerSize(0.6);
  etaspectrumBkg->SetMarkerSize(0.6);
  etaspectrumMass->SetMarkerColor(6);
  etaspectrumSgn->SetMarkerColor(2);
  etaspectrumBkg->SetMarkerColor(4);
  
  TH1F* etaallMass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etaallSgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etaallBkg = (TH1F*)etaspectrumBkg->Clone();
  
  TH1F* etapid3Mass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etapid3Sgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etapid3Bkg = (TH1F*)etaspectrumBkg->Clone();
  
  TH1F* etapid2Mass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etapid2Sgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etapid2Bkg = (TH1F*)etaspectrumBkg->Clone();
  
  TH1F* etapid1Mass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etapid1Sgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etapid1Bkg = (TH1F*)etaspectrumBkg->Clone();
  
  fOutputSpectrum->Add(etaspectrumMass);
  fOutputSpectrum->Add(etaspectrumSgn);
  fOutputSpectrum->Add(etaspectrumBkg);
  
  fOutputAll->Add(etaallMass);
  fOutputAll->Add(etaallSgn);
  fOutputAll->Add(etaallBkg);
  
  fOutputPID3->Add(etapid3Mass);
  fOutputPID3->Add(etapid3Sgn);
  fOutputPID3->Add(etapid3Bkg);
  
  fOutputPID2->Add(etapid2Mass);
  fOutputPID2->Add(etapid2Sgn);
  fOutputPID2->Add(etapid2Bkg);
  
  fOutputPID1->Add(etapid1Mass);
  fOutputPID1->Add(etapid1Sgn);
  fOutputPID1->Add(etapid1Bkg);
  
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, Bool_t PIDon, Int_t nSigma, AliRDHFCutsDStartoKpipi *cuts, TList *listout){
  //
  // Fill histos for D* spectrum
  //

  Int_t ptbin=cuts->PtBin(part->Pt());

  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kCandidate); //selected
  if (!isSelected){
    return;
  }

  Double_t invmassD0   = part->InvMassD0();  
  if (TMath::Abs(invmassD0-1.865)>fD0Window) return;
  
  Double_t pt = part->Pt();
  Double_t eta = part->Eta();
  
  AliAODTrack *softPi = (AliAODTrack*)part->GetBachelor();
  
  //PID of D0 daughters
  AliAODTrack *pos = (AliAODTrack*)part->Get2Prong()->GetDaughter(0);
  AliAODTrack *neg = (AliAODTrack*)part->Get2Prong()->GetDaughter(1);
  
  Bool_t isPID = kTRUE;
  
  //  Double_t cutsN = -1;
    
  //if(part->Charge()<0){
  //  cutsN = (pos->Pt()-neg->Pt())/(part->Get2Prong()->Pt());
  // }
  //if(part->Charge()>0){
  //  cutsN = (neg->Pt()-pos->Pt())/(part->Get2Prong()->Pt());
  // }
 
  //  if(ptbin==1 || ptbin==2){
  //  if(cutsN>0.7) return;
  // }
  //if(ptbin==3){
  //  if(cutsN>0.8) return;
  // }


  if(fAnalysis==1 && ptbin==1){
    if(part->Get2Prong()->DecayLength()>0.1 || part->Get2Prong()->DecayLength()<0.015)return; //(0.05)
   }
  if(fAnalysis==1 && ptbin==2){
    if(part->Get2Prong()->DecayLength()>0.1 || part->Get2Prong()->DecayLength()<0.015)return;
  }
 
  if (PIDon) {
    if(fDebug > 1) printf("AnalysisTaskSEDStar::TPCPIDon \n");
    if(fDebug > 1) printf("AnalysisTaskSEDStar::NSigmaTPC: %d\n", nSigma);
    
    if (part->Charge()>0){
      if(!SelectPID(pos, AliPID::kPion, nSigma)) return;//pion+
      if(!SelectPID(neg, AliPID::kKaon, nSigma)) return;//kaon-
    }else{
      if(!SelectPID(pos, AliPID::kKaon, nSigma)) return;//kaon+
      if(!SelectPID(neg, AliPID::kPion, nSigma)) return;//pion-
    }
    if (ptbin>1){
      isPID =SelectTOFPID(part->Get2Prong(), softPi);
      if(!isPID) return;
    }
  }
  
  Double_t invmassDelta = part->DeltaInvMass();
  Double_t invmassDstar = part->InvMassDstarKpipi();
  
  TString fillthis="";
  Bool_t massInRange=kFALSE;
  if (TMath::Abs(invmassDelta-0.14557)<fPeakWindow) massInRange=kTRUE;
  
  
  if(fUseMCInfo) {
    if(isDStar==1) {
      fillthis="histD0Sgn_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histD0Sgn";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histDstarSgn_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDstarSgn";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDeltaSgn_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      fillthis="histDeltaSgn";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      //if (massInRange) {
      if(ptbin<=1){
	fillthis="ptSgn";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(pt);
	fillthis="etaSgn";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(eta);
      }
    }
    else {//background
      fillthis="histD0Bkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histD0Bkg";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histDstarBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDstarBkg";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDeltaBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      fillthis="histDeltaBkg";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      //if (massInRange) {
      if(ptbin<=1){
	fillthis="ptBkg";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(pt);
	fillthis="etaBkg";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(eta);
      }
    }
  }
  //no MC info, just cut selection
  fillthis="histD0Mass_";
  fillthis+=ptbin;
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
  fillthis="histD0Mass";
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
  fillthis="histDstarMass_";
  fillthis+=ptbin;
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
  fillthis="histDstarMass";
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
  fillthis="histDeltaMass_";
  fillthis+=ptbin;
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
  fillthis="histDeltaMass";
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
  
  if (massInRange) {
    fillthis="ptMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(pt);
    fillthis="etaMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(eta);
  }
 
  return;
}
//______________________________ side band background for D*___________________________________
void AliAnalysisTaskSEDStarSpectra::SideBandBackground(AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma,  AliRDHFCutsDStartoKpipi *cuts, TList *listout){

  //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
  // (expected detector resolution) on the left and right frm the D0 mass. Each band
  //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   

  Int_t ptbin=cuts->PtBin(part->Pt());
  
  Bool_t massInRange=kFALSE;

  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kCandidate); //selected
  if (!isSelected){
    return;
  }
  //  Double_t pt = part->Pt();
  
  AliAODTrack *softPi = (AliAODTrack*)part->GetBachelor();

  // select the side bands intervall
  Double_t invmassD0    = part->InvMassD0();
  if(TMath::Abs(invmassD0-1.865)>4*fD0Window && TMath::Abs(invmassD0-1.865)<8*fD0Window){
    
    // for pt and eta
    Double_t invmassDelta = part->DeltaInvMass();
    if (TMath::Abs(invmassDelta-0.14557)<fPeakWindow) massInRange=kTRUE;
    
    //PID of D0 daughters
    AliAODTrack *pos = (AliAODTrack*)part->Get2Prong()->GetDaughter(0);
    AliAODTrack *neg = (AliAODTrack*)part->Get2Prong()->GetDaughter(1);
    
    //Double_t cutsN = -1;
    

    /* if(part->Charge()<0){
      cutsN = (pos->Pt()-neg->Pt())/(part->Get2Prong()->Pt());
    }
    if(part->Charge()>0){
      cutsN = (neg->Pt()-pos->Pt())/(part->Get2Prong()->Pt());
    }
    
    if(ptbin==1 || ptbin==2){
      if(cutsN>0.5) return;
    }
    if(ptbin==3){
      if(cutsN>0.7) return;
    }
    */
    if(fAnalysis==1 && ptbin==1){
      if(part->Get2Prong()->DecayLength()>0.08 || part->Get2Prong()->DecayLength()<0.022)return; //(0.05)
    }
    if(fAnalysis==1 && ptbin==2){
      if(part->Get2Prong()->DecayLength()>0.08 || part->Get2Prong()->DecayLength()<0.017)return;
    }
    
    Bool_t isPID = kTRUE;

    if (PIDon) {
      if(fDebug > 1) printf("AnalysisTaskSEDStar::TPCPIDon \n");
      if(fDebug > 1) printf("AnalysisTaskSEDStar::NSigmaTPC: %d\n", nSigma);
      
      if (part->Charge()>0){
	if(!SelectPID(pos, AliPID::kPion, nSigma)) return;//pion+
	if(!SelectPID(neg, AliPID::kKaon, nSigma)) return;//kaon-
      }else{
	if(!SelectPID(pos, AliPID::kKaon, nSigma)) return;//kaon+
	if(!SelectPID(neg, AliPID::kPion, nSigma)) return;//pion-
      }
      if (ptbin>2){
	isPID =SelectTOFPID(part->Get2Prong(), softPi);
	if(!isPID) return;
      }
    }
    
    TString fillthis="";
    fillthis="histSideBandMass_";
    fillthis+=ptbin;
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
    fillthis="histSideBandMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
    
  }
}
//________________________________________________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::WrongSignForDStar(AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma,  AliRDHFCutsDStartoKpipi *cuts, TList *listout){
  //
  // assign the wrong charge to the soft pion to create background
  //
  Int_t ptbin=cuts->PtBin(part->Pt());
  
  AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)part->Get2Prong();
  AliAODTrack *theSoftPi = (AliAODTrack*)part->GetBachelor(); 

  Int_t okD0WrongSign,okD0barWrongSign;
  Double_t wrongMassD0=0.;
  
  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kCandidate); //selected
   if (!isSelected){
    return;
  }

  okD0WrongSign =  1;
  okD0barWrongSign = 1;
  
  //if is D*+ than assume D0bar
  if(part->Charge()>0 && (isSelected ==1)) { 
    okD0WrongSign = 0;
  }
  if(part->Charge()<0 && (isSelected ==2)){
    okD0barWrongSign = 0;
  }
  
  // assign the wrong mass in case the cuts return both D0 and D0bar
  if(part->Charge()>0 && (isSelected ==3)) { 
    okD0WrongSign = 0;
  } else if(part->Charge()<0 && (isSelected ==3)){
    okD0barWrongSign = 0;
  }
  
  //wrong D0 inv mass
  if(okD0WrongSign!=0){
    wrongMassD0 = theD0particle->InvMassD0();
  }else if(okD0WrongSign==0){
    wrongMassD0 = theD0particle->InvMassD0bar();
  }
  
  if(TMath::Abs(wrongMassD0-1.865)<fD0Window){
    
    //PID of D0 daughters
    AliAODTrack *pos = (AliAODTrack*)part->Get2Prong()->GetDaughter(0);
    AliAODTrack *neg = (AliAODTrack*)part->Get2Prong()->GetDaughter(1);
    
    Bool_t isPID = kTRUE;

    
    if(fAnalysis==1 && ptbin==1){
      if(part->Get2Prong()->DecayLength()>0.08 || part->Get2Prong()->DecayLength()<0.022)return; //(0.05)
    }
    if(fAnalysis==1 && ptbin==2){
      if(part->Get2Prong()->DecayLength()>0.08 || part->Get2Prong()->DecayLength()<0.017)return;
    }

    if (PIDon) {
      if(fDebug > 1) printf("AnalysisTaskSEDStar::TPCPIDon \n");
      if(fDebug > 1) printf("AnalysisTaskSEDStar::NSigmaTPC: %d\n", nSigma);
      
      if (part->Charge()>0){
	if(!SelectPID(pos, AliPID::kPion, nSigma)) return;//pion+
	if(!SelectPID(neg, AliPID::kKaon, nSigma)) return;//kaon-
      }else{
	if(!SelectPID(pos, AliPID::kKaon, nSigma)) return;//kaon+
	if(!SelectPID(neg, AliPID::kPion, nSigma)) return;//pion-
      }
      if (ptbin>2){
	isPID =SelectTOFPID(theD0particle, theSoftPi);
	if(!isPID) return;
      }
    }
    
    // wrong D* inv mass   
    Double_t e[3];
    if (part->Charge()>0){
      e[0]=theD0particle->EProng(0,321);
      e[1]=theD0particle->EProng(1,211);
    }else{
      e[0]=theD0particle->EProng(0,211);
      e[1]=theD0particle->EProng(1,321);
    }
    e[2]=part->EProng(0,211);
    
    Double_t esum = e[0]+e[1]+e[2];
    Double_t pds = part->P();

    Double_t   wrongMassDstar = TMath::Sqrt(esum*esum-pds*pds);

    TString fillthis="";
    fillthis="histWrongSignMass_";
    fillthis+=ptbin;
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar-wrongMassD0);
    fillthis="histWrongSignMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar-wrongMassD0);
    
  }
}

//_____________________________________________SINGLE TRACK PRE-SELECTION___________________________________________
Bool_t AliAnalysisTaskSEDStarSpectra::SingleTrackSelections(const AliAODRecoDecayHF2Prong* theD0particle, const AliAODTrack *track2){
  
  // Preselection on D0 daughters and the soft pion
  
  // cut in acceptance for the soft pion and for the D0 daughters      
  Bool_t acceptanceProng0 = (TMath::Abs(theD0particle->EtaProng(0))<= 0.9);
  Bool_t acceptanceProng1 = (TMath::Abs(theD0particle->EtaProng(1))<= 0.9);
  // soft pion acceptance ... is it fine 0.9?????
  Bool_t acceptanceProng2 = (TMath::Abs(track2->Eta())<= 1.0);
  
  if (!(acceptanceProng0 && acceptanceProng1 && acceptanceProng2)) return kFALSE;
  AliDebug(2,"D* reco daughters in acceptance");
  
  return kTRUE;
}

//_____________________________ pid _______________________________________-
Bool_t AliAnalysisTaskSEDStarSpectra::SelectPID(const AliAODTrack *track,  AliPID::EParticleType type, Double_t nsig){//type(0-4): {e,mu,pi,K,p}
  //
  // Method to extract the PID for the pion/kaon. The particle type for PID can be set by user
  // At the moment only TPC PID.
  //

  //TPC

  Bool_t isParticle=kTRUE;
  if ((track->GetStatus()&AliESDtrack::kTPCpid )==0) return isParticle;
  AliAODPid *pid = track->GetDetPid();
  static AliTPCPIDResponse theTPCpid;
  Double_t nsigma = theTPCpid.GetNumberOfSigmas(track->P(),pid->GetTPCsignal(),track->GetTPCClusterMap().CountBits(), type);
  if (TMath::Abs(nsigma)>nsig) isParticle=kFALSE;

  return isParticle;
}
//-------------------------------------------------
Bool_t AliAnalysisTaskSEDStarSpectra::SelectTOFPID(const AliAODRecoDecayHF2Prong* d, const AliAODTrack *tracksoft){
  
  
  // ######### SPECIAL PID CUTS #################################
  Int_t isKaon[2]={0,0};
  Bool_t isD0D0barPID[2]={kTRUE,kTRUE};
  for(Int_t daught=0;daught<2;daught++){
    
    
    AliAODTrack *aodtrack=(AliAODTrack*)d->GetDaughter(daught);
    // AliESDtrack *esdtrack=new
    AliESDtrack((AliVTrack*)d->GetDaughter(daught)); 
    if(!(aodtrack->GetStatus()&AliESDtrack::kTOFpid)){
      isKaon[daught]=0;
      //delete esdtrack;
      return kTRUE;
    }
    if(!(aodtrack->GetStatus()&AliESDtrack::kTOFout)){
      isKaon[daught]=0;
      //	delete esdtrack;
      return kTRUE;
    } 
    if(!(aodtrack->GetStatus()&AliESDtrack::kTIME)){
      isKaon[daught]=0;
      //delete esdtrack;
      return kTRUE;
    }
    if(!(aodtrack->GetStatus()&AliESDtrack::kTPCrefit)){
      isKaon[daught]=0;
      //	delete esdtrack;
      return kTRUE;
    } 
    if(!(aodtrack->GetStatus()&AliESDtrack::kITSrefit)){
      isKaon[daught]=0;
      //	delete esdtrack;
      return kTRUE;
    } 
    
    AliAODPid *pid=aodtrack->GetDetPid();
    if(!pid) {
      isKaon[daught]=0;
      //delete esdtrack;
      return kTRUE;
    }
    Double_t tofSig=pid->GetTOFsignal(); 
    
    Double_t times[5];
    //      esdtrack->GetIntegratedTimes(times);
    pid->GetIntegratedTimes(times);
    //fHistCheck->Fill(esdtrack->P(),esdtrack->GetTOFsignal()-times[3]); //
    //3 is the kaon

    //  printf("Test momentum VS time %f, %f \n",aodtrack->P(),tofSig-times[3]);
    if(TMath::Abs(tofSig-times[3])>3.*160.){
      isKaon[daught]=2;
      if(aodtrack->Charge()==-1){
	isD0D0barPID[0]=kFALSE;
      }
      else isD0D0barPID[1]=kFALSE;
    }
    else {
      isKaon[daught]=1;
      
      if(aodtrack->P()<1.5){
	if(aodtrack->Charge()==-1){
	  isD0D0barPID[1]=kFALSE;
	}
	else isD0D0barPID[0]=kFALSE;
	
      }
      //delete esdtrack;
    }
  }
  
  Double_t psCharge = tracksoft->Charge();
  
  if(psCharge>0 &&  !isD0D0barPID[0]) return kFALSE;
  if(psCharge<0 &&  !isD0D0barPID[1]) return kFALSE;
  
  return kTRUE;
}
