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
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//
//                  Base class for DStar Analysis
//
//
//  The D* spectra study is done in pt bins:
//  [0,1] [1,2] [2,3] [3,5] [5,8] [8,14]
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
#include "AliAODJet.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSEDStarSpectra.h"

ClassImp(AliAnalysisTaskSEDStarSpectra)

//__________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra():
  AliAnalysisTaskSE(),
  fEvents(0),
  fAnalysis(0),
  fVHF(0),  
  fVHFloose(0),
  fD0Window(0),
  fPeakWindow(0),
  fMinITSClusters(0),
  fMinITSClustersSoft(0),
  fUseMCInfo(kTRUE), 
  fOutput(0),
  fOutputSpectrum(0),
  fOutputAll(0),
  fOutputPID3(0),
  fOutputPID2(0),
  fOutputPID1(0),
  fNSigma(3),
  fPID(kTRUE),
  fAODTrack(0),
  fCEvents(0),     
  fTrueDiff2(0)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fAnalysis(0),
  fVHF(0),  
  fVHFloose(0),
  fD0Window(0),
  fPeakWindow(0),
  fMinITSClusters(0),
  fMinITSClustersSoft(0),
  fUseMCInfo(kTRUE),
  fOutput(0),
  fOutputSpectrum(0),
  fOutputAll(0),
  fOutputPID3(0),
  fOutputPID2(0),
  fOutputPID1(0),
  fNSigma(3),
  fPID(kTRUE),
  fAODTrack(0),
  fCEvents(0),     
  fTrueDiff2(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarSpectra","Calling Constructor");
  
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());  //Spectrum output
  DefineOutput(3,TList::Class());  //3sigma PID output
  DefineOutput(4,TList::Class());  //2sigma PID output
  DefineOutput(5,TList::Class());  //1sigma PID output
  DefineOutput(6,TList::Class());  //All Entries output
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra& AliAnalysisTaskSEDStarSpectra::operator=(const AliAnalysisTaskSEDStarSpectra& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra(const AliAnalysisTaskSEDStarSpectra& c) :
  AliAnalysisTaskSE(c),
  fEvents(c.fEvents),
  fAnalysis(c.fAnalysis),
  fVHF(c.fVHF),  
  fVHFloose(c.fVHFloose),
  fD0Window(c.fD0Window),
  fPeakWindow(c.fPeakWindow),
  fMinITSClusters(c.fMinITSClusters),
  fMinITSClustersSoft(c.fMinITSClustersSoft),
  fUseMCInfo(c.fUseMCInfo),
  fOutput(c.fOutput),
  fOutputSpectrum(c.fOutputSpectrum),
  fOutputAll(c.fOutputAll),
  fOutputPID3(c.fOutputPID3),
  fOutputPID2(c.fOutputPID2),
  fOutputPID1(c.fOutputPID1),
  fNSigma(c.fNSigma),
  fPID(c.fPID),
  fAODTrack(c.fAODTrack),
  fCEvents(c.fCEvents),     
  fTrueDiff2(c.fTrueDiff2)
{
  //
  // Copy Constructor
  //
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
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  } 
  if (fVHFloose) {
    delete fVHFloose;
    fVHFloose = 0;
  }

}
//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::Init(){
  //
  // Initialization
  //

  if(fDebug > 1) printf("AnalysisTaskSEDStarSpectra::Init() \n");

  //gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fVHFloose = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fVHFloose->SetD0fromDstarCuts(0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.);
  fVHFloose->SetDstarCuts(0.3, 0.1, 0.05, 100000000000.0, 0.5);
  //fVHF->PrintStatus();

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
  
  fCEvents->Fill(1);
  // Load the event
  fEvents++;
  AliInfo(Form("Event %d",fEvents));
  if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayDStartoD0pi=0;
  //  Init();
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
 
 // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();

  // counters for efficiencies
  Int_t icountReco = 0;
  
  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  if (!arrayDStartoD0pi){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast())); 
    
  // loop over the tracks to search for candidates soft pion
  
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {
    
    // D* candidates
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    
    // D0 from the reco cascade
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    Bool_t unsetvtx=kFALSE;
    
    // needed for pointing angle
    if(!theD0particle->GetOwnPrimaryVtx()) {
      theD0particle->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }    

    Int_t isDStar = 0;

    // mc analysis 
    if(fUseMCInfo){
    //MC array need for maching
      TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
      // find associated MC particle for D* ->D0toKpi
      Int_t mcLabel = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      if(mcLabel>=0) isDStar = 1;
    }

    // soft pion
    AliAODTrack *track2 = (AliAODTrack*)dstarD0pi->GetBachelor(); 

    //D0tokpi
    AliAODTrack *track0 = (AliAODTrack*)theD0particle->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)theD0particle->GetDaughter(1);
    
    Double_t pt = dstarD0pi->Pt();
    
       
    // cut in acceptance for the soft pion and for the D0 daughters      

    Bool_t okTracks = SingleTrackSelections(theD0particle, track0, track1, track2);
    if (!okTracks) continue;
 
    // D0 pt needed for the cuts          
    Int_t ptbin =0;
    if (pt>0. && pt<=1.) ptbin =0;
    if (pt>1. && pt<=2.) ptbin =1;
    if (pt>2. && pt<=3.) ptbin =2;
    if (pt>3. && pt<=5.) ptbin =3;
    if (pt>5. && pt<=8.) ptbin =4;
    if (pt>8.) 		 ptbin =5;

    SetSelections(pt);  
    FillSpectrum(ptbin,dstarD0pi,isDStar,1,3,fVHF,fOutputPID3);
    FillSpectrum(ptbin,dstarD0pi,isDStar,1,2,fVHF,fOutputPID2);
    FillSpectrum(ptbin,dstarD0pi,isDStar,1,1,fVHF,fOutputPID1);
    FillSpectrum(ptbin,dstarD0pi,isDStar,fPID,fNSigma,fVHF,fOutputSpectrum);
    FillSpectrum(ptbin,dstarD0pi,isDStar,0,0,fVHFloose,fOutputAll);

    SideBandBackground(ptbin,dstarD0pi,1,3,fVHF,fOutputPID3);
    SideBandBackground(ptbin,dstarD0pi,1,2,fVHF,fOutputPID2);
    SideBandBackground(ptbin,dstarD0pi,1,1,fVHF,fOutputPID1);
    SideBandBackground(ptbin,dstarD0pi,fPID,fNSigma,fVHF,fOutputSpectrum);
    SideBandBackground(ptbin,dstarD0pi,0,0,fVHFloose,fOutputAll);

    WrongSignForDStar(ptbin,dstarD0pi,1,3,fVHF,fOutputPID3);
    WrongSignForDStar(ptbin,dstarD0pi,1,2,fVHF,fOutputPID2);
    WrongSignForDStar(ptbin,dstarD0pi,1,1,fVHF,fOutputPID1);
    WrongSignForDStar(ptbin,dstarD0pi,fPID,fNSigma,fVHF,fOutputSpectrum);
    WrongSignForDStar(ptbin,dstarD0pi,0,0,fVHFloose,fOutputAll);
   
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

}
//___________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserCreateOutputObjects() { 
 // output
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
 
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

  const Int_t nhist=5;
  TString nameMass=" ", nameSgn=" ", nameBkg=" ";

  for(Int_t i=-1;i<nhist;i++){
    nameMass="histDeltaMass_";
    nameMass+=i+1;
    nameSgn="histDeltaSgn_";
    nameSgn+=i+1;
    nameBkg="histDeltaBkg_";
    nameBkg+=i+1; 

    if (i==-1) {
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

    if (i==-1) {
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

    if (i==-1) {
      nameMass="histDstarMass";
      nameSgn="histDstarSgn";
      nameBkg="histDstarBkg";
    }

    TH1F* spectrumDstarMass = new TH1F(nameMass.Data(),"D^{*} invariant mass; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);
    TH1F* spectrumDstarSgn = new TH1F(nameSgn.Data(), "D^{*} Signal invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);
    TH1F* spectrumDstarBkg = new TH1F(nameBkg.Data(), "D^{*} Background invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);

    nameMass="histSideBandMass_";
    nameMass+=i+1;
    if (i==-1) { 
      nameMass="histSideBandMass";
    }
    
    TH1F* spectrumSideBandMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} sideband mass; M(D^{*}) [GeV/c^{2}]; Entries",200,0.1,0.2);

    nameMass="histWrongSignMass_";
    nameMass+=i+1;
    if (i==-1) { 
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
void AliAnalysisTaskSEDStarSpectra::FillSpectrum(Int_t ptbin, AliAODRecoCascadeHF *part, Int_t isDStar, Bool_t PIDon, Int_t nSigma, AliAnalysisVertexingHF *vhf, TList *listout){
  //
  // Fill histos for D* spectrum
  //

  if (ptbin==0) return;

  Double_t invmassDelta = part->DeltaInvMass();
  Double_t invmassD0    = part->InvMassD0();
  Double_t invmassDstar = part->InvMassDstarKpipi();
  if(part->SelectDstar(fVHFloose->GetDstarCuts(),vhf->GetD0toKpiCuts(),kTRUE)) {//selected
    if (TMath::Abs(invmassD0-1.865)>fD0Window) return;
    
    Double_t pt = part->Pt();
    Double_t eta = part->Eta();
    //TVector3 p3Trk0(part->PxProng(0),part->PyProng(0),part->PzProng(0)); // pi_s
    //TVector3 p3Trk1(part->PxProng(1),part->PyProng(1),part->PzProng(1)); // D0
    //Double_t CosOpenAngle = p3Trk0.Dot(p3Trk1)/(p3Trk0.Mag()*p3Trk1.Mag());
    
    //PID of D0 daughters
    AliAODTrack *pos = (AliAODTrack*)part->Get2Prong()->GetDaughter(0);
    AliAODTrack *neg = (AliAODTrack*)part->Get2Prong()->GetDaughter(1);
    
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
      }
    
    TString fillthis="";
    Bool_t massInRange=kFALSE;
    if (TMath::Abs(invmassDelta-0.14557)<fPeakWindow) massInRange=kTRUE;
    
    if(fUseMCInfo) {
      if(isDStar==1) {
	//AliAODMCParticle *partDstar = (AliAODMCParticle*)mcArray->At(labDstar);
	//AliAODMCParticle *partPis = (AliAODMCParticle*)mcArray->At(partDstar->GetDaughter(1));
	//AliAODMCParticle *partD0 = (AliAODMCParticle*)mcArray->At(partDstar->GetDaughter(0));
	//AliAODMCParticle *partD0daughter = (AliAODMCParticle*)mcArray->At(partD0->GetDaughter(0));
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
        if (massInRange) {
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
        if (massInRange) {
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
 
  } //else cout<<"NOT SELECTED"<<endl;
  
  return;
}
//______________________________ side band background for D*___________________________________
void AliAnalysisTaskSEDStarSpectra::SideBandBackground(Int_t ptbin, AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma, AliAnalysisVertexingHF *vhf, TList *listout){

  //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
  // (expected detector resolution) on the left and right frm the D0 mass. Each band
  //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   

  if (ptbin==0) return;

  Double_t invmassDelta = part->DeltaInvMass();
  Double_t invmassD0    = part->InvMassD0();
  
  if(part->SelectDstar(fVHFloose->GetDstarCuts(),vhf->GetD0toKpiCuts(),kTRUE)) {//selected
    if(TMath::Abs(invmassD0-1.865)>2*fD0Window && TMath::Abs(invmassD0-1.865)<4*fD0Window){
      
      //PID of D0 daughters
      AliAODTrack *pos = (AliAODTrack*)part->Get2Prong()->GetDaughter(0);
      AliAODTrack *neg = (AliAODTrack*)part->Get2Prong()->GetDaughter(1);
      
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
	}
      
      TString fillthis="";
      fillthis="histSideBandMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      fillthis="histSideBandMass";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      
    }
  }
}

//________________________________________________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::WrongSignForDStar(Int_t ptbin, AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma, AliAnalysisVertexingHF *vhf, TList *listout){
  //
  // assign the wrong charge to the soft pion to create background
  //
  
  if (ptbin==0) return;
  
  AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)part->Get2Prong();
  Int_t okD0WrongSign,okD0barWrongSign;
  Double_t wrongMassD0=0.,wrongMassDstar=0.;
  theD0particle->SelectD0(vhf->GetD0toKpiCuts(),okD0WrongSign,okD0barWrongSign);
  
  if(part->Charge()>0) { 
    okD0WrongSign = 0;
  } else {
    okD0barWrongSign = 0;
  }
  //wrong D0 inv mass
  if(okD0WrongSign==1){
    wrongMassD0 = theD0particle->InvMassD0();
  }else if(okD0WrongSign==0){
    wrongMassD0 = theD0particle->InvMassD0bar();
  }
  
  if(part->SelectDstar(fVHFloose->GetDstarCuts(),vhf->GetD0toKpiCuts(),kTRUE)) {//selected
    if(TMath::Abs(wrongMassD0-1.865)<0.036){
      
      //PID of D0 daughters
      AliAODTrack *pos = (AliAODTrack*)part->Get2Prong()->GetDaughter(0);
      AliAODTrack *neg = (AliAODTrack*)part->Get2Prong()->GetDaughter(1);
      
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
	}
      
      // wrong D* inv mass    
      Double_t px[3],py[3],pz[3];
      UInt_t pdg[3]={321,211,211};
      pdg[0] = (part->Charge()>0 ? 321 : 211); // positive daughter of D0
      px[0] = theD0particle->PxProng(0);
      py[0] = theD0particle->PyProng(0);
      pz[0] = theD0particle->PzProng(0);
      pdg[1] = (part->Charge()>0 ? 211 : 321); // negative daughter of D0
      px[1] = theD0particle->PxProng(1);
      py[1] = theD0particle->PyProng(1);
      pz[1] = theD0particle->PzProng(1);
      pdg[2] = 211; // soft pion
      px[2] = part->PxProng(0);
      py[2] = part->PyProng(0);
      pz[2] = part->PzProng(0);
      
      Short_t dummycharge=0;
      Double_t dummyd0[3]={0,0,0};
      AliAODRecoDecay *rd = new AliAODRecoDecay(0x0,3,dummycharge,px,py,pz,dummyd0);
      
      wrongMassDstar = rd->InvMass(3,pdg);
      
      delete rd; rd=NULL;
      
      TString fillthis="";
      fillthis="histWrongSignMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar-wrongMassD0);
      fillthis="histWrongSignMass";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar-wrongMassD0);
      
    }
  }
}
//_____________________________________________SINGLE TRACK PRE-SELECTION___________________________________________
Bool_t AliAnalysisTaskSEDStarSpectra::SingleTrackSelections(AliAODRecoDecayHF2Prong* theD0particle, AliAODTrack *track0, AliAODTrack *track1, AliAODTrack *track2){
  
  // Preselection on D0 daughters and the soft pion

  // reft in ITS for soft pion
  //if((!(track2->GetStatus()&AliESDtrack::kITSrefit))) continue;
  
  // cut in acceptance for the soft pion and for the D0 daughters      
  Bool_t acceptanceProng0 = (TMath::Abs(theD0particle->EtaProng(0))<= 0.9 && theD0particle->PtProng(0) >= 0.1);
  Bool_t acceptanceProng1 = (TMath::Abs(theD0particle->EtaProng(1))<= 0.9 && theD0particle->PtProng(1) >= 0.1);
  // soft pion acceptance ... is it fine 0.9?????
  Bool_t acceptanceProng2 = (TMath::Abs(track2->Eta())<= 1.0 && track2->Pt() >= 0.05);
  
  if (!(acceptanceProng0 && acceptanceProng1 && acceptanceProng2)) return kFALSE;
  AliDebug(2,"D* reco daughters in acceptance");
  
  // cut on the min n. of clusters in ITS for the D0 and soft pion
  Int_t ncls0=0,ncls1=0,ncls2=0;
  for(Int_t l=0;l<6;l++) {
    if(TESTBIT(track0->GetITSClusterMap(),l)) ncls0++;
    if(TESTBIT(track1->GetITSClusterMap(),l)) ncls1++;
    if(TESTBIT(track2->GetITSClusterMap(),l)) ncls2++;
  }
  // see AddTask for soft pion and D0 prongs ITS clusters request
  if (!(ncls0 >= fMinITSClusters && ncls1 >= fMinITSClusters && ncls2>=fMinITSClustersSoft)) return  kFALSE;
  
  return kTRUE;
}

//_____________________________________________________________________________________________
Bool_t AliAnalysisTaskSEDStarSpectra::SetSelections(Double_t pt){
  
  //cuts[0] = inv. mass half width [GeV]
  //cuts[1] = dca [cm]
  //cuts[2] = cosThetaStar
  //cuts[3] = pTK [GeV/c]
  //cuts[4] = pTPi [GeV/c]
  //cuts[5] = d0K [cm]   upper limit!
  //cuts[6] = d0Pi [cm]  upper limit!
  //cuts[7] = d0d0 [cm^2]
  //cuts[8] = cosThetaPoint

  if (fAnalysis==0){
    if(pt<=1.){
      fVHF->SetD0toKpiCuts(0.450,0.02,0.7,0.8,0.8,0.1,0.1,0.000002,0.9);
      fD0Window=0.018;
      fPeakWindow=0.0018;
    }
    else if(pt >1. && pt <=2.){
      fVHF->SetD0toKpiCuts(0.450,0.03,0.7,0.8,0.8,0.1,0.1,-0.00002,0.9);
      fD0Window=0.020;
      fPeakWindow=0.0018;
      //fVHF->SetDstarCuts(0.3, 0.0018, 0.05, 100000000000.0, 0.5);
    }
    else if(pt >2. && pt <=3.){
      fVHF->SetD0toKpiCuts(0.450,0.03,0.7,0.8,0.8,0.1,0.1,-0.00002,0.9);
      fD0Window=0.020;
      fPeakWindow=0.0018;
      //fVHF->SetDstarCuts(0.3, 0.0018, 0.05, 100000000000.0, 0.5);
    }
    else if(pt >3. && pt <=5.){
      fVHF->SetD0toKpiCuts(0.450,0.03,0.7,0.9,0.9,0.1,0.1,0.000002,0.8);
      fD0Window=0.022;
      fPeakWindow=0.0016;
        //fVHF->SetDstarCuts(0.3, 0.0016, 0.05, 100000000000.0, 0.5);
    }
    else if(pt >5.){
      fVHF->SetD0toKpiCuts(0.450,0.03,0.7,1.0,1.0,0.1,0.1,0.000002,0.8);
      fD0Window=0.026;
      fPeakWindow=0.0014;
      //fVHF->SetDstarCuts(0.3, 0.0014, 0.05, 100000000000.0, 0.5);
    }
  }

  if(fAnalysis==1){
    if(pt<=1.){
      fVHF->SetD0toKpiCuts(0.450,0.04,0.8,0.21,0.21,0.021,0.021,-0.0002,0.9);
      fD0Window=0.024;
      fPeakWindow=0.0018;
    }
    else if(pt >1. && pt <=2.){
      fVHF->SetD0toKpiCuts(0.450,0.02,0.7,0.8,0.8,0.021,0.021,-0.0002,0.9);
      fD0Window=0.024;
      fPeakWindow=0.0018;
    }  
    else if(pt >2. && pt <=3.){
      fVHF->SetD0toKpiCuts(0.450,0.04,0.8,0.8,0.8,0.035,0.042,-0.000085,0.9);
      fD0Window=0.024;
      fPeakWindow=0.0018;
    } 
    else if(pt >3. && pt <=5.){
      fVHF->SetD0toKpiCuts(0.450,0.016,0.8,1.2,1.2,0.042,0.056,-0.000085,0.9);
      fD0Window=0.024;
      fPeakWindow=0.0016;
    } 
    else if(pt >5.){
      fVHF->SetD0toKpiCuts(0.450,0.08,1.0,1.2,1.2,0.07,0.07,0.0001,0.9);
      fD0Window=0.024;
      fPeakWindow=0.0014;
    }
  }
  return kTRUE;
}

//_____________________________ pid _______________________________________-
Bool_t AliAnalysisTaskSEDStarSpectra::SelectPID(AliAODTrack *track,  AliPID::EParticleType type, Double_t nsig){//type(0-4): {e,mu,pi,K,p}
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
/* //using new AliAODPidHF class
  if ((track->GetStatus()&AliESDtrack::kTPCpid )==0) return kTRUE;
  AliAODPidHF *pidHF = new AliAODPidHF();
  pidHF->SetSigma(nsig);
  if (type == AliPID::kPion) return pidHF->IsPionRaw(track,"TPC");
  if (type == AliPID::kKaon) return pidHF->IsKaonRaw(track,"TPC");
  return kFALSE;
*/
  //ITS
  //
  //
  //
  //TOF
  //
  //
  return isParticle;
}

