/*************************************************************************
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
* Author: Prabi, Feng Fan, and Paolo Bartalini, CCNU 2016                *
**************************************************************************/

/* AliAnaysisTaskMpiUE source code
*
* simple task which can serve as a starting point for building an MPI UE analysis
*
*/

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;


#include <Riostream.h>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
//#include "AliV0vertexer.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"

#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskESDfilter.h"

#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <AliHeader.h>

#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultVariable.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliESDUtils.h"

#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <TTree.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>
#include "TFile.h"
#include "TF1.h"

using std::cout;
using std::endl;

#include "AliAnalysisTaskMpiUE.h"

class AliAnalysisTaskMpiUE;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskMpiUE) // classimp: necessary for root

AliAnalysisTaskMpiUE::AliAnalysisTaskMpiUE() : AliAnalysisTaskSE(),
fESD(0), fStack(0),fMCEvent(0), fOutputList(0),fEventCounter(0),fCuts (0),fCuts2 (0), fUseMC(kFALSE), fEtaCut(0.8), fPtMin(1.0), f1(0),fTrackingEfficiency(0x0),fLeadPtCutMin(7.0), fLeadPtCutMax(40.0),
fHistEvents(0), fHistDCAxy(0), fRecNTracks(0), fGenNTracks(0), fGenLeadPart(0), fGenLeadPhi(0), fRecLeadPhi(0), fHistLeadDPhi(0), fHistLeadDPhiTr(0),fHistLeadDPhiTo(0),fHistLeadDPhiAw(0), fHistMatchParticleDensityTransverse(0), fHistMatchEnergyDensityTransverse(0),
fHistAllTracksDCAxy(0), fHistPrimaryDCAxy(0), fHistStrangeDCAxy(0), fHistPionDecayDCAxy(0), fHistMuonDecayDCAxy(0), fHistMaterialDCAxy(0), fHistOtherDCAxy(0),
fHistAllTracksPt(0), fHistPrimaryPt(0), fHistStrangePt(0), fHistMaterialPt(0), fHistZvtxMult(0), fHistZvtxMultTrue(0), fHistZvtxMultNcon1True(0), fHistZvtxMultNcon1(0), fHistZvtxMultNcon2(0), fHistGenZvtxMult(0), fHistMultResponseMat(0),
fHistMultiplicityStd(0),fHistEta(0), fHistPhi(0), fHistEta_Phi(0),fHistEta_Pt(0), fHistPt(0),fHistLeadPt(0), fHistMatchLeadPt(0), fHistNoMatchLeadPt(0), fHistSumPt(0),fHistAvgPt(0),fHistTracks(0), fHistZvtx(0), fHistZvtxCut(0),
fHistGenEta(0), fHistGenPhi(0), fHistGenPt(0),fHistGenLeadPt(0),fHistGenSumPt(0),fHistGenAvgPt(0), fHistGenTracks(0), fHistGenZvtx(0), fHistGenZvtxCut(0),
fHistParticleDensityTowards(0),fHistParticleDensityTransverse(0),fHistParticleDensityAway(0),fHistEnergyDensityTowards(0),fHistEnergyDensityTransverse(0),fHistEnergyDensityAway(0),
fHistGenParticleDensityTowards(0),fHistGenParticleDensityTransverse(0),fHistGenParticleDensityAway(0),fHistGenEnergyDensityTowards(0),fHistGenEnergyDensityTransverse(0),fHistGenEnergyDensityAway(0),
fHistTrack (0), fHistLeadTrack (0), fHistGenTrack (0), fHistGenLeadTrack (0), fHistDPhiEta(0),fHistGenDPhiEta(0),fHnDelta(0),fHnDeltaTo(0),fHnDeltaTr(0),fHnDeltaAw(0),
fHistGenMultTrans(0), fHistGenMultTowds(0), fHistGenMultAway(0), fHistGenEngyTrans(0), fHistGenEngyTowds(0), fHistGenEngyAway(0),
fHistMultTrans(0), fHistMultTowds(0), fHistMultAway(0), fHistEngyTrans(0), fHistEngyTowds(0), fHistEngyAway(0)
{
  // default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMpiUE::AliAnalysisTaskMpiUE(const char* name) : AliAnalysisTaskSE(name),
fESD(0), fStack(0),fMCEvent(0), fOutputList(0),fEventCounter(0),fCuts (0),fCuts2 (0), fUseMC(kFALSE), fEtaCut(0.8), fPtMin(1.0), f1(0),fTrackingEfficiency(0x0),fLeadPtCutMin(7.0), fLeadPtCutMax(40.0),
fHistEvents(0), fHistDCAxy(0), fRecNTracks(0), fGenNTracks(0), fGenLeadPart(0), fGenLeadPhi(0), fRecLeadPhi(0), fHistLeadDPhi(0), fHistLeadDPhiTr(0),fHistLeadDPhiTo(0),fHistLeadDPhiAw(0), fHistMatchParticleDensityTransverse(0), fHistMatchEnergyDensityTransverse(0),
fHistAllTracksDCAxy(0), fHistPrimaryDCAxy(0), fHistStrangeDCAxy(0), fHistPionDecayDCAxy(0), fHistMuonDecayDCAxy(0), fHistMaterialDCAxy(0), fHistOtherDCAxy(0),
fHistAllTracksPt(0), fHistPrimaryPt(0), fHistStrangePt(0), fHistMaterialPt(0), fHistZvtxMult(0),fHistZvtxMultTrue(0), fHistZvtxMultNcon1True(0), fHistZvtxMultNcon1(0), fHistZvtxMultNcon2(0), fHistGenZvtxMult(0), fHistMultResponseMat(0),
fHistMultiplicityStd(0),fHistEta(0), fHistPhi(0), fHistEta_Phi(0),fHistEta_Pt(0), fHistPt(0),fHistLeadPt(0), fHistMatchLeadPt(0), fHistNoMatchLeadPt(0), fHistSumPt(0),fHistAvgPt(0),fHistTracks(0), fHistZvtx(0), fHistZvtxCut(0),
fHistGenEta(0), fHistGenPhi(0), fHistGenPt(0),fHistGenLeadPt(0),fHistGenSumPt(0),fHistGenAvgPt(0), fHistGenTracks(0), fHistGenZvtx(0), fHistGenZvtxCut(0),
fHistParticleDensityTowards(0),fHistParticleDensityTransverse(0),fHistParticleDensityAway(0),fHistEnergyDensityTowards(0),fHistEnergyDensityTransverse(0),fHistEnergyDensityAway(0),
fHistGenParticleDensityTowards(0),fHistGenParticleDensityTransverse(0),fHistGenParticleDensityAway(0),fHistGenEnergyDensityTowards(0),fHistGenEnergyDensityTransverse(0),fHistGenEnergyDensityAway(0),
fHistTrack (0), fHistLeadTrack (0), fHistGenTrack (0), fHistGenLeadTrack (0), fHistDPhiEta(0),fHistGenDPhiEta(0),fHnDelta(0),fHnDeltaTo(0),fHnDeltaTr(0),fHnDeltaAw(0),
fHistGenMultTrans(0), fHistGenMultTowds(0), fHistGenMultAway(0), fHistGenEngyTrans(0), fHistGenEngyTowds(0), fHistGenEngyAway(0),
fHistMultTrans(0), fHistMultTowds(0), fHistMultAway(0), fHistEngyTrans(0), fHistEngyTowds(0), fHistEngyAway(0)
{

  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it, does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskMpiUE::~AliAnalysisTaskMpiUE()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    fOutputList = 0x0;
  }

  if(fCuts) {
    delete fCuts;
    fCuts = 0x0;
  }
  if(fCuts2) {
    delete fCuts2;
    fCuts2 = 0x0;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskMpiUE::UserCreateOutputObjects()
{

  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
  handler->SetNeedField();


  // fCuts *** is to confirm all tracks quality ***
  if(!fCuts ) {
     //fCuts = new AliESDtrackCuts();
     //fCuts->GetStandardITSTPCTrackCuts2011(kTRUE);  //***  13  ***//

     fCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
     fCuts->SetPtRange(fPtMin);
     fCuts->SetEtaRange(-fEtaCut, fEtaCut);

     // added new below   ** 2 **
          //fCuts->SetMaxFractionSharedTPCClusters(0.4); // track quality
          //fCuts->SetMaxDCAToVertexXY(2.4);   // reject secondaries
  }



   // fCuts2 *** is to get secondaries from MC ***
      fCuts2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
      fCuts2->SetPtRange(fPtMin);
      fCuts2->SetEtaRange(-fEtaCut, fEtaCut);

////*** StandardITSTPCTrackCuts2011() ***  13 ***////

       ///track quality =====
       // TPC
       fCuts2->SetMinNCrossedRowsTPC(70);
       fCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
       fCuts2->SetMaxChi2PerClusterTPC(4);
       fCuts2->SetRequireTPCRefit(kTRUE);
       fCuts2->SetMaxChi2TPCConstrainedGlobal(36);
       //ITS
       fCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
       fCuts2->SetMaxChi2PerClusterITS(25);
       fCuts2->SetRequireITSRefit(kTRUE);
       //
       fCuts2->SetRequireSigmaToVertex(kFALSE);


       /// reject secondaries =====
       fCuts2->SetMaxDCAToVertexZ(2);
       fCuts2->SetDCAToVertex2D(kFALSE);
       fCuts2->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
       fCuts2->SetAcceptKinkDaughters(kFALSE);


       //*** added new below *** 2***
       fCuts2->SetMaxFractionSharedTPCClusters(0.4); // track quality
       fCuts2->SetMaxDCAToVertexXY(2.4);   // reject secondaries


//------------------------------------------------------------------------



  // create output objects

  OpenFile(1);
  fOutputList = new TList();          // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written  to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested

  if(! fEventCounter ){
    fEventCounter = new TH1I( "fEventCounter", ";Evt. Sel. Step;Count",16,0,16);
    fEventCounter->GetXaxis()->SetBinLabel(1, "Total processed Events");
    fEventCounter->GetXaxis()->SetBinLabel(2, "Events selected for Analysis");
    fEventCounter->GetXaxis()->SetBinLabel(3, "Pileup events in SPDInMultBins");
    fEventCounter->GetXaxis()->SetBinLabel(4, "Events not with INEL>0");
    fEventCounter->GetXaxis()->SetBinLabel(5, "Events not in Vertex cut");
    fEventCounter->GetXaxis()->SetBinLabel(6, "Events inconsistent with SPD&TrackVertx");
    fEventCounter->GetXaxis()->SetBinLabel(7, "Events without MinimumBias (kMB)");
    fEventCounter->GetXaxis()->SetBinLabel(8, "ClusterVsTracklet Cut rejected events");
    fEventCounter->GetXaxis()->SetBinLabel(9, "Out of Bunch events in 11 IR");
    fEventCounter->GetXaxis()->SetBinLabel(10, "Incomplete DAQ events");
    fEventCounter->GetXaxis()->SetBinLabel(11, "Events Rejected by events quality criteria");
    fEventCounter->GetXaxis()->SetBinLabel(12, "Events Rejected by V0 decision");
    fEventCounter->GetXaxis()->SetBinLabel(13, "Events Rejected by V0 Assymetry");
    fEventCounter->GetXaxis()->SetBinLabel(14, "Events Rejected by Fastor Online cut");
    fEventCounter->GetXaxis()->SetBinLabel(15, "Events Without BCMod4 == 2 ");
    fEventCounter->GetXaxis()->SetBinLabel(16, "Events with No Good Tracks");




    fOutputList->Add(fEventCounter);
  }


  // Multiplicity and pT histogram

  const Int_t ptNbins = 36;
  const Double_t ptMin = 0.;
  const Double_t ptMax = 50.;

  Double_t ptbins1[ptNbins+1] = {0.0, 0.1, 0.15,  0.2,  0.25,  0.3,   0.35,  0.4,   0.45,  0.5,   0.6,   0.7,   0.8,   0.9,   1.0,   1.5,   2.0,   2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0};


  Int_t binsHisto1[3]   =     {ptNbins,  20,   200};
  Double_t binsminHisto1[3] = {ptMin,   -1.0,  0.0};
  Double_t binsmaxHisto1[3] = {ptMax,    1.0, TMath::TwoPi()};

  Int_t binsHisto2[4]   =     {ptNbins,  ptNbins,200, 100};
  Double_t binsminHisto2[4] = {ptMin,    ptMin,  0.0, -0.5};
  Double_t binsmaxHisto2[4] = {ptMax,    ptMax,  20.0, 99.5};


  fHistEvents = new TH1I("fHistEvents", "; Tracks ; Events", 300,-0.5,299.5);

  fHistAllTracksDCAxy = new TH1D("fHistAllTracksDCAxy", "AllTracksDCAxy; AllTracksDCAxy ; Counts", 200,-2.0, 2.0);
  fHistPrimaryDCAxy = new TH1D("fHistPrimaryDCAxy", "PrimaryDCAxy; PrimaryDCAxy ; Counts", 200,-2.0, 2.0);
  fHistStrangeDCAxy = new TH1D("fHistStrangeDCAxy", "StrangeDCAxy; StrangeDCAxy ; Counts", 200,-2.0, 2.0);
  fHistPionDecayDCAxy = new TH1D("fHistPionDecayDCAxy", "PionDecayDCAxy; PionDecayDCAxy ; Counts", 200,-2.0, 2.0);
  fHistMuonDecayDCAxy = new TH1D("fHistMuonDecayDCAxy", "MuonDecayDCAxy; MuonDecayDCAxy ; Counts", 200,-2.0, 2.0);
  fHistMaterialDCAxy  = new TH1D("fHistMaterialDCAxy", "MaterialDCAxy; MaterialDCAxy ; Counts", 200,-2.0, 2.0);
  fHistOtherDCAxy = new TH1D("fHistOtherDCAxy", "OtherDCAxy; OtherDCAxy ; Counts", 200,-2.0, 2.0);

  fHistAllTracksPt = new TH1D("fHistAllTracksPt", "pT; pT ; Counts",  36, ptbins1);
  fHistPrimaryPt = new TH1D("fHistPrimaryPt", "pT; pT ; Counts",  36, ptbins1);
  fHistStrangePt = new TH1D("fHistStrangePt", "pT; pT ; Counts",  36, ptbins1);
  fHistMaterialPt = new TH1D("fHistMaterialPt", "pT; pT ; Counts",  36, ptbins1);

  fHistDCAxy = new TH1D("fHistDCAxy", "DCAxy; DCAxy ; Counts", 200,-2.0, 2.0);


  fHistMultiplicityStd= new TH1I("fHistMultiplicityStd", "Standard Reference Multiplicity;Std Ref. Multiplicity  ; P(Nch) ",300,-0.5,299.5);       // create Reference Multiplicity histogram
  fHistPt = new TH1D("fHistPt", "pT; pT ; Events",  36, ptbins1);       // create pt histogram
  fHistLeadPt = new TH1D("fHistLeadPt", "LeadpT; LeadpT ; Events", 36, ptbins1);       // create Lead pt hisogram
  fHistMatchLeadPt = new TH1D("fHistMatchLeadPt", "LeadpT; LeadpT ; Events", 36, ptbins1);
  fHistNoMatchLeadPt = new TH1D("fHistNoMatchLeadPt", "LeadpT; LeadpT ; Events", 36, ptbins1);
  fHistSumPt = new TH1D("fHistSumPt", "SumpT; SumpT ; Events", 36, ptbins1);       // create Sum pt histogram
  fHistAvgPt = new TH1D("fHistAvgPt", "AvgpT; AvgpT ; Events", 200, 0, 20);       // create Sum pt histogram
  fHistEta = new TH1D("fHistEta", "#eta; #eta ; Events", 20, -1.0, 1.0);       // create Eta histogram
  fHistPhi = new TH1D("fHistPhi", "#phi; #phi ; Events", 200, 0,TMath::TwoPi());       // create Phi histogram
  fHistEta_Phi = new TH2D("fHistEta_Phi", "Eta_Phi; #eta ; #phi", 20, -1.0, 1.0, 200, 0, TMath::TwoPi());
    fHistEta_Pt = new TH2D("fHistEta_Pt", "Eta_Pt; pt ; #eta", 36,ptbins1,20, -1.0, 1.0);
  fHistZvtx = new TH1D("fHistZvtx", "Zvtx; Zvtx ; Events", 400, -20, 20);                   // Reconstructed Vertex Z
  fHistZvtxCut = new TH1D("fHistZvtxCut", "ZvtxCut; ZvtxCut ; Events", 400, -20, 20);    // Reconstructed Vertex Z with Cut
  fHistTracks = new TH1D("fHistTracks", "Tracks; Tracks ; Events", 300,-0.5,299.5);       // create Tracks histogram
  fHistZvtxMult = new TH2D("fHistZvtxMult", "ZvtxMult; Zvtx ; Mult", 400, -20, 20, 300, -0.5, 299.5);
  fHistZvtxMultTrue = new TH2D("fHistZvtxMultTrue", "ZvtxMultTrue; Zvtx ; MultTrue", 400, -20, 20, 300, -0.5, 299.5);
  fHistZvtxMultNcon1True = new TH2D("fHistZvtxMultNcon1True", "ZvtxMultNcon1True; ZvtxNcon1 ; MultTrue", 400, -20, 20, 300, -0.5, 299.5);
  fHistZvtxMultNcon1 = new TH2D("fHistZvtxMultNcon1", "ZvtxMultNcon1; ZvtxNcon1 ; Mult", 400, -20, 20, 300, -0.5, 299.5);
  fHistZvtxMultNcon2 = new TH2D("fHistZvtxMultNcon2", "ZvtxMultNcon2; ZvtxNcon2 ; Mult", 400, -20, 20, 300, -0.5, 299.5);

  fHistGenPt = new TH1D("fHistGenPt", "GenpT; GenpT ; Events",  36, ptbins1);
  fHistGenLeadPt = new TH1D("fHistGenLeadPt", "GenLeadpT; GenLeadpT ; Events", 36, ptbins1);       // Gen Lead pt histogram
  fHistGenSumPt = new TH1D("fHistGenSumPt", "GenSumpT; GenSumpT ; Events", 36, ptbins1);       // Gen Sum pt histogram
  fHistGenAvgPt = new TH1D("fHistGenAvgPt", "GenAvgpT; GenAvgpT ; Events", 200, 0, 20);       // Gen  Sum pt histogram
  fHistGenEta = new TH1D("fHistGenEta", "#eta; #eta ; Events", 20, -1.0, 1.0);
  fHistGenPhi = new TH1D("fHistGenPhi", "#phi; #phi ; Events", 200, 0,TMath::TwoPi());
  fHistGenZvtx = new TH1D("fHistGenZvtx", "GenZvtx; GenZvtx ; Events", 400, -20, 20);      // Generated Vertex Z
  fHistGenZvtxCut = new TH1D("fHistGenZvtxCut", "GenZvtxCut; GenZvtxCut ; Events", 400, -20, 20);  // Gen Vertex Z with Cut
  fHistGenTracks = new TH1D("fHistGenTracks", "GenTracks; GenTracks ; Events", 300,-0.5,299.5);   // Gen Multiplicity
  fHistGenZvtxMult = new TH2D("fHistGenZvtxMult", "GenZvtxMult; Zvtx ; Mult", 400, -20, 20, 300, -0.5, 299.5);
  fHistMultResponseMat = new TH2D("fHistMultResponseMat", "MultResponseMat; Gen tracks ; Reco tracks", 300,-0.5,299.5, 300,-0.5,299.5);

  fHistLeadDPhi = new TH1D("fHistLeadDPhi", "LeadDPhi; Leading-track #Delta#phi ; Events", 720, -TMath::TwoPi(), TMath::TwoPi());
  fHistLeadDPhiTr = new TH1D("fHistLeadDPhiTr", "LeadDPhiTr; Leading-track #Delta#phi ; Events", 720, -TMath::TwoPi(), TMath::TwoPi());
  fHistLeadDPhiTo = new TH1D("fHistLeadDPhiTo", "LeadDPhiTo; Leading-track #Delta#phi ; Events", 720, -TMath::TwoPi(), TMath::TwoPi());
  fHistLeadDPhiAw = new TH1D("fHistLeadDPhiAw", "LeadDPhiAw; Leading-track #Delta#phi ; Events", 720, -TMath::TwoPi(), TMath::TwoPi());

  fHistTrack = new THnSparseF("fHistTrack", "fHistTrack ; p_{T}(GeV/c) ;  #eta; #phi",3,binsHisto1, binsminHisto1,binsmaxHisto1);       // create pt histogram
  fHistTrack->SetBinEdges(0,ptbins1);
  fHistTrack->GetAxis(0)->SetTitle("p_{T}(GeV/c)");
  fHistTrack->GetAxis(1)->SetTitle(" #eta");
  fHistTrack->GetAxis(2)->SetTitle(" #phi");

  fHistLeadTrack = new THnSparseF("fHistLeadTrack", "fHistLeadTrack ; Lead p_{T}(GeV/c) ; Sum p_{T}(GeV/c) ; Average p_{T}(GeV/c) ; tracks",4,binsHisto2, binsminHisto2,binsmaxHisto2);       // create pt histogram
  fHistLeadTrack->SetBinEdges(0,ptbins1);
  fHistLeadTrack->SetBinEdges(1,ptbins1);

  fHistGenTrack = new THnSparseF("fHistGenTrack", "fHistGenTrack ; Gen p_{T}(GeV/c) ;  #eta; #phi",3,binsHisto1, binsminHisto1,binsmaxHisto1);       // create pt histogram
  fHistGenTrack->SetBinEdges(0,ptbins1);
  fHistGenTrack->GetAxis(0)->SetTitle("p_{T}(GeV/c)");
  fHistGenTrack->GetAxis(1)->SetTitle(" #eta");
  fHistGenTrack->GetAxis(2)->SetTitle(" #phi");

  fHistGenLeadTrack = new THnSparseF("fHistGenLeadTrack", "fHistGenLeadTrack ; Gen Lead p_{T}(GeV/c) ; Sum p_{T}(GeV/c) ; Average p_{T}(GeV/c) ; tracks",4,binsHisto2, binsminHisto2,binsmaxHisto2);       // create pt histogram
  fHistGenLeadTrack->SetBinEdges(0,ptbins1);
  fHistGenLeadTrack->SetBinEdges(1,ptbins1);

  //f3dHistPrimAnalysisPtVsYCMSVsMultMCLambda = new TH3F( "f3dHistPrimAnalysisPtVsYCMSVsMultMCLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; MultMC", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);



   fHistMatchParticleDensityTransverse = new TProfile ("fHistMatchParticleDensityTransverse","Matched Particle Density Transverse region vs P_{T}(leading); LeadpT (Transverse); N_{ch}", 36, ptbins1) ;
   fHistMatchEnergyDensityTransverse = new TProfile ("fHistMatchEnergyDensityTransverse","Matched Energy Density Transverse region vs p_{T}(leading); LeadpT (Transverse) ; SumpT", 36, ptbins1);

   fHistParticleDensityTransverse = new TProfile ("fHistParticleDensityTransverse","Particle Density Transverse region vs P_{T}(leading); LeadpT (Transverse); N_{ch}", 36, ptbins1) ;
   fHistParticleDensityTowards = new TProfile ("fHistParticleDensityTowards","Particle Density Towards region vs p_{T}(leading); LeadpT (Towardss) ; N_{ch}", 36, ptbins1);
   fHistParticleDensityAway = new TProfile ("fHistParticleDensityAway","Particle Density Away region vs p_{T}(leading); LeadpT (Away) ; N_{ch}", 36, ptbins1) ;
   fHistEnergyDensityTransverse = new TProfile ("fHistEnergyDensityTransverse","Energy Density Transverse region vs p_{T}(leading); LeadpT (Transverse) ; SumpT", 36, ptbins1);
   fHistEnergyDensityTowards = new TProfile ("fHistEnergyDensityTowards","Energy Density Towards region vs p_{T}(leading); LeadpT (Towardss) ; SumpT", 36, ptbins1);
   fHistEnergyDensityAway = new TProfile ("fHistEnergyDensityAway","Energy Density Away region vs p_{T}(leading); LeadpT (Away) ; SumpT", 36, ptbins1);

   fHistGenParticleDensityTransverse = new TProfile ("fHistGenParticleDensityTransverse","Particle Density Transverse region vs P_{T}(leading); LeadpT (Transverse); N_{ch}", 36, ptbins1) ;
   fHistGenParticleDensityTowards = new TProfile ("fHistGenParticleDensityTowards","Particle Density Towards region vs p_{T}(leading); LeadpT (Towardss) ; N_{ch}", 36, ptbins1);
   fHistGenParticleDensityAway = new TProfile ("fHistGenParticleDensityAway","Particle Density Away region vs p_{T}(leading); LeadpT (Away) ; N_{ch}", 36, ptbins1) ;
   fHistGenEnergyDensityTransverse = new TProfile ("fHistGenEnergyDensityTransverse","Energy Density Transverse region vs p_{T}(leading); LeadpT (Transverse) ; SumpT", 36, ptbins1);
   fHistGenEnergyDensityTowards = new TProfile ("fHistGenEnergyDensityTowards","Energy Density Towards region vs p_{T}(leading); LeadpT (Towardss) ; SumpT", 36, ptbins1);
   fHistGenEnergyDensityAway = new TProfile ("fHistGenEnergyDensityAway","Energy Density Away region vs p_{T}(leading); LeadpT (Away) ; SumpT", 36, ptbins1);



  Int_t binsHisto[3]   = {40,  31, 50};
  Double_t binsminHisto[3] = {-2.0, 0.0, 0.0};
  Double_t binsmaxHisto[3] = {2.0, TMath::Pi(), 5.0};

  fHnDelta = new THnSparseD("fHnDelta","deltaEta:deltaPhi:deltaR",3,binsHisto,binsminHisto,binsmaxHisto);
  fHnDelta->GetAxis(0)->SetTitle("#Delta#eta");
  fHnDelta->GetAxis(1)->SetTitle("#Delta#phi");
  fHnDelta->GetAxis(2)->SetTitle("#Delta");

  fHnDeltaTo = new THnSparseD("fHnDeltaTo","deltaEta:deltaPhi:deltaR",3,binsHisto,binsminHisto,binsmaxHisto);
  fHnDeltaTo->GetAxis(0)->SetTitle("#Delta#eta (Towardss)");
  fHnDeltaTo->GetAxis(1)->SetTitle("#Delta#phi (Towardss)");
  fHnDeltaTo->GetAxis(2)->SetTitle("#DeltaR (Towardss)");

  fHnDeltaTr = new THnSparseD("fHnDeltaTr","deltaEta:deltaPhi:deltaR",3,binsHisto,binsminHisto,binsmaxHisto);
  fHnDeltaTr->GetAxis(0)->SetTitle("#Delta#eta (Transverse)");
  fHnDeltaTr->GetAxis(1)->SetTitle("#Delta#phi (Transverse)");
  fHnDeltaTr->GetAxis(2)->SetTitle("#DeltaR (Transverse)");

  fHnDeltaAw = new THnSparseD("fHnDeltaAw","deltaEta:deltaPhi:deltaR",3,binsHisto,binsminHisto,binsmaxHisto);
  fHnDeltaAw->GetAxis(0)->SetTitle("#Delta#eta (Away)");
  fHnDeltaAw->GetAxis(1)->SetTitle("#Delta#phi (Away)");
  fHnDeltaAw->GetAxis(2)->SetTitle("#DeltaR (Away)");



  fHistDPhiEta = new TH2D("fHistDPhiEta","Eta Phi distribution (Reconstructed)",124, -(TMath::Pi())/2, (3*TMath::Pi()/2),  40, -2.0, 2.0);
  fHistGenDPhiEta = new TH2D("fHistGenDPhiEta","Eta Phi distribution (Generated)",124, -(TMath::Pi())/2, (3*TMath::Pi()/2),  40, -2.0, 2.0);

    fHistMultTrans = new TH1D("fHistMultTrans", "MultTrans; N_{ch, Trans}; Events", 300,-0.5,299.5);   //  Multiplicity in Trans
    fHistMultTowds = new TH1D("fHistMultTowds", "MultTowds; N_{ch, Towds}; Events", 300,-0.5,299.5);   //  Multiplicity in Towds
    fHistMultAway = new TH1D("fHistMultAway", "MultAway; N_{ch, Away}; Events", 300,-0.5,299.5);       //  Multiplicity in Away
    fHistEngyTrans = new TH1D("fHistEngyTrans", "EngyTrans; Sum p_{T, Trans}(GeV/c); Events", 36, ptbins1);      //  Summed pT in Trans
    fHistEngyTowds = new TH1D("fHistEngyTowds", "EngyTowds; Sum p_{T, Towds}(GeV/c); Events", 36, ptbins1);      //  Summed pT in Towds
    fHistEngyAway = new TH1D("fHistEngyAway", "EngyAway; Sum p_{T, Away}(GeV/c); Events", 36, ptbins1);          //  Summed pT in Away
///............for gen
    fHistGenMultTrans = new TH1D("fHistGenMultTrans", "GenMultTrans; N_{ch, Trans}; Events", 300,-0.5,299.5);   // Gen Multiplicity in Trans
    fHistGenMultTowds = new TH1D("fHistGenMultTowds", "GenMultTowds; N_{ch, Towds}; Events", 300,-0.5,299.5);   // Gen Multiplicity in Towds
    fHistGenMultAway = new TH1D("fHistGenMultAway", "GenMultAway; N_{ch, Away}; Events", 300,-0.5,299.5);       // Gen Multiplicity in Away
    fHistGenEngyTrans = new TH1D("fHistGenEngyTrans", "GenEngyTrans; Sum p_{T, Trans}(GeV/c) ; Events", 36, ptbins1);      // Gen Summed pT in Trans
    fHistGenEngyTowds = new TH1D("fHistGenEngyTowds", "GenEngyTowds; Sum p_{T, Towds}(GeV/c) ; Events", 36, ptbins1);      // Gen Summed pT in Towds
    fHistGenEngyAway = new TH1D("fHistGenEngyAway", "GenEngyAway; Sum p_{T, Away}(GeV/c) ; Events", 36, ptbins1);          // Gen Summed pT in Away
    
 // TF1 *f1 = new TF1( "f1", "([0]+([1]/(x-[2])))", 2., 5.);  //for vertex correction
 //    f1->SetParameter(0,9.999e-01);
 //    f1->SetParameter(1,9.832e-03);
 //    f1->SetParameter(2,-4.789e-02);

  // add Sumw2 for your histogram

  fHistEvents->Sumw2();

  fHistAllTracksDCAxy->Sumw2();
  fHistPrimaryDCAxy->Sumw2();
  fHistStrangeDCAxy->Sumw2();
  fHistPionDecayDCAxy->Sumw2();
  fHistMuonDecayDCAxy->Sumw2();
  fHistMaterialDCAxy->Sumw2();
  fHistOtherDCAxy->Sumw2();

  fHistAllTracksPt->Sumw2();
  fHistPrimaryPt->Sumw2();
  fHistStrangePt->Sumw2();
  fHistMaterialPt->Sumw2();

  fHistDCAxy->Sumw2();


  fHistMultiplicityStd->Sumw2();
  fHistPt->Sumw2();
  fHistLeadPt->Sumw2();
  fHistMatchLeadPt->Sumw2();
  fHistNoMatchLeadPt->Sumw2();
  fHistSumPt->Sumw2();
  fHistAvgPt->Sumw2();
  fHistEta->Sumw2();
  fHistPhi->Sumw2();
  fHistEta_Phi->Sumw2();
    fHistEta_Pt->Sumw2();

  fHistZvtx->Sumw2();
  fHistZvtxCut->Sumw2();
  fHistTracks->Sumw2();
  fHistZvtxMult->Sumw2();
  fHistZvtxMultTrue->Sumw2();
  fHistZvtxMultNcon1->Sumw2();
  fHistZvtxMultNcon1True->Sumw2();
  fHistZvtxMultNcon2->Sumw2();


  fHistGenPt->Sumw2();
  fHistGenLeadPt->Sumw2();
  fHistGenSumPt->Sumw2();
  fHistGenAvgPt->Sumw2();
  fHistGenEta->Sumw2();
  fHistGenPhi->Sumw2();
  fHistGenZvtx->Sumw2();
  fHistGenZvtxCut->Sumw2();
  fHistGenTracks->Sumw2();
  fHistGenZvtxMult->Sumw2();
  fHistMultResponseMat->Sumw2();

  fHistLeadDPhi->Sumw2();
  fHistLeadDPhiTr->Sumw2();
  fHistLeadDPhiTo->Sumw2();
  fHistLeadDPhiAw->Sumw2();

  fHistTrack->Sumw2();
  fHistLeadTrack->Sumw2();
  fHistGenTrack->Sumw2();
  fHistGenLeadTrack->Sumw2();

  fHistMatchParticleDensityTransverse->Sumw2();
  fHistMatchEnergyDensityTransverse->Sumw2();

  fHistParticleDensityTransverse->Sumw2();
  fHistParticleDensityTowards->Sumw2();
  fHistParticleDensityAway->Sumw2();
  fHistEnergyDensityTransverse->Sumw2();
  fHistEnergyDensityTowards->Sumw2();
  fHistEnergyDensityAway->Sumw2();
  fHistGenParticleDensityTransverse->Sumw2();
  fHistGenParticleDensityTowards->Sumw2();
  fHistGenParticleDensityAway->Sumw2();
  fHistGenEnergyDensityTransverse->Sumw2();
  fHistGenEnergyDensityTowards->Sumw2();
  fHistGenEnergyDensityAway->Sumw2();


  fHnDelta->Sumw2();
  fHnDeltaTo->Sumw2();
  fHnDeltaTr->Sumw2();
  fHnDeltaAw->Sumw2();
  fHistDPhiEta->Sumw2();
  fHistGenDPhiEta->Sumw2();

    fHistMultTrans->Sumw2();
    fHistMultTowds->Sumw2();
    fHistMultAway->Sumw2();
    fHistEngyTrans->Sumw2();
    fHistEngyTowds->Sumw2();
    fHistEngyAway->Sumw2();

    fHistGenMultTrans->Sumw2();
    fHistGenMultTowds->Sumw2();
    fHistGenMultAway->Sumw2();
    fHistGenEngyTrans->Sumw2();
    fHistGenEngyTowds->Sumw2();
    fHistGenEngyAway->Sumw2();


  // your histogram in the output file, add it to the list!

  fOutputList->Add(fHistMultiplicityStd);
  fOutputList->Add(fHistEvents);


  fOutputList->Add(fHistAllTracksDCAxy);
  fOutputList->Add(fHistPrimaryDCAxy);
  fOutputList->Add(fHistStrangeDCAxy);
  fOutputList->Add(fHistPionDecayDCAxy);
  fOutputList->Add(fHistMuonDecayDCAxy);
  fOutputList->Add(fHistMaterialDCAxy);
  fOutputList->Add(fHistOtherDCAxy);

  fOutputList->Add(fHistAllTracksPt);
  fOutputList->Add(fHistPrimaryPt);
  fOutputList->Add(fHistStrangePt);
  fOutputList->Add(fHistMaterialPt);

  fOutputList->Add(fHistDCAxy);


  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistPhi);
  fOutputList->Add(fHistEta_Phi);
    fOutputList->Add(fHistEta_Pt);

  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistLeadPt);
  fOutputList->Add(fHistMatchLeadPt);
  fOutputList->Add(fHistNoMatchLeadPt);
  fOutputList->Add(fHistSumPt);
  fOutputList->Add(fHistAvgPt);
  fOutputList->Add(fHistTracks);

  fOutputList->Add(fHistZvtx);
  fOutputList->Add(fHistZvtxCut);
  fOutputList->Add(fHistZvtxMult);
  fOutputList->Add(fHistZvtxMultTrue);
  fOutputList->Add(fHistZvtxMultNcon1);
  fOutputList->Add(fHistZvtxMultNcon1True);
  fOutputList->Add(fHistZvtxMultNcon2);

  fOutputList->Add(fHistGenEta);
  fOutputList->Add(fHistGenPhi);
  fOutputList->Add(fHistGenPt);
  fOutputList->Add(fHistGenLeadPt);
  fOutputList->Add(fHistGenSumPt);
  fOutputList->Add(fHistGenAvgPt);
  fOutputList->Add(fHistGenTracks);

  fOutputList->Add(fHistGenZvtx);
  fOutputList->Add(fHistGenZvtxCut);
  fOutputList->Add(fHistGenZvtxMult);
  fOutputList->Add(fHistMultResponseMat);

  fOutputList->Add(fHistLeadDPhi);
  fOutputList->Add(fHistLeadDPhiTr);
  fOutputList->Add(fHistLeadDPhiTo);
  fOutputList->Add(fHistLeadDPhiAw);

  fOutputList->Add(fHistMatchParticleDensityTransverse);
  fOutputList->Add(fHistMatchEnergyDensityTransverse);

  fOutputList->Add(fHistParticleDensityTransverse);
  fOutputList->Add(fHistParticleDensityTowards);
  fOutputList->Add(fHistParticleDensityAway);
  fOutputList->Add(fHistEnergyDensityTransverse);
  fOutputList->Add(fHistEnergyDensityTowards);
  fOutputList->Add(fHistEnergyDensityAway);
  fOutputList->Add(fHistGenParticleDensityTransverse);
  fOutputList->Add(fHistGenParticleDensityTowards);
  fOutputList->Add(fHistGenParticleDensityAway);
  fOutputList->Add(fHistGenEnergyDensityTransverse);
  fOutputList->Add(fHistGenEnergyDensityTowards);
  fOutputList->Add(fHistGenEnergyDensityAway);


  fOutputList->Add(fHistTrack);
  fOutputList->Add(fHistLeadTrack);
  fOutputList->Add(fHistGenTrack);
  fOutputList->Add(fHistGenLeadTrack);

  fOutputList->Add(fHnDelta);
  fOutputList->Add(fHnDeltaTr);
  fOutputList->Add(fHnDeltaTo);
  fOutputList->Add(fHnDeltaAw);

  fOutputList->Add(fHistDPhiEta);
  fOutputList->Add(fHistGenDPhiEta);

    fOutputList->Add(fHistMultTrans);
    fOutputList->Add(fHistMultTowds);
    fOutputList->Add(fHistMultAway);
    fOutputList->Add(fHistEngyTrans);
    fOutputList->Add(fHistEngyTowds);
    fOutputList->Add(fHistEngyAway);

    
    fOutputList->Add(fHistGenMultTrans);
    fOutputList->Add(fHistGenMultTowds);
    fOutputList->Add(fHistGenMultAway);
    fOutputList->Add(fHistGenEngyTrans);
    fOutputList->Add(fHistGenEngyTowds);
    fOutputList->Add(fHistGenEngyAway);
    

  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output

}
//_____________________________________________________________________________
void AliAnalysisTaskMpiUE::UserExec(Option_t *)
{
    //cout<<"is MC<"<<fUseMC<<">"<<endl;//BUG: in inition, the fUseMC=kFALse, but now <1>


    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD) {
      AliWarning("ERROR: ESD event is not available \n");
      return;                                                 // if the pointer to the event is empty (getting it failed) skip this event
    }

    fEventCounter->Fill(0);

    ///< processed events counter
    if(HasNoPileupSPDInMultBins(fESD) == kFALSE){fEventCounter->Fill(2);}                  ///<  Rejects events tagged with IsPileupFromSPDInMultBins()
    if(IsINELgtZERO( fESD ) == kFALSE){fEventCounter->Fill(3);}                            ///<  INEL > 0 (with SPD tracklets atleast 1)
    if(HasAcceptedVertexPosition( fESD ) == kFALSE){fEventCounter->Fill(4);}               ///<  Checks for accepted vertex position (|eta|<10cm)
    if(HasNoInconsistentSPDandTrackVertices( fESD ) == kFALSE){fEventCounter->Fill(5);}    ///<  Checks for consistent SPD and track vertex (if track vertex exists)
    if(IsMinimumBias( fESD ) == kFALSE){fEventCounter->Fill(6);}                           ///<  kMB minimum bias trigger selection
    if(IsSPDClusterVsTrackletBG( fESD ) == kFALSE){fEventCounter->Fill(7);}                ///<  counts rejected events
    if(IsOutOfBunchPileup( fESD ) == kTRUE){fEventCounter->Fill(8);}                       ///<  counts rejected events
    if(IsInCompleteEvent( fESD ) == kTRUE){fEventCounter->Fill(9);}                        ///<  counts rejected events
    if(IsEventSelected( fESD ) == kFALSE){fEventCounter->Fill(10);}                        ///<  counts rejected events
    if(V0Decision( fESD ) == kFALSE){fEventCounter->Fill(11);}                      ///<  counts rejected events
    if(V0Asymmetry( fESD ) == kFALSE){fEventCounter->Fill(12);}                     ///<  counts rejected events
    if(IsOutOfBunchPileupFastor( fESD ) == kFALSE){fEventCounter->Fill(13);}        ///<  counts rejected events
    if(HasBCMod4( fESD ) == kFALSE){fEventCounter->Fill(14);}                       ///<  counts WITHOUT BCMOD4==2
    if(fCuts->CountAcceptedTracks(fESD)>0 == kFALSE){fEventCounter->Fill(15);}                       ///<  counts WITHOUT good any tracks



    if(!IsEventSelected(fESD) ) {
      PostData(1, fOutputList);   /// Event isn't selected, post output data, done here
      return;
    }

    fEventCounter->Fill(1);

   if (fUseMC){
      FillMCGen();  // Fill all the Generator level information
   }

   fHistMultiplicityStd->Fill(AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE)); // this will give you standard reference multiplicity

   //Int_t Evt_NoGoodTracks =0;

   if (fCuts2->CountAcceptedTracks(fESD)>0) {
         TracksLoop(fESD);
         //cout<<"GenNTracks("<<fGenNTracks<<")***"<<endl;       //$$$  20180108  Test for the different entries between fHistMultResponseMat & fHistTracks (test for the logic of the code)
         //cout<<"RecNTracks("<<fRecNTracks<<")***"<<endl;

         //if (fGenNTracks > 0 && fRecNTracks >0 ) {
         //     fHistMultResponseMat->Fill(fGenNTracks,fRecNTracks);
         //     Double_t fLeadDPhi = deltaPhi(fGenLeadPhi,fRecLeadPhi);
         //     fHistLeadDPhi->Fill(fLeadDPhi);
         //}
         //fGenNTracks = 0; fRecNTracks = 0; //$$$  fGenNTracks & fRecNTracks are two global parameters,
                                             //$$$  if in the current event there is no value for the two paprmeters, then they will be equal to the valus from the previous event.

   }
   else {
      //   AliWarning("<<<< 0 Tracks After All Cuts >>>>>>>");
     return;
   }

   PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file

}


//______________________________________________________________________________
void AliAnalysisTaskMpiUE::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}


//_____________________________________________________________________________
void AliAnalysisTaskMpiUE::FillMCGen() {

  AliAnalysisManager *anmgr=AliAnalysisManager::GetAnalysisManager();
  AliMCEventHandler* eventHandler = 0;
  fMCEvent = 0;
  fStack = 0;
  eventHandler = (AliMCEventHandler*)anmgr->GetMCtruthEventHandler();
  if (!eventHandler) { printf("ERROR: Could not retrieve MC event handler\n"); return; }
  fMCEvent = eventHandler->MCEvent();
  if (!fMCEvent) { printf("ERROR: Could not retrieve MC event\n"); return; }
  fStack = fMCEvent->Stack();
  if (!fStack) { printf("Stack not available\n"); return; }


  AliGenEventHeader* mcGenH = 0;
  AliGenPythiaEventHeader* hPythia=0;

  mcGenH = fMCEvent->GenEventHeader();
  Double_t event_weight = fMCEvent->GenEventHeader()->EventWeight();

  if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
    TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
    TIter next(headers);
    AliGenEventHeader* mcGenH1 = 0;
    while ( (mcGenH1=(AliGenEventHeader*)next()) ) {
      if (mcGenH1->InheritsFrom(AliGenPythiaEventHeader::Class())) {
        hPythia = (AliGenPythiaEventHeader*)mcGenH1;
        break;
      }
    }
  }

  else {} // unknown generator

  TArrayF vtmc(3);
  for (int i = 0; i < 3; i++) vtmc[i] = 0;
  mcGenH->PrimaryVertex(vtmc);
  for (int i=3;i--;) fVertexMC[i] = vtmc[i];
  fHistGenZvtx->Fill(fVertexMC[2],event_weight);

  if (!(fVertexMC[2] >= -10.0 &&  fVertexMC[2] <= 10.0)) return;   //  Vertex Z cut on Gen. particles
  fHistGenZvtxCut->Fill(fVertexMC[2],event_weight);

  Float_t totalNch = 0;
  Double_t fgenlPt = 0;
  Double_t fgensPt = 0;
  Double_t fgenrPt = 0;
  Double_t fgenlPhi = 0;
  Double_t fgenlEta =0;
Double_t Etagencms = 0;
  for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (!fMCEvent->IsPhysicalPrimary(i)) continue;  //  Particles produced including products of strong and electromagnetic decays but excluding feed-down from weak decays of strange particles like Ks,Lambda etc)
    AliMCParticle* particle = (AliMCParticle*)fMCEvent->GetTrack(i);
    if (!particle) continue;
    if (particle->Pt() == 0 || particle->E() <= 0) continue;
    Double_t theta = particle->Theta();
    if (theta<1e-6 || theta>TMath::Pi()-1e-6) continue;
    if (particle->Charge() == 0) continue;
    if (!(TMath::Abs(particle->Eta()) < fEtaCut && particle->Pt()> fPtMin)) continue;
    totalNch++;
Etagencms = (particle->Eta()-0.465);
    fHistGenPt->Fill(particle->Pt(),event_weight);
    //fHistGenEta->Fill(particle->Eta(),event_weight);
      fHistGenEta->Fill(Etagencms,event_weight);
    fHistGenPhi->Fill(particle->Phi(),event_weight);

    Double_t vars_g1a[3] = {particle->Pt(),particle->Eta(),particle->Phi()};
    fHistGenTrack->Fill(vars_g1a,event_weight);

    fgenrPt= particle->Pt();
    fgensPt +=fgenrPt;

    if (fgenrPt > fgenlPt){
      fgenlPt = fgenrPt;
      fgenlPhi = particle->Phi();
      fgenlEta = Etagencms;
      fGenLeadPart = particle;
    }
  }

  fGenNTracks = totalNch;
  //cout<<"in FillMCGen GenN<"<<fGenNTracks<<">-----"<<endl;
  fGenLeadPhi = fgenlPhi;



  if (totalNch < 1) return;
  fHistGenSumPt->Fill(fgensPt,event_weight);
  fHistGenAvgPt->Fill(fgensPt/totalNch,event_weight);
  fHistGenLeadPt->Fill(fgenlPt,event_weight);
  fHistGenTracks->Fill(totalNch,event_weight);
  fHistGenZvtxMult->Fill(fVertexMC[2],totalNch);


  Double_t trg1 = totalNch;
  Double_t varsg1a[4] = {fgenlPt,fgensPt,fgensPt/totalNch,trg1};
  fHistGenLeadTrack->Fill(varsg1a,event_weight);


  Int_t ng1=0;
  Int_t ng2=0;
  Int_t ng3=0;

  Double_t fgsPt1 = 0.0;
  Double_t fgsPt2 = 0.0;
  Double_t fgsPt3 = 0.0;
  Double_t dphi_gen = 0.0;
  Double_t dr_gen = 0.0;


  for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (!fMCEvent->IsPhysicalPrimary(i)) continue; //  Particles produced including products of strong and electromagnetic decays but excluding feed-down from weak decays of strange particles like Ks,Lambda etc)
    AliMCParticle* particle = (AliMCParticle*)fMCEvent->GetTrack(i);
    if (!particle) continue;
    if (particle->Pt() == 0 || particle->E() <= 0) continue;
    Double_t theta = particle->Theta();
    if (theta<1e-6 || theta>TMath::Pi()-1e-6) continue;
    if (particle->Charge() == 0) continue;
    if (!(TMath::Abs(particle->Eta()) < fEtaCut && particle->Pt()> fPtMin)) continue;

    Double_t dphi_gen1 = (particle->Phi()-fgenlPhi);
    dphi_gen = DeltaPhi(particle->Phi(),fgenlPhi);
    Double_t deta_gen = particle->Eta()-fgenlEta;
    Double_t dr_gen = DeltaR(particle->Phi(),particle->Eta(),fgenlPhi,fgenlEta);

    //if (TMath::Abs(dr_gen) < 0.0001) continue;   // Protection to make sure we dont use Leading track itself for the correlation


         if(dphi_gen1>((3./2)*TMath::Pi())) dphi_gen1 = dphi_gen1-2*TMath::Pi();
          else if(dphi_gen1<-TMath::Pi()/2) dphi_gen1=dphi_gen1+2*TMath::Pi();

        fHistGenDPhiEta->Fill(dphi_gen1,deta_gen);
        //cout << " pT =" << particle->Pt() << endl;



        Double_t vars_g[3] = {deta_gen,dphi_gen,dr_gen};

        if (TMath::Abs(dphi_gen)<(1./3)*TMath::Pi())
        {
          fgsPt1 += particle->Pt();
          Double_t vars_g1[3] = {deta_gen,dphi_gen,dr_gen};
          ng1++;
        }

        else if ((1./3)*TMath::Pi()<TMath::Abs(dphi_gen) && TMath::Abs(dphi_gen) <(2./3)*TMath::Pi())
        {
          fgsPt2 += particle->Pt();
          Double_t vars_g2[3] = {deta_gen,dphi_gen,dr_gen};
          ng2++;
        }

        else
        {
          fgsPt3 += particle->Pt();
          Double_t vars_g3[3] = {deta_gen,dphi_gen,dr_gen};
          ng3++;
        }
      }
    if (fgenlPt > fLeadPtCutMin)
    {
        //cout<<"Gen level:     Sum pT(Towds)="<<fgsPt1<<"  Sum pT(Trans)="<<fgsPt2<<"  Sum pT(Away)="<<fgsPt3<<"***************************************"<<endl;
        fHistGenMultTrans->Fill(ng2);
        fHistGenMultTowds->Fill(ng1);
        fHistGenMultAway->Fill(ng3);
        fHistGenEngyTrans->Fill(fgsPt2);
        fHistGenEngyTowds->Fill(fgsPt1);
        fHistGenEngyAway->Fill(fgsPt3);
    }
    
      if ((ng1+ng2+ng3) >= 1)  // Protection just to make sure if less than one track present in an event then do not fill histograms
      {

        fHistGenParticleDensityTowards->Fill(fgenlPt,ng1);
        fHistGenParticleDensityTransverse->Fill(fgenlPt,ng2);
        fHistGenParticleDensityAway->Fill(fgenlPt,ng3);

        fHistGenEnergyDensityTowards->Fill(fgenlPt,fgsPt1);
        fHistGenEnergyDensityTransverse->Fill(fgenlPt,fgsPt2);
        fHistGenEnergyDensityAway->Fill(fgenlPt,fgsPt3);

      }

}

//_____________________________________________________________________________
void AliAnalysisTaskMpiUE::TracksLoop(AliESDEvent *fESD)
{

  const AliVVertex *lPrimaryVtx  = fESD->GetPrimaryVertex();
  //const AliESDVertex *lPrimaryVtx  = fESD->GetPrimaryVertexSPD();
  Double_t Zvertex = lPrimaryVtx->GetZ();
  fHistZvtx->Fill(Zvertex);
  if (!(TMath::Abs(Zvertex) <= 10.0)) return;
  fHistZvtxCut->Fill(Zvertex);


            //correct vertex reconstruction efficiency=========================================================================
            Int_t ntracks = 1;
            if ((TMath::Abs(Zvertex) <= 10.0)) {
              fHistZvtxMult->Fill(Zvertex,AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE));
             // fHistZvtxMultTrue->Fill(Zvertex,(AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE)/0.78)); // pp 5 TeV 0.78
              fHistZvtxMultTrue->Fill(Zvertex,(AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE)/0.79)); // pPb 5 TeV 0.79

         }
         if (TMath::Abs(Zvertex) <= 10.0 && lPrimaryVtx->GetNContributors() > ntracks ) { // (lPrimaryVtx->GetDispersion()<0.04 && lPrimaryVtx->GetZRes()<0.25)
              fHistZvtxMultNcon1->Fill(Zvertex,AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE));
             // fHistZvtxMultNcon1True->Fill(Zvertex,(AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE)/0.78)); // pp 5 TeV 0.78
              fHistZvtxMultNcon1True->Fill(Zvertex,(AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE)/0.79)); // pPb 5 TeV 0.79
            }
            if (TMath::Abs(Zvertex) <= 10.0 && lPrimaryVtx->GetNContributors() > 2 ) {
                 fHistZvtxMultNcon2->Fill(Zvertex,AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE));
            }


  Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
  AliESDtrack* fRecLeadTrack = 0;

  Float_t accepted_tracks = 0;
  Double_t flPt = 0;
  Double_t fsPt = 0;
  Double_t frPt = 0;
  Double_t fphi = 0;
  Double_t flPhi = 0;
  Double_t flEta = 0;
  Double_t trackingEff = 0;
    Double_t Etacms = 0;
  //basic measurements==========================================================================================================================================================

for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;                            // if we failed, skip this track
    if(!fCuts2->AcceptTrack(track)) continue;      // new cuts

         //DCA distribution------------------------------------------------------------------------------------------------------------------------------------------
         Float_t DCAxy[2],Cov_DCAxy[3];
         track->GetImpactParameters(DCAxy,Cov_DCAxy);
         fHistAllTracksDCAxy->Fill(DCAxy[0]);
         fHistAllTracksPt->Fill(track->Pt());

         if (fUseMC){
             AliMCParticle* particle = (AliMCParticle*)fMCEvent->GetTrack(i);
             if (fMCEvent->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))){
                 fHistPrimaryDCAxy->Fill(DCAxy[0]);
                 fHistPrimaryPt->Fill(track->Pt());
             }
             else if(fMCEvent->IsSecondaryFromWeakDecay(TMath::Abs(track->GetLabel()))){
                 fHistStrangeDCAxy->Fill(DCAxy[0]);              //The first mother is strange and it's a decay
                 fHistStrangePt->Fill(track->Pt());
             }
             else if (fMCEvent->IsSecondaryFromMaterial(TMath::Abs(track->GetLabel()))){
                 fHistMaterialDCAxy->Fill(DCAxy[0]);             // If a particle is not a physical primary, check if it comes from material
                 fHistMaterialPt->Fill(track->Pt());
             }
             else {
                    if (TMath::Abs(particle->PdgCode()) ==   11)  printf("Track  Electrons \n");
         	    else if (TMath::Abs(particle->PdgCode()) ==   13) printf("Track  Muon \n");
                    else if (TMath::Abs(particle->PdgCode()) ==   22) printf("Track  Photons \n");
         	    else if (TMath::Abs(particle->PdgCode()) ==  211) printf("Track  Pions \n");
         	    else if (TMath::Abs(particle->PdgCode()) ==  321) printf("Track  K+ \n");
         	    else if (TMath::Abs(particle->PdgCode()) == 2212) printf("Track  Protons \n");
         	    else printf("Track  Unknown \n");
                    fHistOtherDCAxy->Fill(DCAxy[0]);              // If a particle is not a physical primary, check if it comes from other process
             }
         }



    if(!fCuts->AcceptTrack(track)) continue;      // new cuts
    if (track->IsOn(AliESDtrack::kMultInV0)) continue;    // secondary particles from V0s
    if (track->IsOn(AliESDtrack::kMultSec) ) continue;   // secondary particles

    Double_t Sigma = 0.0050+0.0060/TMath::Power(track->Pt(),0.9); // pT dependent cut on track DCAxy to reject secondaries try ("0.0350+0.0420/pt^0.9")
    Double_t D0max = 7.*Sigma;   // can we make the cut even narrower to rejects secondaries based on IP?
    if (TMath::Abs(DCAxy[0]) > D0max) {
      continue;
    }

    accepted_tracks++;

    fHistDCAxy->Fill(DCAxy[0]);

// -------------------------tracking efficiency correction
if (fTrackingEfficiency) {
     trackingEff = fTrackingEfficiency->GetBinContent(fTrackingEfficiency->GetXaxis()->FindFixBin(track->Pt()));
//if (trackingEff == 0) {trackingEff = 0.000001;}

    fHistPt->Fill((track->Pt())/trackingEff);}
    //fHistEta->Fill((track->Eta())/trackingEff);}
else {fHistPt->Fill(track->Pt()); }
//cout << "eff= " << a << " pt=" << track->Pt() <<  " corrected pt = "<< track->Pt()/a << endl;
  Etacms = (track->Eta()-0.465);
    
   // fHistEta->Fill(track->Eta());//}lab frame eta
    
    fHistEta->Fill(Etacms);//
  
    fHistPhi->Fill(track->Phi());                 // plot the pt value of the track in a histogram
    fHistEta_Phi->Fill(track->Eta(), track->Phi());
    fHistEta_Pt->Fill(track->Pt(),track->Eta());

    Double_t vars1[3] = {track->Pt(),track->Eta(),track->Phi()};
    fHistTrack->Fill(vars1);


    fphi = track->Phi();

    //printf("Track  pt = %f\n",  track->Pt());

    frPt= track->Pt();
    fsPt +=frPt;

    if (track->Pt() > flPt){
      flPt = track->Pt();
      flPhi = track->Phi();
      flEta = Etacms;
      fRecLeadTrack = track;
    }

  }  // continue until all the tracks are processed

   fHistEvents->Fill(accepted_tracks);



   fRecNTracks = accepted_tracks;
   //cout<<"in TracksLoop RecN("<<fRecNTracks<<")"<<endl;
   fRecLeadPhi = flPhi;


     //$$$$$$$   The result is same with filling the histogram in UserExec()
     //$$$$$$$   Fill the histogram only when the code and go to here(fRecNTracks = 0 or fRecNTracks > 0 )
     //
     if (fGenNTracks > 0 && fRecNTracks > 0 ) {
          fHistMultResponseMat->Fill(fGenNTracks,fRecNTracks);

          Double_t fLeadDPhi = deltaPhi(fGenLeadPhi,fRecLeadPhi);
          fHistLeadDPhi->Fill(fLeadDPhi);
          if (TMath::Abs(fLeadDPhi)<(1./3)*TMath::Pi()) fHistLeadDPhiTo->Fill(fLeadDPhi);
          else if ((1./3)*TMath::Pi()<TMath::Abs(fLeadDPhi) && TMath::Abs(fLeadDPhi) <(2./3)*TMath::Pi()) fHistLeadDPhiTr->Fill(fLeadDPhi);
          else fHistLeadDPhiAw->Fill(fLeadDPhi);
          //cout<<"GenLeadPhi"<<"( "<<(fGenLeadPhi*(180.0/TMath::Pi()))<<" )---"<<"RecLeadPhi"<<"( "<<(fRecLeadPhi*(180.0/TMath::Pi()))<<" )"<<endl;
          //cout<<"LeadDPhi"<<"( "<<(fLeadDPhi*(180.0/TMath::Pi()))<<" )" <<endl;
     }

   if(accepted_tracks <1) return;
if (fTrackingEfficiency) {
trackingEff = fTrackingEfficiency->GetBinContent(fTrackingEfficiency->GetXaxis()->FindFixBin(flPt));
      fHistLeadPt->Fill(flPt/trackingEff);}
else {      fHistLeadPt->Fill(flPt); }
      fHistSumPt->Fill(fsPt);
      fHistAvgPt->Fill(fsPt/accepted_tracks);
      fHistTracks->Fill(accepted_tracks);

      //$$$$$$  When RecNTrack > 0, will fill the histogram, even though GenNTrack = 0 or don't exit.
      //fRecNTracks = accepted_tracks;
      //fHistMultResponseMat->Fill(fGenNTracks,fRecNTracks);


       if (fUseMC) {
           Int_t label = fRecLeadTrack->GetLabel();
           AliMCParticle* fpart = (AliMCParticle*)fMCEvent->GetTrack(TMath::Abs(label));
           if(fpart == fGenLeadPart) {
               fHistMatchLeadPt->Fill(fRecLeadTrack->Pt());
           }
           else fHistNoMatchLeadPt->Fill(fRecLeadTrack->Pt());
       }

      Double_t tr1 = accepted_tracks;
      Double_t vars1a[4] = {flPt,fsPt,fsPt/accepted_tracks,tr1};
      fHistLeadTrack->Fill(vars1a);

  //printf("Total tracks = %d and Accepted tracks = %d \n",  iTracks, accepted_tracks);
  //printf("Lead pt = %f and Sum Pt = %f and Average Pt = %f \n",  flPt, fsPt,fsPt/accepted_tracks);
  //printf("Lead Pt = %f and Lead Phi = %f and Lead Eta = %f \n",  flPt, flPhi,flEta);

   Int_t n1=0;
   Int_t n2=0;
   Int_t n3=0;

   Double_t fsPt1 = 0.0;
   Double_t fsPt2 = 0.0;
   Double_t fsPt3 = 0.0;
   Double_t dphi = 0.0;

  //particle density and energy density======================================================
  for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;                            // if we failed, skip this track
    if(!fCuts->AcceptTrack(track)) continue;      // new cuts
    if (track->IsOn(AliESDtrack::kMultInV0)) continue;    // secondary particles from V0s
    if (track->IsOn(AliESDtrack::kMultSec) ) continue;   // secondary particles

    Float_t DCAxy[2],Cov_DCAxy[3];
    track->GetImpactParameters(DCAxy,Cov_DCAxy);
    Double_t Sigma = 0.0050+0.0060/TMath::Power(track->Pt(),0.9); // pT dependent cut on track DCAxy to reject secondaries try ("0.0350+0.0420/pt^0.9")
    Double_t D0max = 7.*Sigma;   // can we make the cut even narrower to rejects secondaries based on IP?
    if (TMath::Abs(DCAxy[0]) > D0max) {
      continue;
    }

    Double_t dphi1 = (track->Phi()-flPhi);
    dphi = DeltaPhi(track->Phi(),flPhi);
    Double_t deta = track->Eta()-flEta;
    Double_t dr = DeltaR(track->Phi(),track->Eta(),flPhi,flEta);   // Delta R is the distance between Leading track and the assocaited track (see function definition)

    //printf("dphi = %f \n",  dphi*(180.0/TMath::Pi()));
     //  cout << " lead phi  = "  <<  flPhi << " track phi  = " << track->Phi() << " track eta  = " << track->Eta()  << " DistanceR  = " << dr  << "  DeltaPhi  = " << dphi  << " dPhi = " << dphi1  << endl;



    //if (TMath::Abs(dr) < 0.0001) continue;   // Protection to make sure we dont use Leading track itself for the correlation

     if(dphi1>((3./2)*TMath::Pi())) dphi1 = dphi1-2*TMath::Pi();
     else if(dphi1<-TMath::Pi()/2) dphi1=dphi1+2*TMath::Pi();

     fHistDPhiEta->Fill(dphi1,deta);


    //Double_t vars[3] = {TMath::Abs(deta),dphi,dr};
    Double_t vars[3] = {deta,dphi,dr};

    fHnDelta->Fill(vars,1);


    if (TMath::Abs(dphi)<(1./3)*TMath::Pi())
    {
      fsPt1 += track->Pt();
  //   cout << " Phi_Tow=" << track->Phi() << "(" << (track->Phi()*(180.0/TMath::Pi())) << ")" << " Lead_phi="  <<  flPhi<< "(" << (flPhi*(180.0/TMath::Pi())) << ")"   << " dPhi=" << dphi << " dR=" << dr   << " Nch=" << n1  << " lead_pt="  <<  flPt << " Sum_Pt=" << fsPt1  << endl;


      //Double_t vars1[3] = {TMath::Abs(deta),dphi,dr};
      Double_t vars1[3] = {deta,dphi,dr};
      fHnDeltaTo->Fill(vars1,1);
      n1++;
    }

    else if ((1./3)*TMath::Pi()<TMath::Abs(dphi) && TMath::Abs(dphi) <(2./3)*TMath::Pi())
    {
      fsPt2 += track->Pt();
  //  cout << " Phi_Tra=" << track->Phi() << "(" << (track->Phi()*(180.0/TMath::Pi())) << ")" << " Lead_phi="  <<  flPhi<< "(" << (flPhi*(180.0/TMath::Pi())) << ")"   << " dPhi=" << dphi << " dR=" << dr   << " Nch=" << n2  << " lead_pt="  <<  flPt << " Sum_Pt=" << fsPt2  << endl;


      //Double_t vars2[3] = {TMath::Abs(deta),dphi,dr};
      Double_t vars2[3] = {deta,dphi,dr};
      fHnDeltaTr->Fill(vars2,1);
      n2++;
    }

    else
    {
      fsPt3 +=track->Pt();
    //  cout << " Phi_Awa=" << track->Phi() << "(" << (track->Phi()*(180.0/TMath::Pi())) << ")" << " Lead_phi="  <<  flPhi<< "(" << (flPhi*(180.0/TMath::Pi())) << ")"   << " dPhi=" << dphi << " dR=" << dr   << " Nch=" << n3  << " lead_pt="  <<  flPt << " Sum_Pt=" << fsPt3  << endl;


      //Double_t vars3[3] = {TMath::Abs(deta),dphi,dr};
      Double_t vars3[3] = {deta,dphi,dr};
      fHnDeltaAw->Fill(vars3,1);
      n3++;
    }
  }
    if (flPt > fLeadPtCutMin)
    {
        //cout<<"Rec level:     Sum pT(Towds)="<<fsPt1<<"  Sum pT(Trans) = "<<fsPt2<<"  Sum pT(Away)="<<fsPt3<<"===================================================="<<endl;
        fHistMultTrans->Fill(n2);
        fHistMultTowds->Fill(n1);
        fHistMultAway->Fill(n3);
        fHistEngyTrans->Fill(fsPt2);
        fHistEngyTowds->Fill(fsPt1);
        fHistEngyAway->Fill(fsPt3);
    }

  if(fTrackingEfficiency) {
        trackingEff = fTrackingEfficiency->GetBinContent(fTrackingEfficiency->GetXaxis()->FindFixBin(flPt));

        n1=  n1/trackingEff;
        n2=  n2/trackingEff;
        n3=  n3/trackingEff;
      fsPt1 = fsPt1/trackingEff;
      fsPt2 = fsPt2/trackingEff;
      fsPt3 = fsPt3/trackingEff;

    }
  if ((n1+n2+n3) >= 1)  // Protection just to make sure if less than one track present in an event then do not fill histograms
  //if(flPt > 0)
  {
  //  cout  << "  Lead_pt="  <<  flPt << " dPhi=" << dphi << "(" << (dphi*(180.0/TMath::Pi())) << ")"  << " Nch_Tow=" << n1  << " Nch_Tra=" << n2 << " Nch_Awa=" << n3 << " SumPt_Tow=" << fsPt1 << " SumPt_Tra=" << fsPt2 << " SumPt_Awa=" << fsPt3   << endl;
    //cout  << "  Lead_pt="  <<  flPt << endl;
 // -------------------------applying tracking efficiency on observables ----------------------
if (fTrackingEfficiency) {
trackingEff = fTrackingEfficiency->GetBinContent(fTrackingEfficiency->GetXaxis()->FindFixBin(flPt));
    fHistParticleDensityTowards->Fill(flPt,n1/trackingEff);
    fHistParticleDensityTransverse->Fill(flPt,n2/trackingEff);
    fHistParticleDensityAway->Fill(flPt,n3/trackingEff);

    fHistEnergyDensityTowards->Fill(flPt,fsPt1/trackingEff);
    fHistEnergyDensityTransverse->Fill(flPt,fsPt2/trackingEff);
    fHistEnergyDensityAway->Fill(flPt,fsPt3/trackingEff); }
else {
    fHistParticleDensityTowards->Fill(flPt,n1);
    fHistParticleDensityTransverse->Fill(flPt,n2);
    fHistParticleDensityAway->Fill(flPt,n3);

    fHistEnergyDensityTowards->Fill(flPt,fsPt1);
    fHistEnergyDensityTransverse->Fill(flPt,fsPt2);
    fHistEnergyDensityAway->Fill(flPt,fsPt3);
    }
      if (fUseMC) {
           Int_t label = fRecLeadTrack->GetLabel();
           AliMCParticle* fpart = (AliMCParticle*)fMCEvent->GetTrack(TMath::Abs(label));
           if(fpart == fGenLeadPart) {
               fHistMatchParticleDensityTransverse->Fill(flPt,n2);
               fHistMatchEnergyDensityTransverse->Fill(flPt,fsPt2);
           }
      }

  }

//}//end loop
}// end function

//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsMinimumBias(AliVEvent *fESD)
// Function to check for minimum-bias trigger (AliVEvent::kMB)
{
  //Code to reject events that aren't kMB/INT7/INT1
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t IsTriggered = kFALSE;

  //    IsTriggered = (maskIsSelected & AliVEvent::kHighMultV0) == AliVEvent::kHighMultV0;   // CINT7 trigger
  //        IsTriggered = ((maskIsSelected & AliVEvent::kINT7) || (maskIsSelected & AliVEvent::kHighMultV0));   // CINT7 trigger
  IsTriggered = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;  // CINT7 trigger
  //  IsTriggered = (maskIsSelected & AliVEvent::kINT1) == AliVEvent::kINT1;  // CINT1 trigger
  //    IsTriggered = (maskIsSelected & AliVEvent::kINT10) == AliVEvent::kINT10;  // CINT10 trigger


  return IsTriggered;
}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsINELgtZERO(AliVEvent *fESD)
// Function to check for INEL > 0 condition, need atleast one SPD tracklet
{
  Bool_t returnValue = kFALSE;
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  if ( AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets, 1.0) >= 1 ) returnValue = kTRUE;

  return returnValue;
}

//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::HasNoPileupSPDInMultBins(AliVEvent *fESD)
// Checks if No pileup from SPD (via IsPileupFromSPDInMultBins)
{
  Bool_t returnValue = kTRUE;

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  if ( esdevent->IsPileupFromSPDInMultBins() == kFALSE ) returnValue = kFALSE;
  return returnValue;
}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::HasNoInconsistentSPDandTrackVertices(AliVEvent *fESD)
{
  //Accepts events which has consistent SPD and Global tracks Z Vertices agreementent within 5 mm
  Bool_t returnValue = kTRUE;

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  const AliESDVertex *lPrimaryVtxSPD    = NULL;
  const AliESDVertex *lPrimaryVtxTracks = NULL;

  lPrimaryVtxSPD    = esdevent->GetPrimaryVertexSPD   ();
  lPrimaryVtxTracks = esdevent->GetPrimaryVertexTracks();

  //Only continue if track vertex defined
  if( lPrimaryVtxTracks->GetStatus() && lPrimaryVtxSPD->GetStatus() ){
    //Copy-paste from refmult estimator
    // TODO value of displacement to be studied
    const Float_t maxDisplacement = 0.5;
    //check for displaced vertices
    Double_t displacement = TMath::Abs(lPrimaryVtxSPD->GetZ() - lPrimaryVtxTracks->GetZ());
    if (displacement > maxDisplacement) returnValue = kFALSE;
  }

  return returnValue;
}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::HasAcceptedVertexPosition(AliVEvent *fESD)
// Will accept events only if best primary vertex Z position awailable with |z| < 10cm
{
  Bool_t returnValue = kFALSE;
  const AliVVertex *lPrimaryVtx = NULL;

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  lPrimaryVtx = esdevent->GetPrimaryVertex();

  if ( TMath::Abs( lPrimaryVtx->GetZ() ) <= 10.0 ) returnValue = kTRUE;
  return returnValue;
}
//----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsSPDClusterVsTrackletBG(AliVEvent* fESD)
{
  // rejects BG background based on the cluster vs tracklet correlation
  // returns true if the event is BG slope is 4 and intercept is 65
  const AliVMultiplicity* mult = fESD->GetMultiplicity();
  // if (!mult) { AliFatal("No multiplicity object"); return 0;}
  Int_t ntracklet   = mult->GetNumberOfTracklets();
  Int_t spdClusters = fESD->GetNumberOfITSClusters(0) + fESD->GetNumberOfITSClusters(1);
  //return spdClusters < Float_t(fASPDCvsTCut) + Float_t(ntracklet)*fBSPDCvsTCut;
  return spdClusters < (65.0 + (ntracklet)*4.0);

}
//----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsOutOfBunchPileupFastor(AliVEvent* fESD)
{
  Double_t SPD_Online = (fESD->GetMultiplicity()->GetFastOrFiredChipMap().CountBits(400)+ fESD->GetMultiplicity()->GetFastOrFiredChipMap().CountBits(800));
  Double_t SPD_Offline = (fESD->GetMultiplicity()->GetFiredChipMap().CountBits(400) + fESD->GetMultiplicity()->GetFiredChipMap().CountBits(800));

  return SPD_Online >= (-20.589 + 0.73664*(SPD_Offline));
}
//----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsOutOfBunchPileup(AliVEvent* fESD)
{
  // rejects Outofbunch Pileup events based on the Interaction logic

  TBits fIR11 = fESD->GetHeader()->GetIRInt1InteractionMap();
  //TBits fIR22 = fESD->GetHeader()->GetIRInt2InteractionMap();

  Bool_t isOutOfBunchPileup = 0;
  for (Int_t i=1;i<=3;i++) { isOutOfBunchPileup|=fIR11.TestBitNumber(90-i);}
  for (Int_t i=8;i<=11;i++) { isOutOfBunchPileup|=fIR11.TestBitNumber(90+i);}

  //if(isOutOfBunchPileup == 1) return;

  //for (Int_t i=0;i<=7;i++) isOutOfBunchPileup|=fIR22.TestBitNumber(90-i);
  //for (Int_t i=4;i<=7;i++) isOutOfBunchPileup|=fIR22.TestBitNumber(90+i);

  return isOutOfBunchPileup;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsInCompleteEvent(AliVEvent* fESD)
{
  // rejects Incomplete events in the run
  Bool_t isIncomplete = kFALSE;
  if (fESD->IsIncompleteDAQ()) isIncomplete = kTRUE;

  return isIncomplete;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::V0Decision(AliVEvent* fESD)
{
  // rejects Incomplete events in the run
  AliVVZERO* vzero_decision = fESD->GetVZEROData();

  Bool_t isV0Decision = kFALSE;
  isV0Decision = ((vzero_decision->GetV0ADecision()==1) && (vzero_decision->GetV0CDecision()==1));

  return isV0Decision;
}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::V0Asymmetry(AliVEvent* fESD)
{
  // rejects Incomplete events in the run
  Bool_t isEventSelected_V0Asym = kTRUE;

  AliVVZERO* vzero_asym = fESD->GetVZEROData();

  Double_t v0c012 = vzero_asym->GetMRingV0C(0) + vzero_asym->GetMRingV0C(1) + vzero_asym->GetMRingV0C(2);
  Double_t v0c3   = vzero_asym->GetMRingV0C(3);

  isEventSelected_V0Asym &= vzero_asym->GetMTotV0C() < (330. + 100. * TMath::Power(vzero_asym->GetMTotV0A(), .2));
  isEventSelected_V0Asym &= (v0c012 < 160.) || (v0c3 > 12.*TMath::Power(.01*(v0c012 - 160.), 1.7));

  return isEventSelected_V0Asym;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::HasBCMod4(AliVEvent* fESD)
{
  // Bunch Crossing mod == 4

  Bool_t isBCMod4 = kFALSE;
  Int_t fBC = fESD->GetBunchCrossNumber();
  Int_t bcmod4 = fBC%4;

  if (bcmod4 == 2) isBCMod4 = kTRUE;

  //printf("fBC = %d, bcmod4 = %d, isBCMod4 = %d \n",fBC,bcmod4,isBCMod4);
  return isBCMod4;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskMpiUE::IsEventSelected(AliVEvent *fESD)

{
  Bool_t returnValue = kFALSE;

  if (
   // HasNoPileupSPDInMultBins                 ( fESD ) == kTRUE &&
    //IsINELgtZERO                             ( fESD ) == kTRUE &&
//    HasAcceptedVertexPosition                ( fESD ) == kTRUE  &&
    //HasNoInconsistentSPDandTrackVertices     ( fESD ) == kTRUE &&
    //IsSPDClusterVsTrackletBG                 ( fESD ) == kTRUE &&
    //IsOutOfBunchPileup                       ( fESD ) == kFALSE &&
    //IsOutOfBunchPileupFastor                 ( fESD ) == kTRUE &&
    //IsInCompleteEvent                        ( fESD ) == kFALSE &&
    //V0Decision                               ( fESD ) == kTRUE &&
    //V0Asymmetry                              ( fESD ) == kTRUE &&
    //    HasBCMod4                                 ( fESD ) == kTRUE &&
    IsMinimumBias                             ( fESD ) == kTRUE
  ) returnValue = kTRUE;
  return returnValue;
}

//-----------------------------------------------------------------------------------------
Double_t AliAnalysisTaskMpiUE::DeltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t dphi = TMath::Abs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
  return dphi;
}

//---------------------------------------------------------------------------------------
Double_t AliAnalysisTaskMpiUE::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t dphi = phi1 - phi2;
  if (dphi > TMath::Pi()) dphi = dphi - 2. * TMath::Pi();
  else if (dphi < -TMath::Pi()) dphi = dphi + 2. * TMath::Pi();
  return dphi;
}

//---------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskMpiUE::DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2)
{
  Double_t dphi = DeltaPhi(phi1, phi2);
  Double_t deta = eta1 - eta2;
  return (TMath::Sqrt(dphi * dphi + deta * deta));

}
