 /*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock                                                *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskGammaTriggerQA.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliCaloTrackMatcher.h"
#include <vector>
#include <map>
#include <fstream>

ClassImp(AliAnalysisTaskGammaTriggerQA)

//________________________________________________________________________
AliAnalysisTaskGammaTriggerQA::AliAnalysisTaskGammaTriggerQA(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fOutputContainer(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoCent(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGammaCandidatesBasic(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fHistoNV0Trigger(NULL),
  fHistoNV0TriggerTracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fTreeList(NULL),
  fTreeTriggInfo(NULL),
  fCent(0),
  fT0Trigg(0),
  fV0Mult(0),
  fV0Trigg(0),
  fTPCMult(0),
  fSPDHit(0),
  fSPDTracklet(0),
  fZVertex(0),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fQADetailed(0),
  fIsMC(0),
  fWeightJetJetMC(1),
  fNCurrentClusterBasic(0)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaTriggerQA::AliAnalysisTaskGammaTriggerQA(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fOutputContainer(0),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoCent(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGammaCandidatesBasic(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fHistoNV0Trigger(NULL),
  fHistoNV0TriggerTracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fTreeList(NULL),
  fTreeTriggInfo(NULL),
  fCent(0),
  fT0Trigg(0),
  fV0Mult(0),
  fV0Trigg(0),
  fTPCMult(0),
  fSPDHit(0),
  fSPDTracklet(0),
  fZVertex(0),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fQADetailed(0),
  fIsMC(0),
  fWeightJetJetMC(1),
  fNCurrentClusterBasic(0)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaTriggerQA::~AliAnalysisTaskGammaTriggerQA()
{
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaTriggerQA::UserCreateOutputObjects(){

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  // set common binning in pT for mesons and photons
  Int_t nBinsPt               = 250;
  Float_t minPt               = 0;
  Float_t maxPt               = 25;
  Int_t nBinsQAPt             = 175;
  Float_t maxQAPt             = 25;
  Int_t nBinsClusterPt        = 500;
  Float_t minClusterPt        = 0;
  Float_t maxClusterPt        = 50;
  Double_t *arrPtBinning      = new Double_t[1200];
  Double_t *arrQAPtBinning    = new Double_t[1200];
  Double_t *arrClusPtBinning  = new Double_t[1200];
  for(Int_t i=0; i<nBinsPt+1;i++){
    arrPtBinning[i]         = ((maxPt-minPt)/nBinsPt)*i;
  }
  for(Int_t i=0; i<nBinsClusterPt+1;i++){
    arrClusPtBinning[i]     = ((maxClusterPt-minClusterPt)/nBinsClusterPt)*i;
  }
  for(Int_t i=0; i<nBinsQAPt+1;i++){
    if(i<60) arrQAPtBinning[i]              = 0.05*i;
    else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
    else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
    else if(i<175) arrQAPtBinning[i]        = 20.+1.0*(i-170);
    else arrQAPtBinning[i]                  = maxQAPt;
  }

  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer  = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer  = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  // Array of current cut's gammas
  fClusterCandidates  = new TList();
  fClusterCandidates->SetOwner(kTRUE);

  fCutFolder          = new TList*[fnCuts];
  fESDList            = new TList*[fnCuts];
  fHistoNEvents       = new TH1F*[fnCuts];
  if(fIsMC > 1){
    fHistoNEventsWOWeight   = new TH1F*[fnCuts];
  }
  if(fIsMC == 2){
    fProfileJetJetXSection  = new TProfile*[fnCuts];
    fHistoJetJetNTrials     = new TH1F*[fnCuts];
  }

  if (fQADetailed > 0){
    fTreeList                 = new TList*[fnCuts];
    fTreeTriggInfo            = new TTree*[fnCuts];
  }

  fHistoNGoodESDTracks        = new TH1F*[fnCuts];
  fHistoCent                  = new TH1F*[fnCuts];
  fHistoVertexZ               = new TH1F*[fnCuts];
  fHistoNGammaCandidates      = new TH1F*[fnCuts];
  fHistoNGammaCandidatesBasic = new TH1F*[fnCuts];
  if(!fDoLightOutput){
    fHistoNGoodESDTracksVsNGammaCandidates  = new TH2F*[fnCuts];
    fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
    fHistoNV0Tracks                         = new TH1F*[fnCuts];
    fHistoNV0Trigger                        = new TH1F*[fnCuts];
    fHistoNV0TriggerTracks                  = new TH2F*[fnCuts];
  }
  if(fIsHeavyIon==2) fProfileEtaShift          = new TProfile*[fnCuts];

  fHistoClusGammaPt                 = new TH1F*[fnCuts];
  fHistoClusGammaE                  = new TH1F*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();

    fCutFolder[iCut]        = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s", cutstringEvent.Data(), cutstringCalo.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]          = new TList();
    fESDList[iCut]->SetName(Form("%s_%s ESD histograms", cutstringEvent.Data(), cutstringCalo.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    fHistoNEvents[iCut]     = new TH1F("NEvents", "NEvents", 14, -0.5, 13.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
      TString TriggerNames  = "Not Trigger: ";
      TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"SPD Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut] = new TH1F("NEventsWOWeight", "NEventsWOWeight", 14, -0.5, 13.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames  = "Not Trigger: ";
        TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      } else {
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);
    }
    if (fIsMC == 2){
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 4000, 0, 4000);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 400, 0, 400);
    else
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 200, 0, 200);
    fHistoNGoodESDTracks[iCut]->GetXaxis()->SetTitle("#primary tracks");
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoCent[iCut]             = new TH1F("Centrality", "Centrality; Centrality [%]", 100, 0, 100);
    fESDList[iCut]->Add(fHistoCent[iCut]);

    fHistoVertexZ[iCut]             = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
    fHistoVertexZ[iCut]->GetXaxis()->SetTitle("Z_{vtx} (cm)");
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 600, 0, 600);
    else if(fIsHeavyIon == 2)
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 400, 0, 400);
    else
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 100, 0, 100);
    fHistoNGammaCandidatesBasic[iCut]->GetXaxis()->SetTitle("#cluster candidates basic");
    fESDList[iCut]->Add(fHistoNGammaCandidatesBasic[iCut]);


    if(fIsHeavyIon == 1)
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 200, 0, 200);
    else if(fIsHeavyIon == 2)
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 100, 0, 100);
    else
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 50, 0, 50);
    fHistoNGammaCandidates[iCut]->GetXaxis()->SetTitle("#cluster candidates with current cut");
    fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);

    if(!fDoLightOutput){
      if(fIsHeavyIon == 1)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 4000, 0, 4000, 200, 0, 200);
      else if(fIsHeavyIon == 2)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 400, 0, 400, 100, 0, 100);
      else
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 200, 0, 200, 50, 0, 50);
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetXTitle("#good tracks");
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetYTitle("#cluster candidates");
      fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);

      if(fIsHeavyIon == 1)
        fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 500, 0, 1000, 1000, 0, 4000);
      else if(fIsHeavyIon == 2)
        fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 200, 0, 400, 500, 0, 2000);
      else
        fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 100, 0, 200, 250, 0, 1000);
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 40000, 0, 40000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 2500, 0, 2500);
      else
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 1500, 0, 1500);
      fHistoNV0Tracks[iCut]->SetXTitle("V0 amplitude");
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Trigger[iCut]       = new TH1F("V0 Trigger", "V0 Trigger", 40000, 0, 40000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Trigger[iCut]       = new TH1F("V0 Trigger", "V0 Trigger", 2500, 0, 2500);
      else
        fHistoNV0Trigger[iCut]       = new TH1F("V0 Trigger", "V0 Trigger", 1500, 0, 1500);
      fHistoNV0Trigger[iCut]->SetXTitle("V0 trigger amplitude");
      fESDList[iCut]->Add(fHistoNV0Trigger[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0TriggerTracks[iCut]        = new TH2F("V0 Trigger vs Mult", "V0 Trigger vs Mult", 1000, 0, 40000, 1000, 0, 40000);
      else if(fIsHeavyIon == 2)
        fHistoNV0TriggerTracks[iCut]        = new TH2F("V0 Trigger vs Mult", "V0 Trigger vs Mult", 250, 0, 2500, 250, 0, 2500);
      else
        fHistoNV0TriggerTracks[iCut]        = new TH2F("V0 Trigger vs Mult", "V0 Trigger vs Mult", 150, 0, 1500, 150, 0, 1500);
      fESDList[iCut]->Add(fHistoNV0TriggerTracks[iCut]);


    }

    if(fIsHeavyIon==2) {
      fProfileEtaShift[iCut]        = new TProfile("Eta Shift", "Eta Shift", 1, -0.5, 0.5);
      fProfileEtaShift[iCut]->SetXTitle("#eta shift");
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoNGammaCandidates[iCut]->Sumw2();
      fHistoNGammaCandidatesBasic[iCut]->Sumw2();
      fHistoCent[iCut]->Sumw2();
      if(!fDoLightOutput){
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
        fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
        fHistoNV0Tracks[iCut]->Sumw2();
      }
      if(fIsHeavyIon==2) fProfileEtaShift[iCut]->Sumw2();
    }

    fHistoClusGammaPt[iCut]               = new TH1F("ClusGamma_Pt", "ClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
    fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusGammaE[iCut]               = new TH1F("ClusGamma_E", "ClusGamma_E", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaPt[iCut]->SetXTitle("E_{clus} (GeV/c)");
    fESDList[iCut]->Add(fHistoClusGammaE[iCut]);

    if (fIsMC > 1){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusGammaE[iCut]->Sumw2();
    }

    if (fQADetailed > 0){
      fTreeList[iCut]        = new TList();
      fTreeList[iCut]->SetName(Form("%s_%s TriggerQA tree", cutstringEvent.Data(), cutstringCalo.Data()));
      fTreeList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTreeList[iCut]);

      fTreeTriggInfo[iCut]= new TTree("TriggerInfoTree", "TriggerInfoTree");
      fTreeTriggInfo[iCut]->Branch("Cent",&fCent,"fCent/F");
      fTreeTriggInfo[iCut]->Branch("T0Trigg",&fT0Trigg,"fT0Trigg/s");
      fTreeTriggInfo[iCut]->Branch("V0Mult",&fV0Mult,"fV0Mult/i");
      fTreeTriggInfo[iCut]->Branch("V0Trigg",&fV0Trigg,"fV0Trigg/i");
      fTreeTriggInfo[iCut]->Branch("TPCMult",&fTPCMult,"fTPCMult/i");
      fTreeTriggInfo[iCut]->Branch("SPDTrack",&fSPDHit,"fSPDHit/i");
      fTreeTriggInfo[iCut]->Branch("SPDTracklet",&fSPDTracklet,"fSPDTracklet/i");
      fTreeTriggInfo[iCut]->Branch("ZVertex",&fZVertex,"fZVertex/F");
      fTreeList[iCut]->Add(fTreeTriggInfo[iCut]);
    }
  }


  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
    AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i",iMatcherTask)));
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
    }
  }

  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaTriggerQA::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }
    if(fIsHeavyIon==2) {
      if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        continue; // No Eta Shift requested, continue
      }
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
        continue;
      }
      else{
        printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
            (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      }
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaTriggerQA::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  fInputEvent           = InputEvent();
  if(fIsMC> 0) fMCEvent = MCEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }

  // ------------------- BeginEvent ----------------------------

  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){

    fiCut = iCut;

    fNCurrentClusterBasic       = 0;
    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

    if(fIsMC==2){
      Float_t xsection      = -1.;
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials, fInputEvent );
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}", ntrials);
    }

    if (fIsMC > 0){
      fWeightJetJetMC       = 1;
      Float_t pthard = -1;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent);
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if (!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }

    Bool_t triggered = kTRUE;
    if(eventNotAccepted){
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC>0){
        triggered = kFALSE;
      } else {
        continue;
      }
    }

    if(eventQuality != 0){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
      continue;
    }
    if (triggered == kTRUE) {
      // set variables
      fCent       = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCentrality(fInputEvent);
      if(fInputEvent->IsA()==AliESDEvent::Class()){
        fT0Trigg    = ((AliESDEvent*)fInputEvent)->GetT0Trig();
      } else {
        fT0Trigg    = 0;
      }

      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2){
        fV0Mult     = fInputEvent->GetVZEROData()->GetMTotV0A();
        fV0Trigg    = fInputEvent->GetVZEROData()->GetTriggerChargeA();
      } else{
        fV0Mult     = fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C();
        fV0Trigg    = fInputEvent->GetVZEROData()->GetTriggerChargeA()+fInputEvent->GetVZEROData()->GetTriggerChargeC();
      }
      fTPCMult      = fV0Reader->GetNumberOfPrimaryTracks();
      fSPDTracklet  = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
      fSPDHit       = fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1);
      fZVertex      = fInputEvent->GetPrimaryVertex()->GetZ();

      // fill histograms
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
      fHistoCent[iCut]->Fill(fCent, fWeightJetJetMC);
      fHistoNGoodESDTracks[iCut]->Fill(fTPCMult, fWeightJetJetMC);
      fHistoVertexZ[iCut]->Fill(fZVertex, fWeightJetJetMC);
      if(!fDoLightOutput){
        fHistoSPDClusterTrackletBackground[iCut]->Fill(fSPDTracklet, fSPDHit, fWeightJetJetMC);
        fHistoNV0Tracks[iCut]->Fill(fV0Mult, fWeightJetJetMC);
        fHistoNV0Trigger[iCut]->Fill(fV0Trigg, fWeightJetJetMC);
        fHistoNV0TriggerTracks[iCut]->Fill(fV0Trigg, fV0Mult, fWeightJetJetMC);
      }
      // fill tree
      if (fQADetailed > 0){
          fTreeTriggInfo[iCut]->Fill();
      }
    }

    if(fIsMC> 0){
      // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                                                                 ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                                                                 fMCEvent);
        }
        else if(fInputEvent->IsA()==AliAODEvent::Class()){
          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                                                                 ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                                                                 fInputEvent);
        }
      }
    }


    if (triggered==kFALSE) continue;

    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    ProcessClusters();            // process calo clusters
    fHistoNGammaCandidatesBasic[iCut]->Fill(fNCurrentClusterBasic, fWeightJetJetMC);
    fHistoNGammaCandidates[iCut]->Fill(fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(!fDoLightOutput) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries(), fWeightJetJetMC);

    fClusterCandidates->Clear(); // delete cluster candidates
  }
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaTriggerQA::ProcessClusters()
{

  Int_t nclus                       = 0;
  TClonesArray * arrClustersProcess = NULL;
  fNCurrentClusterBasic             = 0;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaTriggerQA! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }

  if(nclus == 0)  return;
  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersProcess)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }

    if(!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){
      if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsAcceptedForBasicCounting())fNCurrentClusterBasic++;
      delete clus;
      continue;
    }
    fNCurrentClusterBasic++;

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef(i);
    PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent));
    fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
    fHistoClusGammaE[fiCut]->Fill(PhotonCandidate->E(), fWeightJetJetMC);

    fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

    delete clus;
    delete tmpvec;
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaTriggerQA::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}
