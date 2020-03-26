/*************************************************************************
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

/* Class AliAnalysisTaskEmcalJetBtagSV:                               *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/*
 Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it,
 svallero@to.infn.it, s.lapointe@cern.ch
 ycorrale@cern.ch
 */

//--Root--
#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>
#include <THn.h>
#include <THnSparse.h>
#include <TList.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TRandom3.h>

//--AliRoot--
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliGenCocktailEventHeader.h" //FK//
#include "AliGenPythiaEventHeader.h"//FK//
#include "AliMCEvent.h" //FK// 
#include "AliExternalTrackParam.h" //AID//
#include "AliFJWrapper.h"  //EMB_clus

//--AliHFJetsClass--
#include "AliHFJetsTaggingVertex.h"
#include "AliHFJetsContainerVertex.h"

#include "AliAnalysisTaskEmcalJetBtagSV.h"

//_____________________________________________________________________________________

ClassImp(AliAnalysisTaskEmcalJetBtagSV)

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV() :
  AliAnalysisTaskEmcalJet("AnalysisTaskEmcalJetBtagSV", kTRUE),
  fCorrMode(kFALSE),
  fDoBkgRej(kTRUE),
  fDoRndmCone(kFALSE),
  fDoQAVtx(kFALSE),
  fDoFillV0Trks(kFALSE),
  fDoDetRespMtx(kFALSE),
  fDoOnlyMtxAna(kFALSE),
  fUseTriggerData(kFALSE),
  fEmbeddPerpendicular(kFALSE), //EMB_clus
  fRecJetsBranch(),
  fGenJetsBranch(),
  fPtHardName("pthard"),
  fGenNamePattern(""),
  fJetContName(""),
  fTrkContName(""),
  fRhoTaskName(""),
  fMCJetContName(""),
  fMCTrkContName(""),
  fMCRhoTaskName(""),
  fTaggingRadius(0.4),
  fSigmaSVCut(0.04),      //newDeltaPt//
  fMCWeight(1.),
  fMCXsec(0.),
  fMCAvgTrials(0.),
  fZNApercentile(0.),
  fCurrFileName(""),
  fCheckMCCrossSection(kFALSE),
  fSkipWeightInfo(kFALSE),
  fUseWeight(kFALSE),
  fInitialized(kFALSE),
  fOutputList(NULL),
  fhJetVtxSim(NULL),
  fhJetVtxData(NULL),
  fhQaVtx(NULL),
  fhEntries(NULL),
  fhZNApercentQa(NULL),
  fhEvtRej(NULL),
  fhEvtRejBitmap(NULL),
  fhHFjetQa(NULL),
  fhRhoQa(NULL),
  fhMCRhoQa(NULL),
  fhDeltaPt(NULL),
  fhDeltaPtLxy5(NULL), //newDeltaPt//
  fhDeltaPtLxy6(NULL), //newDeltaPt//
  fhDeltaPtLxy7(NULL), //newDeltaPt// 
  fhDeltaPtTrack10(NULL), //newDeltaPt//
  fZVertex(NULL),  //AID//
  fhTrackEta(NULL), //AID//
  fhTrackPhi(NULL), //AID//
  fhJetEta(NULL), //AID//
  fhJetPhi(NULL), //AID//
  fhOneOverPtVsPhiNeg(NULL), //AID//
  fhOneOverPtVsPhiPos(NULL), //AID//
  fhSigmaPtOverPtVsPt(NULL), //AID//
  fhDCAinXVsPt(NULL), //AID//
  fhDCAinYVsPt(NULL), //AID//
  fhDCAinXVsPtPhysPrimary(NULL), //AID//
  fhDCAinYVsPtPhysPrimary(NULL), //AID//
  fhDCAinXVsPtSecondary(NULL), //AID//
  fhDCAinYVsPtSecondary(NULL), //AID//
  fhFractionOfSecInJet(NULL), //AID//
  fhPtTrkTruePrimRec(0x0),
  fhPtTrkTruePrimGen(0x0),
  fhPtTrkSecOrFakeRec(0x0),
  fhnDetRespMtx(NULL),
  fhnGenerated(NULL),
  fhXsec(NULL),
  fhTrials(NULL),
  fEvent(NULL),
  fMCHeader(NULL),
  fMultSelection(NULL),
  fTagger(NULL),
  fCutsHFjets(NULL),
  fAnalysisUtils(NULL),
  fMCTracksCont(NULL),
  fRecJetArray(NULL),
  fRecTrkArray(NULL),
  fMCJetArray(NULL),
  fMCPartArray(NULL),
  fHFvertexing(NULL),
  fV0gTrkMap(NULL),
  fRandom(new TRandom3(0)),
  fGlLogLevel(AliLog::kError),
  fLcDebLevel(1),
  fStartBin(0),
  fMaxFacPtHard(0),  //FK
  fPtCut(0.15),     //AID//
  fEtaCut(0.9),      //AID//
  fDoEmbedding(kFALSE),     //EMB
  fHybridJetContName(""),   //EMB
  fHybridJetCont(0x0),      //EMB
  fhDeltaPtEmbedd(0x0),      //EMB
  fhDeltaPtEmbeddCorrelation(0x0),  //EMB
  fFastJetWrapper(NULL), //EMB_clus
  fTrackGenerator(NULL) //EMB_clus
{
  // default constructor
}

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV(const char* name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCorrMode(kFALSE),
  fDoBkgRej(kTRUE),
  fDoRndmCone(kFALSE),
  fDoQAVtx(kFALSE),
  fDoFillV0Trks(kFALSE),
  fDoDetRespMtx(kFALSE),
  fDoOnlyMtxAna(kFALSE),
  fUseTriggerData(kFALSE),
  fEmbeddPerpendicular(kFALSE), //EMB_clus
  fRecJetsBranch(),
  fGenJetsBranch(),
  fPtHardName("pthard"),
  fGenNamePattern(""),
  fJetContName(""),
  fTrkContName(""),
  fRhoTaskName(""),
  fMCJetContName(""),
  fMCTrkContName(""),
  fMCRhoTaskName(""),
  fTaggingRadius(0.4),
  fSigmaSVCut(0.04),      //newDeltaPt//
  fMCWeight(1.),
  fMCXsec(0.),
  fMCAvgTrials(0.),
  fZNApercentile(0.),
  fCurrFileName(""),
  fCheckMCCrossSection(kFALSE),
  fSkipWeightInfo(kFALSE),
  fUseWeight(kFALSE),
  fInitialized(kFALSE),
  fOutputList(NULL),
  fhJetVtxSim(NULL),
  fhJetVtxData(NULL),
  fhQaVtx(NULL),
  fhEntries(NULL),
  fhZNApercentQa(NULL),
  fhEvtRej(NULL),
  fhEvtRejBitmap(NULL),
  fhHFjetQa(NULL),
  fhRhoQa(NULL),
  fhMCRhoQa(NULL),  
  fhDeltaPt(NULL),
  fhDeltaPtLxy5(NULL), //newDeltaPt//
  fhDeltaPtLxy6(NULL), //newDeltaPt//
  fhDeltaPtLxy7(NULL), //newDeltaPt//
  fhDeltaPtTrack10(NULL), //newDeltaPt//
  fZVertex(NULL),  //AID//
  fhTrackEta(NULL), //AID//
  fhTrackPhi(NULL), //AID//
  fhJetEta(NULL), //AID//
  fhJetPhi(NULL), //AID//
  fhOneOverPtVsPhiNeg(NULL), //AID//
  fhOneOverPtVsPhiPos(NULL), //AID//
  fhSigmaPtOverPtVsPt(NULL), //AID//
  fhDCAinXVsPt(NULL), //AID//
  fhDCAinYVsPt(NULL), //AID//
  fhDCAinXVsPtPhysPrimary(NULL), //AID//
  fhDCAinYVsPtPhysPrimary(NULL), //AID//
  fhDCAinXVsPtSecondary(NULL), //AID//
  fhDCAinYVsPtSecondary(NULL), //AID//
  fhFractionOfSecInJet(NULL), //AID//
  fhPtTrkTruePrimRec(0x0),
  fhPtTrkTruePrimGen(0x0),
  fhPtTrkSecOrFakeRec(0x0),
  fhnDetRespMtx(NULL),
  fhnGenerated(NULL),
  fhXsec(NULL),
  fhTrials(NULL),
  fEvent(NULL),
  fMCHeader(NULL),
  fMultSelection(NULL),
  fTagger(NULL),
  fCutsHFjets(NULL),
  fAnalysisUtils(NULL),
  fMCTracksCont(NULL),
  fRecJetArray(NULL),
  fRecTrkArray(NULL),
  fMCJetArray(NULL),
  fMCPartArray(NULL),
  fHFvertexing(NULL),
  fV0gTrkMap(NULL),
  fRandom(new TRandom3(0)),
  fGlLogLevel(AliLog::kError),
  fLcDebLevel(1),
  fStartBin(0),
  fMaxFacPtHard(0), //FK
  fPtCut(0.15),     //AID//
  fEtaCut(0.9),      //AID//
  fDoEmbedding(kFALSE),     //EMB
  fHybridJetContName(""),   //EMB
  fHybridJetCont(0x0),      //EMB
  fhDeltaPtEmbedd(0x0),      //EMB
  fhDeltaPtEmbeddCorrelation(0x0),  //EMB
  fhDeltaPtEmbeddPerpendicular(0x0),      //EMB_clus
  fhDeltaPtEmbeddCorrelationPerpendicular(0x0),  //EMB_clus
  fFastJetWrapper(NULL), //EMB_clus
  fTrackGenerator(NULL) //EMB_clus
{
  // standard constructor
  AliInfo(MSGINFO("+++ Executing Constructor +++"));

  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::~AliAnalysisTaskEmcalJetBtagSV()
{
  // destructor
  AliInfo(MSGINFO("+++ Executing Destructor +++"));

  // Do not delete outputs in proof mode or merging will fail
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if (fOutputList)  delete fOutputList;
    if (fHFvertexing) delete fHFvertexing;
    if (fV0gTrkMap)   delete fV0gTrkMap;
  }
  
  if (fTagger)     delete fTagger;
  if (fCutsHFjets) delete fCutsHFjets;
  if (fEmbeddPerpendicular) { 		 //EMB_clus
		delete fFastJetWrapper;  //EMB_clus
		delete fTrackGenerator;
  }
  if (fRandom) delete fRandom;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserCreateOutputObjects()
{
  AliInfo(MSGINFO("+++ Executing UserCreateOutputObjects +++"));

  // Initialize output list of containers
  if (fOutputList != NULL) {
    delete fOutputList; fOutputList = NULL;
  }

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  // vertices within the jet - Correction Mode (MC)
  if (fCorrMode) {
    if (!fDoOnlyMtxAna) {
      fhJetVtxSim = new AliHFJetsContainerVertex("kJetVtxSim", AliHFJetsContainerVertex::kJetVtxSim);
      fOutputList->Add(fhJetVtxSim);
    }

    if (fDoDetRespMtx) {
      // detector response matrix for unfolding (from Gyulnara)
      // dimensions: pt_reco, pt_gen, eta_reco, eta_gen, flavor{g=1, L=2, C=3, B=4} BH and BP
      const int kNbins = 6;
      Int_t bins[kNbins]    = {300, 300,  20, 20,    5,   5};
      Double_t xmin[kNbins] = {  0,   0, -1., -1., -.5, -.5};
      Double_t xmax[kNbins] = {300, 300,  1.,  1., 4.5, 4.5};
      fhnDetRespMtx = new THnSparseF("fhnDetRespMtx", "Detector response matrix", kNbins, bins, xmin, xmax);
      fOutputList->Add(fhnDetRespMtx);
      
      // MC generated histogram is needed to calculate efficiency during unfolding
      // dimensions: pt_gen, eta_gen, flavor{g=1, L=2, C=3, B=4} BH and BP
      const Int_t kNhbins = 4;
      Int_t binsh[kNhbins]  =   {300,  20,    5,   5};
      Double_t xminh[kNhbins] = {  0,  -1., -.5, -.5};
      Double_t xmaxh[kNhbins] = {300,   1., 4.5, 4.5};
      fhnGenerated = new THnF("fhnGenerated", "MC Generated histogram", kNhbins, binsh, xminh, xmaxh);
      fOutputList->Add(fhnGenerated);
    }
  } else {
    // vertices within the jet - Data
    fhJetVtxData = new AliHFJetsContainerVertex("kJetVtxData", AliHFJetsContainerVertex::kJetVtxData);
    fOutputList->Add(fhJetVtxData);
  }
  // vertices QA
  if (fDoQAVtx) {
    // vertices within the jet - QA
    fhQaVtx = new AliHFJetsContainerVertex("kQaVtx", AliHFJetsContainerVertex::kQaVtx);
    fOutputList->Add(fhQaVtx);
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++
  // track and event QA histogram //AID//
  fZVertex = new TH1F ("fZVertex","fZVertex",120, -30,30);    // Z vertex distribuition //AID// 
  fZVertex->Sumw2();
  fhTrackEta = new TH2F ("fhTrackEta","fhTrackEta",100, 0,100, 20, -1,1);   // eta track distribuition //AID//  
  fhTrackEta->Sumw2(); 
  fhTrackPhi = new TH2F ("fhTrackPhi","fhTrackPhi",100, 0,100, 60, 0,TMath::TwoPi());// phi track distribuition //AID// 
  fhTrackPhi->Sumw2(); 
  fhJetEta = new TH2F ("fhJetEta","fhJetEta",100, 0,100, 20, -1,1);   // eta track distribuition //AID//    
  fhJetEta->Sumw2(); 
  fhJetPhi = new TH2F ("fhJetPhi","fhJetPhi",100, 0,100, 60, 0,TMath::TwoPi());// phi track distribuition //AID//   
  fhJetPhi->Sumw2(); 

  fhOneOverPtVsPhiNeg = new TH2F("fhOneOverPtVsPhiNeg","1/pt versus track phi negative tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);//AID//
  fhOneOverPtVsPhiNeg->Sumw2();

  fhOneOverPtVsPhiPos = new TH2F("fhOneOverPtVsPhiPos","1/pt versus track phi positive tracks", 36, 0, 2*TMath::Pi(), 40, 0, 0.4);//AID//
  fhOneOverPtVsPhiPos->Sumw2(); 

  fhSigmaPtOverPtVsPt = new TH2F("fhSigmaPtOverPtVsPt",
                                       "track sigma(1/pt)/ 1/pt vs pt", 100, 0, 100, 250, 0, 1);//AID//
  fhSigmaPtOverPtVsPt->Sumw2();

 
  Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
  Int_t nbins = sizeof(bins)/sizeof(Double_t)-1; //pT binning for DCA distribution

  fhDCAinXVsPt = new TH2F("fhDCAinXVsPt","fhDCAinXVsPt",nbins, bins, 200, -10.,10);//AID//
  fhDCAinXVsPt->Sumw2(); 
  fhDCAinYVsPt = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinYVsPt");//AID//
  fhDCAinYVsPt->Sumw2(); 

  if (fCorrMode) {
     fhDCAinXVsPtPhysPrimary = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinXVsPtPhysPrimary"); //AID//
     fhDCAinXVsPtPhysPrimary->Sumw2();
     fhDCAinYVsPtPhysPrimary = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinYVsPtPhysPrimary");//AID//
     fhDCAinXVsPtSecondary   = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinXVsPtSecondary");//AID//
     fhDCAinYVsPtSecondary   = (TH2F*) fhDCAinXVsPt->Clone("fhDCAinYVsPtSecondary");//AID//
     fhFractionOfSecInJet    = new TH2D("fhFractionOfSecInJet", "Frac of jet pT carried by secondary tracks",50,0,50,210,0,1.05);//AID//
     fhFractionOfSecInJet->Sumw2();

     fhPtTrkTruePrimRec = new TH1D("fhPtTrkTruePrimRec","Pt spectrum of reconstructed true generator level phys prim particles", 1000, 0, 100); 
     fhPtTrkTruePrimRec->Sumw2(); 
     fhPtTrkTruePrimGen = new TH1D("fhPtTrkTruePrimGen","Pt spectrum of generator level phys prim particles", 1000, 0, 100);  
     fhPtTrkTruePrimGen->Sumw2();
     fhPtTrkSecOrFakeRec = new TH1D("fhPtTrkSecOrFakeRec","Pt spectrum of reconstructed fake or secondary tracks", 1000, 0, 100);      
     fhPtTrkSecOrFakeRec->Sumw2();
  }

  fOutputList->Add(fZVertex);                                //AID//
  fOutputList->Add(fhTrackEta);                              //AID//  
  fOutputList->Add(fhTrackPhi);                              //AID//
  fOutputList->Add(fhJetEta);                                //AID//  
  fOutputList->Add(fhJetPhi);                                //AID// 
  fOutputList->Add(fhOneOverPtVsPhiNeg);                     //AID//
  fOutputList->Add(fhOneOverPtVsPhiPos);                     //AID//
  fOutputList->Add(fhSigmaPtOverPtVsPt);                     //AID//
  fOutputList->Add(fhDCAinXVsPt);                            //AID//
  fOutputList->Add(fhDCAinYVsPt);                            //AID//
  if (fCorrMode) {
     fOutputList->Add(fhDCAinXVsPtPhysPrimary);              //AID//
     fOutputList->Add(fhDCAinYVsPtPhysPrimary);              //AID//
     fOutputList->Add(fhDCAinXVsPtSecondary);                //AID//
     fOutputList->Add(fhDCAinYVsPtSecondary);                //AID//
     fOutputList->Add(fhFractionOfSecInJet);                 //AID//

     fOutputList->Add(fhPtTrkTruePrimRec); 
     fOutputList->Add(fhPtTrkTruePrimGen);
     fOutputList->Add(fhPtTrkSecOrFakeRec);
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++

    
  // Control histogram
  fhEntries = new TH1F("hEntries", "Analyzed sample properties", 11, -.5, 10.5);
  fhEntries->GetXaxis()->SetBinLabel(1, "nEventsAnal");
  fhEntries->GetXaxis()->SetBinLabel(2, "nEvPhySel");
  fhEntries->GetXaxis()->SetBinLabel(3, "nEvGoodJetArray");
  fhEntries->GetXaxis()->SetBinLabel(4, "nEvPile-up spd");
  fhEntries->GetXaxis()->SetBinLabel(5, "nEvPile-up mv");
  fhEntries->GetXaxis()->SetBinLabel(6, "nTracksEv");
  fhEntries->GetXaxis()->SetBinLabel(7, "nJetsCand");
  fhEntries->GetXaxis()->SetBinLabel(8, "nJetsTagged");
  fhEntries->GetXaxis()->SetBinLabel(9, "nUnexpError");
  fhEntries->GetXaxis()->SetBinLabel(10, "noMCHeader");
  fhEntries->GetXaxis()->SetBinLabel(11, "nEvPtHardOutlier");
  fOutputList->Add(fhEntries);

  fhEvtRej  = new TH1F("fhEvtRej", "Event rejection criteria", 11, -.5, 10.5);
  fOutputList->AddLast(fhEvtRej);
  
  fhEvtRejBitmap  = new TH1F("fhEvtRejBitmap", "Event rejection criteria bitmap", 2049, -.5, 2048.5);
  fOutputList->AddLast(fhEvtRejBitmap);	
	
	
	

  fhZNApercentQa = new TH1F("fhZNApercentQa", "ZNA multiplicity percentile;percent;dN/d(percent)", 100, 0., 100.);
  fOutputList->Add(fhZNApercentQa);

  fhHFjetQa = new TH1F("fhHFjetQa", "Check some QA for HF Jets", 12, 0.5, 12.5);
  fhHFjetQa->GetXaxis()->SetBinLabel(1,  "nWrgProngValue");      //Must be 0;
  fhHFjetQa->GetXaxis()->SetBinLabel(2,  "nJetiRejInFindVtx");   //Must be 0;
  fhHFjetQa->GetXaxis()->SetBinLabel(3,  "nJet_wTrks<2");
  fhHFjetQa->GetXaxis()->SetBinLabel(4,  "nJet_wGoodTrk<2"); fStartBin = 5;
  fhHFjetQa->GetXaxis()->SetBinLabel(5,  "nTrk_wWrongIndex");
  fhHFjetQa->GetXaxis()->SetBinLabel(6,  "nTrk_All");
  fhHFjetQa->GetXaxis()->SetBinLabel(7,  "nTrk_wFB_4");
  fhHFjetQa->GetXaxis()->SetBinLabel(8,  "nTrk_wFB_9");
  fhHFjetQa->GetXaxis()->SetBinLabel(9,  "nTrk_wNegID");
  fhHFjetQa->GetXaxis()->SetBinLabel(10, "nTrk_wFB_4&9");
  fhHFjetQa->GetXaxis()->SetBinLabel(11, "nTrk_wFB_!4!9");
  fhHFjetQa->GetXaxis()->SetBinLabel(12, "nTrkRejected");
  fOutputList->AddLast(fhHFjetQa);

  fhRhoQa = new TH1F("fhRhoQa", "Rho distribution;Probability Density;pt", 600, -10, 50);
  fOutputList->Add(fhRhoQa);
  if (fCorrMode) {
    fhMCRhoQa = new TH1F("fMCRhoQa", "Rho distribution;Probability Density;pt", 600, -10, 50);
    fOutputList->Add(fhMCRhoQa);
  }

  if (fDoRndmCone) {
    fhDeltaPt = new TH1F("fDeltaPt", "DeltaPt distribution", 500, -125, +125);
    fhDeltaPtLxy5 = new TH1F("fhDeltaPtLxy5", "DeltaPt distribution Lxy = 5", 500, -125, +125);  //newDeltaPt//
    fhDeltaPtLxy6 = new TH1F("fhDeltaPtLxy6", "DeltaPt distribution Lxy = 6", 500, -125, +125);  //newDeltaPt//
    fhDeltaPtLxy7 = new TH1F("fhDeltaPtLxy7", "DeltaPt distribution Lxy = 7", 500, -125, +125);  //newDeltaPt//
    fhDeltaPtTrack10 = new TH1F("fhDeltaPtTrack10", "DeltaPt distribution P_{T,track}>10 GeV/c", 500, -125, +125);  //newDeltaPt//
    
    fhDeltaPt ->Sumw2();     //newDeltaPt//
    fhDeltaPtLxy5 ->Sumw2();  //newDeltaPt//
    fhDeltaPtLxy6 ->Sumw2();  //newDeltaPt//
    fhDeltaPtLxy7 ->Sumw2();  //newDeltaPt//
    fhDeltaPtTrack10->Sumw2();  //newDeltaPt//
    
    fOutputList->Add(fhDeltaPtTrack10);  //newDeltaPt//
    fOutputList->Add(fhDeltaPt);
    fOutputList->Add(fhDeltaPtLxy5);  //newDeltaPt//
    fOutputList->Add(fhDeltaPtLxy6);  //newDeltaPt//
    fOutputList->Add(fhDeltaPtLxy7);  //newDeltaPt//  
  }

if(fDoEmbedding)
{   //EMB
     fhDeltaPtEmbedd  = new TH1F("fhDeltaPtEmbedd", "DeltaPt distribution based on track embedding", 500, -125, +125);
     fhDeltaPtEmbedd-> Sumw2();
     fOutputList->Add(fhDeltaPtEmbedd);

     fhDeltaPtEmbeddCorrelation = new TH2F("fhDeltaPtEmbeddCorrelation", "DeltaPt distribution based on track embedding vs pT of Embedded Track", 200, 0, 200, 500, -125, +125);
     fhDeltaPtEmbeddCorrelation->Sumw2();
     fOutputList->Add(fhDeltaPtEmbeddCorrelation);
  }
if(fEmbeddPerpendicular)
{   //EMB_clus
     fhDeltaPtEmbeddPerpendicular  = new TH1F("fhDeltaPtEmbeddPerpendicular", "DeltaPt distribution based on Perpendicular track embedding", 500, -125, +125);
     fhDeltaPtEmbeddPerpendicular-> Sumw2();
     fOutputList->Add(fhDeltaPtEmbeddPerpendicular);

     fhDeltaPtEmbeddCorrelationPerpendicular = new TH2F("fhDeltaPtEmbeddCorrelationPerpendicular", "DeltaPt distribution based on Perpendicular track embedding vs pT of Embedded Track", 200, 0, 200, 500, -125, +125);
     fhDeltaPtEmbeddCorrelationPerpendicular->Sumw2();
     fOutputList->Add(fhDeltaPtEmbeddCorrelationPerpendicular);
  }


 
 
  fhXsec = new TProfile("hXsec", "xsec from pyxsec.root", 1, 0.5, 1.5);
  fhXsec->GetXaxis()->SetBinLabel(1, Form("SelEvent_%s", fPtHardName.Data()));
  fhXsec->GetXaxis()->SetTitle("p_{T} hard bin");
  fhXsec->GetYaxis()->SetTitle("#<sigma>");
  fOutputList->Add(fhXsec);

  fhTrials = new TH1F("hTrials", "trials root file", 1, 0.5, 1.5);
  fhTrials->GetXaxis()->SetBinLabel(1, Form("SelEvent_%s", fPtHardName.Data()));
  fhTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fhTrials->GetYaxis()->SetTitle("#sum{ntrials}");
  fOutputList->Add(fhTrials);


  PostData(1, fOutputList);
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::ExecOnce()
{
  AliInfo(MSGINFO("+++ Executing ExecOnce +++"));
  if (fGlLogLevel) {
    AliLog::SetGlobalLogLevel((AliLog::EType_t)fGlLogLevel);
  }

  if (fEmbeddPerpendicular){
	   fFastJetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper"); //EMB_clus
	   fFastJetWrapper->SetAreaType(fastjet::active_area);		//EMB_clus
	   fFastJetWrapper->SetGhostArea(0.005);
	   fFastJetWrapper->SetR(0.4);					//EMB_clus
	   fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);	//EMB_clus
	   fFastJetWrapper->SetRecombScheme(fastjet::pt_scheme);	//EMB_clus
	   fTrackGenerator = new TRandom(0);
  }

  if ((fGlLogLevel > AliLog::kInfo) && (fLcDebLevel)) {
    AliLog::SetGlobalDebugLevel(0);
    AliLog::SetClassDebugLevel("AliHFJetsTaggingVertex",        fLcDebLevel);
    AliLog::SetClassDebugLevel("AliAnalysisTaskEmcalJetBtagSV", fLcDebLevel);
  }

  fInitialized = kTRUE;
  return;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserExec(Option_t* /*option*/)
{
  if (!fInitialized)
    ExecOnce();

  AliDebug(1, MSGINFO("+++ Executing UserExec +++"));

  // AOD input event
  fEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fEvent) {
    AliWarning(MSGWARNING("Input AOD not available, trying with output handler..."));
    if (AODEvent() && IsStandardAOD()) {
      // In case there is an AOD filter writing a standard AOD, use the AOD
      // event in memory rather than the input event.
      fEvent = dynamic_cast<AliAODEvent*>(AODEvent());
    } else {
      AliError(MSGERROR("No AOD handler found or no standard AOD!"));
      return;
    }
  }

  // load MC header
  fMCHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if (!fMCHeader && fCorrMode) {
    AliError(MSGERROR("MC header branch not found!"));
    return;
  }
  


  if (!fGenNamePattern.IsNull() || !fGenNamePattern.IsWhitespace()) {
    TString title = (fMCHeader) ? fMCHeader->GetGeneratorName() : "";
    if (!title.IsNull() && !title.Contains(fGenNamePattern)) {
      AliDebugF(3, MSGWARNING("Pattern not found in MC header title %s"), title.Data());
      return;
    }
  }



  // get multiplicity and centrality percentile
  if (fEvent) 
    fMultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    
  if(fMultSelection) { 
    fZNApercentile = fMultSelection->GetMultiplicityPercentile("ZNA");
    fhZNApercentQa->Fill(fZNApercentile);
  }
  else {
    fZNApercentile = -1;
    AliWarning("AliMultSelection object not found!");
  }

  // ALL EVENTS
  fhEntries->Fill(0); // EventsAnal
  if (fCheckMCCrossSection && !fSkipWeightInfo) {
    fhXsec->Fill(  1., fMCXsec);
    fhTrials->Fill(1., fMCAvgTrials);
    fSkipWeightInfo = kTRUE;
  }

  if (!fCutsHFjets->IsEventSelected((AliAODEvent*)fEvent)) {
    AliDebug(5, MSGDEBUG("Event did not pass event selection from AliRDHFJetsCuts!"));
    fhEvtRej->Fill(fCutsHFjets->GetWhyRejection(), 1);
	fhEvtRejBitmap->Fill(fCutsHFjets->GetEventRejectionBitMap(),1);
	  
    return;
  }

  fhEntries->Fill(1); // event selected, pileup, trigger, etc...

  if (!GetArrays()) {
    PostData(1, fOutputList);
    return;
  }
  fhEntries->Fill(2); // events with Jet arrays

  fHFvertexing = (!fHFvertexing) ? new TClonesArray("AliAODVertex", 0) : fHFvertexing;
  fHFvertexing->SetOwner(kTRUE);

  if (fDoFillV0Trks && !FillMapOfV0gTrkIDs()) {
    AliError(MSGERROR("Error filling V0 tracks info"));
    PostData(1, fOutputList);
    return;
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++

  if(fEvent){   //AID// QA of tracks and jets 
     Double_t xyz[50];
     Double_t pxpypz[50];
     Double_t cv[21];
     Int_t iTracks =  fEvent->GetNumberOfTracks();

     for(Int_t i = 0; i < iTracks; i++){ 
            
        AliAODTrack *track = static_cast <AliAODTrack*>( fEvent->GetTrack(i));
        if(!track) continue;
    UInt_t trkFilterMap = track->GetFilterMap();  
        if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;
        if (TMath::Abs(track->Eta()) > fEtaCut) continue; 
        if (track->Pt() < fPtCut) continue;

        fhTrackEta->Fill(track->Pt(), track->Eta());
        fhTrackPhi->Fill(track->Pt(), track->Phi());

        //get sigma pT / pT  
        //Taken from AliEMCalTriggerExtraCuts::CalculateTPCTrackLength
        memset(cv, 0, sizeof(Double_t) * 21); //cleanup arrays
        memset(pxpypz, 0, sizeof(Double_t) * 50);
        memset(xyz, 0, sizeof(Double_t) * 50);
        track->GetXYZ(xyz);
        track->GetPxPyPz(pxpypz);
        track->GetCovarianceXYZPxPyPz(cv);
    
        AliExternalTrackParam  par(xyz, pxpypz, cv, track->Charge());
        fhSigmaPtOverPtVsPt->Fill(track->Pt(), TMath::Abs(sqrt(par.GetSigma1Pt2())/par.GetSigned1Pt()));

        if(track->Charge()<0){
           fhOneOverPtVsPhiNeg->Fill(track->Phi(), 1.0/track->Pt());
        }else{
           fhOneOverPtVsPhiPos->Fill(track->Phi(), 1.0/track->Pt());
        }

        //DCA distributions
        fhDCAinXVsPt->Fill(track->Pt(), track->XAtDCA());
        fhDCAinYVsPt->Fill(track->Pt(), track->YAtDCA());

     }
     AliAODVertex* pVtx = (AliAODVertex*) fEvent->GetPrimaryVertex();
     Double_t pvXYZ[3];
     pVtx->GetXYZ(pvXYZ);
     fZVertex->Fill(pvXYZ[2]);

     Double_t ptJetGen_wBkgRej;
     Int_t nJets = fRecJetArray->GetEntries();
     Double_t rho   = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
     AliEmcalJet* jet; //AID// Fraction of non-primary jet constituents 
     for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {
        jet = (AliEmcalJet*) fRecJetArray->UncheckedAt(jetcand);
        if (fCutsHFjets->IsJetSelected(jet)) {
           ptJetGen_wBkgRej = jet->Pt() - (jet->Area() * rho);
           fhJetEta->Fill(jet->Pt(), jet->Eta());
           fhJetPhi->Fill(jet->Pt(), jet->Phi());
        }
     }
  }//AID end
  //++++++++++++++++++++++++++++++++++++++++++++++++
    
  // Execute analysis for current event
  if (fCorrMode)
    if( IsOutlier()){             //FK// Check whether this event is pthard bin outlier 
       fhEntries->Fill(10);       //FK//
       PostData(1, fOutputList); //FK//
       return;                   //FK//
    }else{                       //FK//

        
       //AID//++++++++++++++++++++++++++++++++++++
       Int_t iTracks =  fEvent->GetNumberOfTracks();    //AID
       Int_t label, labelMC;                                     //AID
       Bool_t labelfound=0;
       AliAODMCParticle* particleMC = NULL;             //AID
       AliAODMCParticle* particleMCMother = NULL;       //AID

       for (Int_t it = 0; it < fMCPartArray->GetEntries(); it++) {
           particleMC   = (AliAODMCParticle*) fMCPartArray->At(it);
           if(particleMC->IsPhysicalPrimary()){
               if(particleMC->Pt() < fPtCut) continue;
               if(TMath::Abs(particleMC->Eta()) > fEtaCut) continue; 
               if(! particleMC->Charge()) continue; 
               fhPtTrkTruePrimGen->Fill(particleMC->Pt());
           }
       }
        
       for (Int_t it = 0; it < fRecTrkArray->GetEntries(); it++) {
          AliAODTrack *track = static_cast <AliAODTrack*>( fRecTrkArray->ConstructedAt(it));
          if(!track) continue;
          UInt_t trkFilterMap = track->GetFilterMap();  
          if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;
          if (TMath::Abs(track->Eta()) > fEtaCut) continue; 
          if (track->Pt() < fPtCut) continue;

          label = TMath::Abs(track->GetLabel());        //AID

          particleMC = NULL; 
          labelfound=0;
          for(Int_t it = 0; it < fMCPartArray->GetEntries(); it++) { //find gen level particle with the same label
             particleMC   = (AliAODMCParticle*) fMCPartArray->At(it);
             labelMC = TMath::Abs(particleMC->GetLabel());
             if(labelMC==label && label > -1){
                labelfound=1;
                break;
             }
          }
          if(labelfound && particleMC && particleMC->IsPhysicalPrimary()){
             fhDCAinXVsPtPhysPrimary->Fill(track->Pt(), track->XAtDCA());
             fhDCAinYVsPtPhysPrimary->Fill(track->Pt(), track->YAtDCA());
             fhPtTrkTruePrimRec->Fill(particleMC->Pt());
          }else{
             fhDCAinXVsPtSecondary->Fill(track->Pt(), track->XAtDCA());
             fhDCAinYVsPtSecondary->Fill(track->Pt(), track->YAtDCA());
             fhPtTrkSecOrFakeRec->Fill(track->Pt());
          }//AID
       }//AID
       
       AliAODTrack* constTrackRec = NULL; //AID// jet constituent 
       AliEmcalJet* jet; //AID// Fraction of non-primary jet constituents 
       Double_t sumall = 0.; 
       Double_t sumsec = 0.; 
       Int_t nJets = fRecJetArray->GetEntries();
       Double_t rho   = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
       for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {
          jet = (AliEmcalJet*) fRecJetArray->UncheckedAt(jetcand);
          if(fCutsHFjets->IsJetSelected(jet)) {
             sumall = 0.; 
             sumsec = 0.; 
          
             for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) { //loop over jet constituents
                constTrackRec = ((AliAODTrack*) jet->TrackAt(iq, fRecTrkArray));
                if(!constTrackRec) continue;
                UInt_t trkFilterMap = constTrackRec->GetFilterMap();  
                if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)){  
                   AliError(MSGERROR("Non hybrid tracks in jet")); continue;
                }
 
                label = TMath::Abs(constTrackRec->GetLabel());        //AID

                particleMC = NULL; 
                labelfound=0;
                for(Int_t it = 0; it < fMCPartArray->GetEntries(); it++) { //find gen level particle with the same label
                   particleMC   = (AliAODMCParticle*) fMCPartArray->At(it);
                   labelMC = TMath::Abs(particleMC->GetLabel());
                   if(labelMC==label && label > -1){
                      labelfound=1;
                      break;
                   }
                }
 
                if(!(labelfound && particleMC && particleMC->IsPhysicalPrimary())){
                   sumsec += constTrackRec->Pt();
                }
                sumall += constTrackRec->Pt();
             }
             if(sumall>0){
                Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
                fhFractionOfSecInJet->Fill( ptJet_wBkgRej, sumsec/sumall);
             } 
          }
       }//AID  
       //++++++++++++++++++++++++++++++++++++


       AnalyseCorrectionsMode(); // must be MC, all steps are filled for container kBJets (only)
    }
  else
    AnalyseDataMode();        // can also be MC, only step kCFStepReco is filled also for kBJets

  PostData(1, fOutputList);
  return;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseDataMode()
{

  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO
  AliAODVertex* pVtx = (AliAODVertex*)fEvent->GetPrimaryVertex();

  Double_t pvXYZ[3], pvCov[6];

  pVtx->GetXYZ(pvXYZ);
  pVtx->GetCovarianceMatrix(pvCov);

  AliESDVertex* esdVtx = new AliESDVertex(pvXYZ, pvCov, pVtx->GetChi2(), pVtx->GetNContributors());

  Double_t magzkG = (Double_t)fEvent->GetMagneticField();

  Int_t nJets = fRecJetArray->GetEntries();

  Double_t rho = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
  fhRhoQa->Fill(rho);
  
  Double_t deltapt=99999;   //newDeltaPt
  if(fDoRndmCone){
     deltapt = GetDeltaPtRandomCone(fTaggingRadius, rho);  //newDeltaPt
     if(deltapt<9999){ 
        fhDeltaPt->Fill(deltapt);
        //-------------------fhDeltaPtTrack10-----------------
        Double_t    signalPhi, signalEta;
        Bool_t fillDeltaPt = kFALSE;
        for(Int_t i = 0; i < fRecTrkArray->GetEntries(); i++) {
           AliAODTrack* trk = static_cast<AliAODTrack*>(fRecTrkArray->ConstructedAt(i));
           UInt_t trkFilterMap = trk->GetFilterMap();  
           if(!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;
           if( (fabs(trk->Eta()) < fEtaCut) && (trk->Pt() > 10) ) {
              fillDeltaPt = kTRUE;
              signalPhi =  trk->Phi();
              signalEta =  trk->Eta();
              break;
           }
        }
           
        if(fillDeltaPt)
           fhDeltaPtTrack10->Fill(GetDeltaPtRandomConeWithoutSignalPt(fTaggingRadius,rho, signalEta, signalPhi), fMCWeight);
            //--------------------------------------------------------     
     } 
  }

  vctr_pair_dbl_int aVtxDisp;
  aVtxDisp.reserve(5);   // reserve space for 5 vertex sigma position

  AliEmcalJet* jet;
  Int_t fillDelPtMask = 0; //newDelPt// 
  Int_t fSVreconstucted=0;
  for(Int_t jetcand = 0; jetcand < nJets; ++jetcand) {
     jet = (AliEmcalJet*) fRecJetArray->UncheckedAt(jetcand);
     
     if(!fCutsHFjets->IsJetSelected(jet)) {
        AliDebugF(5, MSGDEBUG("--> Jet with pt=%3.2f and eta=%3.2f not selected in FindVertices"),
                 jet->Pt(), jet->Eta());
        continue;
     }
     Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
     // Run b-tagger
     Int_t  nDauRejCount = 0;
     Int_t nVtx = fTagger->FindVertices(jet,
                                       fRecTrkArray,
                                       (AliAODEvent*)fEvent,
                                       esdVtx,
                                       magzkG,
                                       fHFvertexing,
                                       fV0gTrkMap,
                                       aVtxDisp,
                                       nDauRejCount);
     fhHFjetQa->Fill(12, nDauRejCount);
     if(nVtx < 0){
       fhHFjetQa->Fill(-1 * nVtx);
       continue;
     }
     //------------------------newDeltaPt-------------------------
     if(fDoRndmCone && nVtx > 0 && fillDelPtMask < 7){
     
        fillDelPtMask = FillDeltaPt( rho, nVtx, pVtx, aVtxDisp, jet->Eta(), jet->Phi(), fillDelPtMask);     
        if(((fillDelPtMask & 4) >> 2) == 1) fSVreconstucted=1; 
     }  
     //-------------------------------------------------
    
     fhJetVtxData->FillStepJetVtxData(AliHFJetsContainer::kCFStepReco,
                                     nVtx,
                                     fZNApercentile,
                                     ptJet_wBkgRej,
                                     aVtxDisp,
                                     fHFvertexing,
                                     pVtx,
                                     jet,
                                     fMCWeight);

     fHFvertexing->Clear();
     aVtxDisp.clear();
  }

 


  if(fSVreconstucted == 1 && fDoEmbedding){
     AliVParticle*  hytrk  = NULL;  //track hybrid event jet
     AliEmcalJet*   hyjet  = NULL; //hybrid event jet
     Double_t sumTrkEmbeddedPt=0;
     

     //EMB loop over jets in hybrid event
     for(Int_t i = 0; i < fHybridJetCont->GetEntries(); i++) {
        hyjet = static_cast<AliEmcalJet*>(fHybridJetCont->UncheckedAt(i));
        if(!hyjet)  continue;
        if(!fCutsHFjets->IsJetSelected(hyjet)) continue;

        sumTrkEmbeddedPt=0;
        for(Int_t iq=0; iq < hyjet->GetNumberOfTracks(); iq++) {
           hytrk = static_cast<AliVParticle*> (hyjet->Track(iq));
           if(!hytrk) continue;
           //cout<<"HYTRACK "<<hytrk->Pt()<<"   "<< hytrk->Eta()<<"   "<< hytrk->Phi()<<"  " <<hytrk->Charge()<<endl;
           if(hytrk->Charge()==3) sumTrkEmbeddedPt += hytrk->Pt();
        }

        if(sumTrkEmbeddedPt>0){ 
           Double_t deltaPtEmb = hyjet->Pt() - hyjet->Area() * rho - sumTrkEmbeddedPt;
           fhDeltaPtEmbedd->Fill(deltaPtEmb);
           fhDeltaPtEmbeddCorrelation->Fill(sumTrkEmbeddedPt, deltaPtEmb);
        }
     }
  }//EMB

  delete esdVtx;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseCorrectionsMode()
{
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO
   AliAODVertex* pVtx = (AliAODVertex*)fEvent->GetPrimaryVertex();

  Double_t pvXYZ[3], pvCov[6];

  pVtx->GetXYZ(pvXYZ);
  pVtx->GetCovarianceMatrix(pvCov);

  AliESDVertex* esdVtx = new AliESDVertex(pvXYZ, pvCov, pVtx->GetChi2(), pVtx->GetNContributors());

  Double_t magzkG = (Double_t)fEvent->GetMagneticField();

  // MC primary vertex
  Double_t vtxTrue[3];
  fMCHeader->GetVertex(vtxTrue);

  Int_t nMCJets = fMCJetArray->GetEntries();
  Int_t nJets   = fRecJetArray->GetEntries();

  Double_t rhoMC = (fDoBkgRej) ? GetExternalRho(kTRUE)  : 0.;
  Double_t rho   = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
  fhRhoQa->Fill(rho, fMCWeight);
  fhMCRhoQa->Fill(rhoMC, fMCWeight);

  Double_t deltapt=99999;
  if(fDoRndmCone) {
     deltapt = GetDeltaPtRandomCone(fTaggingRadius, rho); //newDeltaPt
     if(deltapt<9999){ 
        fhDeltaPt->Fill(deltapt, fMCWeight);
        //-------------------fhDeltaPtTrack10-----------------
        Double_t signalPhi, signalEta;    
        Bool_t fillDeltaPt = kFALSE;
        for(Int_t i = 0; i < fRecTrkArray->GetEntries(); i++) {
            AliAODTrack* trk = static_cast<AliAODTrack*>(fRecTrkArray->ConstructedAt(i));
            UInt_t trkFilterMap = trk->GetFilterMap();  
            if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;
            if ( (fabs(trk->Eta()) < fEtaCut) && (trk->Pt() > 10) ) {
               fillDeltaPt = kTRUE;
                signalPhi =  trk->Phi();
                signalEta =  trk->Eta();
               break;
            }
        }       
        if (fillDeltaPt)
         fhDeltaPtTrack10->Fill(GetDeltaPtRandomConeWithoutSignalPt(fTaggingRadius,rho, signalEta, signalPhi), fMCWeight);
       //--------------------------------------------------------     
     } 
  }

  vctr_pair_dbl_int aVtxDisp;
  aVtxDisp.reserve(5);

  // Loop on MC jets
  AliEmcalJet* jetMC;
  Int_t fillDelPtMask=0; //newDeltaPt//
    
  for (Int_t jetcand = 0; jetcand < nMCJets; ++jetcand) {

    jetMC = (AliEmcalJet*)fMCJetArray->UncheckedAt(jetcand);
    if (!jetMC) continue;

    // restrict jet eta and pT ranges
    if (!fCutsHFjets->IsJetSelected(jetMC)) {
      AliDebugF(5, MSGDEBUG("JetMC not selected: pT=%f, eta=%f!"), jetMC->Pt(), jetMC->Eta());
      continue;
    }

    // Get jet flavour from 2 methods
    Double_t partonnatMC[2] = { -1, -1};
    Double_t ptpartMC[2]    = { -1, -1};

    GetFlavour2Methods(jetMC, partonnatMC, ptpartMC, fTaggingRadius);

    Double_t ptJetGen_wBkgRej = jetMC->Pt() - (jetMC->Area() * rhoMC);

    if (!fDoOnlyMtxAna) {
      // Fill container tagger
      // At this point we do not need to fill the secondary vertex QA container
       fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepEventSelected, 0, 0, ptJetGen_wBkgRej, aVtxDisp,
                                     NULL, NULL, jetMC, NULL, partonnatMC, ptpartMC, fMCWeight);

       aVtxDisp.clear(); //just on case
    } 
    
    if (fDoDetRespMtx) {
      Double_t vector[4] = {ptJetGen_wBkgRej, jetMC->Eta(), partonnatMC[0], partonnatMC[1]};
      fhnGenerated->Fill(vector, fMCWeight);
    }
    
  } // end loop on jets
  // Loop on jets (clusterized on RECO particles)
  AliEmcalJet* jet;
  for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {

    jet = (AliEmcalJet*)fRecJetArray->UncheckedAt(jetcand);
    if (!fCutsHFjets->IsJetSelected(jet)) {
      AliDebugF(5, MSGDEBUG("Jet not selected: pT=%f, eta=%f!"), jet->Pt(), jet->Eta());
      continue;
    }

    // Get jet flavour from 3 methods
    Double_t partonnat[2] = { -1, -1};
    Double_t ptpart[2]    = { -1, -1};

    GetFlavour2Methods(jet, partonnat, ptpart, fTaggingRadius);

    Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
    Int_t nVtx = 0;
    if (!fDoOnlyMtxAna) {
      // // Run vertex tagging
      Int_t nDauRejCount = 0;
      nVtx = fTagger->FindVertices(jet,
                                   fRecTrkArray,
                                   (AliAODEvent*)fEvent,
                                   esdVtx,
                                   magzkG,
                                   fHFvertexing,
                                   fV0gTrkMap,
                                   aVtxDisp,
                                   nDauRejCount);
      fhHFjetQa->Fill(12, nDauRejCount);
      if (nVtx < 0) {
        fhHFjetQa->Fill(-1 * nVtx);
        continue;
      }

      //------------------------newDeltaPt-------------------------

      if(fDoRndmCone && nVtx > 0 && fillDelPtMask < 7) {
         fillDelPtMask = FillDeltaPt( rho, nVtx, pVtx,aVtxDisp, jet->Eta(), jet->Phi(), fillDelPtMask);                 
      }  
      //-------------------------------------------------
      // Fill jet-with-vertex container
      fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepReco,
                                     nVtx,
                                     fZNApercentile,
                                     ptJet_wBkgRej,
                                     aVtxDisp,
                                     fHFvertexing,
                                     pVtx,
                                     jet,
                                     fMCPartArray,
                                     partonnat,
                                     ptpart,
                                     fMCWeight);

       if (fDoQAVtx) {
         fhQaVtx->FillStepQaVtx(AliHFJetsContainer::kCFStepReco, nVtx, 0, pVtx, jet,
                                fHFvertexing, fMCPartArray, aVtxDisp, partonnat, fMCWeight);
       }
    }



    AliEmcalJet* matchedjet = jet->ClosestJet();
    if (!matchedjet) continue;

    GetFlavour2Methods(matchedjet, partonnat, ptpart, fTaggingRadius);

    Double_t ptJetMC_wBkgRej = matchedjet->Pt() - (matchedjet->Area() * rhoMC);

    if (!fDoOnlyMtxAna) {
      // step kCFStepMatchedAny
      fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepMatchedAny,
                                     nVtx,
                                     fZNApercentile,
                                     ptJetMC_wBkgRej,
                                     aVtxDisp,
                                     fHFvertexing,
                                     pVtx,
                                     matchedjet,
                                     fMCPartArray,
                                     partonnat,
                                     ptpart,
                                     fMCWeight);

      if (fDoQAVtx) {
        fhQaVtx->FillStepQaVtx(AliHFJetsContainer::kCFStepMatchedAny, nVtx, 0, pVtx, jet,
                               fHFvertexing, fMCPartArray, aVtxDisp, partonnat, fMCWeight);
      }
    }
    fHFvertexing->Clear();
    aVtxDisp.clear();

    // Fill detector response multidimensional histo
    if (fDoDetRespMtx) {
      Double_t vector[6] = {ptJet_wBkgRej, ptJetMC_wBkgRej, jet->Eta(), matchedjet->Eta(), partonnat[0], partonnat[1]};
      fhnDetRespMtx->Fill(vector, fMCWeight);
    }

  }

  delete esdVtx;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::Terminate(const Option_t*)
{
  //TERMINATE METHOD: NOTHING TO DO
  AliInfo(MSGINFO("+++ Executing Terminate +++"));
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::GetPythiaCrossSection()
{
  // Fetch the aod also from the input in,
  // have todo it in notify

  Float_t xsection  = 0;
  Float_t trials    = 1;

  TTree* tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliWarning(MSGWARNING("Analysis Manager could not found the tree...."));
    return;
  }

  TFile* curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliWarning(MSGWARNING("Current file not found"));
    return;
  }

  // Check if file not accessed previously, if so
  // return the previously calculated weight
  if (fCurrFileName == curfile->GetName()) {
    AliDebug(3, MSGDEBUG("File already read. Previous values used"));
    return;
  }
  fCurrFileName = TString(curfile->GetName());

  Bool_t ok = GetPythiaInfoFromFile(fCurrFileName, xsection, trials);
  if (!ok) {
    AliWarning(MSGWARNING("Parameters from file not recovered properly"));
    return;
  }

  fMCXsec = xsection;
  fMCAvgTrials = trials;

  // average number of trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();

  //  if (trials >= nEntries && nEntries > 0.) fMCavgTrials /= nEntries;
  AliDebugF(3, MSGINFO("xs %e, trial %e, avg trials %2.2f, events per file %e"),
            xsection, trials, fMCAvgTrials, nEntries);

  AliDebugF(3, MSGDEBUG("Reading File %s"), curfile->GetName());
  if (fMCAvgTrials > 0.) {
    fMCWeight =  (fUseWeight) ? fMCXsec / fMCAvgTrials : 1.;
  } else {
    AliWarningF(MSGWARNING("Average number of trials is NULL!! Set weight to 1: xs : %e, trials %e, entries %e"),
                xsection, trials, nEntries);
    fMCWeight = 1.;
  }

  AliDebugF(3, MSGINFO("MC Weight: %f"), fMCWeight);
  return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetPythiaInfoFromFile(TString file,
                                                            Float_t& xsec,
                                                            Float_t& trials)
{
  xsec   = 0;
  trials = 1;

  if (file.Contains("root_archive.zip#")) {
    Ssiz_t pos1 = file.Index("root_archive", 12, 0, TString::kExact);
    Ssiz_t pos  = file.Index("#", 1, pos1, TString::kExact);
    Ssiz_t pos2 = file.Index(".root", 5, TString::kExact);
    file.Replace(pos + 1, pos2 - pos1, "");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()), "");
  }

  //Printf("%s",file.Data());

  TFile* fxsec = TFile::Open(Form("%s%s", file.Data(), "pyxsec.root"));
  if (!fxsec || fxsec->IsZombie()) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s", file.Data(), "pyxsec_hists.root"));
    if (!fxsec || fxsec->IsZombie()) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }

      TList* list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }

      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree* xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }

    UInt_t   ntrials  = 0;
    Double_t xsection = 0;
    xtree->SetBranchAddress("xsection", &xsection);
    xtree->SetBranchAddress("ntrials",  &ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
    fxsec->Close();
  }

  return kTRUE;
}


//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::GetFlavour2Methods(AliEmcalJet* jet,
                                                       Double_t (&partonnat)[2],
                                                       Double_t (&ptpart)[2],
                                                       Double_t radius)
{
  // 2 methods to associate jet to mother parton

  /* Nature of the parton (methods 1)      *
   * 1 = gluon (pdg 21)                    *
   * 2 = light quark (pdg < 4)             *
   * 3 = c (pdg 4)                         *
   * 4 = b (pdg 5)                         *
   * Nature of the meson/barion (method 2) *
   * 2 = light                             *
   * 3 = with c                            *
   * 4 = with b                            */

  // Initialize output values
  for (Int_t i = 0; i < 2; ++i) {
    partonnat[i] =   0;
    ptpart[i]    = -1.;
  }

  AliAODMCParticle* parton[2];

  parton[0] = (AliAODMCParticle*) fTagger->IsMCJetParton(fMCPartArray, jet, radius);  // method 2
  parton[1] = (AliAODMCParticle*) fTagger->IsMCJetMeson(fMCPartArray, jet, radius);   // method 3

  if (parton[0]) {
    Int_t pdg = TMath::Abs(parton[0]->PdgCode());
    //if(pdg==4 || pdg==5)
    AliDebugF(6, MSGDEBUG("parton method -> pdg parton: %d"), pdg);

    if      (pdg == 21) partonnat[0] = 1;
    else if (pdg  < 4 ) partonnat[0] = 2;
    else if (pdg == 4 ) partonnat[0] = 3;
    else if (pdg == 5 ) partonnat[0] = 4;

    ptpart[0] = parton[0]->Pt();
  } else {

    AliWarning(MSGWARNING("No parton method output"));
  }

  if (!parton[1])
    partonnat[1] = 2;
  else {
    Int_t pdg = TMath::Abs(parton[1]->PdgCode());

    AliDebugF(6, MSGDEBUG("meson method -> pdg parton: %d"), pdg);

    if      ((pdg >= 400 && pdg <= 500) || (pdg >= 4000 && pdg <= 5000)) partonnat[1] = 3;
    else if ((pdg >= 500 && pdg <= 600) || (pdg >= 5000 && pdg <= 6000)) partonnat[1] = 4;

    ptpart[1] = parton[1]->Pt();
  }
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::CheckTrackQAinJets()
{
  if (fRecJetArray && fRecTrkArray) {
    AliEmcalJet* jet;
    for (Int_t jetcand = 0; jetcand < fRecJetArray->GetEntries(); ++jetcand) {
      jet = (AliEmcalJet*) fRecJetArray->UncheckedAt(jetcand);

      Int_t nTrksInJet = jet->GetNumberOfTracks();
      for (Int_t j = 0; j < nTrksInJet; ++j) {
        //utilize dynamic cast and then check pointer
        AliAODTrack* jTrk   = ((AliAODTrack*)jet->TrackAt(j, fRecTrkArray));
        if (!jTrk) {
          AliWarningF(MSGWARNING("Track in Jet with index %d/%d not found. Total number of AODtracks %d"),
                      j, nTrksInJet, fRecTrkArray->GetEntries());
          fhHFjetQa->Fill(fStartBin);  //fStartBin
          continue;
        }
        fhHFjetQa->Fill(fStartBin + 1);  //fStartBin + 1
        Int_t jTrkID = jTrk->GetID();
        UInt_t trkFilterMap = jTrk->GetFilterMap();
        bool hasFilterBit4 = TESTBIT(trkFilterMap, 4);
        bool hasFilterBit9 = TESTBIT(trkFilterMap, 9);

        if (hasFilterBit4)                    fhHFjetQa->Fill(fStartBin + 2); //bin fStartBin + 2
        if (hasFilterBit9)                    fhHFjetQa->Fill(fStartBin + 3); //bin fStartBin + 3
        if (jTrkID < 0)                       fhHFjetQa->Fill(fStartBin + 4); //bin fStartBin + 4
        if (hasFilterBit4 && hasFilterBit9)   fhHFjetQa->Fill(fStartBin + 5); //bin fStartBin + 5
        if (!hasFilterBit4 && !hasFilterBit9) fhHFjetQa->Fill(fStartBin + 6); //bin fStartBin + 6
      } //loop over tracks
    }   //loop over jets
  } else {
    AliError(MSGERROR("Either Jet or Tracks array found."));
  }//if arrays

  return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetArrays()
{
  // Get jet and track collections
  if (!fJetContName.IsNull()) {
    AliDebugF(4, MSGDEBUG("Retrieve jets %s!"), fJetContName.Data());

    fRecJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fJetContName));
    if (!fRecJetArray) {
      AliErrorF(MSGERROR("%s: Could not retrieve jets %s!"), GetName(), fJetContName.Data());
      return kFALSE;
    }
  }

  if (!fTrkContName.IsNull()) {
    AliDebugF(4, MSGDEBUG("Retrieve tracks %s!"), fTrkContName.Data());

    fRecTrkArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTrkContName));
    if (!fRecTrkArray) {
      AliErrorF(MSGERROR("%s: Could not retrieve tracks %s!"), GetName(), fTrkContName.Data());
      return kFALSE;
    }
  }

  //Get MC jet and particles collections
  if (fCorrMode) {

    if (!fMCJetContName.IsNull()) {
      AliDebugF(4, MSGDEBUG("Retrieve MC jets %s!"), fMCJetContName.Data());

      fMCJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMCJetContName));
      if (!fMCJetArray) {
        AliErrorF(MSGERROR("%s: Could not retrieve MC jets %s!"), GetName(), fMCJetContName.Data());
        return kFALSE;
      }
    }

    TString fMCPartContName(AliAODMCParticle::StdBranchName());
    AliDebugF(4, MSGDEBUG("Retrieve MC particles %s!"), fMCPartContName.Data());

    fMCPartArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMCPartContName));
    if (!fMCPartArray) {
      AliError( MSGERROR("MC particles branch not found!"));
      return kFALSE;
    }
  }

  //EMB
  if(fDoEmbedding){

      fHybridJetCont   = static_cast<TClonesArray*> (fEvent->FindListObject(fHybridJetContName.Data())); 
      if(!fHybridJetCont){
          AliError( MSGERROR("HYBRID JET CONTAINER not found!"));
          return kFALSE;
      } 
  }//EMB


  CheckTrackQAinJets();

  return kTRUE;
}

//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskEmcalJetBtagSV::FillMapOfV0gTrkIDs()
{

  fV0gTrkMap = (!fV0gTrkMap) ? new map_int_bool() : fV0gTrkMap;  //  Map of V0 trks (std::map< Int_t, Bool_t >              map_int_bool;)
  fV0gTrkMap->clear();

  //Fill V0 tracks AODTrack map array

  std::vector<Int_t> TrkIDs;
  FillVecOfV0gTrkIDs(TrkIDs);

  for (Int_t iTrk = 0; iTrk < fEvent->GetNumberOfTracks(); ++iTrk) {

    AliAODTrack* track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrk));

    Int_t trkID = track->GetID();

    //IMPORTANT:
    //Use same good track selection as it's implemented
    //in AliHFJetsTaggingVertex::FindVertex
    if (trkID < 0) continue;
    if (!track->GetFilterMap() && !track->GetTPCNcls()) continue;

    Bool_t trkBelongToV0 = IsAODtrkBelongToV0(TrkIDs, trkID);
    fV0gTrkMap->insert(std::pair<Int_t, Bool_t> (trkID, trkBelongToV0));
  }

  return kTRUE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::FillVecOfV0gTrkIDs(std::vector<Int_t>& vctrTrkIDs)
{

  //Fill V0 tracks AODTrack map array
  if (!fEvent)
    return kFALSE;

  AliAODEvent* aodEvent = (AliAODEvent*)fEvent;

  Int_t nV0s = aodEvent->GetNumberOfV0s();

  vctrTrkIDs.clear();
  vctrTrkIDs.reserve(nV0s);

  for (Int_t iV0 = 0; iV0 < nV0s; ++iV0) {

    AliAODv0* aodV0 = (AliAODv0*)aodEvent->GetV0(iV0);
    if (!aodV0) continue;

    AliAODTrack* pTrack = (AliAODTrack*)aodV0->GetDaughter(0);
    AliAODTrack* nTrack = (AliAODTrack*)aodV0->GetDaughter(1);

    vctrTrkIDs.push_back(pTrack->GetID());
    vctrTrkIDs.push_back(nTrack->GetID());
  }

  return kTRUE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::IsAODtrkBelongToV0(std::vector<Int_t>& vctrTrkIDs, Int_t trkID)
{
  vector<Int_t>::iterator it = std::find(vctrTrkIDs.begin(), vctrTrkIDs.end(), trkID);
  if (it != vctrTrkIDs.end())
    return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetExternalRho(Bool_t isMC)
{
  // Get rho from event using CMS approach
  AliRhoParameter* rho = NULL;
  TString rhoname = (!isMC) ? fRhoTaskName : fMCRhoTaskName;
  if (!rhoname.IsNull()) {
    rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(rhoname.Data()));
    if (!rho) {
      AliWarningF(MSGWARNING("%s: Could not retrieve rho with name %s!"), GetName(), rhoname.Data());
      return 0.;
    }
  } else {
    AliWarningF(MSGWARNING("No %s Rho task name provided"), (!isMC ? "DATA" : "MC"));
    return 0.;
  }

  return rho->GetVal();
}

//_____________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetDeltaPtRandomCone(Double_t jetradius, Double_t rhovalue)
{
    Double_t minConeEta = jetradius - fEtaCut;
    Double_t maxConeEta = fEtaCut - jetradius;

    // throw random cone
    Double_t coneEta = minConeEta + fRandom->Rndm()*(maxConeEta - minConeEta);
    Double_t conePhi = fRandom->Rndm()*TMath::TwoPi();

    // collect track pt within cone
    Double_t conePt = 0.;
    for (Int_t i = 0; i < fEvent->GetNumberOfTracks(); i++) {
        AliAODTrack* trk = static_cast<AliAODTrack*>(fEvent->GetTrack(i));

        // track filter hardwired...
        UInt_t trkFilterMap = trk->GetFilterMap();  
        if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;

        if ( (fabs(trk->Eta()) < fEtaCut) && (trk->Pt() > fPtCut) ) {
            Double_t dphi = TVector2::Phi_mpi_pi((trk->Phi() - conePhi));
            Double_t deta = trk->Eta() - coneEta;
            Double_t dist = sqrt(deta*deta + dphi*dphi);
            if (dist < jetradius) conePt += trk->Pt();
        }
    } // track loop

  return conePt - jetradius*jetradius*TMath::Pi() * rhovalue;

}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::IsOutlier(){ //FK// whole function
   //Checks that this event is pthard bin outlier
   //inspired by Bool_t AliConvEventCuts::IsJetJetMCEventAccepted  

   if(TMath::Abs(fMaxFacPtHard) < 1e-6) return kFALSE; //FK// skip 

   TList *genHeaders         = 0x0;
   AliGenEventHeader* gh     = 0;
   Float_t ptHard;
   AliEmcalJet* jetMC = 0x0;
   Int_t nMCJets = fMCJetArray->GetEntries();
   Bool_t bPythiaHeader = 0; // flag whether pythia header was found

   if(MCEvent()){
      genHeaders = MCEvent()->GetCocktailList(); //get list of MC cocktail headers 
   }

   if(genHeaders){
      for(Int_t i = 0; i<genHeaders->GetEntries(); i++){
         gh = (AliGenEventHeader*)genHeaders->At(i);

         AliGenPythiaEventHeader* pyhead= dynamic_cast<AliGenPythiaEventHeader*>(gh); //identify pythia header

         if(pyhead){
            bPythiaHeader = 1;
            ptHard = pyhead->GetPtHard();

            for(Int_t jetcand = 0; jetcand < nMCJets; ++jetcand) {
               jetMC = (AliEmcalJet*)fMCJetArray->UncheckedAt(jetcand);
               if (!jetMC) continue;
               //Compare jet pT and pt Hard
               if(jetMC->Pt() > fMaxFacPtHard * ptHard){
                  return kTRUE;
               }
            }
         }
      }
      if(!bPythiaHeader){ //ptyhia header was not found
          AliWarning("AliAnalysisTaskEmcalJetBtagSV MC header not found");
          fhEntries->Fill(9);
          return kTRUE; //skip the event
      }
      return kFALSE;  //there was not outlier all jets have pT below fMaxFacPtHard * ptHard
   }else{
      fhEntries->Fill(9);
      AliWarning("AliAnalysisTaskEmcalJetBtagSV MC header not found");
      return kTRUE; //MC header not found
   }
}

//-------------------------------------newDeltaPt-----------------------------------
Int_t AliAnalysisTaskEmcalJetBtagSV::FillDeltaPt(Double_t rho, 
                                                 Int_t nVtx, 
                                                 AliAODVertex* pVtx, 
                                                 vctr_pair_dbl_int aVtxDisp, 
                                                 Double_t signalEta, 
                                                 Double_t signalPhi,  
                                                 Int_t fillMask){
   // fills delta pt for events with SV
   Int_t *idxLxy = new Int_t[nVtx];
   Double_t *sigmavertex     = new Double_t[nVtx];
   Double_t *decLenXY        = new Double_t[nVtx];
   Double_t *sigdecLenXY     = new Double_t[nVtx];
   Double_t sigmaSV, lxy;
   for (Int_t vtxID = 0; vtxID < nVtx; ++vtxID) {
      AliAODVertex *svtx = (AliAODVertex *)fHFvertexing->UncheckedAt(vtxID);  
      decLenXY[vtxID] = pVtx->DistanceXYToVertex(svtx);
      sigdecLenXY[vtxID] = decLenXY[vtxID]/pVtx->ErrorDistanceXYToVertex(svtx);;
      sigmavertex[vtxID] = aVtxDisp[vtxID].first;
   }
   TMath::Sort(nVtx, decLenXY, idxLxy);
   sigmaSV = sigmavertex[idxLxy[0]];
   lxy  = sigdecLenXY[idxLxy[0]];
  
   Double_t fillWeight; // 1 for Data mode, fMCWeight for for correction mode;
   if(fCorrMode){ 
      fillWeight = fMCWeight;
   }else{ 
      fillWeight = 1;
   }
   
   if(sigmaSV < fSigmaSVCut){
      if(lxy > 5 && (fillMask & 1) == 0){
         fhDeltaPtLxy5->Fill(GetDeltaPtRandomConeWithoutSignalPt(fTaggingRadius,rho, signalEta, signalPhi), fillWeight); 
         fillMask = fillMask | 1;
      }
      if(lxy > 6 && ( (fillMask & 2) >> 1) == 0){
         fhDeltaPtLxy6->Fill(GetDeltaPtRandomConeWithoutSignalPt(fTaggingRadius,rho, signalEta, signalPhi), fillWeight);
         fillMask = fillMask | 2;
      }       
      if(lxy > 7 && ((fillMask & 4) >> 2) == 0 ){
         fillMask = fillMask | 4; 
         fhDeltaPtLxy7->Fill(GetDeltaPtRandomConeWithoutSignalPt(fTaggingRadius,rho, signalEta, signalPhi), fillWeight);
//---------------------------------------------------EMB_clus

	if (fEmbeddPerpendicular){
		fFastJetWrapper->Clear();

		
		//----------------Generating NEW perpendicular track
		
		Double_t gen_pt = fTrackGenerator->Uniform(0,100);
		cout<<"!!!!!!!!!!!!!!!!!! gen_pt = "<<gen_pt<<endl;
           	TLorentzVector lVec;
                lVec.SetPtEtaPhiM(gen_pt,signalEta,signalPhi + TMath::Pi()/2,0);
		fFastJetWrapper->AddInputVector(lVec.Px(), lVec.Py(), lVec.Pz(), lVec.E(), -99999);
		//-----Filling old container

		
		 for(Int_t i = 0; i < fRecTrkArray->GetEntries(); i++) {
		    AliAODTrack* trk = static_cast<AliAODTrack*>(fRecTrkArray->ConstructedAt(i));
		    UInt_t trkFilterMap = trk->GetFilterMap();  
		    if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;
		    if ( (fabs(trk->Eta()) < fEtaCut))
			    fFastJetWrapper->AddInputVector(trk->Px(), trk->Py(), trk->Pz(), trk->P(), 1); 

		}  		 
	       


		fFastJetWrapper->Run();
		
		//--------------- DelPt analysis of the new container
		std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets();
		AliVParticle*  hytrk  = NULL;  //track hybrid event jet
		Double_t deltaPtEmb;
		Double_t sumTrkEmbeddedPt=0;
		for (UInt_t ijet = 0; ijet < jets_incl.size(); ++ijet) {
			std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));	

			sumTrkEmbeddedPt=0;
			for (UInt_t ic = 0; ic < constituents.size(); ++ic){ 
	      			if (constituents[ic].user_index() == -99999)  {
					sumTrkEmbeddedPt += constituents[ic].pt();
					break;
				}
	       		}
			
			if(sumTrkEmbeddedPt>0){ 
		   		deltaPtEmb = jets_incl.at(ijet).pt() - jets_incl.at(ijet).area() * rho - sumTrkEmbeddedPt;
		   		fhDeltaPtEmbeddPerpendicular->Fill(deltaPtEmb);
		   		fhDeltaPtEmbeddCorrelationPerpendicular->Fill(sumTrkEmbeddedPt, deltaPtEmb);
				break;
			}

		 }
		cout<<" All is done, default delPt: " << GetDeltaPtRandomConeWithoutSignalPt(fTaggingRadius,rho, signalEta, signalPhi) << " Updated One : " << deltaPtEmb<<endl;
	}

//----------------------------------------------------------

      }       
   }                  

   delete [] sigmavertex;
   delete [] decLenXY;
   delete [] sigdecLenXY;

   return fillMask;

}
//------------------------------------------------------------------------------------
Double_t AliAnalysisTaskEmcalJetBtagSV::GetDeltaPtRandomConeWithoutSignalPt (Double_t jetradius, 
                                                                             Double_t rhovalue, 
                                                                             Double_t signalEta, 
                                                                             Double_t signalPhi )
{
   Double_t minConeEta = jetradius - fEtaCut;
   Double_t maxConeEta = fEtaCut - jetradius;

    // throw random cone
   Double_t coneEta, conePhi, deta, dphi,dist;


   do{
       coneEta = minConeEta + fRandom->Rndm()*(maxConeEta - minConeEta);
       conePhi = fRandom->Rndm()*TMath::TwoPi();
       deta = coneEta - signalEta;
       dphi = TVector2::Phi_mpi_pi((conePhi - signalPhi));
       dist = sqrt(deta*deta + dphi*dphi);
    }
    while(dist<2*jetradius);
    

    // collect track pt within cone
    Double_t conePt = 0.;
    for (Int_t i = 0; i < fEvent->GetNumberOfTracks(); i++) {
        AliAODTrack* trk = static_cast<AliAODTrack*>(fEvent->GetTrack(i));

        // track filter hardwired...
        UInt_t trkFilterMap = trk->GetFilterMap();  
        if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;
            
        
        if ( (fabs(trk->Eta()) < fEtaCut) && (trk->Pt() > fPtCut) ) {
             dphi = TVector2::Phi_mpi_pi((trk->Phi() - conePhi));
             deta = trk->Eta() - coneEta;
             dist = sqrt(deta*deta + dphi*dphi);
            if (dist < jetradius) conePt += trk->Pt();
        }
    } // track loop

  return conePt - jetradius*jetradius*TMath::Pi() * rhovalue;
    
}   


    
