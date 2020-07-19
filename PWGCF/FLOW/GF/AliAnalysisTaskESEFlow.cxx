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

/* AliAnalysisTaskESEFlow
 * Author: Joachim Carlo Kristian Hansen, NBI 2019
 * Flow correlation class
 * Event Shape Engineering: TPC,V0C,V0A
 * Event selection
 * Track selection
 */


#include "AliAnalysisTaskESEFlow.h"

#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TList.h"
#include "TFile.h"
#include "TSpline.h"
#include "TMath.h"
#include "TComplex.h"
#include "TGrid.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliGFWWeights.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "TRandom3.h"

#include "AliUniFlowCorrTask.h"

#include <iostream>

class AliAnalysisTaskESEFlow;

ClassImp(AliAnalysisTaskESEFlow)

AliAnalysisTaskESEFlow::AliAnalysisTaskESEFlow() : AliAnalysisTaskSE(),
    fEventCuts(),
    fFlowRunByRunWeights(kTRUE),
    fV0RunByRunCalibration(kFALSE),
    bUseOwnWeights(0),
    dEtaGap(1),
    bHasGap(kTRUE),
    fSampling(kFALSE),
    fInit(kFALSE),
    fMakeqSelectionRun(kFALSE),
    fMakeRBRweightsRun(kFALSE),

    fAOD(0),
    fOutputList(0),
    fObservables(0),
    fCorrDist(0),
    fpTDiff(0),
    fqnDist(0),
    fpTDiffESETPC(0),
    fcnESETPC(0),
    fpTDiffESEV0C(0),
    fcnESEV0C(0),
    fpTDiffESEV0A(0),
    fcnESEV0A(0),
    SPFlowList(0),
    SPFlowEseList(0),
    fQAEvents(0),

    fFlowWeightsList{nullptr},
    fWeights(0),
    fV0CalibList(0),
    fqSelList(0),

    fHistPhiEtaVz(0),
    fHistPtEtaVz(0),
    fHistPhi(0),
    fHistEta(0),
    fHistPt(0),
    fHistZVertex(0),
    fHistPhiCor(0),
    fHistPhiCor3D(0),
    fHistTPCchi2(0),
    fHistITSchi2(0),
    fhQAEventsfMult32vsCentr(0), 
    fhQAEventsfMult128vsCentr(0),
    fhQAEventsfMult96vsCentr(0),
    fhQAEventsfMultTPCvsTOF(0),
    fhQAEventsfMultTPCvsESD(0),

    fSplq2TPC{0},
    fSplq3TPC{0},
    fSplq2V0C{0},
    fSplq3V0C{0},
    fSplq2V0A{0},
    fSplq3V0A{0},
    fh3Weights{nullptr},
    fhV0Calib{nullptr},
    
    fHistPDG{0},

    fFileTrackEff(0),
    fhTrackNUE(0),

    fq2TPC(0),
    fq3TPC(0),
    fq2V0C(0),
    fq3V0C(0),
    fq2V0A(0),
    fq3V0A(0),
    fQnxV0C{0},
    fQnyV0C{0},
    fQnxV0A{0},
    fQnyV0A{0},
    fQnxTPC{0},
    fQnyTPC{0},
    fQnxV0Cm{0},
    fQnyV0Cm{0},
    fQnxV0Am{0},
    fQnyV0Am{0},
    fQnxTPCm{0},
    fQnyTPCm{0},
    fQnxV0Cs{0},
    fQnyV0Cs{0},
    fQnxV0As{0},
    fQnyV0As{0},

    fQnxV0CEse{0},
    fQnyV0CEse{0},
    fQnxV0AEse{0},
    fQnyV0AEse{0},
    fQnxTPCEse{0},
    fQnyTPCEse{0},

    fQnxV0CCor{0},
    fQnyV0CCor{0},
    fQnxV0ACor{0},
    fQnyV0ACor{0},

    fhEvPlPsi_2V0C{0},
    fhEvPlPsi_2V0A{0},
    fhEvPlPsi_3V0C{0},
    fhEvPlPsi_3V0A{0},

    fProfNPar(0),
    fhV0Multiplicity(0),
    fhV0CorrMult(0),
    fhqnTPCvqnV0C{0},
    fhqnV0CvqnV0A{0},
    fhqnTPCvqnV0A{0},
    fhEventCounter{0},

    fColSystem{kPbPb},
    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fAbsEtaMax(0.8),
    fVtxZCuts(10.0),
    fNPhiBins(120),
    fNEtaBins(32),
    fCentEstimator("V0M"),
    fCutChargedNumTPCclsMin(70),
    fReadMC(kFALSE),
    fMCEvent(0),
    fFlowRFPsPtMin(0.2),
    fFlowRFPsPtMax(5.0),
    fFlowPOIsPtMin(0.0),
    fFlowPOIsPtMax(10.0),
    fRedFlowPtMin(0.0),
    fRedFlowPtMax(20.0),

    fPtAxis(new TAxis()),
    fCentAxis(new TAxis()),
    ESEPercAxis(new TAxis()),
    
    fTPCEse(kTRUE),
    fV0CEse(kTRUE),
    fV0AEse(kFALSE),
    fIndexSampling{0},
    fNumSamples{1},
    fSPAnalysis(kFALSE),

    TPCqnBins(100),
    TPCqnBinMin(0.0),
    TPCqnBinMax(8.0),

    V0qnBins(100),
    V0qnBinMin(0.0),
    V0qnBinMax(15.0),

    fPileupCut(500),
    fCheckChi2TPC(kFALSE),
    fCheckChi2ITS(kFALSE),
    vTPCChi2Bound(4.0),
    vITSChi2Bound(36.0),
    fFillQARej(kFALSE),

    fUseNUEWeights(kFALSE),
    fNUE(1),
    fIs2018Data(kFALSE),
    fBayesUnfoldingInput(kFALSE),
    fActq2Projections(kFALSE),


    fVecCorrTask()
{}
//_____________________________________________________________________________
AliAnalysisTaskESEFlow::AliAnalysisTaskESEFlow(const char* name, ColSystem colSys, Bool_t bUseV0Calibration) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fFlowRunByRunWeights(kTRUE),
    fV0RunByRunCalibration(bUseV0Calibration),
    bUseOwnWeights(0),
    dEtaGap(1),
    bHasGap(kTRUE),
    fSampling(kFALSE),
    fInit(kFALSE),
    fMakeqSelectionRun(kFALSE),
    fMakeRBRweightsRun(kFALSE),

    fAOD(0),
    fOutputList(0),
    fObservables(0),
    fCorrDist(0),
    fpTDiff(0),
    fqnDist(0),
    fpTDiffESETPC(0),
    fcnESETPC(0),
    fpTDiffESEV0C(0),
    fcnESEV0C(0),
    fpTDiffESEV0A(0),
    fcnESEV0A(0),
    SPFlowList(0),
    SPFlowEseList(0),
    fQAEvents(0),

    fFlowWeightsList{nullptr},
    fWeights(0),
    fV0CalibList(0),
    fqSelList(0),

    fHistPhiEtaVz(0),
    fHistPtEtaVz(0),
    fHistPhi(0),
    fHistEta(0),
    fHistPt(0),
    fHistZVertex(0),
    fHistPhiCor(0),
    fHistPhiCor3D(0),
    fHistTPCchi2(0),
    fHistITSchi2(0),
    fhQAEventsfMult32vsCentr(0), 
    fhQAEventsfMult128vsCentr(0),
    fhQAEventsfMult96vsCentr(0),
    fhQAEventsfMultTPCvsTOF(0),
    fhQAEventsfMultTPCvsESD(0),

    fSplq2TPC{0},
    fSplq3TPC{0},
    fSplq2V0C{0},
    fSplq3V0C{0},
    fSplq2V0A{0},
    fSplq3V0A{0},
    fh3Weights{nullptr},
    fhV0Calib{nullptr},

    fHistPDG{0},

    fFileTrackEff(0),
    fhTrackNUE(0),

    fq2TPC(0),
    fq3TPC(0),
    fq2V0C(0),
    fq3V0C(0),
    fq2V0A(0),
    fq3V0A(0),
    fQnxV0C{0},
    fQnyV0C{0},
    fQnxV0A{0},
    fQnyV0A{0},
    fQnxTPC{0},
    fQnyTPC{0},
    fQnxV0Cm{0},
    fQnyV0Cm{0},
    fQnxV0Am{0},
    fQnyV0Am{0},
    fQnxTPCm{0},
    fQnyTPCm{0},
    fQnxV0Cs{0},
    fQnyV0Cs{0},
    fQnxV0As{0},
    fQnyV0As{0},

    fQnxV0CEse{0},
    fQnyV0CEse{0},
    fQnxV0AEse{0},
    fQnyV0AEse{0},
    fQnxTPCEse{0},
    fQnyTPCEse{0},

    fQnxV0CCor{0},
    fQnyV0CCor{0},
    fQnxV0ACor{0},
    fQnyV0ACor{0},

    fhEvPlPsi_2V0C{0},
    fhEvPlPsi_2V0A{0},
    fhEvPlPsi_3V0C{0},
    fhEvPlPsi_3V0A{0},

    fProfNPar(0),
    fhV0Multiplicity(0),
    fhV0CorrMult(0),
    fhqnTPCvqnV0C{0},
    fhqnV0CvqnV0A{0},
    fhqnTPCvqnV0A{0},
    fhEventCounter{0},

    fColSystem{colSys},
    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fAbsEtaMax(0.8),
    fVtxZCuts(10.0),
    fNPhiBins(120),
    fNEtaBins(32),
    fCentEstimator("V0M"),
    fCutChargedNumTPCclsMin(70),
    fReadMC(kFALSE),
    fMCEvent(0),
    fFlowRFPsPtMin(0.2),
    fFlowRFPsPtMax(5.0),
    fFlowPOIsPtMin(0.0),
    fFlowPOIsPtMax(10.0),
    fRedFlowPtMin(0.0),
    fRedFlowPtMax(20.0),

    fPtAxis(new TAxis()),
    fCentAxis(new TAxis()),
    ESEPercAxis(new TAxis()),
    
    fTPCEse(kTRUE),
    fV0CEse(kTRUE),
    fV0AEse(kFALSE),
    fIndexSampling{0},
    fNumSamples{1},
    fSPAnalysis(kFALSE),

    TPCqnBins(100),
    TPCqnBinMin(0.0),
    TPCqnBinMax(8.0),

    V0qnBins(100),
    V0qnBinMin(0.0),
    V0qnBinMax(15.0),

    fPileupCut(500),
    fCheckChi2TPC(kFALSE),
    fCheckChi2ITS(kFALSE),
    vTPCChi2Bound(4.0),
    vITSChi2Bound(36.0),
    fFillQARej(kFALSE),

    fUseNUEWeights(kFALSE),
    fNUE(1),
    fIs2018Data(kFALSE),
    fBayesUnfoldingInput(kFALSE),
    fActq2Projections(kFALSE),

    fVecCorrTask()
{
    DefineInput(0, TChain::Class());
    DefineInput(1, TList::Class());
    DefineInput(2, TList::Class());
    DefineInput(3, TList::Class());

    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
    DefineOutput(4, TList::Class());
    DefineOutput(5, TList::Class());
    DefineOutput(6, TList::Class());
    DefineOutput(7, TList::Class());
    DefineOutput(8, TList::Class());
    DefineOutput(9, TList::Class());
    DefineOutput(10, TList::Class());
    DefineOutput(11, TList::Class());
    DefineOutput(12, TList::Class());
    DefineOutput(13, TList::Class());
    DefineOutput(14, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskESEFlow::~AliAnalysisTaskESEFlow()
{
    if(fOutputList) { delete fOutputList; }
    if(fObservables) { delete fObservables; }
    if(fCorrDist) { delete fCorrDist; }
    if(fpTDiff) { delete fpTDiff; }
    if(fqnDist) { delete fqnDist; }
    if(fpTDiffESETPC) { delete fpTDiffESETPC; }
    if(fcnESETPC) { delete fcnESETPC; }
    if(fpTDiffESEV0C) { delete fpTDiffESEV0C; }
    if(fcnESEV0C) { delete fcnESEV0C; }
    if(fpTDiffESEV0A) { delete fpTDiffESEV0A; }
    if(fcnESEV0A) { delete fcnESEV0A; }
    if(SPFlowList) { delete SPFlowList; }
    if(SPFlowEseList) { delete SPFlowEseList; }
    if(fQAEvents) { delete fQAEvents; }
}
Bool_t AliAnalysisTaskESEFlow::InitializeTask()
{
    if(!fMakeRBRweightsRun){
        if(bUseOwnWeights)
        {
            fFlowWeightsList = static_cast<TList*>(GetInputData(1));
            if(!fFlowWeightsList) { AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n Flow weights list not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }
        }
        else
        {
            fFlowWeightsList = static_cast<TList*>(GetInputData(1));
            if(!fFlowWeightsList) { AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n Flow weights list 2 not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }
        }
    }

    //Load V0 Calibration

    if(fV0RunByRunCalibration){
        fV0CalibList = static_cast<TList*>(GetInputData(2));
        if(!fV0CalibList) { AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n V0 Calibration list not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }
    }

    // load q-splines
    if(!fMakeqSelectionRun){
        fqSelList = static_cast<TList*>(GetInputData(3));
        if(!fqSelList) {AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n q-selection Splines list not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }

        if(!LoadqSelection()) { AliFatal("\n \n \n \n \n \n \n \n \n \n q-Splines not loaded! Terminating! \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }
    }

    if(fSampling && fNumSamples < 2){
        AliFatal("Sampling used, but number of samples < 2! Terminating!");
        return kFALSE;
    }
    
    AliInfo("Initialization of weights and calibration succes");
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskESEFlow::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fObservables = new TList();
    fCorrDist = new TList();
    fpTDiff = new TList();
    fqnDist = new TList();
    fpTDiffESETPC = new TList();
    fcnESETPC = new TList();
    fpTDiffESEV0C = new TList();
    fcnESEV0C = new TList();
    fpTDiffESEV0A = new TList();
    fcnESEV0A = new TList();
    SPFlowList = new TList();
    SPFlowEseList = new TList();
    fQAEvents = new TList();

    fOutputList->SetOwner(kTRUE);
    fObservables->SetOwner(kTRUE);
    fCorrDist->SetOwner(kTRUE);
    fpTDiff->SetOwner(kTRUE);
    fqnDist->SetOwner(kTRUE);
    fpTDiffESETPC->SetOwner(kTRUE);
    fcnESETPC->SetOwner(kTRUE);
    fpTDiffESEV0C->SetOwner(kTRUE);
    fcnESEV0C->SetOwner(kTRUE);
    fpTDiffESEV0A->SetOwner(kTRUE);
    fcnESEV0A->SetOwner(kTRUE);
    SPFlowList->SetOwner(kTRUE);
    SPFlowEseList->SetOwner(kTRUE);
    fQAEvents->SetOwner(kTRUE);

    //RUN INITIALIZE TASK
    fInit = InitializeTask();
    if(!fInit) { return; }

    Int_t nCentBin = fCentAxis->GetNbins();
    Int_t nPtBin = fPtAxis->GetNbins();
    for(Int_t Ptbin(0); Ptbin < nPtBin+1; ++Ptbin){
        PtEdges[Ptbin] = fPtAxis->GetBinLowEdge(Ptbin+1);
    }
    for(Int_t Centbin(0); Centbin < nCentBin+1; ++Centbin){
        CentEdges[Centbin] = fCentAxis->GetBinLowEdge(Centbin+1);
    }

    // Event Shape
    //ESEPercAxis
    //printf("Number ESE percs: ");
    //std::cout << ESEPercAxis->GetNbins() << std::endl;
    for (Int_t Perc(0); Perc < ESEPercAxis->GetNbins()+1; ++Perc){
        EventShapeEdges[Perc] = ESEPercAxis->GetBinLowEdge(Perc+1);
    }
    

    fHistPhiEtaVz = new TH3F("fHistPhiEtaVz", "fHistPhiEtaVz; #phi; #eta; Vz", fNPhiBins, 0.0, TMath::TwoPi(), fNEtaBins, -1.0, 1.0,fVtxZCuts*2,-fVtxZCuts,fVtxZCuts);
    fHistPhiEtaVz->Sumw2();
    fHistPtEtaVz = new TH3F("fHistPtEtaVz","fHistPtEtaVz; pT; #eta; Vz", 100, 0.0, 10.0, fNEtaBins, -1.0, 1.0, fVtxZCuts*2,-fVtxZCuts,fVtxZCuts);
    fHistPtEtaVz->Sumw2();
    fHistPhi = new TH1F("fHistPhi", ";#phi", fNPhiBins, 0.0, TMath::TwoPi());
    fHistEta = new TH1F("fHistEta", ";#eta", fNEtaBins,-1.0, 1.0);
    fHistPt = new TH1F("fHistPt", ";p_{T}", nPtBin,PtEdges);
    fHistZVertex = new TH1F("fHistZVertex", ";Vtx_{Z}", fVtxZCuts,-10,10);
    fHistPhiCor = new TH1F("fHistPhiCor",";#phiCor", fNPhiBins, 0.0, TMath::TwoPi());
    fHistPhiCor->Sumw2();
    fHistPhiCor3D = new TH3F("fHistPhiEtaVz_NUACorr", "fHistPhiEtaVz_NUACorr; #phi; #eta; Vz", fNPhiBins, 0.0, TMath::TwoPi(), fNEtaBins, -1.0, 1.0,fVtxZCuts*2,-fVtxZCuts,fVtxZCuts);
    fHistPhiCor3D->Sumw2();
    fHistTPCchi2 = new TH1F("fHistChi2TPC","fHistChi2TPC",10,0,10);
    fHistITSchi2 = new TH1F("fHistChi2ITS","fHistChi2ITS",40,0,40);
    
    fProfNPar = new TProfile("fProfNparvsCent",";Centrality;N_{Particles}",100,0,100);

    fhV0Multiplicity = new TH2F("fV0Multiplicity","",64,0,64,100,0,1250);
    fhV0Multiplicity->Sumw2();

    fhV0CorrMult = new TH2F("fV0CalibratedMultiplicity","",64,0,64,100,0,1250);
    fhV0CorrMult->Sumw2();

    

    fq2TPC = new TH2D("fq2vCentTPC","",100,0,100,TPCqnBins,TPCqnBinMin,TPCqnBinMax);
    fq2TPC->Sumw2();
    fq3TPC = new TH2D("fq3vCentTPC","",100,0,100,TPCqnBins,TPCqnBinMin,TPCqnBinMax);
    fq3TPC->Sumw2();

    fq2V0C = new TH2D("fq2vCentV0C","",100,0,100,V0qnBins,V0qnBinMin,V0qnBinMax);
    fq2V0C->Sumw2();
    fq3V0C = new TH2D("fq3vCentV0C","",100,0,100,V0qnBins,V0qnBinMin,V0qnBinMax);
    fq3V0C->Sumw2();
    fq2V0A = new TH2D("fq2vCentV0A","",100,0,100,V0qnBins,V0qnBinMin,V0qnBinMax);
    fq2V0A->Sumw2();
    fq3V0A = new TH2D("fq3vCentV0A","",100,0,100,V0qnBins,V0qnBinMin,V0qnBinMax);
    fq3V0A->Sumw2();


    if(fActq2Projections){
      for (Int_t iEse(0);iEse<ESEPercAxis->GetNbins();++iEse){
        TH2* fq2TPC_perc = nullptr;
        TH2* fq2V0C_perc = nullptr;
        TH2* fq2V0A_perc = nullptr;

        fq2TPC_perc = new TH2D(Form("fq2vCentTPC_Perc%i",iEse+1),Form("fq2vCentTPC_Perc%i",iEse+1),100,0,100,TPCqnBins,TPCqnBinMin,TPCqnBinMax);
        fq2TPC_perc->Sumw2();
        fqnDist->Add(fq2TPC_perc);

        fq2V0C_perc = new TH2D(Form("fq2vCentV0C_Perc%i",iEse+1),Form("fq2vCentV0C_Perc%i",iEse+1),100,0,100,V0qnBins,V0qnBinMin,V0qnBinMax);
        fq2V0C_perc->Sumw2();
        fqnDist->Add(fq2V0C_perc);

        fq2V0A_perc = new TH2D(Form("fq2vCentV0A_Perc%i",iEse+1),Form("fq2vCentV0A_Perc%i",iEse+1),100,0,100,V0qnBins,V0qnBinMin,V0qnBinMax);
        fq2V0A_perc->Sumw2();
        fqnDist->Add(fq2V0A_perc);
      }
    }



    for (Int_t qi(0);qi<2;++qi){
        fQnxV0C[qi] = new TH2F(Form("fQ%ixvCentV0C",qi+2),"",100,0,100,100,-1500,1500);
        fQnxV0C[qi]->Sumw2();
        fQnyV0C[qi] = new TH2F(Form("fQ%iyvCentV0C",qi+2),"",100,0,100,100,-1500,1500);
        fQnyV0C[qi]->Sumw2();
        fQnxV0A[qi] = new TH2F(Form("fQ%ixvCentV0A",qi+2),"",100,0,100,100,-1500,1500);
        fQnxV0A[qi]->Sumw2();
        fQnyV0A[qi] = new TH2F(Form("fQ%iyvCentV0A",qi+2),"",100,0,100,100,-1500,1500);
        fQnyV0A[qi]->Sumw2();
        fQnxTPC[qi] = new TH2F(Form("fQ%ixvCentTPC",qi+2),"",100,0,100,100,-1500,1500);
        fQnxTPC[qi]->Sumw2();
        fQnyTPC[qi] = new TH2F(Form("fQ%iyvCentTPC",qi+2),"",100,0,100,100,-1500,1500);
        fQnyTPC[qi]->Sumw2();

        fQnxV0CEse[qi] = new TH2F(Form("fQ%ixvCentV0CEse",qi+2),"",100,0,100,100,-1500,1500);
        fQnxV0CEse[qi]->Sumw2();
        fQnyV0CEse[qi] = new TH2F(Form("fQ%iyvCentV0CEse",qi+2),"",100,0,100,100,-1500,1500);
        fQnyV0CEse[qi]->Sumw2();
        fQnxV0AEse[qi] = new TH2F(Form("fQ%ixvCentV0AEse",qi+2),"",100,0,100,100,-1500,1500);
        fQnxV0AEse[qi]->Sumw2();
        fQnyV0AEse[qi] = new TH2F(Form("fQ%iyvCentV0AEse",qi+2),"",100,0,100,100,-1500,1500);
        fQnyV0AEse[qi]->Sumw2();
        fQnxTPCEse[qi] = new TH2F(Form("fQ%ixvCentTPCEse",qi+2),"",100,0,100,100,-1500,1500);
        fQnxTPCEse[qi]->Sumw2();
        fQnyTPCEse[qi] = new TH2F(Form("fQ%iyvCentTPCEse",qi+2),"",100,0,100,100,-1500,1500);
        fQnyTPCEse[qi]->Sumw2();        

        fQnxV0CCor[qi] = new TH2F(Form("fQ%ixvCentV0CCor",qi+2),"",100,0,100,100,-10,10);
        fQnxV0CCor[qi]->Sumw2();
        fQnyV0CCor[qi] = new TH2F(Form("fQ%iyvCentV0CCor",qi+2),"",100,0,100,100,-10,10);
        fQnyV0CCor[qi]->Sumw2();
        fQnxV0ACor[qi] = new TH2F(Form("fQ%ixvCentV0ACor",qi+2),"",100,0,100,100,-10,10);
        fQnxV0ACor[qi]->Sumw2();
        fQnyV0ACor[qi] = new TH2F(Form("fQ%iyvCentV0ACor",qi+2),"",100,0,100,100,-10,10);
        fQnyV0ACor[qi]->Sumw2();



        fhqnTPCvqnV0C[qi] = new TH2F(Form("fq%iTPCvq%iV0C",qi+2,qi+2),"",100,0,16,100,0,16);
        fhqnTPCvqnV0C[qi]->Sumw2();
        fhqnV0CvqnV0A[qi] = new TH2F(Form("fq%iV0Cvq%iV0A",qi+2,qi+2),"",100,0,16,100,0,16);
        fhqnV0CvqnV0A[qi]->Sumw2();
        fhqnTPCvqnV0A[qi] = new TH2F(Form("fq%iTPCvq%iV0A",qi+2,qi+2),"",100,0,16,100,0,16);
        fhqnTPCvqnV0A[qi]->Sumw2();
    }

    fhEvPlPsi_2V0C = new TH2F("EvPlPsi_2_V0C","EvPlPsi_2_V0C",100,0,100,100,-2,2);
    fhEvPlPsi_2V0C->Sumw2();
    fqnDist->Add(fhEvPlPsi_2V0C);
    fhEvPlPsi_2V0A = new TH2F("EvPlPsi_2_V0A","EvPlPsi_2_V0A",100,0,100,100,-2,2);
    fhEvPlPsi_2V0A->Sumw2();
    fqnDist->Add(fhEvPlPsi_2V0A);

    fhEvPlPsi_3V0C = new TH2F("EvPlPsi_3_V0C","EvPlPsi_3_V0C",100,0,100,100,-2,2);
    fhEvPlPsi_3V0C->Sumw2();
    fqnDist->Add(fhEvPlPsi_3V0C);
    fhEvPlPsi_3V0A = new TH2F("EvPlPsi_3_V0A","EvPlPsi_3_V0A",100,0,100,100,-2,2);
    fhEvPlPsi_3V0A->Sumw2();
    fqnDist->Add(fhEvPlPsi_3V0A);


    for (Int_t iQnR(0);iQnR<2;++iQnR){
        TH1* tempAvgQn_TPC = nullptr;
        TH1* tempAvgQn_V0C = nullptr;
        TH1* tempAvgQn_V0A = nullptr;

        tempAvgQn_TPC = new TProfile(Form("fPq_%i_TPC",iQnR+2),"avgqn",100,0,100);
        tempAvgQn_TPC->Sumw2();
        fOutputList->Add(tempAvgQn_TPC);

        tempAvgQn_V0C = new TProfile(Form("fPq_%i_V0C",iQnR+2),"avgqn",100,0,100);
        tempAvgQn_V0C->Sumw2();
        fOutputList->Add(tempAvgQn_V0C);

        tempAvgQn_V0A = new TProfile(Form("fPq_%i_V0A",iQnR+2),"avgqn",100,0,100);
        tempAvgQn_V0A->Sumw2();
        fOutputList->Add(tempAvgQn_V0A);

    }


    Int_t iSizeTask = fVecCorrTask.size();

    for (Int_t iTask(0); iTask<iSizeTask; ++iTask){
        AliUniFlowCorrTask* task = fVecCorrTask.at(iTask);
        if(!task) { AliError(Form("AliUniFlowCorrTask%d does not exist",iTask)); return; }

        const char* CorrName = task->fsName.Data();
        const char* CorrLabel = task->fsLabel.Data();
        Int_t CorrOrder = task->fiNumHarm;

        for (Int_t iSample(0); iSample < fNumSamples; ++iSample){

            if(iSample > 0 && !fSampling) { break; }

            TH1* cn = nullptr;
            TH1* dn = nullptr;

            TH1* cnESETPC = nullptr;
            TH1* dnESETPC = nullptr;
            TH1* cnESEV0C = nullptr;
            TH1* dnESEV0C = nullptr;
            TH1* cnESEV0A = nullptr;
            TH1* dnESEV0A = nullptr;

            

            cn = new TProfile(Form("%s_sample%d",CorrName,iSample),Form("%s",CorrLabel),nCentBin,CentEdges);

            if(!cn) { AliError("Centrality profile not created"); task->PrintTask(); return; }
            if(fCorrDist->FindObject(cn->GetName())) {
                AliError(Form("Task %d: Profile '%s' already exists",iTask,cn->GetName()));
                task->PrintTask();
                delete cn;
                return;
            }

            cn->Sumw2();
            fCorrDist->Add(cn);



            for (Int_t qi(0);qi<2;++qi){
                for (Int_t iEse(0);iEse<ESEPercAxis->GetNbins();++iEse){ //10
                    if(fTPCEse){
                        cnESETPC = new TProfile(Form("%s_q%iTPC_PerCode%i_sample%d",CorrName,qi+2,iEse+1,iSample),Form("%s_q%iTPCPerCode%i",CorrLabel,qi+2,iEse+1),nCentBin,CentEdges);

                        if(!cnESETPC) { AliError("ESETPC Centrality profile not created"); task->PrintTask(); return; }
                        if(fcnESETPC->FindObject(cnESETPC->GetName())) {
                            AliError(Form("Task %d: Profile '%s' already exists",iTask,cnESETPC->GetName()));
                            task->PrintTask();
                            delete cnESETPC;
                            return;
                        }

                        cnESETPC->Sumw2();
                        fcnESETPC->Add(cnESETPC);
                    }

                    if(fV0CEse){
                        cnESEV0C = new TProfile(Form("%s_q%iV0C_PerCode%i_sample%d",CorrName,qi+2,iEse+1,iSample),Form("%s_q%iV0CPerCode%i",CorrLabel,qi+2,iEse+1),nCentBin,CentEdges);

                        if(!cnESEV0C) { AliError("ESEV0C Centrality profile not created"); task->PrintTask(); return; }
                        if(fcnESEV0C->FindObject(cnESEV0C->GetName())) {
                            AliError(Form("Task %d: Profile '%s' already exists",iTask,cnESEV0C->GetName()));
                            task->PrintTask();
                            delete cnESEV0C;
                            return;
                        }

                        cnESEV0C->Sumw2();
                        fcnESEV0C->Add(cnESEV0C);
                    }

                    if(fV0AEse){
                        cnESEV0A = new TProfile(Form("%s_q%iV0A_PerCode%i_sample%d",CorrName,qi+2,iEse+1,iSample),Form("%s_q%iV0APerCode%i",CorrLabel,qi+2,iEse+1),nCentBin,CentEdges);

                        if(!cnESEV0A) { AliError("ESEV0A Centrality profile not created"); task->PrintTask(); return; }
                        if(fcnESEV0A->FindObject(cnESEV0A->GetName())) {
                            AliError(Form("Task %d: Profile '%s' already exists",iTask,cnESEV0A->GetName()));
                            task->PrintTask();
                            delete cnESEV0A;
                            return;
                        }

                        cnESEV0A->Sumw2();
                        fcnESEV0A->Add(cnESEV0A);
                    }
                }
            }


            //pt differentials
            if (CorrOrder < 6) {
                for(Int_t fCentNum(0) ; fCentNum<nCentBin; ++fCentNum){


                    dn = new TProfile(Form("%s_diff_%.0f_%.0f_sample%d",CorrName,CentEdges[fCentNum],CentEdges[fCentNum+1],iSample),Form("%s_diff_%.0f_%.0f",CorrLabel,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);

                    if(!dn) { AliError("pt differential profile not created"); task->PrintTask(); return; }
                    if(fpTDiff->FindObject(dn->GetName())) {
                        AliError(Form("Task %d: Profile '%s' already exists",iTask,dn->GetName()));
                        task->PrintTask();
                        delete dn;
                        return;
                    }

                    dn->Sumw2();
                    fpTDiff->Add(dn);


                    for (Int_t qi(0);qi<2;++qi){
                        for (Int_t iEse(0);iEse<ESEPercAxis->GetNbins();++iEse){

                            if(fTPCEse){
                                dnESETPC = new TProfile(Form("%s_diff_q%iTPC_PerCode%i_%.0f_%.0f_sample%d",CorrName,qi+2,iEse+1,CentEdges[fCentNum],CentEdges[fCentNum+1],iSample),Form("%s_q%iTPCPerCode%i_%.0f_%.0f",CorrLabel,qi+2,iEse+1,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);

                                if(!dnESETPC) { AliError("ESETPC pt diff profile not created"); task->PrintTask(); return; }
                                if(fpTDiffESETPC->FindObject(dnESETPC->GetName())) {
                                    AliError(Form("Task %d: Profile '%s' already exists",iTask,dnESETPC->GetName()));
                                    task->PrintTask();
                                    delete dnESETPC;
                                    return;
                                }

                                dnESETPC->Sumw2();
                                fpTDiffESETPC->Add(dnESETPC);
                            }

                            if(fV0CEse){
                                dnESEV0C = new TProfile(Form("%s_diff_q%iV0C_PerCode%i_%.0f_%.0f_sample%d",CorrName,qi+2,iEse+1,CentEdges[fCentNum],CentEdges[fCentNum+1],iSample),Form("%s_q%iV0CPerCode%i_%.0f_%.0f",CorrLabel,qi+2,iEse+1,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);

                                if(!dnESEV0C) { AliError("ESEV0C pt diff profile not created"); task->PrintTask(); return; }
                                if(fpTDiffESEV0C->FindObject(dnESEV0C->GetName())) {
                                    AliError(Form("Task %d: Profile '%s' already exists",iTask,dnESEV0C->GetName()));
                                    task->PrintTask();
                                    delete dnESEV0C;
                                    return;
                                }

                                dnESEV0C->Sumw2();
                                fpTDiffESEV0C->Add(dnESEV0C);
                            }

                            if(fV0AEse){
                                dnESEV0A = new TProfile(Form("%s_diff_q%iV0A_PerCode%i_%.0f_%.0f_sample%d",CorrName,qi+2,iEse+1,CentEdges[fCentNum],CentEdges[fCentNum+1],iSample),Form("%s_q%iV0APerCode%i_%.0f_%.0f",CorrLabel,qi+2,iEse+1,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);

                                if(!dnESEV0A) { AliError("ESEV0A pt diff profile not created"); task->PrintTask(); return; }
                                if(fpTDiffESEV0A->FindObject(dnESEV0A->GetName())) {
                                    AliError(Form("Task %d: Profile '%s' already exists",iTask,dnESEV0A->GetName()));
                                    task->PrintTask();
                                    delete dnESEV0A;
                                    return;
                                }

                                dnESEV0A->Sumw2();
                                fpTDiffESEV0A->Add(dnESEV0A);
                            }
                        }
                    }
                }
            }
        }
    }

    /// SCALAR PRODUCT ---- WORK IN PROGRESS
    if(fSPAnalysis){
        for (Int_t iVN(0);iVN<2;++iVN){
            TH1* vnSPV0A = nullptr;
            TH1* vnSPV0C = nullptr;
            TH1* QQSPV0ATPC = nullptr;
            TH1* QQSPV0AV0C = nullptr;
            TH1* QQSPTPCV0C = nullptr;  

            TH1* vnSPV0A_q2V0C = nullptr;
            TH1* vnSPV0C_q2V0C = nullptr;
            TH1* QQSPV0ATPC_q2V0C = nullptr;
            TH1* QQSPV0AV0C_q2V0C = nullptr;
            TH1* QQSPTPCV0C_q2V0C = nullptr;        

            
            for(Int_t fCentNum(0) ; fCentNum<nCentBin; ++fCentNum){
                vnSPV0A = new TProfile(Form("v_{%i}{SP V0A}_%.0f_%.0f",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1]),Form("v_{%i}{SP V0A}_%.0f_%.0f",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);
                if(!vnSPV0A) { AliError("Centrality profile not created"); return; }
                if(SPFlowList->FindObject(vnSPV0A->GetName())) {
                    AliError(Form("Profile '%s' already exists",vnSPV0A->GetName()));
                    delete vnSPV0A;
                    return;
                }
                vnSPV0A->Sumw2();
                SPFlowList->Add(vnSPV0A); //make new folders for SP method

                vnSPV0C = new TProfile(Form("v_{%i}{SP V0C}_%.0f_%.0f",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1]),Form("v_{%i}{SP V0C}_%.0f_%.0f",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);
                if(!vnSPV0C) { AliError("Centrality profile not created"); return; }
                if(SPFlowList->FindObject(vnSPV0C->GetName())) {
                    AliError(Form("Profile '%s' already exists",vnSPV0C->GetName()));
                    delete vnSPV0C;
                    return;
                }
                vnSPV0C->Sumw2();
                SPFlowList->Add(vnSPV0C); //make new folders for SP method

                for (Int_t iEse(0);iEse<ESEPercAxis->GetNbins();++iEse){
                    vnSPV0A_q2V0C = new TProfile(Form("v_{%i}{SP V0A}_%.0f_%.0f_PerCode%i",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1],iEse+1),Form("v_{%i}{SP V0A}_%.0f_%.0f",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);
                    if(!vnSPV0A_q2V0C) { AliError("Centrality profile not created"); return; }
                    if(SPFlowEseList->FindObject(vnSPV0A_q2V0C->GetName())) {
                        AliError(Form("Profile '%s' already exists",vnSPV0A_q2V0C->GetName()));
                        delete vnSPV0A_q2V0C;
                        return;
                    }
                    vnSPV0A_q2V0C->Sumw2();
                    SPFlowEseList->Add(vnSPV0A_q2V0C); //make new folders for SP method

                    vnSPV0C_q2V0C = new TProfile(Form("v_{%i}{SP V0C}_%.0f_%.0f_PerCode%i",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1],iEse+1),Form("v_{%i}{SP V0C}_%.0f_%.0f",iVN+2,CentEdges[fCentNum],CentEdges[fCentNum+1]),nPtBin,PtEdges);
                    if(!vnSPV0C_q2V0C) { AliError("Centrality profile not created"); return; }
                    if(SPFlowEseList->FindObject(vnSPV0C_q2V0C->GetName())) {
                        AliError(Form("Profile '%s' already exists",vnSPV0C_q2V0C->GetName()));
                        delete vnSPV0C_q2V0C;
                        return;
                    }
                    vnSPV0C_q2V0C->Sumw2();
                    SPFlowEseList->Add(vnSPV0C_q2V0C); //make new folders for SP method
                }
            }

            QQSPV0ATPC = new TProfile(Form("<QQ*>_%i{SP_V0A_TPC}",iVN+2),Form("<QQ*>_%i{SP_V0A_TPC}",iVN+2),nCentBin,CentEdges);
            if(!QQSPV0ATPC) { AliError("<QQ*> V0A TPC profile not created"); return; }
            if(SPFlowList->FindObject(QQSPV0ATPC->GetName())) {
                AliError(Form("Profile '%s' already exists",QQSPV0ATPC->GetName()));
                delete QQSPV0ATPC;
                return;
            }
            QQSPV0ATPC->Sumw2();
            SPFlowList->Add(QQSPV0ATPC);

            QQSPV0AV0C = new TProfile(Form("<QQ*>_%i{SP_V0A_V0C}",iVN+2),Form("<QQ*>_%i{SP_V0A_V0C}",iVN+2),nCentBin,CentEdges);
            if(!QQSPV0AV0C) { AliError("<QQ*> V0A V0C profile not created"); return; }
            if(SPFlowList->FindObject(QQSPV0AV0C->GetName())) {
                AliError(Form("Profile '%s' already exists",QQSPV0AV0C->GetName()));
                delete QQSPV0AV0C;
                return;
            }
            QQSPV0AV0C->Sumw2();
            SPFlowList->Add(QQSPV0AV0C);

            QQSPTPCV0C = new TProfile(Form("<QQ*>_%i{SP_V0C_TPC}",iVN+2),Form("<QQ*>_%i{SP_V0C_TPC}",iVN+2),nCentBin,CentEdges);
            if(!QQSPTPCV0C) { AliError("<QQ*> V0C TPC profile not created"); return; }
            if(SPFlowList->FindObject(QQSPTPCV0C->GetName())) {
                AliError(Form("Profile '%s' already exists",QQSPTPCV0C->GetName()));
                delete QQSPTPCV0C;
                return;
            }
            QQSPTPCV0C->Sumw2();
            SPFlowList->Add(QQSPTPCV0C);

            /// SP ESE ///
            for (Int_t iEse(0);iEse<ESEPercAxis->GetNbins();++iEse){ //10
                QQSPV0ATPC_q2V0C = new TProfile(Form("<QQ*>_%i{SP_V0A_TPC}_PerCode%i",iVN+2,iEse+1),Form("<QQ*>_%i{SP_V0A_TPC}",iVN+2),nCentBin,CentEdges);
                if(!QQSPV0ATPC_q2V0C) { AliError("<QQ*> V0A TPC profile not created"); return; }
                if(SPFlowEseList->FindObject(QQSPV0ATPC_q2V0C->GetName())) {
                    AliError(Form("Profile '%s' already exists",QQSPV0ATPC_q2V0C->GetName()));
                    delete QQSPV0ATPC_q2V0C;
                    return;
                }
                QQSPV0ATPC_q2V0C->Sumw2();
                SPFlowEseList->Add(QQSPV0ATPC_q2V0C);

                QQSPV0AV0C_q2V0C = new TProfile(Form("<QQ*>_%i{SP_V0A_V0C}_PerCode%i",iVN+2,iEse+1),Form("<QQ*>_%i{SP_V0A_V0C}",iVN+2),nCentBin,CentEdges);
                if(!QQSPV0AV0C_q2V0C) { AliError("<QQ*> V0A V0C profile not created"); return; }
                if(SPFlowEseList->FindObject(QQSPV0AV0C_q2V0C->GetName())) {
                    AliError(Form("Profile '%s' already exists",QQSPV0AV0C_q2V0C->GetName()));
                    delete QQSPV0AV0C_q2V0C;
                    return;
                }
                QQSPV0AV0C_q2V0C->Sumw2();
                SPFlowEseList->Add(QQSPV0AV0C_q2V0C);

                QQSPTPCV0C_q2V0C = new TProfile(Form("<QQ*>_%i{SP_V0C_TPC}_PerCode%i",iVN+2,iEse+1),Form("<QQ*>_%i{SP_V0C_TPC}",iVN+2),nCentBin,CentEdges);
                if(!QQSPTPCV0C_q2V0C) { AliError("<QQ*> V0C TPC profile not created"); return; }
                if(SPFlowEseList->FindObject(QQSPTPCV0C_q2V0C->GetName())) {
                    AliError(Form("Profile '%s' already exists",QQSPTPCV0C_q2V0C->GetName()));
                    delete QQSPTPCV0C_q2V0C;
                    return;
                }
                QQSPTPCV0C_q2V0C->Sumw2();
                SPFlowEseList->Add(QQSPTPCV0C_q2V0C);
            }

        }
    }
    

    if(fReadMC)
    {
        fHistPDG = new TH1F("fHistPDG","",500,0,500);
        fOutputList->Add(fHistPDG);
    }

    
    //load q-selection - testing connected splines via addtask macro currently
    /*if(!gGrid) { TGrid::Connect("alien://"); }

    if(fTPCEse){
        TFile* fFileSpq2TPC = TFile::Open("alien:///alice/cern.ch/user/j/joachimh/q2TPCSpRun15oPbPb.root");
        if(!fFileSpq2TPC) { printf("q_2 TPC Spline file cannot be opened \n"); return; }
        TFile* fFileSpq3TPC = TFile::Open("alien:///alice/cern.ch/user/j/joachimh/q3TPCSpRun15oPbPb.root");
        if(!fFileSpq3TPC) { printf("q_3 TPC Spline file cannot be opened \n"); return; }

        for (Int_t iSp(0); iSp<90; ++iSp){
            fSplq2TPC[iSp] = (TSpline3*)fFileSpq2TPC->Get(Form("sp_q2TPC_%i",iSp));
            fSplq3TPC[iSp] = (TSpline3*)fFileSpq3TPC->Get(Form("sp_q3TPC_%i",iSp));
        }
    }
    if(fV0CEse){
        TFile* fFileSpq2V0C = TFile::Open("alien:///alice/cern.ch/user/j/joachimh/q2V0CSpRun15oPbPb.root");
        if(!fFileSpq2V0C) { printf("q_2 V0C Spline file cannot be opened \n"); return; }
        TFile* fFileSpq3V0C = TFile::Open("alien:///alice/cern.ch/user/j/joachimh/q3V0CSpRun15oPbPb.root");
        if(!fFileSpq3V0C) { printf("q_3 V0C Spline file cannot be opened \n"); return; }

        for (Int_t iSp(0); iSp<90; ++iSp){
            fSplq2V0C[iSp] = (TSpline3*)fFileSpq2V0C->Get(Form("sp_q2V0C_%i",iSp));
            fSplq3V0C[iSp] = (TSpline3*)fFileSpq3V0C->Get(Form("sp_q3V0C_%i",iSp));
        }
    }

    if(fV0AEse){
        TFile* fFileSpq2V0A = TFile::Open("alien:///alice/cern.ch/user/j/joachimh/q2V0ASpRun15oPbPb.root");
        if(!fFileSpq2V0A) { printf("q_2 V0A Spline file cannot be opened \n"); return; }
        TFile* fFileSpq3V0A = TFile::Open("alien:///alice/cern.ch/user/j/joachimh/q3V0ASpRun15oPbPb.root");
        if(!fFileSpq3V0A) { printf("q_3 V0A Spline file cannot be opened \n"); return; }

        for (Int_t iSp(0); iSp<90; ++iSp){
            fSplq2V0A[iSp] = (TSpline3*)fFileSpq2V0A->Get(Form("sp_q2V0A_%i",iSp));
            fSplq3V0A[iSp] = (TSpline3*)fFileSpq3V0A->Get(Form("sp_q3V0A_%i",iSp));
        }
    }*/
    // Load NUE
    if(fUseNUEWeights){
        if(!gGrid) { TGrid::Connect("alien://"); }
        if(!fIs2018Data){
          if(fNUE==1){
              fFileTrackEff = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_AMPT.root");
          }
          if(fNUE==2){
              fFileTrackEff = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR.root");
          }

          if(!fFileTrackEff) { printf("Problem loading NUE file \n"); }
        }
        
    }

    if(fFillQARej){
    fhQAEventsfMult32vsCentr = new TH2D("fhQAEventsfMult32vsCentr", "; centr: V0M; TPC Mult (FB32)", 100, 0, 100, 100, 0, 3000);
    fObservables->Add(fhQAEventsfMult32vsCentr);
    fhQAEventsfMult96vsCentr = new TH2D("fhQAEventsfMult96vsCentr", "; centr: V0M; TPC Mult (FB96)", 100,0,100,100,0,4000);
    fObservables->Add(fhQAEventsfMult96vsCentr);
    fhQAEventsfMult128vsCentr = new TH2D("fhQAEventsfMult128vsCentr", "; centr: V0M; TPC Mult (FB128)", 100, 0, 100, 100, 0, 5000);
    fObservables->Add(fhQAEventsfMult128vsCentr);
    fhQAEventsfMultTPCvsTOF = new TH2D("fhQAEventsfMultTPCvsTOF", "; TPC FB32 Mult; TOF Mult", 200, 0, 4000, 200, 0, 2000);
    fObservables->Add(fhQAEventsfMultTPCvsTOF);
    fhQAEventsfMultTPCvsESD = new TH2D("fhQAEventsfMultTPCvsESD", "; TPC FB128 Mult; TOF Mult)", 200, 0, 7000, 300, -1000, 35000);
    fObservables->Add(fhQAEventsfMultTPCvsESD);
    }

    TString sEventCounterLabel[] = {"Input"};
    const Int_t iEventCounterBins = sizeof(sEventCounterLabel)/sizeof(sEventCounterLabel[0]);
    fhEventCounter = new TH1D("fhEventCounter","Event Counter", iEventCounterBins,0,iEventCounterBins);
    for(Int_t i(0); i < iEventCounterBins; ++i) { fhEventCounter->GetXaxis()->SetBinLabel(i+1,sEventCounterLabel[i].Data() ); }
    fObservables->Add(fhEventCounter);

    if(fBayesUnfoldingInput){
      TH2* v2Raw = nullptr;

      TH2* v2XRaw = nullptr;
      TH2* v2YRaw = nullptr;

      TH2* v2RawA = nullptr;
      TH2* v2RawB = nullptr;

      TH2* v2Raw_q2V0C = nullptr;

      v2Raw = new TH2D("v2ObsRaw","v2ObsRaw",100,0,100,500,0,0.3);

      if(!v2Raw) { AliError("Centrality 2Dhist not created"); return; }
      if(fOutputList->FindObject(v2Raw->GetName())) {
          AliError(Form("Hist '%s' already exists",v2Raw->GetName()));
          delete v2Raw;
          return;
      }

      v2Raw->Sumw2();
      fOutputList->Add(v2Raw);

      v2RawA = new TH2D("v2ObsRawA","v2ObsRawA",100,0,100,500,0,0.3);

      if(!v2RawA) { AliError("Centrality 2Dhist not created"); return; }
      if(fOutputList->FindObject(v2RawA->GetName())) {
          AliError(Form("Hist '%s' already exists",v2RawA->GetName()));
          delete v2RawA;
          return;
      }

      v2RawA->Sumw2();
      fOutputList->Add(v2RawA);

      v2RawB = new TH2D("v2ObsRawB","v2ObsRawB",100,0,100,500,0,0.3);

      if(!v2RawB) { AliError("Centrality 2Dhist not created"); return; }
      if(fOutputList->FindObject(v2RawB->GetName())) {
          AliError(Form("Hist '%s' already exists",v2RawB->GetName()));
          delete v2RawB;
          return;
      }

      v2RawB->Sumw2();
      fOutputList->Add(v2RawB);

      v2XRaw = new TH2D("v2XObsRaw","v2XObsRaw",100,0,100,100,-0.5,0.5);

      if(!v2XRaw) { AliError("Centrality 2Dhist not created"); return; }
      if(fOutputList->FindObject(v2XRaw->GetName())) {
          AliError(Form("Hist '%s' already exists",v2XRaw->GetName()));
          delete v2XRaw;
          return;
      }

      v2XRaw->Sumw2();
      fOutputList->Add(v2XRaw);

      v2YRaw = new TH2D("v2YObsRaw","v2YObsRaw",100,0,100,100,-0.5,0.5);

      if(!v2YRaw) { AliError("Centrality 2Dhist not created"); return; }
      if(fOutputList->FindObject(v2YRaw->GetName())) {
          AliError(Form("Hist '%s' already exists",v2YRaw->GetName()));
          delete v2YRaw;
          return;
      }

      v2YRaw->Sumw2();
      fOutputList->Add(v2YRaw);

      for (Int_t iEse(0);iEse<ESEPercAxis->GetNbins();++iEse){
        v2Raw_q2V0C = new TH2D(Form("v2ObsRaw_q2V0C_Perc%i",iEse+1),Form("v2ObsRaw_q2V0C_Perc%i",iEse+1),100,0,100,500,0,0.3);

        if(!v2Raw_q2V0C) { AliError("Centrality 2Dhist not created"); return; }
        if(fOutputList->FindObject(v2Raw_q2V0C->GetName())) {
            AliError(Form("Hist '%s' already exists",v2Raw_q2V0C->GetName()));
            delete v2Raw_q2V0C;
            return;
        }

        v2Raw_q2V0C->Sumw2();
        fOutputList->Add(v2Raw_q2V0C);
      }
    }
    

    fEventCuts.AddQAplotsToList(fQAEvents); //QA plots
    
    fObservables->Add(fHistPhiEtaVz);
    fObservables->Add(fHistPhi);
    fObservables->Add(fHistEta);
    fObservables->Add(fHistPt);
    fObservables->Add(fHistZVertex);
    fObservables->Add(fHistPhiCor);
    fObservables->Add(fHistPhiCor3D);
    fObservables->Add(fHistTPCchi2);
    fObservables->Add(fHistITSchi2);
    fObservables->Add(fProfNPar);
    fObservables->Add(fhV0Multiplicity);
    fObservables->Add(fHistPtEtaVz);

    fObservables->Add(fhV0CorrMult);

    fqnDist->Add(fq2TPC);
    fqnDist->Add(fq3TPC);
    fqnDist->Add(fq2V0C);
    fqnDist->Add(fq3V0C);
    fqnDist->Add(fq2V0A);
    fqnDist->Add(fq3V0A);



    for (Int_t qi(0);qi<2;++qi){
        fqnDist->Add(fQnxV0C[qi]);
        fqnDist->Add(fQnyV0C[qi]);
        fqnDist->Add(fQnxV0A[qi]);
        fqnDist->Add(fQnyV0A[qi]);
        fqnDist->Add(fQnxTPC[qi]);
        fqnDist->Add(fQnyTPC[qi]);

        fqnDist->Add(fQnxV0CEse[qi]);
        fqnDist->Add(fQnyV0CEse[qi]);
        fqnDist->Add(fQnxV0AEse[qi]);
        fqnDist->Add(fQnyV0AEse[qi]);
        fqnDist->Add(fQnxTPCEse[qi]);
        fqnDist->Add(fQnyTPCEse[qi]);

        fqnDist->Add(fQnxV0CCor[qi]);
        fqnDist->Add(fQnyV0CCor[qi]);
        fqnDist->Add(fQnxV0ACor[qi]);
        fqnDist->Add(fQnyV0ACor[qi]);
    

        fqnDist->Add(fhqnTPCvqnV0C[qi]);
        fqnDist->Add(fhqnV0CvqnV0A[qi]);
        fqnDist->Add(fhqnTPCvqnV0A[qi]);
    }



    PostData(1, fOutputList);
    PostData(2, fObservables);
    PostData(3, fCorrDist);
    PostData(4, fpTDiff);
    PostData(5, fqnDist);
    PostData(6, fpTDiffESETPC);
    PostData(7, fcnESETPC);
    PostData(8, fpTDiffESEV0C);
    PostData(9, fcnESEV0C);
    PostData(10, fpTDiffESEV0A);
    PostData(11, fcnESEV0A);
    PostData(12, SPFlowList);
    PostData(13, SPFlowEseList);
    PostData(14, fQAEvents);
}
//_____________________________________________________________________________
void AliAnalysisTaskESEFlow::UserExec(Option_t *)
{

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) { return; }

  fhEventCounter->Fill("Input",1);

  if(!IsEventSelected()) { return; }
  
  fhEventCounter->Fill("Event OK",1);

  //centrality
  Float_t centrality(0);
  AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
  if(multSelection) centrality = multSelection->GetMultiplicityPercentile(fCentEstimator);

  if(!fMakeRBRweightsRun){
      if(fFlowRunByRunWeights){
          if(!LoadWeights()) { AliFatal("\n \n \n \n \n \n \n \n \n \n Weights not loaded, terminating! \n \n \n \n \n \n \n \n \n \n "); return; }
      }
  }

  if(fV0RunByRunCalibration){
      if(!LoadV0Calibration()) { AliFatal("\n \n \n \n \n \n \n \n \n \n V0 Calibration not loaded, terminating! \n \n \n \n \n \n \n \n \n \n "); return; }
  }

  if(fUseNUEWeights){
      if(!LoadNUE()) { AliFatal("\n \n \n NUE Weights not loaded properly, terminating! \n \n \n"); return; }
  }
  

  fhEventCounter->Fill("Selected OK",1);

  FillObsDistributions(centrality);

  Int_t fSpCent = (Int_t)centrality;

  CorrelationTask(centrality, fSpCent);

  if(fReadMC) 
  { 
      fMCEvent = MCEvent();
      if(fMCEvent) ProcessMCParticles();
  }

  
  PostData(1, fOutputList);
  PostData(2, fObservables);
  PostData(3, fCorrDist);
  PostData(4, fpTDiff);
  PostData(5, fqnDist);
  PostData(6, fpTDiffESETPC);
  PostData(7, fcnESETPC);
  PostData(8, fpTDiffESEV0C);
  PostData(9, fcnESEV0C);
  PostData(10, fpTDiffESEV0A);
  PostData(11, fcnESEV0A);
  PostData(12, SPFlowList);
  PostData(13, SPFlowEseList);
  PostData(14, fQAEvents);
}
//_____________________________________________________________________________
void AliAnalysisTaskESEFlow::Terminate(Option_t *)
{
}
void AliAnalysisTaskESEFlow::CorrelationTask(const Float_t centrality, Int_t fSpCent)
{
    Int_t CentrCode = GetCentralityCode(centrality);

    //if( (CentrCode < 0) || (CentrCode > 9)) { return; }
    if ( (centrality < 0) || (centrality > 90)) { return; }

    ReducedqVectorsTPC(centrality, fSpCent);

    //if( (CentrCode < 0) || (CentrCode > 8)) { return; }
    if(fV0RunByRunCalibration){
      ReducedqVectorsV0(centrality, fSpCent);
      
      for (Int_t i1(0);i1<2;++i1){
        fhqnTPCvqnV0C[i1]->Fill(qnTPC[i1],qnV0C[i1]);
        fhqnV0CvqnV0A[i1]->Fill(qnV0C[i1],qnV0A[i1]);
        fhqnTPCvqnV0A[i1]->Fill(qnTPC[i1],qnV0A[i1]);

        TProfile* prof_avg_TPC = (TProfile*)fOutputList->FindObject(Form("fPq_%i_TPC",i1+2));
        if(!prof_avg_TPC) { AliError("avg qn TPC not found"); return; }
        prof_avg_TPC->Fill(centrality, qnTPC[i1]);

        TProfile* prof_avg_V0C = (TProfile*)fOutputList->FindObject(Form("fPq_%i_V0C",i1+2));
        if(!prof_avg_V0C) { AliError("avg qn V0C not found"); return; }
        prof_avg_V0C->Fill(centrality, qnV0C[i1]);

        TProfile* prof_avg_V0A = (TProfile*)fOutputList->FindObject(Form("fPq_%i_V0A",i1+2));
        if(!prof_avg_V0A) { AliError("avg qn V0A not found"); return; }
        prof_avg_V0A->Fill(centrality, qnV0A[i1]);
      }
    }

    if( (CentrCode < 0) || (CentrCode > 8)) { return; }

    fIndexSampling = GetSamplingIndex();


    if((!fMakeqSelectionRun && fV0RunByRunCalibration)){

        
        Double_t q2TPCInp = 100.*fSplq2TPC[fSpCent]->Eval(qnTPC[0]);
        Double_t q3TPCInp = 100.*fSplq3TPC[fSpCent]->Eval(qnTPC[1]);

        Int_t q2ESECodeTPC = GetEsePercentileCode(q2TPCInp);
        if (q2ESECodeTPC<0) { printf("Problem with q_2 TPC percentile: negative percentile \n"); return; } 
        Int_t q3ESECodeTPC = GetEsePercentileCode(q3TPCInp);
        if (q2ESECodeTPC<0) { printf("Problem with q_3 TPC percentile: negative percentile \n"); return; } 
        
        Double_t q2V0CInp = 100.*fSplq2V0C[fSpCent]->Eval(qnV0C[0]);
        Double_t q3V0CInp = 100.*fSplq3V0C[fSpCent]->Eval(qnV0C[1]);

        Int_t q2ESECodeV0C = GetEsePercentileCode(q2V0CInp);
        if (q2ESECodeV0C<0) { printf("Problem with q_2 V0C percentile: negative percentile \n"); return; } 
        Int_t q3ESECodeV0C = GetEsePercentileCode(q3V0CInp);
        if (q3ESECodeV0C<0) { printf("Problem with q_3 V0C percentile: negative percentile \n"); return; } 
        


        Double_t q2V0AInp = 100.*fSplq2V0A[fSpCent]->Eval(qnV0A[0]); // do q-selection for V0A
        Double_t q3V0AInp = 100.*fSplq3V0A[fSpCent]->Eval(qnV0A[1]);

        Int_t q2ESECodeV0A = GetEsePercentileCode(q2V0AInp);
        if (q2ESECodeV0A<0) { printf("Problem with q_2 V0A percentile: negative percentile \n"); return; } 
        Int_t q3ESECodeV0A = GetEsePercentileCode(q3V0AInp);
        if (q3ESECodeV0A<0) { printf("Problem with q_3 V0A percentile: negative percentile \n"); return; } 
        

        RFPVectors(centrality);

        
        FillCorrelation(centrality, -1, q2ESECodeTPC, q3ESECodeTPC, q2ESECodeV0C, q3ESECodeV0C, q2ESECodeV0A, q3ESECodeV0A, 1, 0);

        POIVectors(centrality, q2ESECodeTPC, q3ESECodeTPC, q2ESECodeV0C, q3ESECodeV0C, q2ESECodeV0A, q3ESECodeV0A);

        // SCALAR PROD ANALYZER
        if(fSPAnalysis){
          SPVienna(centrality, q2ESECodeV0C);
        }

        if(fActq2Projections){
          TH2D* histq2TPC_perc = (TH2D*)fqnDist->FindObject(Form("fq2vCentTPC_Perc%i",q2ESECodeTPC+1));
          if(!histq2TPC_perc) { AliError(Form("fq2vCentTPC_Perc%i not found",q2ESECodeTPC+1)); return; }
          histq2TPC_perc->Fill(centrality, qnTPC[0]);

          TH2D* histq2V0C_perc = (TH2D*)fqnDist->FindObject(Form("fq2vCentV0C_Perc%i",q2ESECodeV0C+1));
          if(!histq2V0C_perc) { AliError(Form("fq2vCentV0C_Perc%i not found",q2ESECodeV0C+1)); return; }
          histq2V0C_perc->Fill(centrality, qnV0C[0]);

          TH2D* histq2V0A_perc = (TH2D*)fqnDist->FindObject(Form("fq2vCentV0A_Perc%i",q2ESECodeV0A+1));
          if(!histq2V0A_perc) { AliError(Form("fq2vCentV0A_Perc%i not found",q2ESECodeV0A+1)); return; }
          histq2V0A_perc->Fill(centrality, qnV0A[0]);
        }


        if(fBayesUnfoldingInput){
          /// Input for Unfolding
          BayesianUnfolding(centrality, q2ESECodeV0C);
        }
        
    }

    return;
}
void AliAnalysisTaskESEFlow::BayesianUnfolding(Float_t centrality, Int_t q2ESECodeV0C)
{
  Int_t iTracks(fAOD->GetNumberOfTracks());
  //Double_t dVz = fAOD->GetPrimaryVertex()->GetZ();

  if(iTracks < 1 ) { return; }

  Double_t V2xA = 0.0;
  Double_t V2yA = 0.0;
  Double_t V2xB = 0.0;
  Double_t V2yB = 0.0;

  Double_t MA = 0.0;
  Double_t MB = 0.0;

  for(Int_t i(0); i < iTracks; ++i) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track || !IsTrackSelected(track)) { continue; }

    double dEta = track->Eta();
    double dPhi = track->Phi();
    //double dPt = track->Pt();
    
    Double_t dWeight = 1.0;
    if((0.) < dEta && dEta < (0.8)){
      MA += dWeight;
      Double_t dCos = dWeight * TMath::Cos( 2 * dPhi);
      Double_t dSin = dWeight * TMath::Sin( 2 * dPhi);
      V2xA += dCos;
      V2yA += dSin;
    }

    if(- (0.8) < dEta && dEta < (0.)){
      MB += dWeight;
      Double_t dCos = dWeight * TMath::Cos( 2 * dPhi);
      Double_t dSin = dWeight * TMath::Sin( 2 * dPhi);
      V2xB += dCos;
      V2yB += dSin;
    }

  } // ending track loop

  Double_t Ev2xA = V2xA;
  Double_t Ev2yA = V2yA;

  Double_t Ev2xB = V2xB;
  Double_t Ev2yB = V2yB;

  Double_t dq2A = Ev2xA*Ev2xA + Ev2yA*Ev2yA;
  Double_t dq2B = Ev2xB*Ev2xB + Ev2yB*Ev2yB;


  if(MA>0 && MB>0){
    Double_t v2RawIn = TMath::Sqrt(dq2A*dq2B)/(MA*MB);

    Double_t v22PC = TMath::Sqrt(v2RawIn);

    TH2D* v2raw = (TH2D*)fOutputList->FindObject("v2ObsRaw");
    if(!v2raw) { AliError("v2ObsRaw 2Dhist not found"); return; }
    v2raw->Fill(centrality, v22PC);

    TH2D* v2Rawq2V0C = (TH2D*)fOutputList->FindObject(Form("v2ObsRaw_q2V0C_Perc%i",q2ESECodeV0C+1));
    if(!v2Rawq2V0C) {AliError(Form("v2ObsRaw_q2V0C_Perc%i not found in list",q2ESECodeV0C+1)); return;}
    v2Rawq2V0C->Fill(centrality, v22PC);
  }
  
  if(MA>0){
    Double_t v2A = TMath::Sqrt(dq2A)/MA;

    TH2D* v2rawA = (TH2D*)fOutputList->FindObject("v2ObsRawA");
    if(!v2rawA) { AliError("v2ObsRawA 2Dhist not found"); return; }
    v2rawA->Fill(centrality, v2A);
  }


  if(MB>0){
    Double_t v2B = TMath::Sqrt(dq2B)/MB;

    TH2D* v2rawB = (TH2D*)fOutputList->FindObject("v2ObsRawB");
    if(!v2rawB) { AliError("v2ObsRawB 2Dhist not found"); return; }
    v2rawB->Fill(centrality, v2B);
  }

  TH2D* v2Xraw = (TH2D*)fOutputList->FindObject("v2XObsRaw");
  if(!v2Xraw) { AliError("v2XObsRaw 2Dhist not found"); return; }
  v2Xraw->Fill(centrality, Ev2xA);

  TH2D* v2Yraw = (TH2D*)fOutputList->FindObject("v2YObsRaw");
  if(!v2Yraw) { AliError("v2YObsRaw 2Dhist not found"); return; }
  v2Yraw->Fill(centrality, Ev2yA);

}
void AliAnalysisTaskESEFlow::FillObsDistributions(const Float_t centrality)
{

    Int_t iTracks(fAOD->GetNumberOfTracks());
    Double_t dVz = fAOD->GetPrimaryVertex()->GetZ();

    if(iTracks < 1 ) { return; }
    fHistZVertex->Fill(dVz);
    fProfNPar->Fill(centrality,iTracks);

    

    for(Int_t i(0); i < iTracks; ++i) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }

        double dEta = track->Eta();
        double dPhi = track->Phi();
        double dPt = track->Pt();
        
        fHistPhiEtaVz->Fill(dPhi, dEta, dVz);
        fHistPtEtaVz->Fill(dPt, dEta, dVz);
        fHistPhi->Fill(dPhi);
        fHistEta->Fill(dEta);
        fHistPt->Fill(dPt);

        Float_t chi2tpc = track->GetTPCchi2();
        Float_t ntpccls = track->GetTPCNcls();

        Float_t chi2its = track->GetITSchi2();
        Float_t nitscls = track->GetITSNcls();

        Float_t tpcchi2percls = (ntpccls==0)?0.0:chi2tpc/ntpccls;
        Float_t itschi2percls = (nitscls==0)?0.0:chi2its/nitscls;
        fHistTPCchi2->Fill(tpcchi2percls);
        fHistITSchi2->Fill(itschi2percls);
        

        Double_t dWeight = GetFlowWeight(track, dVz);
        fHistPhiCor->Fill(dPhi,dWeight);
        fHistPhiCor3D->Fill(dPhi, dEta, dVz, dWeight);
        //fHistPtCor3D->Fill(dPt,dEta,dVz,dPtWeight);

    } // ending track loop

    AliAODVZERO* aodV0 = fAOD->GetVZEROData();

    for (Int_t iV0=0; iV0 < 64; ++iV0){
        Float_t multV0 = aodV0->GetMultiplicity(iV0);
        fhV0Multiplicity->Fill(iV0,multV0);

        if(fV0RunByRunCalibration){
            Double_t multCor = -10;

            //V0 calibration applied to each channel
            if (iV0 < 8)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(1);
            else if (iV0 >= 8 && iV0 < 16)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(9);
            else if (iV0 >= 16 && iV0 < 24)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(17);
            else if (iV0 >= 24 && iV0 < 32)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(25);
            else if (iV0 >= 32 && iV0 < 40)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(33);
            else if (iV0 >= 40 && iV0 < 48)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(41);
            else if (iV0 >= 48 && iV0 < 56)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(49);
            else if (iV0 >= 56 && iV0 < 64)
                multCor = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(57);            

            if (multCor < 0){
                printf("Problem with multiplicity in V0 in observation func");
                continue;
            }

            fhV0CorrMult->Fill(iV0,multCor);
        }

    } // ending V0 calibration loop

    if(fFillQARej){
        QAMultFiller(centrality);
    }

    return;
}
void AliAnalysisTaskESEFlow::RFPVectors(const Float_t centrality)
{
    ResetFlowVector(Qvector);
    ResetFlowVector(Qvector10P);
    ResetFlowVector(Qvector10M);
    Double_t dEtaLimit = 0.5*dEtaGap;

    Int_t iTracks(fAOD->GetNumberOfTracks());
    Double_t dVz = fAOD->GetPrimaryVertex()->GetZ();

    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; ++i) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        if(!WithinRFP(track)){continue;} // check if also within reference particle

        double dEta = track->Eta();
        double dPhi = track->Phi();
        double dPt = track->Pt();

        Double_t dWeight = GetFlowWeight(track,dVz);

        Double_t dPtWeight = 1.0;
        if(fUseNUEWeights){
            if(dPt < 3.0){
                dPtWeight = GetNUEPtWeight(dPt, dEta, dVz);
            }
            if(dPt > 3.0){
                dPtWeight = GetNUEPtWeight(2.5, dEta, dVz);
            }
        }
        
        // no eta gap
        for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
        {
            for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
            {
                Double_t dCos = TMath::Power(dWeight*dPtWeight,iPower) * TMath::Cos(iHarm * dPhi);
                Double_t dSin = TMath::Power(dWeight*dPtWeight,iPower) * TMath::Sin(iHarm * dPhi);
                Qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            }
        }
    
        // eta gap
        if(dEta > dEtaLimit)
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
            {
                for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                {
                    Double_t dCos = TMath::Power(dWeight*dPtWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight*dPtWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    Qvector10P[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                }
            }
        } 
        // RFP in negative eta acceptance
        else if(dEta < - dEtaLimit)
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
            { 
                for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                {
                    Double_t dCos = TMath::Power(dWeight*dPtWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight*dPtWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    Qvector10M[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                }
            }
        }            
    } // ending track loop

    return;
}
void AliAnalysisTaskESEFlow::FillPOI(const Double_t dPtLow, const Double_t dPtHigh)
{

    ResetFlowVector(qvector);
    ResetFlowVector(pvector);
    ResetFlowVector(pvector10P);
    ResetFlowVector(pvector10M);

    Double_t dEtaLimit = 0.5*dEtaGap;

    Int_t iPart(fAOD->GetNumberOfTracks());
    Double_t dVz = fAOD->GetPrimaryVertex()->GetZ();

    if(iPart < 1 ) { return; }
    for(Int_t i(0); i < iPart; ++i)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }

        double dEta = track->Eta();
        double dPhi = track->Phi();
        double dPt = track->Pt();

        Bool_t bIsWithinPOI = WithinPOI(track);

        Bool_t bIsWithinRFP = WithinRFP(track);

        if(!bIsWithinPOI) {continue;}

        Double_t dWeight = GetFlowWeight(track,dVz);
        
        if(dPt > dPtLow && dPt <= dPtHigh)
        {
            // q vector
            if(bIsWithinRFP) // check if also within reference particle
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
                { 
                    for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                    {   
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                }
            }

            // NO eta gap for p-vector
            
            for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
            {
                for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    pvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                }
            }

            // Eta gap for p-vector
            if(dEta > dEtaLimit)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
                {
                    for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                    {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    pvector10P[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                }
            } 
            if(dEta < -dEtaLimit)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
                { 
                    for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                    {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    pvector10M[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                }
            }     
        }
    } // ending track loop

    return;
}
void AliAnalysisTaskESEFlow::POIVectors(const Float_t centrality, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C, Int_t q2ESECodeV0A, Int_t q3ESECodeV0A)
{
    // DIFFERENTIAL FLOW
    // loop over p_T histogram
    Int_t iNumPtBins = fHistPt->GetXaxis()->GetNbins();

    for(Int_t iPt(1); iPt< iNumPtBins+1; ++iPt)
    {
        // reset flow vector for each pt bin

        Double_t dPtLow = fHistPt->GetXaxis()->GetBinLowEdge(iPt);
        Double_t dPtHigh = fHistPt->GetXaxis()->GetBinUpEdge(iPt);
        Double_t dPtL = fHistPt->GetXaxis()->GetBinCenter(iPt);

        FillPOI(dPtLow, dPtHigh);

        FillCorrelation(centrality, dPtL, q2ESECodeTPC, q3ESECodeTPC, q2ESECodeV0C, q3ESECodeV0C, q2ESECodeV0A, q3ESECodeV0A, 0, 1);

    }

    return;
}
void AliAnalysisTaskESEFlow::ReducedqVectorsTPC(const Float_t centrality, const Int_t SPCode)
{
    ResetReducedqVector(qnTPC);
    ResetReducedqVector(QxnTPC);
    ResetReducedqVector(QynTPC);
    ResetReducedqVector(QxnTPCEse);
    ResetReducedqVector(QynTPCEse);

    ResetReducedqVector(qnTPC_SP);
    ResetReducedqVector(QxnTPCSP);
    ResetReducedqVector(QynTPCSP);

    if(!bHasGap) { AliWarning("Warning, no eta gap. Overlap between reduced q-vectors and flow vectors! \n");}

    Int_t iTracks(fAOD->GetNumberOfTracks());
    Double_t dVz = fAOD->GetPrimaryVertex()->GetZ();

    Double_t M = 0;
    Double_t M_SP = 0;

    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; ++i) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        if(!WithinRedRP(track)){continue;} // 

        double dEta = track->Eta();
        double dPhi = track->Phi();
        //double dPt = track->Pt();

        Double_t dWeight = GetFlowWeight(track,dVz);

        if(- (0.4) < dEta && dEta < (0.4))
        {
            M += dWeight;
            for(Int_t iHarm(0); iHarm < 2; ++iHarm)
            {                
                Double_t dCos = dWeight * TMath::Cos( (iHarm+2) * dPhi);
                Double_t dSin = dWeight * TMath::Sin( (iHarm+2) * dPhi);
                QxnTPC[iHarm] += dCos;
                QynTPC[iHarm] += dSin;
            }
        }

        if( TMath::Abs(dEta) < fAbsEtaMax){
            M_SP += 1;
            for(Int_t iHarm(0); iHarm < 2; ++iHarm)
            {                
                Double_t dCos = 1 * TMath::Cos( (iHarm+2) * dPhi);
                Double_t dSin = 1 * TMath::Sin( (iHarm+2) * dPhi);
                QxnTPCSP[iHarm] += dCos;
                QynTPCSP[iHarm] += dSin;
            }
        }
    }


    if(fV0RunByRunCalibration){
    for (Int_t iQn(0); iQn < 2; ++iQn){
        QxnTPCEse[iQn] = QxnTPC[iQn] - fQnxTPCm[iQn]->GetBinContent(SPCode+1);
        QynTPCEse[iQn] = QynTPC[iQn] - fQnyTPCm[iQn]->GetBinContent(SPCode+1);
    }
    }

    for (Int_t nQ(0); nQ<2;++nQ){
        fQnxTPC[nQ]->Fill(centrality,QxnTPC[nQ]);
        fQnyTPC[nQ]->Fill(centrality,QynTPC[nQ]);


        if(fV0RunByRunCalibration){
        fQnxTPCEse[nQ]->Fill(centrality,QxnTPCEse[nQ]);
        fQnyTPCEse[nQ]->Fill(centrality,QynTPCEse[nQ]);
        }
    } //used for recentering of Qnx/Qny, using normal for ESE

    if(M>0)
    {
        for(Int_t iHarm(0); iHarm < 2; ++iHarm)
        { 
            Double_t dqnEse = QxnTPC[iHarm]*QxnTPC[iHarm]+QynTPC[iHarm]*QynTPC[iHarm];
            qnTPC[iHarm] = TMath::Sqrt(dqnEse)/TMath::Sqrt(M);
        }
        FillqnRedTPC(centrality);
    }
    if(M_SP>0){
        for(Int_t iHarm(0); iHarm < 2; ++iHarm)
        { 
            Double_t dqn_SP = QxnTPCSP[iHarm]*QxnTPCSP[iHarm]+QynTPCSP[iHarm]*QynTPCSP[iHarm];
            qnTPC_SP[iHarm] = TMath::Sqrt(dqn_SP)/TMath::Sqrt(M_SP);
        }
    }

    return;
}
void AliAnalysisTaskESEFlow::ReducedqVectorsV0(const Float_t centrality,const Int_t SPCode)
{
    ResetReducedqVector(qnV0C);
    ResetReducedqVector(QxnV0C);
    ResetReducedqVector(QynV0C);
    ResetReducedqVector(qnV0A);
    ResetReducedqVector(QxnV0A);
    ResetReducedqVector(QynV0A);
    ResetReducedqVector(QxnV0CEse);
    ResetReducedqVector(QynV0CEse);
    ResetReducedqVector(QxnV0AEse);
    ResetReducedqVector(QynV0AEse);

    ResetReducedqVector(QxnV0CCorr);
    ResetReducedqVector(QynV0CCorr);
    ResetReducedqVector(QxnV0ACorr);
    ResetReducedqVector(QynV0ACorr);

    AliAODVZERO* aodV0 = fAOD->GetVZEROData();

    Double_t MC=0;

    Double_t MA=0;

    for (Int_t iV0=0; iV0 < 64; ++iV0){

        Double_t PhiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);

        Float_t multV0 = aodV0->GetMultiplicity(iV0);

        if(iV0<32){

            Double_t multCorC = -10;

            //V0 calibration
            if (iV0 < 8)
                multCorC = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(1);
            else if (iV0 >= 8 && iV0 < 16)
                multCorC = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(9);
            else if (iV0 >= 16 && iV0 < 24)
                multCorC = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(17);
            else if (iV0 >= 24 && iV0 < 32)
                multCorC = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(25);
            

            if (multCorC < 0){
                printf("Problem with multiplicity in V0C");
                continue;
            }

            for(Int_t iHarm(0); iHarm < 2; ++iHarm)
            {                
                Double_t dCosC = multCorC * TMath::Cos( (iHarm+2) * PhiV0);
                Double_t dSinC = multCorC * TMath::Sin( (iHarm+2) * PhiV0);
                QxnV0C[iHarm] += dCosC;
                QynV0C[iHarm] += dSinC;
            }

            MC += multCorC;
        }
        else {
            Double_t multCorA = -10;

            if (iV0 >= 32 && iV0 < 40)
                multCorA = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(33);
            else if (iV0 >= 40 && iV0 < 48)
                multCorA = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(41);
            else if (iV0 >= 48 && iV0 < 56)
                multCorA = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(49);
            else if (iV0 >= 56 && iV0 < 64)
                multCorA = multV0/fhV0Calib->GetBinContent(iV0+1)*fhV0Calib->GetBinContent(57);            

            if (multCorA < 0){
                printf("Problem with multiplicity in V0A");
                continue;
            }

            for(Int_t iHarm(0); iHarm < 2; ++iHarm)
            {                
                Double_t dCosA = multCorA * TMath::Cos( (iHarm+2) * PhiV0);
                Double_t dSinA = multCorA * TMath::Sin( (iHarm+2) * PhiV0);
                QxnV0A[iHarm] += dCosA;
                QynV0A[iHarm] += dSinA;
            }

            MA += multCorA;

        }
    }


    for (Int_t iQn(0); iQn < 2; ++iQn){
        QxnV0CEse[iQn] = QxnV0C[iQn] - fQnxV0Cm[iQn]->GetBinContent(SPCode+1);
        QynV0CEse[iQn] = QynV0C[iQn] - fQnyV0Cm[iQn]->GetBinContent(SPCode+1);

        QxnV0AEse[iQn] = QxnV0A[iQn] - fQnxV0Am[iQn]->GetBinContent(SPCode+1);
        QynV0AEse[iQn] = QynV0A[iQn] - fQnyV0Am[iQn]->GetBinContent(SPCode+1);

        // recentering
        if (SPCode >88) { continue; }
        QxnV0CCorr[iQn] = (QxnV0C[iQn] - fQnxV0Cm[iQn]->GetBinContent(SPCode+1))/fQnxV0Cs[iQn]->GetBinContent(SPCode+1);
        QynV0CCorr[iQn] = (QynV0C[iQn] - fQnyV0Cm[iQn]->GetBinContent(SPCode+1))/fQnyV0Cs[iQn]->GetBinContent(SPCode+1);

        QxnV0ACorr[iQn] = (QxnV0A[iQn] - fQnxV0Am[iQn]->GetBinContent(SPCode+1))/fQnxV0As[iQn]->GetBinContent(SPCode+1);
        QynV0ACorr[iQn] = (QynV0A[iQn] - fQnyV0Am[iQn]->GetBinContent(SPCode+1))/fQnyV0As[iQn]->GetBinContent(SPCode+1);
    }


    // QA PLOTS
    for (Int_t nQ(0);nQ<2;++nQ){
        fQnxV0C[nQ]->Fill(centrality,QxnV0C[nQ]); // fill for recentering
        fQnyV0C[nQ]->Fill(centrality,QynV0C[nQ]);
        fQnxV0A[nQ]->Fill(centrality,QxnV0A[nQ]);
        fQnyV0A[nQ]->Fill(centrality,QynV0A[nQ]); // End recenter fill

        fQnxV0CEse[nQ]->Fill(centrality,QxnV0CEse[nQ]); // fill after ese corr
        fQnyV0CEse[nQ]->Fill(centrality,QynV0CEse[nQ]);
        fQnxV0AEse[nQ]->Fill(centrality,QxnV0AEse[nQ]);
        fQnyV0AEse[nQ]->Fill(centrality,QynV0AEse[nQ]); // end after 

        if (SPCode >88) { continue; }
        fQnxV0CCor[nQ]->Fill(centrality,QxnV0CCorr[nQ]);
        fQnyV0CCor[nQ]->Fill(centrality,QynV0CCorr[nQ]);
        fQnxV0ACor[nQ]->Fill(centrality,QxnV0ACorr[nQ]);
        fQnyV0ACor[nQ]->Fill(centrality,QynV0ACorr[nQ]);
    } 

    if(MC>0)
    {
        for(Int_t iHarm(0); iHarm < 2; ++iHarm)
        {
            Double_t dqnV0CEse = QxnV0CEse[iHarm]*QxnV0CEse[iHarm]+QynV0CEse[iHarm]*QynV0CEse[iHarm];
            qnV0C[iHarm] = TMath::Sqrt(dqnV0CEse)/TMath::Sqrt(MC);
        }
        FillqnRedV0(centrality,"V0C");
    }
    if(MA>0)
    {
        for(Int_t iHarm(0); iHarm < 2; ++iHarm)
        {
            Double_t dqnV0AEse = QxnV0AEse[iHarm]*QxnV0AEse[iHarm]+QynV0AEse[iHarm]*QynV0AEse[iHarm];
            qnV0A[iHarm] = TMath::Sqrt(dqnV0AEse)/TMath::Sqrt(MA);
        }
        FillqnRedV0(centrality,"V0A");
    }

    return;
}
void AliAnalysisTaskESEFlow::SPVienna(const Float_t centrality, Int_t q2ESECodeV0C)
{
    Int_t CentrCode = GetCentralityCode(centrality);

    // event plane angle for the second harmonic
    Double_t Psi_V0C2 = TMath::ATan2(QynV0CCorr[0],QxnV0CCorr[0])/2.0;
    Double_t Psi_V0A2 = TMath::ATan2(QynV0ACorr[0],QxnV0ACorr[0])/2.0;

    Double_t Psi_V0C3 = TMath::ATan2(QynV0CCorr[1],QxnV0CCorr[1])/3.0;
    Double_t Psi_V0A3 = TMath::ATan2(QynV0ACorr[1],QxnV0ACorr[1])/3.0;

    fhEvPlPsi_2V0C->Fill(centrality,Psi_V0C2);
    fhEvPlPsi_2V0A->Fill(centrality,Psi_V0A2);

    fhEvPlPsi_3V0C->Fill(centrality,Psi_V0C3);
    fhEvPlPsi_3V0A->Fill(centrality,Psi_V0A3);

    Int_t iTracks(fAOD->GetNumberOfTracks());
    Double_t dVz = fAOD->GetPrimaryVertex()->GetZ();

    //vnSPV0A = 0;
    //vnSPV0C = 0;
    ResetReducedqVector(vnSPV0A);
    ResetReducedqVector(vnSPV0C);

    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; ++i) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        if(!WithinPOI(track)){continue;}
        //if(dPt > dPtLow && dPt <= dPtHigh){

        Double_t dPt = track->Pt();
        Double_t dWeight = GetFlowWeight(track,dVz);

        for (Int_t iVienna(0); iVienna<2; ++iVienna){
            Int_t nHarmv = iVienna+2;
            vnSPV0A[iVienna] =  (TMath::Cos(nHarmv * track->Phi())*QxnV0ACorr[iVienna] + TMath::Sin(nHarmv * track->Phi())*QynV0ACorr[iVienna]); 
            vnSPV0C[iVienna] =  (TMath::Cos(nHarmv * track->Phi())*QxnV0CCorr[iVienna] + TMath::Sin(nHarmv * track->Phi())*QynV0CCorr[iVienna]);

            TProfile* SPvProfV0A = (TProfile*)SPFlowList->FindObject(Form("v_{%i}{SP V0A}_%.0f_%.0f",iVienna+2,CentEdges[CentrCode],CentEdges[CentrCode+1]));
            if(!SPvProfV0A) { AliError(Form("Profile '%s' not found","<v>_2{SP_V0A}")); return; }
            SPvProfV0A->Fill(dPt, vnSPV0A[iVienna]);

            TProfile* SPvProfV0C = (TProfile*)SPFlowList->FindObject(Form("v_{%i}{SP V0C}_%.0f_%.0f",iVienna+2,CentEdges[CentrCode],CentEdges[CentrCode+1]));
            if(!SPvProfV0C) { AliError(Form("Profile '%s' not found","<v>_2{SP_V0C}")); return; }
            SPvProfV0C->Fill(dPt, vnSPV0C[iVienna]);

            TProfile* SPvProfV0A_q2V0C = (TProfile*)SPFlowEseList->FindObject(Form("v_{%i}{SP V0A}_%.0f_%.0f_PerCode%i",iVienna+2,CentEdges[CentrCode],CentEdges[CentrCode+1],q2ESECodeV0C+1));
            if(!SPvProfV0A_q2V0C) { AliError(Form("Profile '%s' not found","_q2V0C <v>_2{SP_V0A}")); return; }
            SPvProfV0A_q2V0C->Fill(dPt, vnSPV0A[iVienna]);

            TProfile* SPvProfV0C_q2V0C = (TProfile*)SPFlowEseList->FindObject(Form("v_{%i}{SP V0C}_%.0f_%.0f_PerCode%i",iVienna+2,CentEdges[CentrCode],CentEdges[CentrCode+1],q2ESECodeV0C+1));
            if(!SPvProfV0C_q2V0C) { AliError(Form("Profile '%s' not found","_q2V0C <v>_2{SP_V0C}")); return; }
            SPvProfV0C_q2V0C->Fill(dPt, vnSPV0C[iVienna]);
        }

        //}
    }

    for (Int_t iVienna(0); iVienna<2; ++iVienna){
        Double_t QQsV0AV0C = QxnV0CCorr[iVienna]*QxnV0ACorr[iVienna] + QynV0CCorr[iVienna]*QynV0ACorr[iVienna];
        Double_t QQsV0ATPC = QxnV0ACorr[iVienna]*QxnTPCSP[iVienna] + QynV0ACorr[iVienna]*QynTPCSP[iVienna];
        Double_t QQsV0CTPC = QxnV0CCorr[iVienna]*QxnTPCSP[iVienna] + QynV0CCorr[iVienna]*QynTPCSP[iVienna];

        TProfile* prof1 = (TProfile*)SPFlowList->FindObject(Form("<QQ*>_%i{SP_V0A_TPC}",iVienna+2));
        if(!prof1) { AliError(Form("Profile '%s' not found","<QQ*>_2{SP_V0A_TPC}")); return; }
        prof1->Fill(centrality, QQsV0ATPC);

        TProfile* prof2 = (TProfile*)SPFlowList->FindObject(Form("<QQ*>_%i{SP_V0A_V0C}",iVienna+2));
        if(!prof2) { AliError(Form("Profile '%s' not found","<QQ*>_2{SP_V0A_V0C}")); return; }
        prof2->Fill(centrality, QQsV0AV0C);

        TProfile* prof3 = (TProfile*)SPFlowList->FindObject(Form("<QQ*>_%i{SP_V0C_TPC}",iVienna+2));
        if(!prof3) { AliError(Form("Profile '%s' not found","<QQ*>_2{SP_V0C_TPC}")); return; }
        prof3->Fill(centrality, QQsV0CTPC);

        TProfile* prof1V0C = (TProfile*)SPFlowEseList->FindObject(Form("<QQ*>_%i{SP_V0A_TPC}_PerCode%i",iVienna+2,q2ESECodeV0C+1));
        if(!prof1V0C) { AliError(Form("Profile '%s' not found","ESE V0C <QQ*>_2{SP_V0A_TPC}")); return; }
        prof1V0C->Fill(centrality, QQsV0ATPC);

        TProfile* prof2V0C = (TProfile*)SPFlowEseList->FindObject(Form("<QQ*>_%i{SP_V0A_V0C}_PerCode%i",iVienna+2,q2ESECodeV0C+1));
        if(!prof2V0C) { AliError(Form("Profile '%s' not found","ESE V0C <QQ*>_2{SP_V0A_V0C}")); return; }
        prof2V0C->Fill(centrality, QQsV0AV0C);

        TProfile* prof3V0C = (TProfile*)SPFlowEseList->FindObject(Form("<QQ*>_%i{SP_V0C_TPC}_PerCode%i",iVienna+2,q2ESECodeV0C+1));
        if(!prof3V0C) { AliError(Form("Profile '%s' not found","<QQ*>_2{SP_V0C_TPC}")); return; }
        prof3V0C->Fill(centrality, QQsV0CTPC);
    }
    
    return;
}
void AliAnalysisTaskESEFlow::FillCorrelation(Float_t centrality, Double_t dPt, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C, Int_t q2ESECodeV0A, Int_t q3ESECodeV0A, Bool_t doRef, Bool_t doDiff)
{
    Int_t iNumTask = fVecCorrTask.size();
    for(Int_t iTask(0); iTask < iNumTask; ++iTask){

        const AliUniFlowCorrTask* const task = fVecCorrTask.at(iTask);
        if(!task) { AliError("AliUniFlowCorrTask does not exist"); return; }

        Int_t CentrCode = GetCentralityCode(centrality);

        Int_t corrOrder= task->fiNumHarm;

        TComplex cNum = TComplex(0.0,0.0,kFALSE);
        TComplex cDn = TComplex(0.0,0.0,kFALSE);
        TComplex cNumDiff = TComplex(0.0,0.0,kFALSE);
        TComplex cDnDiff = TComplex(0.0,0.0,kFALSE);

        switch (corrOrder)
        {
        case 2 :
            if(!bHasGap) { 
                if(doDiff) {
                    cDnDiff = TwoDiff(0,0);
                    cNumDiff = TwoDiff(task->fiHarm[0],task->fiHarm[1]);
                }
                if(doRef) { 
                    cDn = Two(0,0);
                    cNum = Two(task->fiHarm[0],task->fiHarm[1]);
                }
            }
            else {
                if(doDiff) {
                    cDnDiff = TwoDiffGap10M(0,0);
                    cNumDiff = TwoDiffGap10M(task->fiHarm[0],task->fiHarm[1]);
                }
                if(doRef) {
                    cDn = TwoGap10(0,0);
                    cNum = TwoGap10(task->fiHarm[0],task->fiHarm[1]);
                }
            }
            break;
        case 3 :
            if(!bHasGap){
                if(doDiff){
                    cDnDiff = ThreeDiff(0,0,0);
                    cNumDiff = ThreeDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
                }
                if(doRef){
                    cDn = Three(0,0,0);
                    cNum = Three(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
                }
            }
            else{
                if(doDiff){
                    cDnDiff = ThreeDiffGapM(0,0,0);
                    cNumDiff = ThreeDiffGapM(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
                }
                if(doRef){
                    cDn = ThreeGapM(0,0,0);
                    cNum = ThreeGapM(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
                }
            }
            break;
        case 4 :
            if(!bHasGap){
                if(doDiff){
                    cDnDiff = FourDiff(0,0,0,0);
                    cNumDiff = FourDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
                }
                if(doRef) {
                    cDn = Four(0,0,0,0);
                    cNum = Four(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
                }
            }
            else {
                if(doDiff){
                    cDnDiff = FourDiffGap10M(0,0,0,0);
                    cNumDiff = FourDiffGap10M(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
                }
                if(doRef) {
                    cDn = FourGap10(0,0,0,0);
                    cNum = FourGap10(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
                }
            }
            break;
        case 6 :
            if(!bHasGap){
                if(doDiff){
                    cDnDiff = SixDiff(0,0,0,0,0,0);
                    cNumDiff = SixDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
                }
                if(doRef){
                    cDn = Six(0,0,0,0,0,0);
                    cNum = Six(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
                }
            }
            else{
                if(doDiff){
                    cDnDiff = SixDiffGap10M(0,0,0,0,0,0);
                    cNumDiff = SixDiffGap10M(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
                }
                if(doRef){
                    cDn = SixGap10(0,0,0,0,0,0);
                    cNum = SixGap10(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
                }
            }
            break;
        case 8 :
            if(!bHasGap){
                if(doDiff){
                    return;
                }
                if(doRef){
                    cDn = Eight(0,0,0,0,0,0,0,0);
                    cNum = Eight(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
                }
            }
            else{
                return;
            }
            break;
        default:
            return;
        }

        if(doRef){
            Double_t dDn = cDn.Re();
            Double_t dNum = cNum.Re();
            Double_t dValue = 0.0;
            Bool_t bFillPos = kFALSE;

            if(dDn > 0.0) {bFillPos = kTRUE; dValue = dNum/dDn; }
            if(bFillPos && TMath::Abs(dValue) > 1.0) { bFillPos = kFALSE; }
            if(!bFillPos) { return; }

            TProfile* prof = (TProfile*)fCorrDist->FindObject(Form("%s_sample%d",task->fsName.Data(),fIndexSampling));
            if(!prof) { AliError(Form("Profile %s_sample%d not found",task->fsName.Data(),fIndexSampling)); return; }
            prof->Fill(centrality, dValue, dDn);

            

            if(fTPCEse){
                TProfile* profESETPCq2 = (TProfile*)fcnESETPC->FindObject(Form("%s_q2TPC_PerCode%i_sample%d",task->fsName.Data(),q2ESECodeTPC+1,fIndexSampling));
                if(!profESETPCq2) { AliError(Form("Profile %s_q2TPC_PerCode%i_sample%d not found",task->fsName.Data(),q2ESECodeTPC+1,fIndexSampling)); return; }
                profESETPCq2->Fill(centrality, dValue, dDn);
        
                TProfile* profESETPCq3 = (TProfile*)fcnESETPC->FindObject(Form("%s_q3TPC_PerCode%i_sample%d",task->fsName.Data(),q3ESECodeTPC+1,fIndexSampling));
                if(!profESETPCq3) { AliError(Form("Profile %s_q3TPC_PerCode%i_sample%d not found",task->fsName.Data(),q3ESECodeTPC+1,fIndexSampling)); return; }
                profESETPCq3->Fill(centrality, dValue, dDn);        
            }

            if(fV0CEse){
                TProfile* profESEV0Cq2 = (TProfile*)fcnESEV0C->FindObject(Form("%s_q2V0C_PerCode%i_sample%d",task->fsName.Data(),q2ESECodeV0C+1,fIndexSampling));
                if(!profESEV0Cq2) { AliError(Form("Profile %s_q2V0C_PerCode%i_sample%d not found",task->fsName.Data(),q2ESECodeV0C+1,fIndexSampling)); return; }
                profESEV0Cq2->Fill(centrality, dValue, dDn);
            
                TProfile* profESEV0Cq3 = (TProfile*)fcnESEV0C->FindObject(Form("%s_q3V0C_PerCode%i_sample%d",task->fsName.Data(),q3ESECodeV0C+1,fIndexSampling));
                if(!profESEV0Cq3) { AliError(Form("Profile %s_q3V0C_PerCode%i_sample%d not found",task->fsName.Data(),q3ESECodeV0C+1,fIndexSampling)); return; }
                profESEV0Cq3->Fill(centrality, dValue, dDn);
            }

            if(fV0AEse){
                TProfile* profESEV0Aq2 = (TProfile*)fcnESEV0A->FindObject(Form("%s_q2V0A_PerCode%i_sample%d",task->fsName.Data(),q2ESECodeV0A+1,fIndexSampling));
                if(!profESEV0Aq2) { AliError(Form("Profile %s_q2V0A_PerCode%i_sample%d not found",task->fsName.Data(),q2ESECodeV0A+1,fIndexSampling)); return; }
                profESEV0Aq2->Fill(centrality, dValue, dDn);
            
                TProfile* profESEV0Aq3 = (TProfile*)fcnESEV0A->FindObject(Form("%s_q3V0A_PerCode%i_sample%d",task->fsName.Data(),q3ESECodeV0A+1,fIndexSampling));
                if(!profESEV0Aq3) { AliError(Form("Profile %s_q3V0C_PerCode%i_sample%d not found",task->fsName.Data(),q3ESECodeV0A+1,fIndexSampling)); return; }
                profESEV0Aq3->Fill(centrality, dValue, dDn);
            }
        }
        

        if(doDiff && corrOrder < 6){
            Double_t dDnDiff = cDnDiff.Re();
            Double_t dNumDiff = cNumDiff.Re();
            Double_t dValueDiff = 0.0;
            Bool_t bFillDiff = kFALSE;
            if(dDnDiff > 0.0) { bFillDiff = kTRUE; dValueDiff = dNumDiff/dDnDiff; }
            if(bFillDiff && TMath::Abs(dValueDiff) > 1.0) { bFillDiff = kFALSE; }

            if(!bFillDiff) { return; }
            
            TProfile* dn = (TProfile*)fpTDiff->FindObject(Form("%s_diff_%.0f_%.0f_sample%d",task->fsName.Data(),CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));

            if(!dn) { AliError(Form("Profile %s_diff_sample%d not found in centrality range: %.0f-%.0f",task->fsName.Data(),fIndexSampling,CentEdges[CentrCode],CentEdges[CentrCode+1])); return;}
            dn->Fill(dPt,dValueDiff,dDnDiff);

            if(fTPCEse){
                TProfile* profESETPCptdiffq2 = (TProfile*)fpTDiffESETPC->FindObject(Form("%s_diff_q2TPC_PerCode%i_%.0f_%.0f_sample%d",task->fsName.Data(),q2ESECodeTPC+1,CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));
                if(!profESETPCptdiffq2) { AliError(Form("Profile %s_q2TPC_PerCode%i_sample%d pt differential not found",task->fsName.Data(),q2ESECodeTPC+1,fIndexSampling)); return; }
                profESETPCptdiffq2->Fill(dPt,dValueDiff,dDnDiff);
            
                TProfile* profESETPCptdiffq3 = (TProfile*)fpTDiffESETPC->FindObject(Form("%s_diff_q3TPC_PerCode%i_%.0f_%.0f_sample%d",task->fsName.Data(),q3ESECodeTPC+1,CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));
                if(!profESETPCptdiffq3) { AliError(Form("Profile %s_q3TPC_PerCode%i_sample%d pt differential not found",task->fsName.Data(),q3ESECodeTPC+1,fIndexSampling)); return; }
                profESETPCptdiffq3->Fill(dPt,dValueDiff,dDnDiff);
            }

            if(fV0CEse){
                TProfile* profESEV0Cptdiffq2 = (TProfile*)fpTDiffESEV0C->FindObject(Form("%s_diff_q2V0C_PerCode%i_%.0f_%.0f_sample%d",task->fsName.Data(),q2ESECodeV0C+1,CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));
                if(!profESEV0Cptdiffq2) { AliError(Form("Profile %s_q2V0C_PerCode%i_sample%d pt differential not found",task->fsName.Data(),q2ESECodeV0C+1,fIndexSampling)); return; }
                profESEV0Cptdiffq2->Fill(dPt,dValueDiff,dDnDiff);
            
                TProfile* profESEV0Cptdiffq3 = (TProfile*)fpTDiffESEV0C->FindObject(Form("%s_diff_q3V0C_PerCode%i_%.0f_%.0f_sample%d",task->fsName.Data(),q3ESECodeV0C+1,CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));
                if(!profESEV0Cptdiffq3) { AliError(Form("Profile %s_q3V0C_PerCode%i_sample%d pt differential not found",task->fsName.Data(),q3ESECodeV0C+1,fIndexSampling)); return; }
                profESEV0Cptdiffq3->Fill(dPt,dValueDiff,dDnDiff);
            }

            if(fV0AEse){
                TProfile* profESEV0Aptdiffq2 = (TProfile*)fpTDiffESEV0A->FindObject(Form("%s_diff_q2V0A_PerCode%i_%.0f_%.0f_sample%d",task->fsName.Data(),q2ESECodeV0A+1,CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));
                if(!profESEV0Aptdiffq2) { AliError(Form("Profile %s_q2V0A_PerCode%i_sample%d pt differential not found",task->fsName.Data(),q2ESECodeV0A+1,fIndexSampling)); return; }
                profESEV0Aptdiffq2->Fill(dPt,dValueDiff,dDnDiff);
            
                TProfile* profESEV0Aptdiffq3 = (TProfile*)fpTDiffESEV0A->FindObject(Form("%s_diff_q3V0A_PerCode%i_%.0f_%.0f_sample%d",task->fsName.Data(),q3ESECodeV0A+1,CentEdges[CentrCode],CentEdges[CentrCode+1],fIndexSampling));
                if(!profESEV0Aptdiffq3) { AliError(Form("Profile %s_q3V0A_PerCode%i_sample%d pt differential not found",task->fsName.Data(),q3ESECodeV0A+1,fIndexSampling)); return; }
                profESEV0Aptdiffq3->Fill(dPt,dValueDiff,dDnDiff);
            }
        }
    }

    return;
}
void AliAnalysisTaskESEFlow::FillqnRedTPC(const Float_t centrality)
{
    fq2TPC->Fill(centrality,qnTPC[0]);
    fq3TPC->Fill(centrality,qnTPC[1]);

    return;
}
void AliAnalysisTaskESEFlow::FillqnRedV0(const Float_t centrality, TString V0type)
{
    if(V0type=="V0C"){
    fq2V0C->Fill(centrality,qnV0C[0]);
    fq3V0C->Fill(centrality,qnV0C[1]);
    }
    if(V0type=="V0A"){
    fq2V0A->Fill(centrality,qnV0A[0]);
    fq3V0A->Fill(centrality,qnV0A[1]);
    }
    else { return; } 

    return;
}
Bool_t AliAnalysisTaskESEFlow::WithinRFP(const AliVParticle* track) const
{
    if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fFlowRFPsPtMin > 0.0 && track->Pt() < fFlowRFPsPtMin) { return kFALSE; }
    if(fFlowRFPsPtMax > 0.0 && track->Pt() > fFlowRFPsPtMax) { return kFALSE; }
    
    return kTRUE;
}
Bool_t AliAnalysisTaskESEFlow::WithinPOI(const AliVParticle* track) const
{
    if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fFlowPOIsPtMin > 0.0 && track->Pt() < fFlowPOIsPtMin) { return kFALSE; }
    if(fFlowPOIsPtMax > 0.0 && track->Pt() > fFlowPOIsPtMax) { return kFALSE; }
    
    return kTRUE;
}
Bool_t AliAnalysisTaskESEFlow::WithinRedRP(const AliVParticle* track) const
{
    if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fRedFlowPtMin > 0.0 && track->Pt() < 0.0) { return kFALSE; }
    if(fRedFlowPtMax > 0.0 && track->Pt() > 20.0) { return kFALSE; }
    
    return kTRUE;
}
Bool_t AliAnalysisTaskESEFlow::LoadWeights()
{
    if(bUseOwnWeights)
    {
        TList* listFlowWeights = nullptr;
        
        Int_t runno = fAOD->GetRunNumber();
        listFlowWeights = (TList*) fFlowWeightsList->FindObject(Form("%i", runno));

        if(!listFlowWeights) 
        {
            // run-specific weights not found for this run; loading run-averaged instead
            AliWarning(Form("TList with flow weights (run %i) not found. returning", runno));
            return kFALSE;
        }
        fh3Weights = (TH3F*) listFlowWeights->FindObject("3DRefs");
    }
    else
    {
        Int_t runno = fAOD->GetRunNumber();
        fWeights = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",runno));
        if(!fWeights)
        {
            printf("Weights could not be found in list!\n");
            return kFALSE;
        }
        fWeights->CreateNUA();
    }
    
    return kTRUE;
}
Bool_t AliAnalysisTaskESEFlow::LoadNUE()
{
    if(!fFileTrackEff) { return kFALSE; }

    //fFileTrackEff->ls();

    Int_t runno = fAOD->GetRunNumber();

    if(!fIs2018Data){
      if (fNUE == 1){
          fhTrackNUE = (TH3F*)fFileTrackEff->Get(Form("eff_LHC15o_AMPT_%i",runno));
      }
      if (fNUE == 2){
          fhTrackNUE = (TH3F*)fFileTrackEff->Get(Form("eff_LHC15o_HIJING_%i",runno));
      }
    }
    if(fIs2018Data){
      printf("There is currently no NUE correction for 18 data! \n");
      return kFALSE; //find 2018 NUE
    }
    
    if(!fhTrackNUE) { printf("Problems loading NUE in LoadNUE()! \n"); return kFALSE; }

    return kTRUE;
}
Bool_t AliAnalysisTaskESEFlow::LoadV0Calibration()
{
    TList* listV0CalibRbr = nullptr;

    Int_t runno = fAOD->GetRunNumber();

    listV0CalibRbr = (TList*) fV0CalibList->FindObject(Form("%i",runno));

    if(!listV0CalibRbr){ AliWarning(Form("V0 Calibration for run %i not found. \n",runno)); return kFALSE; }

    fhV0Calib = (TH1F*) listV0CalibRbr->FindObject("V0MultCorr");

    for (Int_t i(0);i<2;++i){
        fQnxV0Cm[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQx%iV0Cm",i+2));
        fQnyV0Cm[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQy%iV0Cm",i+2));
        fQnxV0Am[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQx%iV0Am",i+2));
        fQnyV0Am[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQy%iV0Am",i+2));
        fQnxTPCm[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQx%iTPCm",i+2));
        fQnyTPCm[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQy%iTPCm",i+2));


        fQnxV0Cs[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQx%iV0Cs",i+2));
        fQnyV0Cs[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQy%iV0Cs",i+2));
        fQnxV0As[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQx%iV0As",i+2));
        fQnyV0As[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQy%iV0As",i+2));
    }

    return kTRUE;
}
Bool_t AliAnalysisTaskESEFlow::LoadqSelection()
{
    if(!fqSelList){ printf("q-selection not found when loading \n"); return kFALSE; }

    for (Int_t iSp(0); iSp<90; ++iSp){
        if(fTPCEse){
            fSplq2TPC[iSp] = (TSpline3*)fqSelList->FindObject(Form("sp_q2TPC_%i",iSp));
            fSplq3TPC[iSp] = (TSpline3*)fqSelList->FindObject(Form("sp_q3TPC_%i",iSp));
        }
        if(fV0CEse){
            fSplq2V0C[iSp] = (TSpline3*)fqSelList->FindObject(Form("sp_q2V0C_%i",iSp));
            fSplq3V0C[iSp] = (TSpline3*)fqSelList->FindObject(Form("sp_q3V0C_%i",iSp));
        }
        if(fV0AEse){
            fSplq2V0A[iSp] = (TSpline3*)fqSelList->FindObject(Form("sp_q2V0A_%i",iSp));
            fSplq3V0A[iSp] = (TSpline3*)fqSelList->FindObject(Form("sp_q3V0A_%i",iSp));
        }
    }

    return kTRUE;
}
Double_t AliAnalysisTaskESEFlow::GetFlowWeight(const AliAODTrack* track, const float dVz) const
{    
    Double_t dWeight = 1.0;

    if(fFlowRunByRunWeights){
        if(bUseOwnWeights)
        {
        Int_t iBin = fh3Weights->FindFixBin(track->Phi(),track->Eta(),dVz); //fh3Weights
        dWeight = fh3Weights->GetBinContent(iBin);
        }
        else
        {
            dWeight = fWeights->GetNUA(track->Phi(),track->Eta(),dVz);
        }

        if (dWeight <= 0.0) { dWeight = 1.0; }
    }

    return dWeight;
}
Double_t AliAnalysisTaskESEFlow::GetNUEPtWeight(Double_t pt, Double_t eta, const float dVz) const
{
    Double_t dPtWeight = 1;

    if(!fFileTrackEff) { return dPtWeight; }
    if(!fhTrackNUE) { return dPtWeight; }

    Double_t binPt = fhTrackNUE->GetXaxis()->FindBin(pt);
    Double_t binEta = fhTrackNUE->GetYaxis()->FindBin(eta);
    Double_t binVz = fhTrackNUE->GetZaxis()->FindBin(dVz);

    Double_t eff = fhTrackNUE->GetBinContent(binPt, binEta, binVz);
    Double_t error = fhTrackNUE->GetBinError(binPt, binEta, binVz);

    if((eff < 0.03) || ((error/eff) > 0.1)) dPtWeight = 1.0;
    else{
        TRandom3 r(0);
        Double_t efficiency = 0.0;
        efficiency = r.Gaus(eff, error);
        dPtWeight = 1./efficiency;
    }

    return dPtWeight;
}
Int_t AliAnalysisTaskESEFlow::GetSamplingIndex() const
{
    if(!fSampling) { return 0; }

    TRandom3 r(0);

    Double_t RandNum = r.Rndm();
    Double_t RandGen = RandNum*fNumSamples;

    Int_t index = 0;
    for(Int_t i(0); i < fNumSamples; ++i)
    {
        if(RandGen < (i+1)) { index = i; break; }
    }

    return index;
}
Int_t AliAnalysisTaskESEFlow::GetCentralityCode(const Float_t centrality)
{
    Int_t centrcode = -1;

    /*if ((centrality >= 0) && (centrality <= 5.)){
        centrcode = 0;
    }
    else if ((centrality > 5.) && (centrality <= 10.)){
        centrcode = 1;
    }
    else if ((centrality > 10.) && (centrality <= 20.)){
        centrcode = 2;
    }
    else if ((centrality > 20.) && (centrality <= 30.)){
        centrcode = 3;
    }
    else if ((centrality > 30.) && (centrality <= 40.)){
        centrcode = 4;
    }
    else if ((centrality > 40.) && (centrality <= 50.)){
        centrcode = 5;
    }
    else if ((centrality > 50.) && (centrality <= 60.)){
        centrcode = 6;
    }
    else if ((centrality > 60.) && (centrality <= 70.)){
        centrcode = 7;
    }
    else if ((centrality > 70.) && (centrality <= 80.)){
        centrcode = 8;
    }
    else if ((centrality > 80.) && (centrality <= 90.)){
        centrcode = 9;
    }
    else if (centrality > 90.){
        centrcode = 10;
    }*/
    
    for (Int_t centr(0); centr<fCentAxis->GetNbins(); ++centr){
        if ((centrality > CentEdges[centr]) && (centrality < CentEdges[centr+1])){
            centrcode = centr;
            return centrcode;
        }
    }
    
    return centrcode;
}
Int_t AliAnalysisTaskESEFlow::GetEsePercentileCode(Double_t qPerc) const
{
    Int_t qPerccode = -1;

    /*if ((qPerc >= 0) && (qPerc <= 10.)){
        qPerccode = 0;
    }
    else if ((qPerc > 10.) && (qPerc <= 20.)){
        qPerccode = 1;
    }
    else if ((qPerc > 20.) && (qPerc <= 30.)){
        qPerccode = 2;
    }
    else if ((qPerc > 30.) && (qPerc <= 40.)){
        qPerccode = 3;
    }
    else if ((qPerc > 40.) && (qPerc <= 50.)){
        qPerccode = 4;
    }
    else if ((qPerc > 50.) && (qPerc <= 60.)){
        qPerccode = 5;
    }
    else if ((qPerc > 60.) && (qPerc <= 70.)){
        qPerccode = 6;
    }
    else if ((qPerc > 70.) && (qPerc <= 80.)){
        qPerccode = 7;
    }
    else if ((qPerc > 80.) && (qPerc <= 90.)){
        qPerccode = 8;
    }
    else if (qPerc > 90.){
        qPerccode = 9;
    }*/
    
    for (Int_t code(0); code<ESEPercAxis->GetNbins(); ++code){
        if ((qPerc > EventShapeEdges[code]) && (qPerc < EventShapeEdges[code+1])){
            qPerccode = code;
            return qPerccode;
        }
    }
    
    return qPerccode;
}
// ######################### Generic FW #########################
void AliAnalysisTaskESEFlow::ResetFlowVector(TComplex (&array)[fNumHarms][fNumPowers])
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm) {
    for(Int_t iPower(0); iPower < fNumPowers; ++iPower) {
      array[iHarm][iPower](0.0,0.0);
    }
  }
  return;
}
void AliAnalysisTaskESEFlow::ResetReducedqVector(double (&array)[2])
{
  // RESET Reduced q vector
  // *************************************************************
  for(Int_t iHarm(0); iHarm < 2; ++iHarm) {
      array[iHarm]=(0.0);
  }
  return;
}
//____________________________________________________________________
Bool_t AliAnalysisTaskESEFlow::ProcessMCParticles()
{
    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!AODMCTrackArray) return kFALSE;

    for(Int_t iPart(0); iPart<AODMCTrackArray->GetEntriesFast(); ++iPart)
    {
        AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(iPart));
        if (!particle) continue;

        fHistPDG->Fill(particle->GetPdgCode());
    }

    return kTRUE;
}
//____________________________________________________________________
Bool_t AliAnalysisTaskESEFlow::IsEventSelected()
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  fhEventCounter->Fill("Loaded OK",1);
  
  AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
  if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
  Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator); 
  

  if(fIs2018Data){
    fEventCuts.SetupPbPb2018();
    fhEventCounter->Fill("2018 OK",1);
  }

  if(!fIs2018Data){
    if(!(fSelectMask & fTrigger)) { return kFALSE; }
  }
  else{
    if((dPercentile<10) || (dPercentile>30 && dPercentile<50)){
      if(!(fSelectMask & (AliVEvent::kCentral|AliVEvent::kSemiCentral|fTrigger))) { return kFALSE; }
    }
    else{
      if(!(fSelectMask & fTrigger)) { return kFALSE;}
    }
  }
  fhEventCounter->Fill("Triggers OK",1);



  
  if(!fEventCuts.AcceptEvent(fAOD)) { return kFALSE; }
  fhEventCounter->Fill("EventCuts OK",1);

  if(dPercentile > 100 || dPercentile < 0) { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1; }
  fhEventCounter->Fill("Centrality Cut OK",1);

  if( (fColSystem == kPbPb || fColSystem == kXeXe) && fEventRejectAddPileUp && dPercentile > 0 && dPercentile < 10 && IsEventRejectedAddPileUp()) { return kFALSE; }
  fhEventCounter->Fill("Pileup Cut OK",1);

  if(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ()) > fVtxZCuts) { return kFALSE; }
  fhEventCounter->Fill("VtxZ Cut OK",1);

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskESEFlow::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < fCutChargedNumTPCclsMin && fFilterBit != 2) { return kFALSE; }
  //if(fPtMin > 0 && track->Pt() < fPtMin) { return kFALSE; }
  //if(fPtMax > 0 && track->Pt() > fPtMax) { return kFALSE; }
  if(fAbsEtaMax > 0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }

  if(fCheckChi2TPC){
    Float_t chi2tpc = track->GetTPCchi2();
    Float_t ntpccls = track->GetTPCNcls();

    Float_t tpcchi2percls = (ntpccls==0)?0.0:chi2tpc/ntpccls;

    if (tpcchi2percls > vTPCChi2Bound) { return kFALSE; }
  } 
  if(fCheckChi2ITS){
    Float_t chi2its = track->GetITSchi2();
    Float_t nitscls = track->GetITSNcls();

    Float_t itschi2percls = (nitscls==0)?0.0:chi2its/nitscls;

    if (itschi2percls > vITSChi2Bound) { return kFALSE; }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskESEFlow::IsEventRejectedAddPileUp() const
{
  // Check for additional pile-up rejection in Run 2 Pb-Pb collisions (15o, 17n)
  // based on multiplicity correlations
  // ***************************************************************************

  Bool_t bIs17n = kFALSE;
  Bool_t bIs15o = kFALSE;

  Int_t iRunNumber = fAOD->GetRunNumber();
  if(iRunNumber >= 244824 && iRunNumber <= 246994) { bIs15o = kTRUE; }
  else if(iRunNumber == 280235 || iRunNumber == 20234) { bIs17n = kTRUE; }
  else { return kFALSE; }

  // recounting multiplcities
  const Int_t multESD = ((AliAODHeader*) fAOD->GetHeader())->GetNumberOfESDTracks();
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTPC32 = 0;
  Int_t multTPC128 = 0;
  Int_t multTPC96 = 0;
  Int_t multTOF = 0;
  Int_t multTrk = 0;
  Double_t multESDTPCdif = 0.0;
  Double_t v0Centr = 0.0;

  for(Int_t it(0); it < nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*) fAOD->GetTrack(it);
    if(!track) { continue; }

    if(track->TestFilterBit(32))
    {
      multTPC32++;
      if(TMath::Abs(track->GetTOFsignalDz()) <= 10.0 && track->GetTOFsignal() >= 12000.0 && track->GetTOFsignal() <= 25000.0) { multTOF++; }
      if((TMath::Abs(track->Eta())) < fAbsEtaMax && (track->GetTPCNcls() >= fCutChargedNumTPCclsMin) && (track->Pt() >= fFlowRFPsPtMin) && (track->Pt() < fFlowRFPsPtMax)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }

    if(track->TestFilterBit(96)) { multTPC96++; }
  }

  if(bIs17n)
  {
    multESDTPCdif = multESD - (6.6164 + 3.64583*multTPC128 + 0.000126397*multTPC128*multTPC128);
    if(multESDTPCdif > 1000) { return kTRUE; }
    if( ((AliAODHeader*) fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) { return kTRUE; }
  }

  if(bIs15o)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > fPileupCut) { return kTRUE; }

    TF1 fMultTOFLowCut = TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFLowCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) < fMultTOFLowCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    TF1 fMultTOFHighCut = TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFHighCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) > fMultTOFHighCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
    v0Centr = multSelection->GetMultiplicityPercentile("V0M");

    TF1 fMultCentLowCut = TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCentLowCut.SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
    if(Double_t(multTrk) < fMultCentLowCut.Eval(v0Centr)) { return kTRUE; }
  }

  



  return kFALSE;
}
void AliAnalysisTaskESEFlow::QAMultFiller(Float_t v0Centr)
{
    // recounting multiplcities
  const Int_t multESD = ((AliAODHeader*) fAOD->GetHeader())->GetNumberOfESDTracks();
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTPC32 = 0;
  Int_t multTPC128 = 0;
  Int_t multTPC96 = 0;
  Int_t multTOF = 0;
  Int_t multTrk = 0;
  //Double_t multESDTPCdif = 0.0;
  //Double_t v0Centr = 0.0;


  for(Int_t it(0); it < nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*) fAOD->GetTrack(it);
    if(!track) { continue; }

    if(track->TestFilterBit(32))
    {
      multTPC32++;
      if(TMath::Abs(track->GetTOFsignalDz()) <= 10.0 && track->GetTOFsignal() >= 12000.0 && track->GetTOFsignal() <= 25000.0) { multTOF++; }
      if((TMath::Abs(track->Eta())) < fAbsEtaMax && (track->GetTPCNcls() >= fCutChargedNumTPCclsMin) && (track->Pt() >= fFlowRFPsPtMin) && (track->Pt() < fFlowRFPsPtMax)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }

    if(track->TestFilterBit(96)) { multTPC96++; }
  }


  fhQAEventsfMult32vsCentr->Fill(v0Centr, multTrk);
  fhQAEventsfMult128vsCentr->Fill(v0Centr, multTPC128);
  fhQAEventsfMult96vsCentr->Fill(v0Centr, multTPC96);
  fhQAEventsfMultTPCvsTOF->Fill(multTPC32, multTOF);
  fhQAEventsfMultTPCvsESD->Fill(multTPC128, multESD);

    
    

  return;
}
//_____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::QGap10M(int n, int p)
{

	if(n>=0) return Qvector10M[n][p];
  else return TComplex::Conjugate(Qvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::QGap10P(int n, int p)
{

	if(n>=0) return Qvector10P[n][p];
  else return TComplex::Conjugate(Qvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::pGap10M(int n, int p)
{

	if(n>=0) return pvector10M[n][p];
	else return TComplex::Conjugate(pvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::pGap10P(int n, int p)
{

	if(n>=0) return pvector10P[n][p];
	else return TComplex::Conjugate(pvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::qGap10M(int n, int p)
{

    if(n>=0) return pvector10M[n][p];
    else return TComplex::Conjugate(pvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::qGap10P(int n, int p)
{

    if(n>=0) return pvector10P[n][p];
    else return TComplex::Conjugate(pvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::pPtA(int n, int p)
{

    if(n>=0) return pvector[n][p];
    else return TComplex::Conjugate(pvector[-n][p]);

}

//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::pPtB(int n, int p)
{

    if(n>=0) return pvectorPtB[n][p];
    else return TComplex::Conjugate(pvectorPtB[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::qPtA(int n, int p)
{

    if(n>=0) return qvector[n][p];
    else return TComplex::Conjugate(qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::qPtB(int n, int p)
{

    if(n>=0) return qvectorPtB[n][p];
    else return TComplex::Conjugate(qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::pPtBGap10M(int n, int p)
{

    if(n>=0) return pvectorPtB10M[n][p];
    else return TComplex::Conjugate(pvectorPtB10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::pPtBGap10P(int n, int p)
{

    if(n>=0) return pvectorPtB10P[n][p];
    else return TComplex::Conjugate(pvectorPtB10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Two(int n1, int n2)
{
	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoGap10(int n1, int n2)
{
	TComplex formula = QGap10M(n1,1)*QGap10P(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiff(int n1, int n2)
{
	TComplex formula = p(n1,1)*Q(n2,1) - q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiffGap10P(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*QGap10P(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiffGap10M(int n1, int n2)
{
    TComplex formula = pGap10P(n1,1)*QGap10M(n2,1);
    return formula;
}
//____________________________________________________________________
// 2-particles from the same pt
TComplex AliAnalysisTaskESEFlow::TwoDiff_Pt(int n1, int n2)
{
    TComplex formula = pPtA(n1,1)*pPtA(n2,1) - pPtA(n1+n2,2);
    return formula;
}
//____________________________________________________________________
// 2-particles from the same pt but two different eta regions
TComplex AliAnalysisTaskESEFlow::TwoDiffGap10_Pt(int n1, int n2)
{
    TComplex formula = pGap10P(n1,1)*pGap10M(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiff_PtA(int n1, int n2)
{
    TComplex formula = pPtA(n1,1)*pPtA(n2,1) - qPtA(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiffGap10M_PtA(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*pGap10M(n2,1) - pGap10M(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiffGap10P_PtB(int n1, int n2)
{
    TComplex formula = pPtBGap10P(n1,1)*pPtBGap10P(n2,1) - pPtBGap10P(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiff_PtB(int n1, int n2)
{
    TComplex formula = pPtB(n1,1)*pPtB(n2,1) - qPtB(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiff_PtA_PtB(int n1, int n2)
{
    TComplex formula = pPtA(n1,1)*pPtB(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::TwoDiffGap10_PtA_PtB(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*pPtBGap10P(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Three(int n1, int n2, int n3)
{
    
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
    - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
    return formula;
    
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::ThreeGapM(int n1, int n2, int n3)
{
    TComplex formula = QGap10M(n1,1)*QGap10P(n2,1)*QGap10P(n3,1) - QGap10M(n1,1)*QGap10P(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::ThreeGapP(int n1, int n2, int n3)
{
    TComplex formula = QGap10P(n1,1)*QGap10M(n2,1)*QGap10M(n3,1) - QGap10P(n1,1)*QGap10M(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::ThreeDiff(int n1, int n2, int n3)
{

    TComplex formula = p(n1,1)*Q(n2,1)*Q(n3,1)-q(n1+n2,2)*Q(n3,1)-q(n1+n3,2)*Q(n2,1)
    - p(n1,1)*Q(n2+n3,2)+2.*q(n1+n2+n3,3);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::ThreeDiffGapP(int n1, int n2, int n3)
{
    TComplex formula = pGap10P(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)- pGap10P(n1,1)*QGap10M(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::ThreeDiffGapM(int n1, int n2, int n3)
{
    TComplex formula = pGap10M(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)- pGap10M(n1,1)*QGap10P(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Four(int n1, int n2, int n3, int n4)
{
	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                 		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                 		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                 		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                 		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourGap10(int n1, int n2, int n3, int n4)
{
    TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3,1)*QGap10M(n4,1)-QGap10P(n1+n2,2)*QGap10M(n3,1)*QGap10M(n4,1)
    -QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3+n4,2)+QGap10P(n1+n2,2)*QGap10M(n3+n4,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourDiff(int n1, int n2, int n3, int n4)
{

	TComplex formula = p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*q(n1+n3,2)*Q(n4,1)
                 		- p(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*q(n1+n4,2)
                 		+ Q(n2+n3,2)*q(n1+n4,2)-p(n1,1)*Q(n3,1)*Q(n2+n4,2)+q(n1+n3,2)*Q(n2+n4,2)
                 		+ 2.*Q(n3,1)*q(n1+n2+n4,3)-p(n1,1)*Q(n2,1)*Q(n3+n4,2)+q(n1+n2,2)*Q(n3+n4,2)
                 		+ 2.*Q(n2,1)*q(n1+n3+n4,3)+2.*p(n1,1)*Q(n2+n3+n4,3)-6.*q(n1+n2+n3+n4,4);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourDiff_PtA_PtA(int n1, int n2, int n3, int n4)
{
    TComplex formula = pPtA(n1,1)*pPtA(n2,1)*pPtA(n3,1)*pPtA(n4,1)-pPtA(n1+n2,2)*pPtA(n3,1)*pPtA(n4,1)-pPtA(n2,1)*pPtA(n1+n3,2)*pPtA(n4,1)
    - pPtA(n1,1)*pPtA(n2+n3,2)*pPtA(n4,1)+2.*pPtA(n1+n2+n3,3)*pPtA(n4,1)-pPtA(n2,1)*pPtA(n3,1)*pPtA(n1+n4,2)
    + pPtA(n2+n3,2)*pPtA(n1+n4,2)-pPtA(n1,1)*pPtA(n3,1)*pPtA(n2+n4,2)+pPtA(n1+n3,2)*pPtA(n2+n4,2)
    + 2.*pPtA(n3,1)*pPtA(n1+n2+n4,3)-pPtA(n1,1)*pPtA(n2,1)*pPtA(n3+n4,2)+pPtA(n1+n2,2)*pPtA(n3+n4,2)
    + 2.*pPtA(n2,1)*pPtA(n1+n3+n4,3)+2.*pPtA(n1,1)*pPtA(n2+n3+n4,3)-6.*pPtA(n1+n2+n3+n4,4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourDiff_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_PtA(n1, n2)*TwoDiff_PtB(n3, n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourDiffGap10_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiffGap10M_PtA(n1, n2)*TwoDiffGap10P_PtB(n3, n4);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourDiffGap10P(int n1, int n2, int n3, int n4)
{
    TComplex formula = pGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3,1)*QGap10M(n4,1)-qGap10P(n1+n2,2)*QGap10M(n3,1)*QGap10M(n4,1)
    -pGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3+n4,2)+qGap10P(n1+n2,2)*QGap10M(n3+n4,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskESEFlow::FourDiffGap10M(int n1, int n2, int n3, int n4)
{
    TComplex formula = pGap10M(n1,1)*QGap10M(n2,1)*QGap10P(n3,1)*QGap10P(n4,1)-qGap10M(n1+n2,2)*QGap10P(n3,1)*QGap10P(n4,1)
    -pGap10M(n1,1)*QGap10M(n2,1)*QGap10P(n3+n4,2)+qGap10M(n1+n2,2)*QGap10P(n3+n4,2);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Five(int n1, int n2, int n3, int n4, int n5)
{
    
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
    return formula;
    
}
//___________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
              - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
              - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
              - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
              + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
              - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
              - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
              - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
              - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
              - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
              - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
              - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
              + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
              - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
              - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
              - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
              - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
              + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
              + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
              - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
              - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
              - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
              + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
              - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
              + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
              + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
              - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
              - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
              - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
              - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
              + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
              - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
              - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
              - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
              - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
              + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
              - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
              - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
              - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
              + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
              - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
              + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
              - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
              - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
              - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
              + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
              + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
              - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
              - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
              - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
              - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
              + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
              + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
              - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
              - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
              - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
              - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
              + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
              - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
              - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
              - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
              + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
              + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
              + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
              + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
              + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
              + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
              - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
              + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
              - 120.*Q(n1+n2+n3+n4+n5+n6,6);
    return formula;
}
//_____________________________________________________________________________
TComplex AliAnalysisTaskESEFlow::SixDiff(int n1, int n2, int n3, int n4, int n5, int n6)
{
    TComplex formula = p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    - Q(n2,1)*q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-p(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    + 2.*q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*q(n1+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n2+n3,2)*q(n1+n4,2)*Q(n5,1)*Q(n6,1)-p(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
    + q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
    - p(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n2,1)*q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*p(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
    - 6.*q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*q(n1+n5,2)*Q(n6,1)
    + Q(n2+n3,2)*Q(n4,1)*q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*q(n1+n5,2)*Q(n6,1)
    + Q(n2,1)*Q(n3+n4,2)*q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*q(n1+n5,2)*Q(n6,1)
    - p(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
    + Q(n3,1)*q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+p(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
    - 2.*q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*q(n1+n2+n5,3)*Q(n6,1)
    - 2.*Q(n3+n4,2)*q(n1+n2+n5,3)*Q(n6,1)-p(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
    + q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
    + p(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n4,1)*q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*q(n1+n3+n5,3)*Q(n6,1)
    + 2.*p(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
    - 6.*Q(n4,1)*q(n1+n2+n3+n5,4)*Q(n6,1)-p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
    + q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
    + p(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n3,1)*q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*q(n1+n4+n5,3)*Q(n6,1)
    + 2.*p(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
    - 6.*Q(n3,1)*q(n1+n2+n4+n5,4)*Q(n6,1)+2.*p(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
    - 2.*q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*q(n1+n3+n4+n5,4)*Q(n6,1)
    - 6.*p(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*q(n1+n2+n3+n4+n5,5)*Q(n6,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*q(n1+n6,2)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*q(n1+n6,2)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*q(n1+n6,2)
    - Q(n3+n4,2)*Q(n2+n5,2)*q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*q(n1+n6,2)
    - Q(n2+n4,2)*Q(n3+n5,2)*q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*q(n1+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*q(n1+n6,2)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*q(n1+n6,2)
    + 6.*Q(n2+n3+n4+n5,4)*q(n1+n6,2)-p(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
    + q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
    + p(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
    + Q(n3,1)*Q(n4,1)*q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*q(n1+n5,2)*Q(n2+n6,2)
    + p(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
    - 2.*Q(n4,1)*q(n1+n3+n5,3)*Q(n2+n6,2)+p(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
    - q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*q(n1+n4+n5,3)*Q(n2+n6,2)
    - 2.*p(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*q(n1+n3+n4+n5,4)*Q(n2+n6,2)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*q(n1+n2+n6,3)
    - 2.*Q(n4,1)*Q(n3+n5,2)*q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*q(n1+n2+n6,3)
    + 4.*Q(n3+n4+n5,3)*q(n1+n2+n6,3)-p(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
    + q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
    + p(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
    + Q(n2,1)*Q(n4,1)*q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*q(n1+n5,2)*Q(n3+n6,2)
    + p(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
    - 2.*Q(n4,1)*q(n1+n2+n5,3)*Q(n3+n6,2)+p(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
    - q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*q(n1+n4+n5,3)*Q(n3+n6,2)
    - 2.*p(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*q(n1+n2+n4+n5,4)*Q(n3+n6,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*q(n1+n3+n6,3)
    - 2.*Q(n4,1)*Q(n2+n5,2)*q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*q(n1+n3+n6,3)
    + 4.*Q(n2+n4+n5,3)*q(n1+n3+n6,3)+2.*p(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
    - 2.*q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*q(n1+n5,2)*Q(n2+n3+n6,3)
    - 2.*p(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*q(n1+n4+n5,3)*Q(n2+n3+n6,3)
    - 6.*Q(n4,1)*Q(n5,1)*q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*q(n1+n2+n3+n6,4)
    - p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
    + Q(n2,1)*q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+p(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
    - 2.*q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*q(n1+n5,2)*Q(n4+n6,2)
    - Q(n2+n3,2)*q(n1+n5,2)*Q(n4+n6,2)+p(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
    - q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*q(n1+n2+n5,3)*Q(n4+n6,2)
    + p(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
    - 2.*Q(n2,1)*q(n1+n3+n5,3)*Q(n4+n6,2)-2.*p(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
    + 6.*q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*q(n1+n4+n6,3)
    - 2.*Q(n2+n3,2)*Q(n5,1)*q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*q(n1+n4+n6,3)
    - 2.*Q(n2,1)*Q(n3+n5,2)*q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*q(n1+n4+n6,3)
    + 2.*p(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
    - 2.*Q(n3,1)*q(n1+n5,2)*Q(n2+n4+n6,3)-2.*p(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
    + 4.*q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*q(n1+n2+n4+n6,4)
    + 6.*Q(n3+n5,2)*q(n1+n2+n4+n6,4)+2.*p(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
    - 2.*q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*q(n1+n5,2)*Q(n3+n4+n6,3)
    - 2.*p(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*q(n1+n2+n5,3)*Q(n3+n4+n6,3)
    - 6.*Q(n2,1)*Q(n5,1)*q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*q(n1+n3+n4+n6,4)
    - 6.*p(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*q(n1+n5,2)*Q(n2+n3+n4+n6,4)
    + 24.*Q(n5,1)*q(n1+n2+n3+n4+n6,5)-p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
    + q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
    + p(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
    + Q(n2,1)*Q(n3,1)*q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*q(n1+n4,2)*Q(n5+n6,2)
    + p(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
    - 2.*Q(n3,1)*q(n1+n2+n4,3)*Q(n5+n6,2)+p(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
    - q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*q(n1+n3+n4,3)*Q(n5+n6,2)
    - 2.*p(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*q(n1+n2+n3+n4,4)*Q(n5+n6,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*q(n1+n5+n6,3)
    - 2.*Q(n3,1)*Q(n2+n4,2)*q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*q(n1+n5+n6,3)
    + 4.*Q(n2+n3+n4,3)*q(n1+n5+n6,3)+2.*p(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
    - 2.*q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*q(n1+n4,2)*Q(n2+n5+n6,3)
    - 2.*p(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*q(n1+n3+n4,3)*Q(n2+n5+n6,3)
    - 6.*Q(n3,1)*Q(n4,1)*q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*q(n1+n2+n5+n6,4)
    + 2.*p(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
    - 2.*Q(n2,1)*q(n1+n4,2)*Q(n3+n5+n6,3)-2.*p(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
    + 4.*q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*q(n1+n3+n5+n6,4)
    + 6.*Q(n2+n4,2)*q(n1+n3+n5+n6,4)-6.*p(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
    + 6.*q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*q(n1+n2+n3+n5+n6,5)
    + 2.*p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
    - 2.*Q(n2,1)*q(n1+n3,2)*Q(n4+n5+n6,3)-2.*p(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
    + 4.*q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*q(n1+n4+n5+n6,4)
    + 6.*Q(n2+n3,2)*q(n1+n4+n5+n6,4)-6.*p(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
    + 6.*q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*q(n1+n2+n4+n5+n6,5)
    - 6.*p(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*q(n1+n2,2)*Q(n3+n4+n5+n6,4)
    + 24.*Q(n2,1)*q(n1+n3+n4+n5+n6,5)+24.*p(n1,1)*Q(n2+n3+n4+n5+n6,5)
    - 120.*q(n1+n2+n3+n4+n5+n6,6);
    return formula;
}

//_________________________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{
    
    TComplex Correlation = {0, 0};
    int Narray[] = {n1, n2, n3, n4, n5, n6};
    
    for(int k=7; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4
        
        int array[6] = {0,1,2,3,4,5};
        int iPerm = 0;
        //int argument = 0;
        int count = 0;
        
        // k==6: there is just one combination, we can add it manually
        if(k==6){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
            Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
        }// k==6
        
        else if(k==5){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                    Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
                         Narray[int(array[3])], Narray[int(array[4])])*
                    Q(Narray[int(array[5])]+n7, 7-k);
                }
            }while(std::next_permutation(array, array+6));
        }// k==5
        
        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
                             Narray[int(array[3])])*
                        Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==4
        
        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
                        Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==3
        
        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Two(Narray[int(array[0])], Narray[int(array[1])])*
                        Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
                          +Narray[int(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==2
        
        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
        }// k==1
        
        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
        }// k==0
        
        else{
            std::cout<<"invalid range of k"<< std::endl;
            return {0,0};
        }
        
    }// loop over k
    
    return Correlation;
    
}
//_____________________________________________________________________________
TComplex AliAnalysisTaskESEFlow::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
    
    TComplex Correlation = {0, 0};
    int Narray[] = {n1, n2, n3, n4, n5, n6, n7};
    
    for(int k=8; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4
        
        int array[7] = {0,1,2,3,4,5,6};
        int iPerm = 0;
        //int argument = 0;
        int count = 0;
        
        // k==7: there is just one combination, we can add it manually
        if(k==7){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
            Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
        }// k==7
        
        else if(k==6){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                    Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
                        Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
                    Q(Narray[int(array[6])]+n8, 8-k);
                }
            }while(std::next_permutation(array, array+7));
        }// k==6
        
        else if(k==5){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
                             Narray[int(array[3])], Narray[int(array[4])])*
                        Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==5
        
        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
                        Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==4
        
        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
                        Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==3
        
        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%120 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Two(Narray[int(array[0])], Narray[int(array[1])])*
                        Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
                          +Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==2
        
        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
        }// k==1
        
        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
        }// k==0
        
        else{
            std::cout<<"invalid range of k"<<std::endl;
            return {0,0};
        }
        
    }// loop over k
    
    return Correlation;
}
TComplex AliAnalysisTaskESEFlow::SixGap10(int n1, int n2, int n3, int n4, int n5, int n6)
{
  TComplex formula = ThreePos(n1,n2,n3)*ThreeNeg(n4,n5,n6);
	return formula;
}
TComplex AliAnalysisTaskESEFlow::ThreePos(int n1, int n2, int n3)
{
  TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)-QGap10P(n1+n2,2)*QGap10P(n3,1)-QGap10P(n2,1)*QGap10P(n1+n3,2)
 		                 - QGap10P(n1,1)*QGap10P(n2+n3,2)+2.*QGap10P(n1+n2+n3,3);
  return formula;
}
TComplex AliAnalysisTaskESEFlow::ThreeNeg(int n1, int n2, int n3)
{
  TComplex formula = QGap10M(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)-QGap10M(n1+n2,2)*QGap10M(n3,1)-QGap10M(n2,1)*QGap10M(n1+n3,2)
 		                 - QGap10M(n1,1)*QGap10M(n2+n3,2)+2.*QGap10M(n1+n2+n3,3);
  return formula;
}
TComplex AliAnalysisTaskESEFlow::SixDiffGap10M(int n1, int n2, int n3, int n4, int n5, int n6)
{
  TComplex formula = ThreeDiffNeg(n1,n2,n3)*ThreePos(n4,n5,n6);
  return formula;
}
TComplex AliAnalysisTaskESEFlow::ThreeDiffNeg(int n1, int n2, int n3)
{
  TComplex formula = pGap10M(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)-qGap10M(n1+n2,2)*QGap10M(n3,1)-qGap10M(n1+n3,2)*QGap10M(n2,1)
 		                 - pGap10M(n1,1)*QGap10M(n2+n3,2)+2.0*qGap10M(n1+n2+n3,3);
  return formula;
}
TComplex AliAnalysisTaskESEFlow::ThreeDiffPos(int n1, int n2, int n3)
{
  TComplex formula = pGap10P(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)-qGap10P(n1+n2,2)*QGap10P(n3,1)-qGap10P(n1+n3,2)*QGap10P(n2,1)
 		                 - pGap10P(n1,1)*QGap10P(n2+n3,2)+2.0*qGap10P(n1+n2+n3,3);
  return formula;
}
//_____________________________________________________________________________
