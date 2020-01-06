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
 * flow class
 * Event Shape Engineering
 * event selection
 * track selection
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

#include <iostream>

class AliAnalysisTaskESEFlow;

ClassImp(AliAnalysisTaskESEFlow)

AliAnalysisTaskESEFlow::AliAnalysisTaskESEFlow() : AliAnalysisTaskSE(),
    fEventCuts(),
    fFlowRunByRunWeights(kTRUE),
    bUseOwnWeights(0),
    dGap(0.5),
    fInit(kFALSE),
    fqRun(kFALSE),
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
    fQAEvents(0),
    fFlowWeightsList{nullptr},
    fWeights(0),
    fV0CalibList(0),
    fqSelList(0),
    fHistPhiEtaVz(0),
    fHistPhi(0),
    fHistEta(0),
    fHistPt(0),
    fHistZVertex(0),
    fSplq2TPC{0},
    fSplq3TPC{0},
    fSplq2V0C{0},
    fSplq3V0C{0},
    fSplq2V0A{0},
    fSplq3V0A{0},
    fcn2Gap{0},
    fcn2GapInclusive{0},
    fdn2GapPt{0},
    fdn2GapPtB{0},
    fh2Weights{nullptr},
    fhV0Calib{nullptr},
    fHistPDG{0},
    fcn2GapESETPC{0},
    fdn2GapESETPC{0},
    fcn4Gap(0),
    fdn4GapPt{0},
    fcn4GapESETPC{0},
    fdn4GapESETPC{0},
    fcn2GapESEV0C{0},
    fdn2GapESEV0C{0},
    fcn4GapESEV0C{0},
    fdn4GapESEV0C{0},
    fcn2GapESEV0A{0},
    fdn2GapESEV0A{0},
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
    fQnxV0CEse{0},
    fQnyV0CEse{0},
    fQnxV0AEse{0},
    fQnyV0AEse{0},
    fQnxTPCEse{0},
    fQnyTPCEse{0},
    fCentq2TPCvsv22(0),
    fCentq2V0Cvsv22(0),
    fProfNPar(0),
    fhV0Multiplicity(0),
    fhV0CorrMult(0),
    fhq2TPCvq2V0C(0),
    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fAbsEtaMax(0.8),
    fCentEstimator(),
    fReadMC(kFALSE),
    fMCEvent(0),
    fFlowRFPsPtMin(0.2),
    fFlowRFPsPtMax(5.0),
    fFlowPOIsPtMin(0.0),
    fFlowPOIsPtMax(10.0),
    fnTwoCorr(kFALSE),
    fnFourCorr(kFALSE),
    fTPCEse(kTRUE),
    fV0CEse(kTRUE),
    fV0AEse(kFALSE)
{}
//_____________________________________________________________________________
AliAnalysisTaskESEFlow::AliAnalysisTaskESEFlow(const char* name) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fFlowRunByRunWeights(kTRUE),
    bUseOwnWeights(0),
    dGap(0.5),
    fInit(kFALSE),
    fqRun(kFALSE),
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
    fQAEvents(0),
    fFlowWeightsList{nullptr},
    fWeights(0),
    fV0CalibList(0),
    fqSelList(0),
    fHistPhiEtaVz(0),
    fHistPhi(0),
    fHistEta(0),
    fHistPt(0),
    fHistZVertex(0),
    fSplq2TPC{0},
    fSplq3TPC{0},
    fSplq2V0C{0},
    fSplq3V0C{0},
    fSplq2V0A{0},
    fSplq3V0A{0},
    fcn2Gap{0},
    fcn2GapInclusive{0},
    fdn2GapPt{0},
    fdn2GapPtB{0},
    fh2Weights{nullptr},
    fhV0Calib{nullptr},
    fHistPDG{0},
    fcn2GapESETPC{0},
    fdn2GapESETPC{0},
    fcn4Gap(0),
    fdn4GapPt{0},
    fcn4GapESETPC{0},
    fdn4GapESETPC{0},
    fcn2GapESEV0C{0},
    fdn2GapESEV0C{0},
    fcn4GapESEV0C{0},
    fdn4GapESEV0C{0},
    fcn2GapESEV0A{0},
    fdn2GapESEV0A{0},
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
    fQnxV0CEse{0},
    fQnyV0CEse{0},
    fQnxV0AEse{0},
    fQnyV0AEse{0},
    fQnxTPCEse{0},
    fQnyTPCEse{0},
    fCentq2TPCvsv22(0),
    fCentq2V0Cvsv22(0),
    fProfNPar(0),
    fhV0Multiplicity(0),
    fhV0CorrMult(0),
    fhq2TPCvq2V0C(0),
    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fAbsEtaMax(0.8),
    fCentEstimator(),
    fReadMC(kFALSE),
    fMCEvent(0),
    fFlowRFPsPtMin(0.2),
    fFlowRFPsPtMax(5.0),
    fFlowPOIsPtMin(0.0),
    fFlowPOIsPtMax(10.0),
    fnTwoCorr(kFALSE),
    fnFourCorr(kFALSE),
    fTPCEse(kTRUE),
    fV0CEse(kTRUE),
    fV0AEse(kFALSE)
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
    if(fQAEvents) { delete fQAEvents; }
}
Bool_t AliAnalysisTaskESEFlow::InitializeTask()
{
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

    //Load V0 Calibration
    fV0CalibList = static_cast<TList*>(GetInputData(2));
    if(!fV0CalibList) { AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n V0 Calibration list not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }

    if(!fqRun){
        fqSelList = static_cast<TList*>(GetInputData(3));
        if(!fqSelList) {AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n q-selection Splines list not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }

        if(!LoadqSelection()) { AliFatal("\n \n \n \n \n \n \n \n \n \n q-Splines not loaded! Terminating! \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }
    }
    
    AliInfo("Initialization succes");
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
    fQAEvents->SetOwner(kTRUE);

    //RUN INITIALIZE TASK
    fInit = InitializeTask();
    if(!fInit) { return; }

    const int NvnPtBin = 28;
    double PtEdgesvn[NvnPtBin+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,7.0,8.0,9.0,10.0};

    fHistPhiEtaVz = new TH3F("fHistPhiEtaVz", "fHistPhiEtaVz; #phi; #eta; Vz", 120, 0.0, TMath::TwoPi(), 120, -1.0, 1.0,100,-15,15);
    fHistPhiEtaVz->Sumw2();
    fHistPhi = new TH1F("fHistPhi", ";#phi", 120, 0.0, TMath::TwoPi());
    fHistEta = new TH1F("fHistEta", ";#eta", 120,-1.0, 1.0);
    fHistPt = new TH1F("fHistPt", ";p_{T}", NvnPtBin,PtEdgesvn);
    fHistZVertex = new TH1F("fHistZVertex", ";Vtx_{Z}", 100,-15,15);
    fProfNPar = new TProfile("fProfNparvsCent",";Centrality;N_{Particles}",100,0,100);

    fhV0Multiplicity = new TH2F("fV0Multiplicity","",64,0,64,100,0,1250);
    fhV0Multiplicity->Sumw2();
    fhV0CorrMult = new TH2F("fV0CalibratedMultiplicity","",64,0,64,100,0,1250);
    fhV0CorrMult->Sumw2();

    fhq2TPCvq2V0C = new TH2F("fq2TPCvq2V0C","",100,0,8,100,0,15);
    fhq2TPCvq2V0C->Sumw2();

    fq2TPC = new TH2D("fq2vCentTPC","",100,0,100,100,0,8);
    fq2TPC->Sumw2();
    fq3TPC = new TH2D("fq3vCentTPC","",100,0,100,100,0,8);
    fq3TPC->Sumw2();

    fq2V0C = new TH2D("fq2vCentV0C","",100,0,100,100,0,15);
    fq2V0C->Sumw2();
    fq3V0C = new TH2D("fq3vCentV0C","",100,0,100,100,0,15);
    fq3V0C->Sumw2();
    fq2V0A = new TH2D("fq2vCentV0A","",100,0,100,100,0,15);
    fq2V0A->Sumw2();
    fq3V0A = new TH2D("fq3vCentV0A","",100,0,100,100,0,15);
    fq3V0A->Sumw2();

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
    }

    const int nBins = 11;
    double binedge[nBins+1] = {0, 5., 10., 20., 30., 40., 50., 60., 70., 80.,90.,100.};

    //Differential PT v's
    if(fnTwoCorr){
        for(Int_t fHarmNum(0); fHarmNum<fNumHarmHists; ++fHarmNum)
        {
            fcn2Gap[fHarmNum] = new TProfile(Form("c%i2GapCent",fHarmNum+2),"",nBins,binedge);
            fcn2Gap[fHarmNum]->Sumw2();
            fCorrDist->Add(fcn2Gap[fHarmNum]);

            fcn2GapInclusive[fHarmNum] = new TProfile(Form("c%i2GapnPart",fHarmNum+2),"",14,0,400);
            fcn2GapInclusive[fHarmNum]->Sumw2();
            fOutputList->Add(fcn2GapInclusive[fHarmNum]);

            for(Int_t qn(0); qn<2; ++qn){
                for(Int_t qBin(0); qBin<10; ++qBin){
                    if(fTPCEse){
                    fcn2GapESETPC[fHarmNum][qn][qBin] = new TProfile(Form("c%i2GapESETPCq%iPercCode%i",fHarmNum+2,qn+2,qBin+1),"",nBins, binedge);
                    fcn2GapESETPC[fHarmNum][qn][qBin]->Sumw2();
                    fcnESETPC->Add(fcn2GapESETPC[fHarmNum][qn][qBin]);
                    }
                    if(fV0CEse){
                    fcn2GapESEV0C[fHarmNum][qn][qBin] = new TProfile(Form("c%i2GapESEV0Cq%iPercCode%i",fHarmNum+2,qn+2,qBin+1),"",nBins, binedge);
                    fcn2GapESEV0C[fHarmNum][qn][qBin]->Sumw2();
                    fcnESEV0C->Add(fcn2GapESEV0C[fHarmNum][qn][qBin]);
                    }
                    if(fV0AEse){
                    fcn2GapESEV0A[fHarmNum][qn][qBin] = new TProfile(Form("c%i2GapESEV0Aq%iPercCode%i",fHarmNum+2,qn+2,qBin+1),"",nBins, binedge);
                    fcn2GapESEV0A[fHarmNum][qn][qBin]->Sumw2();
                    fcnESEV0C->Add(fcn2GapESEV0A[fHarmNum][qn][qBin]);
                    }
                }
            }

            for(Int_t fCentNum(0) ; fCentNum<fNumCentHists; ++fCentNum)
            {
                fdn2GapPt[fHarmNum][fCentNum] = new TProfile(Form("d%i2Gap_%.0f_%.0f",fHarmNum+2,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin,PtEdgesvn);
                fdn2GapPt[fHarmNum][fCentNum]->Sumw2();
                fpTDiff->Add(fdn2GapPt[fHarmNum][fCentNum]);
            
                fdn2GapPtB[fHarmNum][fCentNum] = new TProfile(Form("d%i2GapB_%.0f_%.0f",fHarmNum+2,binedge[fCentNum],binedge[fCentNum+1]),"", NvnPtBin, PtEdgesvn);
                fdn2GapPtB[fHarmNum][fCentNum]->Sumw2();
                fpTDiff->Add(fdn2GapPtB[fHarmNum][fCentNum]);

                //ESE
                for(Int_t qn(0); qn<2; ++qn){
                    for(Int_t qBin(0); qBin<10; ++qBin){
                        if(fTPCEse){
                        fdn2GapESETPC[fHarmNum][qn][fCentNum][qBin] = new TProfile(Form("d%i2GapESETPCq%iPerCode%i_%.0f_%.0f",fHarmNum+2,qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                        fdn2GapESETPC[fHarmNum][qn][fCentNum][qBin]->Sumw2();
                        fpTDiffESETPC->Add(fdn2GapESETPC[fHarmNum][qn][fCentNum][qBin]);
                        }
                        if(fV0CEse){
                        fdn2GapESEV0C[fHarmNum][qn][fCentNum][qBin] = new TProfile(Form("d%i2GapESEV0Cq%iPerCode%i_%.0f_%.0f",fHarmNum+2,qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                        fdn2GapESEV0C[fHarmNum][qn][fCentNum][qBin]->Sumw2();
                        fpTDiffESEV0C->Add(fdn2GapESEV0C[fHarmNum][qn][fCentNum][qBin]);
                        }
                        if(fV0AEse){
                        fdn2GapESEV0A[fHarmNum][qn][fCentNum][qBin] = new TProfile(Form("d%i2GapESEV0Aq%iPerCode%i_%.0f_%.0f",fHarmNum+2,qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                        fdn2GapESEV0A[fHarmNum][qn][fCentNum][qBin]->Sumw2();
                        fpTDiffESEV0C->Add(fdn2GapESEV0A[fHarmNum][qn][fCentNum][qBin]);
                        }
                    }
                }
            }
        }
    }

    if(fnFourCorr){
        fcn4Gap = new TProfile("c24GapCent","",nBins,binedge);
        fcn4Gap->Sumw2();
        fCorrDist->Add(fcn4Gap);

        for(Int_t qn(0); qn<2; ++qn){
            for(Int_t qBin(0); qBin<10; ++qBin){
                if(fTPCEse){
                fcn4GapESETPC[qn][qBin] = new TProfile(Form("c24GapESETPCq%iPercCode%i",qn+2,qBin+1),"",nBins, binedge);
                fcn4GapESETPC[qn][qBin]->Sumw2();
                fcnESETPC->Add(fcn4GapESETPC[qn][qBin]);
                }
                if(fV0CEse){
                fcn4GapESEV0C[qn][qBin] = new TProfile(Form("c24GapESEV0Cq%iPercCode%i",qn+2,qBin+1),"",nBins, binedge);
                fcn4GapESEV0C[qn][qBin]->Sumw2();
                fcnESEV0C->Add(fcn4GapESEV0C[qn][qBin]);
                }
            }
        }

        for(int fCentNum(0);fCentNum<fNumCentHists;++fCentNum){
            fdn4GapPt[fCentNum] = new TProfile(Form("d24_%.0f_%.0f",binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin,PtEdgesvn);
            fdn4GapPt[fCentNum]->Sumw2();
            fpTDiff->Add(fdn4GapPt[fCentNum]);

            for(Int_t qn(0); qn<2; ++qn){
                for(Int_t qBin(0); qBin<10; ++qBin){
                    if(fTPCEse){
                    fdn4GapESETPC[qn][fCentNum][qBin] = new TProfile(Form("d24GapESETPCq%iPerCode%i_%.0f_%.0f",qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                    fdn4GapESETPC[qn][fCentNum][qBin]->Sumw2();
                    fpTDiffESETPC->Add(fdn4GapESETPC[qn][fCentNum][qBin]);
                    }
                    if(fV0CEse){
                    fdn4GapESEV0C[qn][fCentNum][qBin] = new TProfile(Form("d24GapESEV0Cq%iPerCode%i_%.0f_%.0f",qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                    fdn4GapESEV0C[qn][fCentNum][qBin]->Sumw2();
                    fpTDiffESEV0C->Add(fdn4GapESEV0C[qn][fCentNum][qBin]);
                    }
                }
            }
        }
    }

    fCentq2TPCvsv22 = new TH3F("Centq2TPCv22","",100,0,100,100,0,10,100,0,0.02);
    fCentq2TPCvsv22->Sumw2();
    fOutputList->Add(fCentq2TPCvsv22);
    fCentq2V0Cvsv22 = new TH3F("Centq2V0Cv22","",100,0,100,100,0,15,100,0,0.02);
    fCentq2V0Cvsv22->Sumw2();
    fOutputList->Add(fCentq2V0Cvsv22);

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
    

    fEventCuts.AddQAplotsToList(fQAEvents); //QA plots
    
    fObservables->Add(fHistPhiEtaVz);
    fObservables->Add(fHistPhi);
    fObservables->Add(fHistEta);
    fObservables->Add(fHistPt);
    fObservables->Add(fHistZVertex);
    fObservables->Add(fProfNPar);
    fObservables->Add(fhV0Multiplicity);
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
    }
    fqnDist->Add(fhq2TPCvq2V0C);


    PostData(1, fOutputList);
    PostData(2, fObservables);
    PostData(3, fCorrDist);
    PostData(4, fpTDiff);
    PostData(5, fqnDist);
    PostData(6, fpTDiffESETPC);
    PostData(7, fcnESETPC);
    PostData(8, fpTDiffESEV0C);
    PostData(9, fcnESEV0C);
    PostData(10, fQAEvents);
}
//_____________________________________________________________________________
void AliAnalysisTaskESEFlow::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { return; }

    if(!IsEventSelected()) { return; }

    
    Int_t iTracks(fAOD->GetNumberOfTracks());
    //VERTEX
    float dVz = fAOD->GetPrimaryVertex()->GetZ();

    //centrality
    Float_t centrality(0);
    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(multSelection) centrality = multSelection->GetMultiplicityPercentile(fCentEstimator);

    if(fFlowRunByRunWeights){
        if(!LoadWeights()) { AliFatal("\n \n \n \n \n \n \n \n \n \n Weights not loaded! \n \n \n \n \n \n \n \n \n \n "); return; }
    }

    if(!LoadV0Calibration()) { AliFatal("\n \n \n \n \n \n \n \n \n \n V0 Calibration not loaded! \n \n \n \n \n \n \n \n \n \n "); return; }
    

    FillObsDistributions(iTracks, fAOD, dVz, centrality);

    Int_t fSpCent = (Int_t)centrality;

    CorrelationTask(centrality, iTracks, fAOD, dVz, fSpCent);

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
    PostData(10, fQAEvents);
}
//_____________________________________________________________________________
void AliAnalysisTaskESEFlow::Terminate(Option_t *)
{
}
void AliAnalysisTaskESEFlow::CorrelationTask(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, Int_t fSpCent)
{
    ReducedqVectorsTPC(centrality, iTracks, fAOD, dVz);
    ReducedqVectorsV0(centrality,fAOD, fSpCent);

    fhq2TPCvq2V0C->Fill(qnTPC[0],qnV0C[0]);

    if(!fqRun){
    Int_t CenterCode = GetCentrCode(centrality);

    if( (CenterCode < 0) || (CenterCode > 9)) { return; }

    Double_t q2TPCInp = 100.*fSplq2TPC[fSpCent]->Eval(qnTPC[0]);
    Double_t q3TPCInp = 100.*fSplq3TPC[fSpCent]->Eval(qnTPC[1]);

    Int_t q2ESECodeTPC = GetPercCode(q2TPCInp);
    if (q2ESECodeTPC<0) { printf("Problem with q_2 TPC percentile: negative percentile \n"); return; } 
    Int_t q3ESECodeTPC = GetPercCode(q3TPCInp);
    if (q2ESECodeTPC<0) { printf("Problem with q_3 TPC percentile: negative percentile \n"); return; } 

    Double_t q2V0CInp = 100.*fSplq2V0C[fSpCent]->Eval(qnV0C[0]);
    Double_t q3V0CInp = 100.*fSplq3V0C[fSpCent]->Eval(qnV0C[1]);

    Int_t q2ESECodeV0C = GetPercCode(q2V0CInp);
    if (q2ESECodeV0C<0) { printf("Problem with q_2 V0C percentile: negative percentile \n"); return; } 
    Int_t q3ESECodeV0C = GetPercCode(q3V0CInp);
    if (q3ESECodeV0C<0) { printf("Problem with q_3 V0C percentile: negative percentile \n"); return; } 

    Double_t q2V0AInp = 2.;//100.*fSplq2V0A[fSpCent]->Eval(qnV0A[0]); do q-selection for V0A
    Double_t q3V0AInp = 2.;//100.*fSplq3V0A[fSpCent]->Eval(qnV0A[1]);

    Int_t q2ESECodeV0A = GetPercCode(q2V0AInp);
    if (q2ESECodeV0A<0) { printf("Problem with q_2 V0A percentile: negative percentile \n"); return; } 
    Int_t q3ESECodeV0A = GetPercCode(q3V0AInp);
    if (q3ESECodeV0A<0) { printf("Problem with q_3 V0A percentile: negative percentile \n"); return; } 

    RFPVectors(centrality, iTracks, fAOD, dVz,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
    POIVectors(CenterCode, iTracks, fAOD, dVz,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
    }

    return;
}
void AliAnalysisTaskESEFlow::FillObsDistributions(const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, const Float_t centrality)
{
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
        fHistPhi->Fill(dPhi);
        fHistEta->Fill(dEta);
        fHistPt->Fill(dPt);
    } // ending track loop

    AliAODVZERO* aodV0 = fAOD->GetVZEROData();

    for (Int_t iV0=0; iV0 < 64; ++iV0){
        Float_t multV0 = aodV0->GetMultiplicity(iV0);
        fhV0Multiplicity->Fill(iV0,multV0);



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

    } // ending V0 calibration loop

    return;
}
void AliAnalysisTaskESEFlow::RFPVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C)
{
    ResetFlowVector(Qvector);
    ResetFlowVector(Qvector10P);
    ResetFlowVector(Qvector10M);


    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; ++i) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        if(!WithinRFP(track)){continue;} // check if also within reference particle

        double dEta = track->Eta();
        double dPhi = track->Phi();
        //double dPt = track->Pt();

        Double_t dWeight = GetFlowWeight(track,dVz);
        
        // no eta gap
        for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
        {
            for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
            {
                Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                Qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            }
        }
    
        // eta gap
        if(dEta > dGap)
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
            {
                for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                {
                Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                Qvector10P[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                }
            }
        } 
        // RFP in negative eta acceptance
        if(dEta < -dGap)
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
            { 
                for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                {
                Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                Qvector10M[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                }
            }
        }            
    } // ending track loop
    if(fnTwoCorr){
    FillRFP(centrality,iTracks,2,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);  //with gap for 2-particle correlation nHarm=2
    FillRFP(centrality,iTracks,3,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);  // nHarm=3
    FillRFP(centrality,iTracks,4,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);  // nHarm=4
    FillRFP(centrality,iTracks,5,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
    FillRFP(centrality,iTracks,6,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
    }
    if(fnFourCorr){
    FillRFP(centrality,iTracks,2,4,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);  //with gap for 4-particle correlation nHarm=2
    }

    return;
}
void AliAnalysisTaskESEFlow::FillPOI(const Double_t dPtL, const Double_t dPtLow, const Double_t dPtHigh, const float dVz, const Int_t iTracks)
{
    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; ++i)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }

        double dEta = track->Eta();
        double dPhi = track->Phi();
        double dPt = track->Pt();

        Bool_t bIsWithinPOI = WithinPOI(track);
        if(!bIsWithinPOI) {continue;}

        Double_t dWeight = GetFlowWeight(track,dVz);
        
        if(dPt > dPtLow && dPt <= dPtHigh)
        {
            // q vector
            if(WithinRFP(track)) // check if also within reference particle
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
            if(dEta > dGap)
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
            if(dEta < -dGap)
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
void AliAnalysisTaskESEFlow::POIVectors(const Int_t CenterCode, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C)
{
    // DIFFERENTIAL FLOW
    // loop over p_T histogram
    Int_t iNumPtBins = fHistPt->GetXaxis()->GetNbins();

    for(Int_t iPt(1); iPt< iNumPtBins+1; ++iPt)
    {
        // reset flow vector for each pt bin
        ResetFlowVector(qvector);
        ResetFlowVector(pvector);
        ResetFlowVector(pvector10P);
        ResetFlowVector(pvector10M);

        Double_t dPtLow = fHistPt->GetXaxis()->GetBinLowEdge(iPt);
        Double_t dPtHigh = fHistPt->GetXaxis()->GetBinUpEdge(iPt);
        Double_t dPtL = fHistPt->GetXaxis()->GetBinCenter(iPt);

        FillPOI(dPtL, dPtLow, dPtHigh, dVz, iTracks);
        
        if(fnTwoCorr){
        Filldn(CenterCode,dPtL,2,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); //Fill d2_2 <<2'>>  fnParCorr-particle correlation (std 2)
        Filldn(CenterCode,dPtL,3,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); 
        Filldn(CenterCode,dPtL,4,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); 
        Filldn(CenterCode,dPtL,5,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); 
        Filldn(CenterCode,dPtL,6,2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); 
        }
        if(fnFourCorr){
        Filldn(CenterCode,dPtL,2,4,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
        }
    }

    return;
}
void AliAnalysisTaskESEFlow::ReducedqVectorsTPC(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz)
{
    ResetReducedqVector(qnTPC);
    ResetReducedqVector(QxnTPC);
    ResetReducedqVector(QynTPC);

    Double_t M = 0;
    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; ++i) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        if(!WithinRFP(track)){continue;} // check if also within reference particle

        double dEta = track->Eta();
        double dPhi = track->Phi();
        //double dPt = track->Pt();

        Double_t dWeight = GetFlowWeight(track,dVz);

        if(- (0.4) < dEta && dEta < (0.4))
        {
            M += 1;
            for(Int_t iHarm(0); iHarm < 2; ++iHarm)
            {                
                Double_t dCos = dWeight * TMath::Cos( (iHarm+2) * dPhi);
                Double_t dSin = dWeight * TMath::Sin( (iHarm+2) * dPhi);
                QxnTPC[iHarm] += dCos;
                QynTPC[iHarm] += dSin;
            }
        }
    }

    for (Int_t nQ(0); nQ<2;++nQ){
        fQnxTPC[nQ]->Fill(centrality,QxnTPC[nQ]);
        fQnyTPC[nQ]->Fill(centrality,QynTPC[nQ]);
    } //used for recentering of Qnx/Qny

    if(M>0)
    {
        for(Int_t iHarm(0); iHarm < 2; ++iHarm)
        { 
            Double_t dqn = QxnTPC[iHarm]*QxnTPC[iHarm]+QynTPC[iHarm]*QynTPC[iHarm];
            qnTPC[iHarm] = TMath::Sqrt(dqn)/TMath::Sqrt(M);
        }
        FillqnRedTPC(centrality);
    }

    return;
}
void AliAnalysisTaskESEFlow::ReducedqVectorsV0(const Float_t centrality, const AliAODEvent* fAOD,const Int_t SPCode)
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
    }

    for (Int_t nQ(0);nQ<2;++nQ){
        fQnxV0C[nQ]->Fill(centrality,QxnV0C[nQ]); // fill for recentering
        fQnyV0C[nQ]->Fill(centrality,QynV0C[nQ]);
        fQnxV0A[nQ]->Fill(centrality,QxnV0A[nQ]);
        fQnyV0A[nQ]->Fill(centrality,QynV0A[nQ]); // End recenter fill

        fQnxV0CEse[nQ]->Fill(centrality,QxnV0CEse[nQ]); // fill after recenter
        fQnyV0CEse[nQ]->Fill(centrality,QynV0CEse[nQ]);
        fQnxV0AEse[nQ]->Fill(centrality,QxnV0AEse[nQ]);
        fQnyV0AEse[nQ]->Fill(centrality,QynV0AEse[nQ]); // end after rec
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
void AliAnalysisTaskESEFlow::FillRFP(const Float_t centrality,const Int_t iTracks,const int nHarm, const int nCorr, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C)
{
    //                          Calculate particle correlations 
    int bhistN = nHarm-2;

    double cNum=0.0;
    double cDenom=0.0;

    switch(nCorr){
        case 2: {
            cNum = TwoGap10(nHarm,-nHarm).Re();
            cDenom = TwoGap10(0,0).Re();
            break;
        }
        case 4: {
            cNum = FourGap10(nHarm,nHarm,-nHarm,-nHarm).Re();
            cDenom = FourGap10(0,0,0,0).Re();
            break;
        }
    }

    if(cDenom>0.0){
        double cn_fill = cNum/cDenom;
        if(TMath::Abs(cn_fill < 1.0)){
        if(nCorr==2) {
            fcn2Gap[bhistN]->Fill(centrality,cn_fill,cDenom);

            FillESEcn(centrality, nHarm, cn_fill, cDenom, 2,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); //fill integrated cqncut here

            if(nHarm==2){
                fCentq2TPCvsv22->Fill(centrality,qnTPC[0],cn_fill,cDenom);
                fCentq2V0Cvsv22->Fill(centrality,qnV0C[0],cn_fill,cDenom);
            }
        }// end 2 corr
        if(nCorr==4 && nHarm==2){
            fcn4Gap->Fill(centrality,cn_fill,cDenom);

            FillESEcn(centrality, nHarm, cn_fill, cDenom, 4,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C); //fill integrated cqncut here
        } // end 4-corr
    }
    }

    return;
}
void AliAnalysisTaskESEFlow::FillESEcn(const Float_t centrality, const int nHarm, const double c, const double c_weight,const int nCorr, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C)
{
    int nHist = nHarm-2;
    
    if(fTPCEse){
    if(nCorr==2) {
        fcn2GapESETPC[nHist][0][q2ESECodeTPC]->Fill(centrality,c,c_weight);
        fcn2GapESETPC[nHist][1][q3ESECodeTPC]->Fill(centrality,c,c_weight);
    }
    if(nCorr==4) {
        fcn4GapESETPC[0][q2ESECodeTPC]->Fill(centrality,c,c_weight);
        fcn4GapESETPC[1][q3ESECodeTPC]->Fill(centrality,c,c_weight);
    }
    }

    if(fV0CEse){
    if(nCorr==2) {
        fcn2GapESEV0C[nHist][0][q2ESECodeV0C]->Fill(centrality,c,c_weight);
        fcn2GapESEV0C[nHist][1][q3ESECodeV0C]->Fill(centrality,c,c_weight);
    }
    if(nCorr==4) {
        fcn4GapESEV0C[0][q2ESECodeV0C]->Fill(centrality,c,c_weight);
        fcn4GapESEV0C[1][q3ESECodeV0C]->Fill(centrality,c,c_weight);
    }
    }

    if(fV0AEse){
    if(nCorr==2) {
        fcn2GapESEV0A[nHist][0][q2ESECodeV0C]->Fill(centrality,c,c_weight);
        fcn2GapESEV0A[nHist][1][q3ESECodeV0C]->Fill(centrality,c,c_weight);
    }
    }

    return;
}
void AliAnalysisTaskESEFlow::FillESEdnPt(const Int_t CenterCode, const int nHarm, const Double_t dPt, const double d, const double d_weight,const int nCorr, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C)
{
    int nHist = nHarm-2;

    if(fTPCEse){
    if(nCorr==2) {
        fdn2GapESETPC[nHist][0][CenterCode][q2ESECodeTPC]->Fill(dPt,d,d_weight);
        fdn2GapESETPC[nHist][1][CenterCode][q3ESECodeTPC]->Fill(dPt,d,d_weight);
    }
    if(nCorr==4) {
        fdn4GapESETPC[0][CenterCode][q2ESECodeTPC]->Fill(dPt,d,d_weight);
        fdn4GapESETPC[1][CenterCode][q3ESECodeTPC]->Fill(dPt,d,d_weight);
    }
    }

    if(fV0CEse){
    if(nCorr==2) {
        fdn2GapESEV0C[nHist][0][CenterCode][q2ESECodeV0C]->Fill(dPt,d,d_weight);
        fdn2GapESEV0C[nHist][1][CenterCode][q3ESECodeV0C]->Fill(dPt,d,d_weight);
    }
    if(nCorr==4) {
        fdn4GapESEV0C[0][CenterCode][q2ESECodeV0C]->Fill(dPt,d,d_weight);
        fdn4GapESEV0C[1][CenterCode][q3ESECodeV0C]->Fill(dPt,d,d_weight);
    }
    }

    if(fV0AEse){
    if(nCorr==2) {
        fdn2GapESEV0A[nHist][0][CenterCode][q2ESECodeV0C]->Fill(dPt,d,d_weight);
        fdn2GapESEV0A[nHist][1][CenterCode][q3ESECodeV0C]->Fill(dPt,d,d_weight);
    }
    }

    return;
}
void AliAnalysisTaskESEFlow::Filldn(const Int_t CenterCode, const double dPt, const int nHarm, const int nCorr, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C)
{
    int bhistnumb = nHarm-2;
    double dNum=0.0;
    double dDenom=0.0;
    switch(nCorr){
        case 2: {
            dNum = TwoDiffGap10M(nHarm,-nHarm).Re();
            dDenom = TwoDiffGap10M(0,0).Re();
            break;
        }
        case 4: {
            dNum = FourDiffGap10M(nHarm,nHarm,-nHarm,-nHarm).Re();
            dDenom = FourDiffGap10M(0,0,0,0).Re();
            break;
        }
    }


    if(dDenom>0.0)
    {
        double dn_pt = dNum/dDenom;
        if(TMath::Abs(dn_pt)<1)
        {
            
            if(nCorr==2){
            fdn2GapPt[bhistnumb][CenterCode]->Fill(dPt,dn_pt,dDenom);
            }
            if(nCorr==4){
                if(nHarm==2){
                    fdn4GapPt[CenterCode]->Fill(dPt,dn_pt,dDenom);
                }
            }

            if(nCorr==2){
            FillESEdnPt(CenterCode, nHarm, dPt, dn_pt, dDenom, nCorr,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
            }
            if(nCorr==4){
            FillESEdnPt(CenterCode, nHarm, dPt, dn_pt, dDenom, nCorr,q2ESECodeTPC,q3ESECodeTPC,q2ESECodeV0C,q3ESECodeV0C);
            }
        }
    }

    if(nCorr==2){
    double dn2ptB = TwoDiffGap10P(nHarm,-nHarm).Re();
    double dn2pt_wB = TwoDiffGap10P(0,0).Re();
    if(dn2pt_wB!=0)
    {
        double dn_2ptB = dn2ptB/dn2pt_wB;
        if(TMath::Abs(dn_2ptB)<1)
        {
            fdn2GapPtB[bhistnumb][CenterCode]->Fill(dPt,dn_2ptB,dn2pt_wB);
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
Bool_t AliAnalysisTaskESEFlow::LoadWeights()
{
    if(bUseOwnWeights)
    {
        TList* listFlowWeights = nullptr;
        if(!fFlowRunByRunWeights)
        {
            listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
        }
        else
        {
            listFlowWeights = (TList*) fFlowWeightsList->FindObject(Form("%d", fAOD->GetRunNumber()));

            if(!listFlowWeights) 
            {
                // run-specific weights not found for this run; loading run-averaged instead
                AliWarning(Form("TList with flow weights (run %d) not found. Using run-averaged weights instead (as a back-up)", fAOD->GetRunNumber()));
                listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
                if(!listFlowWeights) { AliError("Loading run-averaged weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
            }
        }
        fh2Weights = (TH2D*) listFlowWeights->FindObject("Refs");
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
        //fQnxTPCm[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQx%iTPCm",i+2));
        //fQnyTPCm[i] = (TH1F*) listV0CalibRbr->FindObject(Form("hQy%iTPCm",i+2));
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
        Int_t iBin = fh2Weights->FindFixBin(track->Phi(),track->Eta()); //fh2Weights
        dWeight = fh2Weights->GetBinContent(iBin);
        }
        else
        {
            dWeight = fWeights->GetNUA(track->Phi(),track->Eta(),dVz);
        }

        if (dWeight <= 0.0) { dWeight = 1.0; }
    }

    return dWeight;
}
Int_t AliAnalysisTaskESEFlow::GetCentrCode(const Float_t centrality)
{
    Int_t centrcode = -1;

    if ((centrality >= 0) && (centrality <= 5.)){
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
    }
    
    return centrcode;
}
Int_t AliAnalysisTaskESEFlow::GetPercCode(Double_t qPerc) const
{
    Int_t qPerccode = -1;

    if ((qPerc >= 0) && (qPerc <= 10.)){
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
  if(!(fSelectMask & fTrigger)) { return kFALSE; }
  if(!fEventCuts.AcceptEvent(fAOD)) { return kFALSE; }
  AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
  if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
  Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
  if(dPercentile > 100 || dPercentile < 0) { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1; }
  if(fEventRejectAddPileUp && dPercentile > 0 && dPercentile < 10 && IsEventRejectedAddPileUp()) { return kFALSE; }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskESEFlow::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < 70 && fFilterBit != 2) { return kFALSE; }
  //if(fPtMin > 0 && track->Pt() < fPtMin) { return kFALSE; }
  //if(fPtMax > 0 && track->Pt() > fPtMax) { return kFALSE; }
  if(fAbsEtaMax > 0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
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
      if((TMath::Abs(track->Eta())) < fAbsEtaMax && (track->GetTPCNcls() >= 70) && (track->Pt() >= fFlowRFPsPtMin) && (track->Pt() < fFlowRFPsPtMax)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }
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
    if(multESDTPCdif > 500) { return kTRUE; }

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
//_____________________________________________________________________________
