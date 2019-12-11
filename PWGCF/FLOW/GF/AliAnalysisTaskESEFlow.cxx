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

#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TList.h"
#include "TMath.h"
#include "TComplex.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliGFWWeights.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskESEFlow.h"

#include <iostream>

class AliAnalysisTaskESEFlow;

ClassImp(AliAnalysisTaskESEFlow)

AliAnalysisTaskESEFlow::AliAnalysisTaskESEFlow() : AliAnalysisTaskSE(),
    fEventCuts(),
    fFlowRunByRunWeights(kTRUE),
    bUseOwnWeights(0),
    dGap(0.5),
    fInit{0},
    fqRun{0},
    fAOD(0),
    fOutputList(0),
    fObservables(0),
    fCorrDist(0),
    fpTDiff(0),
    fqnDist(0),
    fpTDiffqselec(0),
    fcnqselec(0),
    fQAEvents(0),
    fFlowWeightsList{nullptr},
    fqCutsList{nullptr},
    fHistqCuts{0},
    fWeights(0),
    fHistPhiEtaVz(0),
    fHistPhi(0),
    fHistEta(0),
    fHistPt(0),
    fHistZVertex(0),
    fProfcn_2gap{0},
    fProfcn_2gapINC{0},
    fProfPTdn_2gap{0},
    fHistqn_reduced{0},
    fProfPTdn_2gap_B{0},
    fh2Weights{nullptr},
    fHistPDG{0},
    fProfcn_2gap_qn{0},
    fProfdn_2gap_qn{0},
    fProfcn_2gap_q2{0},
    fProfcn_4gap_qn{0},
    fProfdn_4gap_qn{0},    
    fProfcn_4gap(0),
    fProfPTdn_4gap{0},
    fHistq2Cent(0),
    fHistq3Cent(0),
    fProfNPar{0},
    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fAbsEtaMax(0.8),
    fCentEstimator(),
    fReadMC{kFALSE},
    fMCEvent{0},
    fCentInterval{0,5,10,20,30,40,50,60,70,80,90,100},
    fFlowRFPsPtMin(0.2),
    fFlowRFPsPtMax(5.0),
    fFlowPOIsPtMin(0.0),
    fFlowPOIsPtMax(10.0),
    fnTwoCorr(kFALSE),
    fnFourCorr(kFALSE)
{}
//_____________________________________________________________________________
AliAnalysisTaskESEFlow::AliAnalysisTaskESEFlow(const char* name) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fFlowRunByRunWeights(kTRUE),
    bUseOwnWeights(0),
    dGap(0.5),
    fInit{0},
    fqRun{0},
    fAOD(0),
    fOutputList(0),
    fObservables(0),
    fCorrDist(0),
    fpTDiff(0),
    fqnDist(0),
    fpTDiffqselec(0),
    fcnqselec(0),
    fQAEvents(0),
    fFlowWeightsList{nullptr},
    fqCutsList{nullptr},
    fHistqCuts{0},
    fWeights(0),
    fHistPhiEtaVz(0),
    fHistPhi(0),
    fHistEta(0),
    fHistPt(0),
    fHistZVertex(0),
    fProfcn_2gap{0},
    fProfcn_2gapINC{0},
    fProfPTdn_2gap{0},
    fHistqn_reduced{0},
    fProfPTdn_2gap_B{0},
    fh2Weights{nullptr},
    fHistPDG{0},
    fProfcn_2gap_qn{0},
    fProfdn_2gap_qn{0},
    fProfcn_2gap_q2{0},
    fProfcn_4gap_qn{0},
    fProfdn_4gap_qn{0},    
    fProfcn_4gap(0),
    fProfPTdn_4gap{0},
    fHistq2Cent(0),
    fHistq3Cent(0),
    fProfNPar{0},
    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fAbsEtaMax(0.8),
    fCentEstimator(),
    fReadMC{kFALSE},
    fMCEvent{0},
    fCentInterval{0,5,10,20,30,40,50,60,70,80,90,100},
    fFlowRFPsPtMin(0.2),
    fFlowRFPsPtMax(5.0),
    fFlowPOIsPtMin(0.0),
    fFlowPOIsPtMax(10.0),
    fnTwoCorr(kFALSE),
    fnFourCorr(kFALSE)
{
    //define input and output
    DefineInput(0, TChain::Class());
    DefineInput(1, TList::Class());
    DefineInput(2, TList::Class());

    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
    DefineOutput(4, TList::Class());
    DefineOutput(5, TList::Class());
    DefineOutput(6, TList::Class());
    DefineOutput(7, TList::Class());
    DefineOutput(8, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskESEFlow::~AliAnalysisTaskESEFlow()
{
    if(fOutputList) { delete fOutputList; }
    if(fObservables) { delete fObservables; }
    if(fCorrDist) { delete fCorrDist; }
    if(fpTDiff) { delete fpTDiff; }
    if(fqnDist) { delete fqnDist; }
    if(fpTDiffqselec) { delete fpTDiffqselec; }
    if(fcnqselec) { delete fcnqselec; }
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

    fqCutsList = static_cast<TList*>(GetInputData(2));
    if(!fqCutsList) { AliFatal("\n \n \n \n \n \n \n \n \n \n \n \n q-selection cuts file not found! Terminating! \n \n \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }

    if(!LoadqCuts()) { AliFatal("\n \n \n \n \n \n \n \n \n \n q selection cuts not loaded! \n \n \n \n \n \n \n \n \n \n "); return kFALSE; }

    //printf("\n \n \n \n");
    //for(Int_t i(0); i<11;++i){
    //    std::cout << GetqSelectionCut(2,0,i) << std::endl;
    //}
    //printf("\n \n \n \n");
    
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
    fpTDiffqselec = new TList();
    fcnqselec = new TList();
    fQAEvents = new TList();

    fOutputList->SetOwner(kTRUE);
    fObservables->SetOwner(kTRUE);
    fCorrDist->SetOwner(kTRUE);
    fpTDiff->SetOwner(kTRUE);
    fqnDist->SetOwner(kTRUE);
    fpTDiffqselec->SetOwner(kTRUE);
    fcnqselec->SetOwner(kTRUE);
    fQAEvents->SetOwner(kTRUE);

    //RUN INITIALIZE TASK
    fInit = InitializeTask();
    if(!fInit) { return; }


    const int NvnPtBin = 28;
    double PtEdgesvn[NvnPtBin+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,7.0,8.0,9.0,10.0};

    

    fHistPhiEtaVz = new TH3F("fHistPhiEtaVz", "fHistPhiEtaVz; phi; eta; Vz", 120, 0.0, TMath::TwoPi(), 120, -1.0, 1.0,100,-15,15);
    fHistPhi = new TH1F("fHistPhi", "fHistPhi", 120, 0.0, TMath::TwoPi());
    fHistEta = new TH1F("fHistEta", "fHistEta", 120,-1.0, 1.0);
    fHistPt = new TH1F("fHistPt", "fHistPt", NvnPtBin,PtEdgesvn);
    fHistZVertex = new TH1F("fHistZVertex", "fHistZVertex", 100,-15,15);
    fProfNPar = new TProfile("fProfNparvsCent",";Centrality;N_{particles}",100,0,100);

    fHistq2Cent = new TH2D("Histq2vCent","",100,0,100,100,0,10);
    fHistq2Cent->Sumw2();
    fHistq3Cent = new TH2D("Histq3vCent","",100,0,100,100,0,10);
    fHistq3Cent->Sumw2();

    const int nBins = 9;
    double binedge[nBins+1] = {0, 5., 10., 20., 30., 40., 50., 60., 70., 80.};


    if(fnFourCorr){
        fProfcn_4gap = new TProfile("c_{2}{4} gap vs cent","",nBins,binedge);
        fProfcn_4gap->Sumw2();
        fCorrDist->Add(fProfcn_4gap);

        for(int fCentNum(0);fCentNum<fNumCentHists;++fCentNum){
            fProfPTdn_4gap[fCentNum] = new TProfile(Form("d_{2}{4} in centrality percentile: %.0f-%.0f",binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin,PtEdgesvn);
            fProfPTdn_4gap[fCentNum]->Sumw2();
            fpTDiff->Add(fProfPTdn_4gap[fCentNum]);

            for(Int_t qn(0); qn<2; ++qn){
            for(Int_t qBin(0); qBin<10; ++qBin){
                fProfcn_4gap_qn[qn][fCentNum][qBin] = new TProfile(Form("c2_4gap_q%isel_point_%i_in_%.0f_%.0f",qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",nBins, binedge);
                fProfcn_4gap_qn[qn][fCentNum][qBin]->Sumw2();
                fcnqselec->Add(fProfcn_4gap_qn[qn][fCentNum][qBin]);
                
                fProfdn_4gap_qn[qn][fCentNum][qBin] = new TProfile(Form("d2_4gap_q%isel_point_%i_in_%.0f_%.0f",qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                fProfdn_4gap_qn[qn][fCentNum][qBin]->Sumw2();
                fpTDiffqselec->Add(fProfdn_4gap_qn[qn][fCentNum][qBin]);
            }
            }
        }
    }

    //Differential PT v's
    if(fnTwoCorr){
    for(Int_t fHarmNum(0); fHarmNum<fNumHarmHists; ++fHarmNum)
    {
        fProfcn_2gap[fHarmNum] = new TProfile(Form("c_{%i}{2} gap vs cent",fHarmNum+2),"",nBins,binedge);
        fProfcn_2gap[fHarmNum]->Sumw2();
        fCorrDist->Add(fProfcn_2gap[fHarmNum]);

        fProfcn_2gapINC[fHarmNum] = new TProfile(Form("c_{%i}{2,gap}_vs_npart",fHarmNum+2),"",50,0,400);
        fProfcn_2gapINC[fHarmNum]->Sumw2();
        fOutputList->Add(fProfcn_2gapINC[fHarmNum]);

        // test alternative q-sel method
        fProfcn_2gap_q2[fHarmNum] = new TProfile2D(Form("c_%i_2gap_q2sel",fHarmNum+2),"",nBins,binedge,10,0,10);
        fProfcn_2gap_q2[fHarmNum]->Sumw2();
        fOutputList->Add(fProfcn_2gap_q2[fHarmNum]);

        for(Int_t fCentNum(0) ; fCentNum<fNumCentHists; ++fCentNum)
        {
            fProfPTdn_2gap[fHarmNum][fCentNum] = new TProfile(Form("d_{%i}{2} in centrality percentile: %.0f-%.0f",fHarmNum+2,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin,PtEdgesvn);
            fProfPTdn_2gap[fHarmNum][fCentNum]->Sumw2();
            fpTDiff->Add(fProfPTdn_2gap[fHarmNum][fCentNum]);
        
            fProfPTdn_2gap_B[fHarmNum][fCentNum] = new TProfile(Form("d_{%i}{2} B in centrality percentile: %.0f-%.0f",fHarmNum+2,binedge[fCentNum],binedge[fCentNum+1]),"", NvnPtBin, PtEdgesvn);
            fProfPTdn_2gap_B[fHarmNum][fCentNum]->Sumw2();
            fpTDiff->Add(fProfPTdn_2gap_B[fHarmNum][fCentNum]);

            //ESE
            for(Int_t qn(0); qn<2; ++qn){
            for(Int_t qBin(0); qBin<10; ++qBin){
                fProfcn_2gap_qn[qn][fHarmNum][fCentNum][qBin] = new TProfile(Form("c%i_2gap_q%isel_point_%i_in_%.0f_%.0f",fHarmNum+2,qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",nBins, binedge);
                fProfcn_2gap_qn[qn][fHarmNum][fCentNum][qBin]->Sumw2();
                fcnqselec->Add(fProfcn_2gap_qn[qn][fHarmNum][fCentNum][qBin]);

                fProfdn_2gap_qn[qn][fHarmNum][fCentNum][qBin] = new TProfile(Form("d%i_2gap_q%isel_point_%i_in_%.0f_%.0f",fHarmNum+2,qn+2,qBin+1,binedge[fCentNum],binedge[fCentNum+1]),"",NvnPtBin, PtEdgesvn);
                fProfdn_2gap_qn[qn][fHarmNum][fCentNum][qBin]->Sumw2();
                fpTDiffqselec->Add(fProfdn_2gap_qn[qn][fHarmNum][fCentNum][qBin]);
            }
            }
        }
    }
    }

    for(Int_t fHarmNum(0); fHarmNum<fNumHarmHists; ++fHarmNum){
        for(Int_t fCentNum(0); fCentNum<fNumCentHists; ++fCentNum){
            fHistqn_reduced[fHarmNum][fCentNum] = new TH1F(Form("q_{%i} red in centrality percentile: %.0f-%.0f",fHarmNum+2,binedge[fCentNum],binedge[fCentNum+1]),"",100,0,10);
            fqnDist->Add(fHistqn_reduced[fHarmNum][fCentNum]);
        }
    }

    if(fReadMC)
    {
        fHistPDG = new TH1F("fHistPDG","",500,0,500);
        fOutputList->Add(fHistPDG);
    }

    fEventCuts.AddQAplotsToList(fQAEvents); //QA plots
    
    fObservables->Add(fHistPhiEtaVz);
    fObservables->Add(fHistPhi);
    fObservables->Add(fHistEta);
    fObservables->Add(fHistPt);
    fObservables->Add(fHistZVertex);
    fObservables->Add(fProfNPar);
    fqnDist->Add(fHistq2Cent);
    fqnDist->Add(fHistq3Cent);

    PostData(1, fOutputList);
    PostData(2, fObservables);
    PostData(3, fCorrDist);
    PostData(4, fpTDiff);
    PostData(5, fqnDist);
    PostData(6, fpTDiffqselec);
    PostData(7, fcnqselec);
    PostData(8, fQAEvents);
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

    FillObsDistributions(iTracks, fAOD, dVz);

    CorrelationTask(centrality, iTracks, fAOD, dVz);
    

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
    PostData(6, fpTDiffqselec);
    PostData(7, fcnqselec);
    PostData(8, fQAEvents);
}
//_____________________________________________________________________________
void AliAnalysisTaskESEFlow::Terminate(Option_t *)
{
}
void AliAnalysisTaskESEFlow::CorrelationTask(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz)
{
    ReducedRFPVectors(centrality, iTracks, fAOD, dVz);
    fProfNPar->Fill(centrality,iTracks);
    if(!fqRun){
    RFPVectors(centrality, iTracks, fAOD, dVz);
    POIVectors(centrality, iTracks, fAOD, dVz);
    }
}
void AliAnalysisTaskESEFlow::FillObsDistributions(const Int_t iTracks, const AliAODEvent* fAOD, const float dVz)
{
    if(iTracks < 1 ) { return; }
    fHistZVertex->Fill(dVz);
            for(Int_t i(0); i < iTracks; i++) 
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
            }
}
void AliAnalysisTaskESEFlow::RFPVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz)
{
    ResetFlowVector(Qvector);
    ResetFlowVector(Qvector10P);
    ResetFlowVector(Qvector10M);


    if(iTracks < 1 ) { return; }
            for(Int_t i(0); i < iTracks; i++) 
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
                        Qvector[iHarm][iPower] += TComplex(dCos,dSin);
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
                        Qvector10P[iHarm][iPower] += TComplex(dCos,dSin);
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
                        Qvector10M[iHarm][iPower] += TComplex(dCos,dSin);
                        }
                    }
                }
            
            if(fnTwoCorr){
            FillRFP(centrality,iTracks,2,2);  //with gap for 2-particle correlation nHarm=2
            FillRFP(centrality,iTracks,3,2);  // nHarm=3
            FillRFP(centrality,iTracks,4,2);  // nHarm=4
            FillRFP(centrality,iTracks,5,2);
            FillRFP(centrality,iTracks,6,2);
            }
            if(fnFourCorr){
            FillRFP(centrality,iTracks,2,4);  //with gap for 4-particle correlation nHarm=2
            }
            
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
                            qvector[iHarm][iPower] += TComplex(dCos,dSin);
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
                        pvector[iHarm][iPower] += TComplex(dCos,dSin);
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
                        pvector10P[iHarm][iPower] += TComplex(dCos,dSin);
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
                        pvector10M[iHarm][iPower] += TComplex(dCos,dSin);
                        }
                    }
                }     
            }
        }
}
void AliAnalysisTaskESEFlow::POIVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz)
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
        Filldn(centrality,dPtL,2,2); //Fill d2_2 <<2'>>  fnParCorr-particle correlation (std 2)
        Filldn(centrality,dPtL,3,2); 
        Filldn(centrality,dPtL,4,2); 
        Filldn(centrality,dPtL,5,2); 
        Filldn(centrality,dPtL,6,2); 
        }
        if(fnFourCorr){
        Filldn(centrality,dPtL,2,4);
        }
    }
}
void AliAnalysisTaskESEFlow::ReducedRFPVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz)
{
    ResetReducedFlowVector(qvector_red);
    ResetReducedFlowVector(sumCos);
    ResetReducedFlowVector(sumSin);

    Double_t M = 0;
    if(iTracks < 1 ) { return; }
    for(Int_t i(0); i < iTracks; i++) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        if(!WithinRFP(track)){continue;} // check if also within reference particle

        double dEta = track->Eta();
        double dPhi = track->Phi();
        //double dPt = track->Pt();

        Double_t dWeight = GetFlowWeight(track,dVz);

        if(- (dGap-0.1) < dEta && dEta < (dGap-0.1))
        {
            M += 1;
            for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
            { 
                for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
                {
                
                Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                sumCos[iHarm][iPower] += dCos;
                sumSin[iHarm][iPower] += dSin;
                }
            }
        }
    }
    if(M!=0)
    {
        for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm)
        { 
            for(Int_t iPower(0); iPower < fNumPowers; ++iPower)
            {
            Double_t dqn = sumCos[iHarm][iPower]*sumCos[iHarm][iPower]+sumSin[iHarm][iPower]*sumSin[iHarm][iPower];
            qvector_red[iHarm][iPower] = TMath::Sqrt(dqn)/TMath::Sqrt(M);
            }
        }
        Fillqnreduced(centrality);
    }
}
void AliAnalysisTaskESEFlow::FillRFP(const Float_t centrality,const Int_t iTracks,const int nHarm, const int nCorr)
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
            fProfcn_2gap[bhistN]->Fill(centrality,cn_fill,cDenom);
            fProfcn_2gapINC[bhistN]->Fill(iTracks,cn_fill,cDenom);

            FillIntegratecqnCut(centrality, nHarm, cn_fill, cDenom, 2,2); //fill integrated cqncut here
            FillIntegratecqnCut(centrality, nHarm, cn_fill, cDenom, 3,2); //fill integrated cqncut here
        }
        if(nCorr==4 && nHarm==2){
            fProfcn_4gap->Fill(centrality,cn_fill,cDenom);

            FillIntegratecqnCut(centrality, nHarm, cn_fill, cDenom, 2,4); //fill integrated cqncut here
            FillIntegratecqnCut(centrality, nHarm, cn_fill, cDenom, 3,4); //fill integrated cqncut here
        }
    }
    }
}
void AliAnalysisTaskESEFlow::FillIntegratecqnCut(const Float_t centrality, const int nHarm, const double c, const double c_weight, const int q_i,const int nCorr)
{
    int nHist = nHarm-2;
    int q_iHist = q_i-2;

    for(Int_t nCentHist(0); nCentHist<fNumCentHists; ++nCentHist){
        for(Int_t iCut(0); iCut<10; ++iCut){
            if(qvector_red[q_i][1] > GetqSelectionCut(q_i,nCentHist,iCut) && qvector_red[q_i][1] < GetqSelectionCut(q_i,nCentHist,iCut+1)){
                if(nCorr==2) {fProfcn_2gap_qn[q_iHist][nHist][nCentHist][iCut]->Fill(centrality,c,c_weight);}
                if(nCorr==4) {fProfcn_4gap_qn[q_iHist][nCentHist][iCut]->Fill(centrality,c,c_weight);}
            }
        }
    }

    if(q_i==2){
    for(Int_t nCentHist(0); nCentHist<fNumCentHists;++nCentHist){
        if(centrality > fCentInterval[nCentHist] && centrality <= fCentInterval[nCentHist+1]){
            Double_t CentBin = fProfcn_2gap_q2[nHist]->GetXaxis()->GetBinCenter(nCentHist);
            for(Int_t nBin(1); nBin<11; ++nBin){
                if(qvector_red[q_i][1] > GetqSelectionCut(q_i,nCentHist,nBin-1) && qvector_red[q_i][1] <= GetqSelectionCut(q_i,nCentHist,nBin)){
                    fProfcn_2gap_q2[nHist]->Fill(CentBin,nBin,c,c_weight);
                }
            }
        }
    }
    }
}
void AliAnalysisTaskESEFlow::FillDifferentialdqnCut(const Float_t centrality, const int nHarm, const Double_t dPt, const double d, const double d_weight, const int q_i,const int nCorr)
{
    int nHist = nHarm-2;
    int q_iHist = q_i-2;

    for(Int_t nCentHist(0); nCentHist<fNumCentHists; ++nCentHist){
        if(centrality > fCentInterval[nCentHist] && centrality < fCentInterval[nCentHist+1]){
            for(Int_t iCut(0); iCut<10; ++iCut){
                if(qvector_red[q_i][1] > GetqSelectionCut(q_i,nCentHist,iCut) && qvector_red[q_i][1] < GetqSelectionCut(q_i,nCentHist,iCut+1)){
                    if(nCorr==2) {fProfdn_2gap_qn[q_iHist][nHist][nCentHist][iCut]->Fill(dPt,d,d_weight);}
                    if(nCorr==4) {fProfdn_4gap_qn[q_iHist][nCentHist][iCut]->Fill(dPt,d,d_weight);}
                }
            }
        }
    }
}
void AliAnalysisTaskESEFlow::Filldn(const Float_t centrality, const double dPt, const int nHarm, const int nCorr)
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
            for(int i(0); i<fNumCentHists; ++i)
            {
                if(centrality > fCentInterval[i] && centrality < fCentInterval[i+1])
                {
                    if(nCorr==2){
                    fProfPTdn_2gap[bhistnumb][i]->Fill(dPt,dn_pt,dDenom);
                    }
                    if(nCorr==4){
                        if(nHarm==2){
                            fProfPTdn_4gap[i]->Fill(dPt,dn_pt,dDenom);
                        }
                    }
                }
            }

            if(nCorr==2){
            FillDifferentialdqnCut(centrality, nHarm, dPt, dn_pt, dDenom, 2,nCorr);
            FillDifferentialdqnCut(centrality, nHarm, dPt, dn_pt, dDenom, 3,nCorr);
            }
            if(nCorr==4){
            FillDifferentialdqnCut(centrality, nHarm, dPt, dn_pt, dDenom, 2,nCorr);
            FillDifferentialdqnCut(centrality, nHarm, dPt, dn_pt, dDenom, 3,nCorr);
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
            for(int i(0); i<fNumCentHists; ++i)
            {
                if(centrality > fCentInterval[i] && centrality < fCentInterval[i+1])
                {
                    fProfPTdn_2gap_B[bhistnumb][i]->Fill(dPt,dn_2ptB,dn2pt_wB);
                }
            }
        }
    }
    }
}
void AliAnalysisTaskESEFlow::Fillqnreduced(const Float_t centrality)
{
    for(int j(0); j<fNumHarmHists; ++j)
    {
        for(int i(0); i<fNumCentHists; ++i)
        {
            if(centrality > fCentInterval[i] && centrality < fCentInterval[i+1])
            {
                fHistqn_reduced[j][i]->Fill(qvector_red[j+2][1]);              // filling n=j+2 since j starting from 0
            }
        }
    }
    fHistq2Cent->Fill(centrality,qvector_red[2][1]);
    fHistq3Cent->Fill(centrality,qvector_red[3][1]);
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
Bool_t AliAnalysisTaskESEFlow::LoadqCuts()
{
    for (int i(0); i<2;++i){
        for (int j(0); j < fNumCentHists; ++j){
            fHistqCuts[i][j] = static_cast<TH1F*>(fqCutsList->FindObject(Form("q%i_%.0f_%.0f",i+2,fCentInterval[j],fCentInterval[j+1])));
        }    
    }




    return kTRUE;
}
Double_t AliAnalysisTaskESEFlow::GetqSelectionCut(Int_t nHarm, Int_t CentRange, Int_t Entry)
{
    // CentRange: 
    // 0: 0-5
    // 1: 5-10
    // 2: 10-20
    // 3: 20-30
    // 4: 30-40
    // 5: 40-50
    // 6: 50-60
    Double_t temp;
    temp = fHistqCuts[nHarm-2][CentRange]->GetBinContent(Entry);

    return temp;
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
void AliAnalysisTaskESEFlow::ResetReducedFlowVector(double (&array)[fNumHarms][fNumPowers])
{
  // RESET Reduced RFP vector
  // *************************************************************
  for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm) {
    for(Int_t iPower(0); iPower < fNumPowers; ++iPower) {
      array[iHarm][iPower]=(0.0);
    }
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