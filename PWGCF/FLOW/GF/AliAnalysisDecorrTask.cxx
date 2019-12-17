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



#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliGFWWeights.h"
#include <iostream>
#include "TComplex.h"
#include "TRandom3.h"
#include "AliUniFlowCorrTask.h"

#include "AliAnalysisDecorrTask.h"


class AliAnalysisDecorrTask;


using namespace std;
using namespace TMath;


ClassImp(AliAnalysisDecorrTask)

AliAnalysisDecorrTask::AliAnalysisDecorrTask() : AliAnalysisTaskSE(),
    fEventCuts(),
    fFlowList{nullptr},
    fFlowWeights{nullptr},
    fWeights(0),
    fWeightList{nullptr},
    fh2Weights(nullptr),
    fh3Weights(nullptr),

    fIndexSampling{0},
    fAOD(nullptr),
    fInitTask{kFALSE},   
    fVecCorrTask{},

    fSampling{kFALSE},
    fFillQA(kFALSE),

    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fCentEstimator("V0M"), 
    fFilterBit(96),
    fPtAxis(new TAxis()),
    fCentAxis(new TAxis()),
    fCentMin{0.0},
    fCentMax{100.0},
    fPVtxCutZ{10.0},

    fCutChargedTrackFilterBit{96},
    fCutNumTPCclsMin{70},
    fCutDCAzMax{0.0},
    fCutDCAxyMax{0.0},
    bUseLikeSign(kFALSE),
    iSign(0),

    fAbsEtaMax(0.8),
    dEtaGap(1.0),
    fEtaBinNum{0},
    fPhiBinNum{60},
    fUseWeights3D(kTRUE),
    fFillWeights(kFALSE),
    fNumSamples{1},

    bHasGap(kTRUE),
    bDiff(kFALSE),
    bRef(kTRUE),
    bPtB(kFALSE),
    bHigherOrder(kFALSE),

    fPOIsPtmax(10.0),
    fPOIsPtmin(0.2),
    fRFPsPtMax(5.0),
    fRFPsPtMin(0.2)
{}
//_____________________________________________________________________________
AliAnalysisDecorrTask::AliAnalysisDecorrTask(const char* name) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fFlowList{nullptr},
    fFlowWeights{nullptr},
    fWeights(0),
    fWeightList{nullptr},
    fh2Weights(nullptr),
    fh3Weights(nullptr),

    fIndexSampling{0},
    fAOD(nullptr),
    fInitTask{kFALSE},   
    fVecCorrTask{},

    fSampling{kFALSE},
    fFillQA(kFALSE),

    fTrigger(AliVEvent::kINT7),
    fEventRejectAddPileUp(kFALSE),
    fCentEstimator("V0M"), 
    fFilterBit(96),
    fPtAxis(new TAxis()),
    fCentAxis(new TAxis()),
    fCentMin{0.0},
    fCentMax{100.0},
    fPVtxCutZ{10.0},

    fCutChargedTrackFilterBit{96},
    fCutNumTPCclsMin{70},
    fCutDCAzMax{0.0},
    fCutDCAxyMax{0.0},
    bUseLikeSign(kFALSE),
    iSign(0),

    fAbsEtaMax(0.8),
    dEtaGap(1.0),
    fEtaBinNum{0},
    fPhiBinNum{60},
    fUseWeights3D(kTRUE),
    fFillWeights(kFALSE),
    fNumSamples{1},

    bHasGap(kTRUE),
    bDiff(kFALSE),
    bRef(kTRUE),
    bPtB(kFALSE),
    bHigherOrder(kFALSE),

    fPOIsPtmax(10.0),
    fPOIsPtmin(0.2),
    fRFPsPtMax(5.0),
    fRFPsPtMin(0.2)
{
    DefineInput(0, TChain::Class());
    DefineInput(1, TList::Class()); 
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisDecorrTask::~AliAnalysisDecorrTask()
{
    
    if(fFlowList) delete fFlowList;
    if(fFlowWeights) delete fFlowWeights;
}

Bool_t AliAnalysisDecorrTask::InitTask()
{
    if(fFillWeights) { AliWarning("\n\n\n\n\n Filling weights. No analysis will be produced!\n\n\n\n\n"); return kTRUE; }

    if(fEtaBinNum < 1)
    {
        AliWarning("fEtaBinNum not set. Setting automatically");
        fEtaBinNum = 2.0*fAbsEtaMax/0.05;
    }
    if(fPhiBinNum < 1)
    {
        AliFatal("fPhiBin wrong! Terminating!");
        return kFALSE;
    }

    if (!fUseWeights3D) { 
        fWeightList = (TList*) GetInputData(1);
        if(!fWeightList) { AliFatal("\n\n\n\n\n\n\n\n\n Weight List not found! \n\n\n\n\n\n\n\n"); return kFALSE; }
    } else {
        fWeightList = (TList*) GetInputData(1);
        if(!fWeightList) { AliFatal("\n\n\n\n\n\n\n\n\n Weight List not found! \n\n\n\n\n\n\n\n"); return kFALSE; }
    }

    AliInfo("Weight List loaded");

    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisDecorrTask::UserCreateOutputObjects()
{

    fFlowList = new TList();
    fFlowList->SetOwner(kTRUE);
    fFlowWeights = new TList();
    fFlowWeights->SetOwner(kTRUE);

    fInitTask = InitTask();
    if(!fInitTask) { return; }
    
    NcentBin = fCentAxis->GetNbins();
    NPtBin = fPtAxis->GetNbins();
    for(Int_t bin(0); bin < NPtBin+1; bin++)
    {
        PtEdges[bin] = fPtAxis->GetBinLowEdge(bin+1);
    }
    for(Int_t bin(0); bin < NcentBin+1; bin++)
    {
        centEdges[bin] = fCentAxis->GetBinLowEdge(bin+1);
    }

    Int_t iNumTasks = fVecCorrTask.size();

    if(iNumTasks < 1) { return; }

    for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
    {
        AliUniFlowCorrTask* task = fVecCorrTask.at(iTask);
        if(!task) { fInitTask = kFALSE; AliError(Form("AliUniFlowCorrTask%d does not exist",iTask)); return; }

        Bool_t bHasGap = task->HasGap();
        Int_t CorrOrder = task->fiNumHarm;
        Bool_t bHarmSign = kFALSE;
        if(task->fiHarm[0] == task->fiHarm[1]) { bHarmSign = kTRUE; }
        const char* CorrName = task->fsName.Data();
        const char* CorrLabel = task->fsLabel.Data();

        for(Int_t iSample(0); iSample < fNumSamples; ++iSample)
        {
            if(iSample > 0 && !fSampling) { break; }
            TH1* profile = nullptr;
            TH1* profDiff = nullptr;
            TH1* profPtA = nullptr;
            TH1* profPtAPtB = nullptr;

            if(bRef && CorrOrder < 4)
            {
                profile = new TProfile(Form("%s_sample%d",CorrName,iSample),Form("%s",CorrLabel),NcentBin,centEdges);

                if(!profile) { fInitTask = kFALSE; AliError("Centrality profile not created"); task->PrintTask(); return; }
                if(fFlowList->FindObject(profile->GetName())) {
                    AliError(Form("Task %d: Profile '%s' already exists",iTask,profile->GetName()));
                    fInitTask=kFALSE;
                    task->PrintTask();
                    delete profile;
                    return;
                }

                profile->Sumw2();
                fFlowList->Add(profile);
            }

            if(bDiff) 
            {
                if( (bHigherOrder && bHarmSign) || CorrOrder < 4)
                {
                    profDiff = new TProfile2D(Form("%s_diff_sample%d",CorrName,iSample),Form("%s_diff",CorrLabel),NcentBin,centEdges,NPtBin,PtEdges);
                    if(!profDiff) { fInitTask = kFALSE; AliError("Differential profile not created"); task->PrintTask(); return; }
                    if(fFlowList->FindObject(profDiff->GetName())) {
                        AliError(Form("Task %d: Profile '%s' already exists",iTask,profDiff->GetName()));
                        fInitTask=kFALSE;
                        task->PrintTask();
                        delete profDiff;
                        return;
                    }

                profDiff->Sumw2();
                fFlowList->Add(profDiff);
                }

                if(CorrOrder < 4)
                {
                    profPtA = new TProfile2D(Form("%s_PtA_sample%d",CorrName,iSample),Form("%s_PtA",CorrLabel),NcentBin,centEdges,NPtBin,PtEdges);
                    if(!profPtA) { fInitTask = kFALSE; AliError("\n\n\nPtA profile not created\n\n\n"); task->PrintTask(); return; }
                    if(fFlowList->FindObject(profPtA->GetName())) {
                        AliError(Form("Task %d: Profile '%s' already exists",iTask,profPtA->GetName()));
                        fInitTask=kFALSE;
                        task->PrintTask();
                        delete profPtA;
                        return;
                    }

                    profPtA->Sumw2();
                    fFlowList->Add(profPtA);
                }
            }
            if(bPtB)
            { 
                profPtAPtB = new TProfile3D(Form("%s_PtAPtB_sample%d",CorrName,iSample),Form("%s_PtAPtB",CorrLabel),NcentBin,centEdges, NPtBin,PtEdges, NPtBin, PtEdges);
                if(!profPtAPtB) { fInitTask = kFALSE; AliError("PtAPtB profile not created"); task->PrintTask(); return; }
                if(fFlowList->FindObject(profPtAPtB->GetName())) {
                    AliError(Form("Task %d: Profile '%s' already exists",iTask,profPtAPtB->GetName()));
                    fInitTask=kFALSE;
                    task->PrintTask();
                    delete profPtAPtB;
                    return;
                }

                profPtAPtB->Sumw2();
                fFlowList->Add(profPtAPtB);
            }

        } //End for iSample
    } //End for iTask

    if(fFillWeights)
    {
        if(fUseWeights3D)
        {
            const char* weightName = Form("fh3Weights");
            const char* weightLabel = Form("3D weights; #varphi; #eta; PVtxZ");
            fh3Weights = new TH3D(weightName, weightLabel, fPhiBinNum, 0, TMath::TwoPi(), fEtaBinNum, -fAbsEtaMax, fAbsEtaMax, 2*fPVtxCutZ, -fPVtxCutZ, fPVtxCutZ);
            fh3Weights->Sumw2();
            fFlowWeights->Add(fh3Weights);
        }
        else
        {
            const char* weightName = Form("fh2Weights");
            const char* weightLabel = Form("2D weights; #varphi; #eta");
            fh2Weights = new TH2D(weightName, weightLabel, fPhiBinNum, 0, TMath::TwoPi(), fEtaBinNum, -fAbsEtaMax, fAbsEtaMax);    
            fh2Weights->Sumw2();
            fFlowWeights->Add(fh2Weights);
        } 
    }

    

    PostData(1, fFlowList);
    PostData(2, fFlowWeights);

}

Bool_t AliAnalysisDecorrTask::LoadWeights()
{
    if(!fUseWeights3D)
    {
        TList* listFlowWeights = nullptr;
       
        listFlowWeights = (TList*) fWeightList->FindObject(Form("%d", fAOD->GetRunNumber()));

        if(!listFlowWeights) 
        {
            // run-specific weights not found for this run; loading run-averaged instead
            AliWarning(Form("TList with flow weights (run %d) not found. Using run-averaged weights instead (as a back-up)", fAOD->GetRunNumber()));
            listFlowWeights = (TList*) fWeightList->FindObject("averaged");
            if(!listFlowWeights) { AliError("Loading run-averaged weights failed!"); fWeightList->ls(); return kFALSE; }
        }
        
        fh2Weights = (TH2D*) listFlowWeights->FindObject("Refs");
        return kTRUE;
    }
    else 
    {
        fWeights = (AliGFWWeights*)fWeightList->FindObject(Form("w%i",fAOD->GetRunNumber()));
        if(!fWeights)
        {
            printf("Weights could not be found in list!\n");
            return kFALSE;
        }
        fWeights->CreateNUA();
        return kTRUE;
    }
}

Bool_t AliAnalysisDecorrTask::FillWeights()
{
    int iPart(fAOD->GetNumberOfTracks());
    if(iPart < 1) { return kFALSE; }
    double dVz = fAOD->GetPrimaryVertex()->GetZ();

    for(int index(0); index < iPart; ++index)
    {   
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(index));
        if(!track || !IsTrackSelected(track)) { continue; }

        //double dPt = POItrack->Pt(); Keep for efficiency corrections one day
        double dPhi = track->Phi();
        double dEta = track->Eta(); 

        if(fUseWeights3D) { fh3Weights->Fill(dPhi,dEta,dVz); }
        else { fh2Weights->Fill(dPhi,dEta); }
    }
    

    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisDecorrTask::UserExec(Option_t *)
{


    if(!fInitTask) { AliFatal("Something went wrong! Task not initialized"); return; }
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { return; }
    
    //Event selection
    if(!IsEventSelected()) { return; }

    if(fFillWeights) 
    { 
        FillWeights();
        PostData(1, fFlowList);
        PostData(2, fFlowWeights);
        return; 
    }

    if(!LoadWeights())
    {
        AliFatal("\n\n\n\n\n\n\n\n Weights could not be loaded \n\n\n\n\n\n\n\n");
        return;
    }
    //Get centrality of event
    Float_t centrality(0);
    AliMultSelection *multSelect =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(multSelect) centrality = multSelect->GetMultiplicityPercentile(fCentEstimator);

    fIndexSampling = GetSamplingIndex();

    Int_t iNumTask = fVecCorrTask.size();
    for(Int_t iTask(0); iTask < iNumTask; ++iTask)
    {
        //Multiplicity
        Int_t iTracks(fAOD->GetNumberOfTracks());
        //Fill RP vectors
        FillRPvectors(dEtaGap);
        Bool_t doRef = kFALSE;
        Bool_t doDiff = kFALSE;
        Bool_t doPtB = kFALSE;
        if(bRef) { doRef = kTRUE; }
        CalculateCorrelations(fVecCorrTask.at(iTask), centrality, -1.0, -1.0, doRef, doDiff, doPtB);
        doRef = kFALSE;
        //Fill POI vectors
        int iNumPtBins = fPtAxis->GetNbins();


        if (bDiff || bPtB) {
        //Loop over Pt bins
            for(int iPtA(1); iPtA < iNumPtBins+1; ++iPtA)
            {
                double dPt = fPtAxis->GetBinCenter(iPtA);
                double dPtLow = fPtAxis->GetBinLowEdge(iPtA);
                double dPtHigh = fPtAxis->GetBinUpEdge(iPtA);

                FillPOIvectors(dEtaGap, dPtLow, dPtHigh);
                if(bDiff) { doDiff = kTRUE; }
                CalculateCorrelations(fVecCorrTask.at(iTask), centrality, dPt, -1.0, doRef, doDiff, doPtB);
                doDiff = kFALSE;
            
                if(bPtB && dPt < 5.0)
                {
                    // Too slow  -- reimplement maybe
                    for(int iPtB(1); iPtB < iNumPtBins+1; ++iPtB)
                    { 
                        double dPtB = fPtAxis->GetBinCenter(iPtB);
                        double dPtBLow = fPtAxis->GetBinLowEdge(iPtB);
                        double dPtBHigh = fPtAxis->GetBinUpEdge(iPtB);
                        FillPtBvectors(dEtaGap, dPtBLow, dPtBHigh);
                        if(bPtB) { doPtB = kTRUE;}
                        CalculateCorrelations(fVecCorrTask.at(iTask), centrality, dPt, dPtB, doRef, doDiff, doPtB); 
                        doPtB = kFALSE;   
                    }
                }
            } 
        }
    }

    
    PostData(1, fFlowList);
    PostData(2, fFlowWeights);

}

void AliAnalysisDecorrTask::CalculateCorrelations(const AliUniFlowCorrTask* const task, double centrality, double dPtA, double dPtB, Bool_t doRef, Bool_t doDiff, Bool_t doPtB)
{
    if(!task) { AliError("AliUniFlowCorrTask does not exist"); return; }

    //Bool_t bHasGap = task->HasGap();
    Int_t corrOrder= task->fiNumHarm;
    Bool_t bHarmSign = kFALSE;
    if(task->fiHarm[0] == task->fiHarm[1]) { bHarmSign = kTRUE; }

    TComplex cNum = TComplex(0.0,0.0,kFALSE);
    TComplex cDn = TComplex(0.0,0.0,kFALSE);
    TComplex cNumDiff = TComplex(0.0,0.0,kFALSE);
    TComplex cDnDiff = TComplex(0.0,0.0,kFALSE);
    TComplex cNumPtB = TComplex(0.0,0.0,kFALSE);
    TComplex cDnPtB = TComplex(0.0,0.0,kFALSE);
    TComplex cNumPtA = TComplex(0.0,0.0,kFALSE);
    TComplex cDnPtA = TComplex(0.0,0.0,kFALSE);


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
                cDnPtA = TwoDiffGap10_Pt(0,0);
                cNumPtA = TwoDiffGap10_Pt(task->fiHarm[0],task->fiHarm[1]);
            }
            if(doPtB) {
                cDnPtB = TwoDiffGap10_PtA_PtB(0,0);
                cNumPtB = TwoDiffGap10_PtA_PtB(task->fiHarm[0],task->fiHarm[1]);
            }
            if(doRef) {
                cDn = TwoGap10(0,0);
                cNum = TwoGap10(task->fiHarm[0],task->fiHarm[1]);
            }
        }
        break;
    case 4 :
        if(!bHasGap){
            if(doDiff){
                cDnDiff = FourDiff(0,0,0,0);
                cNumDiff = FourDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
            }
            if(doPtB) {
                cDnPtB = FourDiff_PtA_PtB(0,0,0,0);
                cNumPtB = FourDiff_PtA_PtB(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
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
            if(doPtB) {
                if(task->fiHarm[1] > 0)         //if associate particle have same sign take associate from eta regions: M:AA and P:TT    (M = negative, P = positive, A = associate, T = trigger)
                {
                    cDnPtB = FourDiffGap10_PtA_PtB(0,0,0,0);
                    cNumPtB = FourDiffGap10_PtA_PtB(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
                }
                else                            //if associate particle have opposite sign take from eta regions: M:AT and P:AT
                {
                    cDnPtB = FourDiffGap10_OS_PtA_PtB(0,0,0,0);
                    cNumPtB = FourDiffGap10_OS_PtA_PtB(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
                }
            }
            if(doRef) {
                cDn = FourGap10(0,0,0,0);
                cNum = FourGap10(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
            }
        }
        break;
    default:
        return;
    }

    if(doRef && corrOrder < 4)
    {
        Double_t dDn = cDn.Re();
        Double_t dNum = cNum.Re();
        Double_t dValue = 0.0;
        Bool_t bFillPos = kFALSE;

        if(dDn > 0.0) {bFillPos = kTRUE; dValue = dNum/dDn; }
        if(bFillPos && TMath::Abs(dValue > 1.0)) { bFillPos = kFALSE; }
        if(!bFillPos) { return; }
        TProfile* prof = (TProfile*)fFlowList->FindObject(Form("%s_sample0",task->fsName.Data()));
        if(!prof) { AliError(Form("Profile %s_sample0 not found",task->fsName.Data())); return; }
        prof->Fill(centrality, dValue, dDn);

    }
    if(doDiff)
    {
        Double_t dDnDiff = cDnDiff.Re();
        Double_t dNumDiff = cNumDiff.Re();
        Double_t dValueDiff = 0.0;
        Bool_t bFillDiff = kFALSE;

        Double_t dDnPtA = cDnPtA.Re();
        Double_t dNumPtA = cNumPtA.Re();
        Double_t dValuePtA = 0.0;
        Bool_t bFillPtA = kFALSE;

        if(dDnDiff > 0.0) { bFillDiff = kTRUE; dValueDiff = dNumDiff/dDnDiff; }
        if(bFillDiff && TMath::Abs(dValueDiff > 1.0)) { bFillDiff = kFALSE; }
        if(dDnPtA > 0.0) { bFillPtA = kTRUE; dValuePtA = dNumPtA/dDnPtA; }
        if(bFillPtA && TMath::Abs(dValuePtA) > 1.0) { bFillPtA = kFALSE; }

        if(!bFillDiff && !bFillPtA) { return; }
        
        if(corrOrder < 4)
        {
            TProfile2D* profDiff = (TProfile2D*)fFlowList->FindObject(Form("%s_diff_sample0",task->fsName.Data()));
            TProfile2D* profPtA = (TProfile2D*)fFlowList->FindObject(Form("%s_PtA_sample0",task->fsName.Data()));

            if(!profDiff) { AliError(Form("Profile %s_diff_sample0 not found",task->fsName.Data())); return; }
            if(!profPtA) { AliError(Form("Profile_%s_PtA_sample0 not found",task->fsName.Data())); return; }
            profDiff->Fill(centrality, dPtA, dValueDiff, dDnDiff);
            profPtA->Fill(centrality, dPtA, dValuePtA, dDnPtA);
        }
        if(bHigherOrder && bHarmSign)
        {
            TProfile2D* profDiff = (TProfile2D*)fFlowList->FindObject(Form("%s_diff_sample%d",task->fsName.Data(),fIndexSampling));
            if(!profDiff) { AliError(Form("Profile %s_diff_sample%d not found",task->fsName.Data(),fIndexSampling)); return; }
            profDiff->Fill(centrality, dPtA, dValueDiff, dDnDiff);
        }        
    }
    if(doPtB)
    {
        Double_t dDnPtB = cDnPtB.Re();
        Double_t dNumPtB = cNumPtB.Re();
        Double_t dValuePtB = 0.0;
        Bool_t bFillPtB = kFALSE;

        if(dDnPtB > 0.0) { bFillPtB = kTRUE; dValuePtB = dNumPtB/dDnPtB; }
        if(bFillPtB && TMath::Abs(dValuePtB > 1.0)) { bFillPtB = kFALSE; }
        if(!bFillPtB) { return; }
        if(corrOrder < 4)
        {
            TProfile3D* profPtAPtB = (TProfile3D*)fFlowList->FindObject(Form("%s_PtAPtB_sample0",task->fsName.Data()));
            if(!profPtAPtB) { AliError(Form("Profile %s_PtAPtB_sample0 not found",task->fsName.Data())); }
            profPtAPtB->Fill(centrality,dPtA,dPtB, dValuePtB, dDnPtB);
        }
        else
        {
            TProfile3D* profPtAPtB = (TProfile3D*)fFlowList->FindObject(Form("%s_PtAPtB_sample%d",task->fsName.Data(),fIndexSampling));
            if(!profPtAPtB) { AliError(Form("Profile %s_PtAPtB_sample%d not found",task->fsName.Data(),fIndexSampling)); }
            profPtAPtB->Fill(centrality,dPtA,dPtB, dValuePtB, dDnPtB);
        }
    }

    return;
}

//Method to fill RP vectors and stat histograms
void AliAnalysisDecorrTask::FillRPvectors(double dEtaGap)
{

    ResetFlowVector(Qvector);
    ResetFlowVector(Qvector10P);
    ResetFlowVector(Qvector10M);

    int iTracks(fAOD->GetNumberOfTracks());
    if(iTracks < 1 ) { return; }

    bool bIsRP; 
    double dEtaLimit = 0.5*dEtaGap;
    double dVz = fAOD->GetPrimaryVertex()->GetZ();
    
    for(Int_t i(0); i < iTracks; i++) 
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        bIsRP = IsWithinRP(track);
        if (!bIsRP) { continue; }

        double dPhi = track->Phi();
        double dEta = track->Eta();
        double dPt = track->Pt();

        //Calculating weights    
        double dWeight = GetWeights(dPhi, dEta, dVz);
        if(dWeight <= 0.0) { dWeight = 1.0; }
        
        //Filling Q-vectors for RPs
            
        for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++) 
        {
            for(Int_t iPower(0); iPower < fNumPowers; iPower++)
            {
                Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                Qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            } //End for iPower
        }  //End for iHarm
        // RFP in positive and negative eta acceptance
        if(dEta > dEtaLimit)
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++)
            {
                for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    Qvector10P[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                }  //End for iPower
            }  //End for iHarm
        }
        else if(dEta < -dEtaLimit)
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++)
            {
                for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    Qvector10M[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                } //End for iPower
            } //End for iHarm
        }  //end if eta gap
    }  //end track loop

    return;
}

//Method to fill POIs into pvectors
void AliAnalysisDecorrTask::FillPOIvectors(const double dEtaGap, const double dPtLow, const double dPtHigh)
{

    ResetFlowVector(pvector);
    ResetFlowVector(pvector10M);
    ResetFlowVector(pvector10P);
    ResetFlowVector(qvector);


    int iPart(fAOD->GetNumberOfTracks());
    if(iPart < 1) { return; }
    double dEtaLimit = 0.5*dEtaGap;
    double dVz = fAOD->GetPrimaryVertex()->GetZ();

    for(int index(0); index < iPart; ++index)
    {   
        AliAODTrack* POItrack = static_cast<AliAODTrack*>(fAOD->GetTrack(index));
        if(!POItrack || !IsTrackSelected(POItrack)) { continue; }

        double dPt = POItrack->Pt();
        double dPhi = POItrack->Phi();
        double dEta = POItrack->Eta();

        //if (Abs(dEta) < dEtaLimit) { continue; } Why the hell is this here?

        //Check for overlap with RPs
        bool bIsWithinRP = IsWithinRP(POItrack);

        //Check that particle is in POI range
        bool bIsWithinPOI = IsWithinPOI(POItrack);

        if(!bIsWithinPOI) { continue; }

        //Load weights
        double dWeight = GetWeights(dPhi,dEta,dVz);
        if(dWeight <= 0.0) { dWeight = 1.0; }

        //POI with no eta gap
        if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++) 
            {
                for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    pvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    
                    //Check if there is overlap of POI and RP
                    if(bIsWithinRP)
                    {
                        qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                }  //End for iPower
            }  //End for iHarm

            //POI with eta gap
            if(dEta > dEtaLimit)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++)
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        pvector10P[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                        
                    } //End for iPower
                }  //End for iHarm
            }
            if (dEta < -dEtaLimit)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++)
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        pvector10M[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                        
                    } //End for iPower
                } //End for iHarm
            }  // end if eta gap
        } //end if dPtLow < dPt < dPtHigh
    } //end for track

    return;
}

//Method to fill POIs into pvectors
void AliAnalysisDecorrTask::FillPtBvectors(const double dEtaGap, const double dPtLow, const double dPtHigh)
{

    ResetFlowVector(pvectorPtB);
    ResetFlowVector(pvectorPtB10M);
    ResetFlowVector(pvectorPtB10P);
    ResetFlowVector(qvectorPtB);

    int iPart(fAOD->GetNumberOfTracks());
    if(iPart < 1) { return; }
    double dEtaLimit = 0.5*dEtaGap;
    double dVz = fAOD->GetPrimaryVertex()->GetZ();

    for(int index(0); index < iPart; ++index)
    {   
        AliAODTrack* POItrack = static_cast<AliAODTrack*>(fAOD->GetTrack(index));
        if(!POItrack || !IsTrackSelected(POItrack)) { continue; }

        double dPt = POItrack->Pt();
        double dPhi = POItrack->Phi();
        double dEta = POItrack->Eta();

        //if (Abs(dEta) < dEtaLimit) { continue; } Why the hell is this here?

        //Check for overlap with RPs
        bool bIsWithinRP = IsWithinRP(POItrack);

        //Check that particle is in POI range
        bool bIsWithinPOI = IsWithinPOI(POItrack);

        if(!bIsWithinPOI) { continue; }

        //Load weights
        double dWeight = GetWeights(dPhi,dEta,dVz);
        if(dWeight <= 0.0) { dWeight = 1.0; }

        //POI with no eta gap
        if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
        {
            for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++) 
            {
                for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                {
                    Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                    Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                    pvectorPtB[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

                    //Check if there is overlap of POI and RP
                    if(bIsWithinRP)
                    {
                        qvectorPtB[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                }  //End for iPower
            }  //End for iHarm

            //POI with eta gap
            if(dEta > dEtaLimit)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++)
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        pvectorPtB10P[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                       
                    } //End for iPower
                }  //End for iHarm
            }
            if (dEta < -dEtaLimit)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++)
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        pvectorPtB10M[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                        
                    } //End for iPower
                } //End for iHarm
            }  // end if eta gap
        } //end if dPtLow < dPt < dPtHigh
    } //end for track

    return;
}

bool AliAnalysisDecorrTask::IsWithinRP(const AliAODTrack* track) const
{
    if(fAbsEtaMax > 0.0 && Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fRFPsPtMax > 0.0 && track->Pt() < fRFPsPtMin) { return kFALSE; }
    if(fRFPsPtMax > 0.0 && track->Pt() > fRFPsPtMax) { return kFALSE; }

    return kTRUE;
}

bool AliAnalysisDecorrTask::IsWithinPOI(const AliAODTrack* track) const
{
    if(fAbsEtaMax > 0.0 && Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fPOIsPtmax > 0.0 && track->Pt() < fPOIsPtmin) { return kFALSE; }
    if(fPOIsPtmin > 0.0 && track->Pt() > fPOIsPtmax) { return kFALSE; }

    return kTRUE;
}

Int_t AliAnalysisDecorrTask::GetSamplingIndex() const
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

//_____________________________________________________________________________
Bool_t AliAnalysisDecorrTask::IsEventSelected()
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
  if(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ()) > fPVtxCutZ) { return kFALSE; }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisDecorrTask::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < fCutNumTPCclsMin && fFilterBit != 2) { return kFALSE; }
  if(fAbsEtaMax > 0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }

  Double_t dTrackXYZ[3] = {0.0,0.0,0.0};
  Double_t dVtxXYZ[3] = {0.0,0.0,0.0};
  Double_t dDCAXYZ[3] = {0.0,0.0,0.0};
  if(fCutDCAzMax > 0.0 || fCutDCAxyMax > 0.0)
  {
      const AliAODVertex* vtx = fAOD->GetPrimaryVertex();
      if(!vtx) { return kFALSE; }

      track->GetXYZ(dTrackXYZ);
      vtx->GetXYZ(dVtxXYZ);
      for(Int_t i(0); i < 3; ++i) { dDCAXYZ[i] = dTrackXYZ[i]-dVtxXYZ[i]; }

      if(fCutDCAzMax > 0.0 && TMath::Abs(dDCAXYZ[2]) > fCutDCAzMax) { return kFALSE; }
      if(fCutDCAxyMax > 0.0 && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0]+dDCAXYZ[1]*dDCAXYZ[1]) > fCutDCAxyMax) { return kFALSE; }
  }
  if(bUseLikeSign)
  {
      if(!(track->Charge() == iSign)) { return kFALSE; }
  }


  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisDecorrTask::IsEventRejectedAddPileUp() const
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
      if((TMath::Abs(track->Eta())) < fAbsEtaMax && (track->GetTPCNcls() >= fCutNumTPCclsMin) && (track->Pt() >= fRFPsPtMin) && (track->Pt() < fRFPsPtMax)) { multTrk++; }
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

double AliAnalysisDecorrTask::GetWeights(double dPhi, double dEta, double dVz)
{
    double dWeight = 1.0;
    if(!fUseWeights3D)
    {
        Int_t iBin = fh2Weights->FindFixBin(dPhi,dEta);
        dWeight = fh2Weights->GetBinContent(iBin);
        return dWeight;
    }
    else
    {
        dWeight = fWeights->GetNUA(dPhi, dEta, dVz);
        return dWeight;
    }
    
}


//Implemented from You's Generic Framework

//_____________________________________________________________________
TComplex AliAnalysisDecorrTask::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);

}

TComplex AliAnalysisDecorrTask::QGap10M(int n, int p)
{

	if(n>=0) return Qvector10M[n][p];
  else return TComplex::Conjugate(Qvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::QGap10P(int n, int p)
{

	if(n>=0) return Qvector10P[n][p];
  else return TComplex::Conjugate(Qvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::pGap10M(int n, int p)
{

	if(n>=0) return pvector10M[n][p];
	else return TComplex::Conjugate(pvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::pGap10P(int n, int p)
{

	if(n>=0) return pvector10P[n][p];
	else return TComplex::Conjugate(pvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::qGap10M(int n, int p)
{

    if(n>=0) return pvector10M[n][p];
    else return TComplex::Conjugate(pvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::qGap10P(int n, int p)
{

    if(n>=0) return pvector10P[n][p];
    else return TComplex::Conjugate(pvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::pPtA(int n, int p)
{

    if(n>=0) return pvector[n][p];
    else return TComplex::Conjugate(pvector[-n][p]);

}

//____________________________________________________________________
TComplex AliAnalysisDecorrTask::pPtB(int n, int p)
{

    if(n>=0) return pvectorPtB[n][p];
    else return TComplex::Conjugate(pvectorPtB[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::qPtA(int n, int p)
{

    if(n>=0) return qvector[n][p];
    else return TComplex::Conjugate(qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::qPtB(int n, int p)
{

    if(n>=0) return qvectorPtB[n][p];
    else return TComplex::Conjugate(qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::pPtBGap10M(int n, int p)
{

    if(n>=0) return pvectorPtB10M[n][p];
    else return TComplex::Conjugate(pvectorPtB10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::pPtBGap10P(int n, int p)
{

    if(n>=0) return pvectorPtB10P[n][p];
    else return TComplex::Conjugate(pvectorPtB10P[-n][p]);

}
//____________________________________________________________________
void AliAnalysisDecorrTask::ResetFlowVector(TComplex (&array)[fNumHarms][fNumPowers])
{
  for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm) {
    for(Int_t iPower(0); iPower < fNumPowers; ++iPower) {
      array[iHarm][iPower](0.0,0.0);
    }
  }
  return;

}

//Correlations
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::Two(int n1, int n2)
{
	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoGap10(int n1, int n2)
{
	TComplex formula = QGap10M(n1,1)*QGap10P(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff(int n1, int n2)
{
	TComplex formula = p(n1,1)*Q(n2,1) - q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiffGap10P(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*QGap10P(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiffGap10M(int n1, int n2)
{
    TComplex formula = pGap10P(n1,1)*QGap10M(n2,1);
    return formula;
}
//____________________________________________________________________
// 2-particles from the same pt
TComplex AliAnalysisDecorrTask::TwoDiff_Pt(int n1, int n2)
{
    TComplex formula = pPtA(n1,1)*pPtA(n2,1) - pPtA(n1+n2,2);
    return formula;
}
//____________________________________________________________________
// 2-particles from the same pt but two different eta regions
TComplex AliAnalysisDecorrTask::TwoDiffGap10_Pt(int n1, int n2)
{
    TComplex formula = pGap10P(n1,1)*pGap10M(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_PtA(int n1, int n2)
{
    TComplex formula = pPtA(n1,1)*pPtA(n2,1) - qPtA(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiffGap10M_PtA(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*pGap10M(n2,1) - pGap10M(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiffGap10P_PtB(int n1, int n2)
{
    TComplex formula = pPtBGap10P(n1,1)*pPtBGap10P(n2,1) - pPtBGap10P(n1+n2,2);
    return formula;
}
TComplex AliAnalysisDecorrTask::TwoDiffGap10_PtB(int n1, int n2)
{
    TComplex formula = pPtBGap10P(n1,1)*pPtBGap10M(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_PtB(int n1, int n2)
{
    TComplex formula = pPtB(n1,1)*pPtB(n2,1) - qPtB(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_PtA_PtB(int n1, int n2)
{
    TComplex formula = pPtA(n1,1)*pPtB(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiffGap10_PtA_PtB(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*pPtBGap10P(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::Three(int n1, int n2, int n3)
{
    
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
    - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
    return formula;
    
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::ThreeGapM(int n1, int n2, int n3)
{
    TComplex formula = QGap10M(n1,1)*QGap10P(n2,1)*QGap10P(n3,1) - QGap10M(n1,1)*QGap10P(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::ThreeGapP(int n1, int n2, int n3)
{
    TComplex formula = QGap10P(n1,1)*QGap10M(n2,1)*QGap10M(n3,1) - QGap10P(n1,1)*QGap10M(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::ThreeDiff(int n1, int n2, int n3)
{

    TComplex formula = p(n1,1)*Q(n2,1)*Q(n3,1)-q(n1+n2,2)*Q(n3,1)-q(n1+n3,2)*Q(n2,1)
    - p(n1,1)*Q(n2+n3,2)+2.*q(n1+n2+n3,3);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::ThreeDiffGapP(int n1, int n2, int n3)
{
    TComplex formula = pGap10P(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)- pGap10P(n1,1)*QGap10M(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::ThreeDiffGapM(int n1, int n2, int n3)
{
    TComplex formula = pGap10M(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)- pGap10M(n1,1)*QGap10P(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::Four(int n1, int n2, int n3, int n4)
{
	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                 		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                 		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                 		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                 		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourGap10(int n1, int n2, int n3, int n4)
{
    TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3,1)*QGap10M(n4,1)-QGap10P(n1+n2,2)*QGap10M(n3,1)*QGap10M(n4,1)
    -QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3+n4,2)+QGap10P(n1+n2,2)*QGap10M(n3+n4,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiff(int n1, int n2, int n3, int n4)
{

	TComplex formula = p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*q(n1+n3,2)*Q(n4,1)
                 		- p(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*q(n1+n4,2)
                 		+ Q(n2+n3,2)*q(n1+n4,2)-p(n1,1)*Q(n3,1)*Q(n2+n4,2)+q(n1+n3,2)*Q(n2+n4,2)
                 		+ 2.*Q(n3,1)*q(n1+n2+n4,3)-p(n1,1)*Q(n2,1)*Q(n3+n4,2)+q(n1+n2,2)*Q(n3+n4,2)
                 		+ 2.*Q(n2,1)*q(n1+n3+n4,3)+2.*p(n1,1)*Q(n2+n3+n4,3)-6.*q(n1+n2+n3+n4,4);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiff_PtA_PtA(int n1, int n2, int n3, int n4)  
{
    TComplex formula = pPtA(n1,1)*pPtA(n2,1)*pPtA(n3,1)*pPtA(n4,1)-pPtA(n1+n2,2)*pPtA(n3,1)*pPtA(n4,1)-pPtA(n2,1)*pPtA(n1+n3,2)*pPtA(n4,1)
    - pPtA(n1,1)*pPtA(n2+n3,2)*pPtA(n4,1)+2.*pPtA(n1+n2+n3,3)*pPtA(n4,1)-pPtA(n2,1)*pPtA(n3,1)*pPtA(n1+n4,2)
    + pPtA(n2+n3,2)*pPtA(n1+n4,2)-pPtA(n1,1)*pPtA(n3,1)*pPtA(n2+n4,2)+pPtA(n1+n3,2)*pPtA(n2+n4,2)
    + 2.*pPtA(n3,1)*pPtA(n1+n2+n4,3)-pPtA(n1,1)*pPtA(n2,1)*pPtA(n3+n4,2)+pPtA(n1+n2,2)*pPtA(n3+n4,2)
    + 2.*pPtA(n2,1)*pPtA(n1+n3+n4,3)+2.*pPtA(n1,1)*pPtA(n2+n3+n4,3)-6.*pPtA(n1+n2+n3+n4,4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiff_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_PtA(n1, n2)*TwoDiff_PtB(n3, n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiffGap10M_PtA(n1, n2)*TwoDiffGap10P_PtB(n3, n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10_OS_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiffGap10_Pt(n1, n2)*TwoDiffGap10_PtB(n3, n4);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10P(int n1, int n2, int n3, int n4)
{
    TComplex formula = pGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3,1)*QGap10M(n4,1)-qGap10P(n1+n2,2)*QGap10M(n3,1)*QGap10M(n4,1)
    -pGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3+n4,2)+qGap10P(n1+n2,2)*QGap10M(n3+n4,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10M(int n1, int n2, int n3, int n4)
{
    TComplex formula = pGap10M(n1,1)*QGap10M(n2,1)*QGap10P(n3,1)*QGap10P(n4,1)-qGap10M(n1+n2,2)*QGap10P(n3,1)*QGap10P(n4,1)
    -pGap10M(n1,1)*QGap10M(n2,1)*QGap10P(n3+n4,2)+qGap10M(n1+n2,2)*QGap10P(n3+n4,2);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::Five(int n1, int n2, int n3, int n4, int n5)
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
TComplex AliAnalysisDecorrTask::Six(int n1, int n2, int n3, int n4, int n5, int n6)
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
TComplex AliAnalysisDecorrTask::SixDiff(int n1, int n2, int n3, int n4, int n5, int n6)
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
TComplex AliAnalysisDecorrTask::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{
    
    TComplex Correlation = {0, 0};
    int Narray[] = {n1, n2, n3, n4, n5, n6};
    
    for(int k=7; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4
        
        int array[6] = {0,1,2,3,4,5};
        int iPerm = 0;
        int argument = 0;
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
            cout<<"invalid range of k"<<endl;
            return {0,0};
        }
        
    }// loop over k
    
    return Correlation;
    
}
//_____________________________________________________________________________
TComplex AliAnalysisDecorrTask::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
    
    TComplex Correlation = {0, 0};
    int Narray[] = {n1, n2, n3, n4, n5, n6, n7};
    
    for(int k=8; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4
        
        int array[7] = {0,1,2,3,4,5,6};
        int iPerm = 0;
        int argument = 0;
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
            cout<<"invalid range of k"<<endl;
            return {0,0};
        }
        
    }// loop over k
    
    return Correlation;
    
}
//_____________________________________________________________________________
void AliAnalysisDecorrTask::Terminate(Option_t *)
{

}