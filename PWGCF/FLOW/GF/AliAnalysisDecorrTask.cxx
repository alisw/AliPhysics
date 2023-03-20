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
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliDecorrFlowCorrTask.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisDecorrTask.h"


class AliAnalysisDecorrTask;


using namespace std;
using namespace TMath;


ClassImp(AliAnalysisDecorrTask)

AliAnalysisDecorrTask::AliAnalysisDecorrTask() : AliAnalysisTaskSE(),
    fEventCuts(),
    fFlowList{nullptr},
    //fFlowWeights{nullptr},
    fQA{nullptr},
    lSampleList{nullptr},
    fWeights(0),
    fWeightList{nullptr},
    fV0MMulti(0),
    fMulti(0),
    hITSclsB{nullptr},
    hTPCclsB{nullptr},
    hTPCchi2B{nullptr},
    hDCAxyB{nullptr},   
    hITSclsA{nullptr},
    hTPCclsA{nullptr},
    hTPCchi2A{nullptr},
    hDCAxyA{nullptr}, 
    hPtPhiEtaB{nullptr},  
    hPtPhiEtaA{nullptr},  
    hNumTracksB{nullptr},
    hNumTracksA{nullptr},
    hNumHighPtTracksA{nullptr},
    hCharge{nullptr},
    fhEventSel{nullptr},

    fIndexSampling{0},
    fTaskCounter(0),
    fAOD(nullptr),
    fIsMC(kFALSE),
    fMCEvent(nullptr),
    fEvent(nullptr),
    fInitTask{kFALSE},  
    fOnTheFly(kFALSE),
    fImpactParameterMC{-1.0},
    fVecCorrTask{},
    fCorrName{""},
    fCorrLabel{""},
    fCorrOrder(0),

    fSampling{kFALSE},
    fRedTracks{kFALSE},
    fTrackprevent{0.5},
    fFillQA(kFALSE),

    fTrigger(AliVEvent::kMB+AliVEvent::kINT7),
    fEventRejectAddPileUp(kTRUE),
    fCentralPileupCut(500),
    fDefaultPileupCut(15000),
    fCentEstimator("V0M"), 
    fFilterBit(96),
    fPtAxis(new TAxis()),
    fCentAxis(new TAxis()),
    fCentMin{0.0},
    fCentMax{50.0},
    fOverrideCentrality(-1.0),
    fPVtxCutX{3.},
    fPVtxCutY{3.},
    fPVtxCutZ{10.0},
    fRequireHighPtTracks{kFALSE},
    fNHighPtTracks{1},
    fHighPtCut{7.0},
    bOnlyPrimariesAndCh(false),

    fCutChargedTrackFilterBit{96},
    fCutNumTPCclsMin{70},
    fCutDCAzMax{0.0},
    fCutDCAxyMax{0.0},
    fChi2Cut{0.0},
    bUseLikeSign(kFALSE),
    iSign(0),

    fAbsEtaMax(0.8),
    //dEtaGap(1.0),
    fEtaBinNum{0},
    fPhiBinNum{60},
    fUseWeightsOne(kFALSE),
    fCurrSystFlag(0),
    fFillWeights(kFALSE),
    fNumSamples{1},

    //bHasGap(kTRUE),
    bDiff(kFALSE),
    bRef(kTRUE),
    bPtA(kFALSE),
    bPtRef(kFALSE),
    bPtB(kFALSE),

    fPOIsPtmax(10.0),
    fPOIsPtmin(0.2),
    fRFPsPtMax(5.0),
    fRFPsPtMin(0.2),
    fRequireTwoPart(kFALSE),
    bEqualPt(kFALSE)
{}
//_____________________________________________________________________________
AliAnalysisDecorrTask::AliAnalysisDecorrTask(const char* name, bool IsMC, bool isOnTheFly) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fFlowList{nullptr},
    //fFlowWeights{nullptr},
    fQA{nullptr},
    lSampleList{nullptr},
    fWeights(0),
    fWeightList{nullptr},
    fV0MMulti(0),
    fMulti(0),
    hITSclsB{nullptr},
    hTPCclsB{nullptr},
    hTPCchi2B{nullptr},
    hDCAxyB{nullptr},   
    hITSclsA{nullptr},
    hTPCclsA{nullptr},
    hTPCchi2A{nullptr},
    hDCAxyA{nullptr}, 
    hPtPhiEtaB{nullptr},
    hPtPhiEtaA{nullptr},  
    hNumTracksB{nullptr},
    hNumTracksA{nullptr},
    hNumHighPtTracksA{nullptr},
    hCharge{nullptr},
    fhEventSel{nullptr},

    fIndexSampling{0},
    fTaskCounter(0),
    fAOD(nullptr),
    fIsMC(IsMC),
    fMCEvent(nullptr),
    fEvent(nullptr),
    fInitTask{kFALSE}, 
    fOnTheFly(isOnTheFly),
    fImpactParameterMC{-1.0},
    fVecCorrTask{},
    fCorrName{""},
    fCorrLabel{""},
    fCorrOrder(0),

    fSampling{kFALSE},
    fRedTracks{kFALSE},
    fTrackprevent{0.5},
    fFillQA(kFALSE),

    fTrigger(AliVEvent::kMB+AliVEvent::kINT7),
    fEventRejectAddPileUp(kTRUE),
    fCentralPileupCut(500),
    fDefaultPileupCut(15000),
    fCentEstimator("V0M"), 
    fFilterBit(96),
    fPtAxis(new TAxis()),
    fCentAxis(new TAxis()),
    fCentMin{0.0},
    fCentMax{50.0},
    fOverrideCentrality(-1.0),
    fPVtxCutX{3.},
    fPVtxCutY{3.},
    fPVtxCutZ{10.0},
    fRequireHighPtTracks{kFALSE},
    fNHighPtTracks{0},
    fHighPtCut{7.0},
    bOnlyPrimariesAndCh(false),

    fCutChargedTrackFilterBit{96},
    fCutNumTPCclsMin{70},
    fCutDCAzMax{0.0},
    fCutDCAxyMax{0.0},
    fChi2Cut{0.0},
    bUseLikeSign(kFALSE),
    iSign(0),

    fAbsEtaMax(0.8),
    //dEtaGap(1.0),
    fEtaBinNum{0},
    fPhiBinNum{60},
    fUseWeightsOne(kFALSE),
    fCurrSystFlag(0),
    fFillWeights(kFALSE),
    fNumSamples{1},

    //bHasGap(kTRUE),
    bDiff(kFALSE),
    bRef(kTRUE),
    bPtA(kFALSE),
    bPtRef(kFALSE),
    bPtB(kFALSE),

    fPOIsPtmax(10.0),
    fPOIsPtmin(0.2),
    fRFPsPtMax(5.0),
    fRFPsPtMin(0.2),
    fRequireTwoPart(kFALSE),
    bEqualPt(kFALSE)
{

    if(fOnTheFly)
    {
      vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
      vector<double> cent = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0};
      for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i];
    }
    
    DefineInput(0, TChain::Class());
    if(!fIsMC&&!fOnTheFly) DefineInput(1, TList::Class()); 
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisDecorrTask::~AliAnalysisDecorrTask()
{
    
    if(fFlowList) delete fFlowList;
    //if(fFlowWeights) delete fFlowWeights;
    if(fQA) delete fQA;
}
bool AliAnalysisDecorrTask::InitTask()
{
    if(fFillWeights) { AliWarning("\n Filling weights. \n"); }

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
    if(!fIsMC&&!fOnTheFly) 
    {
        fWeightList = (TList*) GetInputData(1);
        if(!fWeightList) { AliFatal("Weight List not found!\n"); return kFALSE; }
        AliInfo("Weight List loaded");
    }
    if(fIsMC) AliInfo("Running over MC data!");
    if(fOnTheFly) AliInfo("Running over generated MC data!");
    AliInfo("Creating output objects");
    return kTRUE;
}
void AliAnalysisDecorrTask::NotifyRun() {
    if(fOnTheFly) return;
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    fEventCuts.fESDvsTPConlyLinearCut[0] = fCentralPileupCut;
}
void AliAnalysisDecorrTask::CreateProfiles(AliDecorrFlowCorrTask* task, double* lCentEdges, double* lPtEdges)
{
    TH1* profile = nullptr;
    TH1* profDiff = nullptr;
    TH1* profPtA = nullptr;
    TH1* profPtRef = nullptr;
    TH1* profPtRefP = nullptr;
    TH1* profPtAPtB = nullptr;
    TH1* profPtAPtBP = nullptr;

    lSampleList = new TList();
    lSampleList->SetOwner(kTRUE);
    lSampleList->SetName(fCorrName.Data());
    if(bRef)
    {
        profile = new TProfile(Form("%s",fCorrName.Data()),Form("%s",fCorrLabel.Data()),NcentBin,lCentEdges);
        if(!profile) { fInitTask = kFALSE; AliError("Centrality profile not created"); task->PrintTask(); return; }
        profile->Sumw2();
        lSampleList->Add(profile);
    }
    if(bDiff) 
    {      
        profDiff = new TProfile2D(Form("%s_diff",fCorrName.Data()),Form("%s_diff",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
        if(!profDiff) { fInitTask = kFALSE; AliError("Differential profile not created"); task->PrintTask(); return; }
        profDiff->Sumw2();
        lSampleList->Add(profDiff);
    }
    if(bPtA)
    {
        profPtA = new TProfile2D(Form("%s_PtA",fCorrName.Data()),Form("%s_PtA",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
        if(!profPtA) { fInitTask = kFALSE; AliError("\n\n\nPtA profile not created\n\n\n"); task->PrintTask(); return; }
        profPtA->Sumw2();
        lSampleList->Add(profPtA);      
    }
    if(bPtRef)
    {
        profPtRef = new TProfile2D(Form("%s_PtRef",fCorrName.Data()),Form("%s_PtRef",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
        if(!profPtRef) { fInitTask = kFALSE; AliError("\n\n\nPtRef profile not created\n\n\n"); task->PrintTask(); return; }
        profPtRef->Sumw2();
        lSampleList->Add(profPtRef); 
        if(fCorrOrder==4 && task->fiHarm[1] > 0)
        {
            profPtRefP = new TProfile2D(Form("%s_PtRefP",fCorrName.Data()),Form("%s_PtRefP",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
            if(!profPtRefP) { fInitTask = kFALSE; AliError("\n\n\nPtRefP profile not created\n\n\n"); task->PrintTask(); return; }
            profPtRefP->Sumw2();
            lSampleList->Add(profPtRefP); 
        }        
    }
    if(bPtB)
    { 
        profPtAPtB = new TProfile3D(Form("%s_PtAPtB",fCorrName.Data()),Form("%s_PtAPtB",fCorrLabel.Data()),NcentBin,lCentEdges, NPtBin,lPtEdges, NPtBin, lPtEdges);
        if(!profPtAPtB) { fInitTask = kFALSE; AliError("PtAPtB profile not created"); task->PrintTask(); return; }
        profPtAPtB->Sumw2();   
        lSampleList->Add(profPtAPtB);
        if(fCorrOrder==4 && task->fiHarm[1] > 0)
        {
            profPtAPtBP = new TProfile3D(Form("%s_PtAPtBP",fCorrName.Data()),Form("%s_PtAPtBP",fCorrLabel.Data()),NcentBin,lCentEdges, NPtBin,lPtEdges, NPtBin, lPtEdges);
            if(!profPtAPtBP) { fInitTask = kFALSE; AliError("PtAPtBP profile not created"); task->PrintTask(); return; }
            profPtAPtBP->Sumw2();
            lSampleList->Add(profPtAPtBP);
        }
    }

    for(Int_t iSample(0); iSample < fNumSamples; ++iSample)
    {
        if(!fSampling) { break; }
        if(bRef)
        {
            profile = new TProfile(Form("%s_sample%d",fCorrName.Data(),iSample),Form("%s",fCorrLabel.Data()),NcentBin,lCentEdges);
            if(!profile) { fInitTask = kFALSE; AliError("Centrality profile not created"); task->PrintTask(); return; }

            profile->Sumw2();
            lSampleList->Add(profile);
        }
        if(bDiff) 
        {
            
            profDiff = new TProfile2D(Form("%s_diff_sample%d",fCorrName.Data(),iSample),Form("%s_diff",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
            if(!profDiff) { fInitTask = kFALSE; AliError("Differential profile not created"); task->PrintTask(); return; }

            profDiff->Sumw2();
            lSampleList->Add(profDiff);
        }
        if(bPtA)
        {
            profPtA = new TProfile2D(Form("%s_PtA_sample%d",fCorrName.Data(),iSample),Form("%s_PtA",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
            if(!profPtA) { fInitTask = kFALSE; AliError("\n\n\nPtA profile not created\n\n\n"); task->PrintTask(); return; }


            profPtA->Sumw2();
            lSampleList->Add(profPtA);
            
        }
        if(bPtRef)
        {
            profPtRef = new TProfile2D(Form("%s_PtRef_sample%d",fCorrName.Data(),iSample),Form("%s_PtRef",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
            if(!profPtRef) { fInitTask = kFALSE; AliError("\n\n\nPtRef profile not created\n\n\n"); task->PrintTask(); return; }
            profPtRef->Sumw2();
            lSampleList->Add(profPtRef); 
            if(fCorrOrder==4 && task->fiHarm[1] > 0)
            {
                profPtRefP = new TProfile2D(Form("%s_PtRefP_sample%d",fCorrName.Data(),iSample),Form("%s_PtRefP",fCorrLabel.Data()),NcentBin,lCentEdges,NPtBin,lPtEdges);
                if(!profPtRefP) { fInitTask = kFALSE; AliError("\n\n\nPtRefP profile not created\n\n\n"); task->PrintTask(); return; }            
                profPtRefP->Sumw2();
                lSampleList->Add(profPtRefP); 
            }

        }
        if(bPtB)
        { 
            profPtAPtB = new TProfile3D(Form("%s_PtAPtB_sample%d",fCorrName.Data(),iSample),Form("%s_PtAPtB",fCorrLabel.Data()),NcentBin,lCentEdges, NPtBin,lPtEdges, NPtBin, lPtEdges);
            if(!profPtAPtB) { fInitTask = kFALSE; AliError("PtAPtB profile not created"); task->PrintTask(); return; }
            profPtAPtB->Sumw2();
            lSampleList->Add(profPtAPtB);
            if(fCorrOrder==4 && task->fiHarm[1] > 0) {
                profPtAPtBP = new TProfile3D(Form("%s_PtAPtBP_sample%d",fCorrName.Data(),iSample),Form("%s_PtAPtBP",fCorrLabel.Data()),NcentBin,lCentEdges, NPtBin,lPtEdges, NPtBin, lPtEdges);
                if(!profPtAPtBP) { fInitTask = kFALSE; AliError("PtAPtBP profile not created"); task->PrintTask(); return; }
                profPtAPtBP->Sumw2();
                lSampleList->Add(profPtAPtBP);
            }            
        }
    } 
    fFlowList->Add(lSampleList);
}
void AliAnalysisDecorrTask::UserCreateOutputObjects()
{
    fFlowList = new TList();
    fFlowList->SetOwner(kTRUE);
    fQA = new TList();
    fQA->SetOwner(kTRUE);

    fInitTask = InitTask();
    if(!fInitTask) { return; }
    
    int NV0M_default = 10;
    double lV0M_default[] = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
    NcentBin = (fCentAxis)?fCentAxis->GetNbins():NV0M_default;
    double* centEdges = (fCentAxis)?GetBinsFromAxis(fCentAxis):lV0M_default;
    NPtBin = fPtAxis->GetNbins();
    double* PtEdges = GetBinsFromAxis(fPtAxis);
    Int_t iNumTasks = fVecCorrTask.size();
    if(iNumTasks < 1) { AliFatal("UserCreateOutputObjects: At least on correlation needs to be specified"); return; }
    for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
    {
        AliDecorrFlowCorrTask* task = fVecCorrTask.at(iTask);
        if(!task) { fInitTask = kFALSE; AliError(Form("AliDecorrFlowCorrTask%d does not exist",iTask)); return; }
        InitTypes(task);
        fCorrName = task->fsName;
        fCorrLabel = task->fsLabel;
        fCorrOrder = task->fiNumHarm;
        CreateProfiles(task,centEdges,PtEdges);
    } 

    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",NcentBin,centEdges);
    fQA->Add(fV0MMulti);
    fMulti = new TH1D("Multi","Multi",NV0M_default,lV0M_default);
    fQA->Add(fMulti);

    TString sEventSelLabel[] = {"Input","AliEventCuts","Trigger","Centrality","Pileup","V_{z}","High p_{T} tracks"};
    const Int_t iNBinsEventSelCounter = sizeof(sEventSelLabel)/sizeof(sEventSelLabel[0]);
    fhEventSel = new TH1D("fhEventSel","Event selection",iNBinsEventSelCounter,0,iNBinsEventSelCounter);
    for(Int_t i(0); i < iNBinsEventSelCounter; i++) fhEventSel->GetXaxis()->SetBinLabel(i+1,sEventSelLabel[i].Data());
    fQA->Add(fhEventSel);

    if(fFillQA)
    {
        fEventCuts.AddQAplotsToList(fQA,true);

        hITSclsB = new TH1I("ITS_clusters_on_trackB","ITS clusters on track",8,0,8);
        hITSclsA = (TH1I*)hITSclsB->Clone("ITS_clusters_on_trackA");
        fQA->Add(hITSclsB); fQA->Add(hITSclsA);    
        hTPCclsB = new TH1I("TPC_clusters_on_trackB","TPC clusters on track",159,1,160);
        hTPCclsA = (TH1I*)hTPCclsB->Clone("TPC_clusters_on_trackA");
        fQA->Add(hTPCclsB); fQA->Add(hTPCclsA); 
        hTPCchi2B = new TH1D("TPC_chi2_pr_clusterB","TPC #chi^{2}/clusters",100,0.0,5.0);
        hTPCchi2A = (TH1D*)hTPCchi2B->Clone("TPC_chi2_pr_clusterA");
        fQA->Add(hTPCchi2B); fQA->Add(hTPCchi2A); 
        hDCAxyB = new TH2D("DCAxyB","DCAxy before",50,0.2,5.0,200,-5.0, 5.0);
        fQA->Add(hDCAxyB); 
        hDCAxyA = new TH2D("DCAxyA","DCAxy after",50,0.2,5.0,200,-5.0, 5.0);
        fQA->Add(hDCAxyA); 
        hPtPhiEtaB = new TH3D("pt_phi_etaB","Pt, phi, eta distribution",50,0.2,20.0,60,0.0,TMath::TwoPi(),50,-1.0,1.0);
        fQA->Add(hPtPhiEtaB);
        hPtPhiEtaA = new TH3D("pt_phi_etaA","Pt, phi, eta distribution",50,0.2,20.0,60,0.0,TMath::TwoPi(),50,-1.0,1.0);
        fQA->Add(hPtPhiEtaA);
        hNumTracksB = new TH1D("hNumTracksB","Number of tracks before",100,0.0,2500.0);
        fQA->Add(hNumTracksB);
        hNumTracksA = new TH1D("hNumTracksA","Number of tracks after",100,0.0,2500.0);
        fQA->Add(hNumTracksA);
        hNumHighPtTracksA = new TH1D("hNumHighPtTracksA",Form("Number of track above %.1f GeV",fHighPtCut),50,0.0,50.0);
        fQA->Add(hNumHighPtTracksA);
        int nChBins = 3;
        double chBins[] = {-1.5,-0.5,0.5,1.5};
        hCharge = new TH1I("hCharge","Charge distribution",nChBins,chBins);
        fQA->Add(hCharge);

    }    
    /*
    if(fEventRejectAddPileUp)
    {
        fhQAEventsfMult32vsCentr = new TH2D("fhQAEventsfMult32vsCentr", "; centrality V0M; TPC multiplicity (FB32)", 100, 0, 100, 100, 0, 3000);
        fQA->Add(fhQAEventsfMult32vsCentr);
        fhQAEventsMult128vsCentr = new TH2D("fhQAEventsfMult128vsCentr", "; centrality V0M; TPC multiplicity (FB128)", 100, 0, 100, 100, 0, 5000);
        fQA->Add(fhQAEventsMult128vsCentr);
        fhQAEventsfMultTPCvsTOF = new TH2D("fhQAEventsfMultTPCvsTOF", "; TPC FB32 multiplicity; TOF multiplicity", 200, 0, 4000, 200, 0, 2000);
        fQA->Add(fhQAEventsfMultTPCvsTOF);
        fhQAEventsfMultTPCvsESD = new TH2D("fhQAEventsfMultTPCvsESD", "; TPC FB128 multiplicity; ESD multiplicity", 200, 0, 7000, 300, -1000, 35000);
        fQA->Add(fhQAEventsfMultTPCvsESD);
    }
    */
    
    if(!fInitTask) AliFatal("UserCreateOutputObjects: Task not initialized!");
    AliInfo("Output objects created");
    PostData(1, fFlowList);
    PostData(2, fQA);

}
bool AliAnalysisDecorrTask::LoadWeights()
{
    if(fCurrSystFlag == 0) fWeights = (AliGFWWeights*)fWeightList->FindObject(Form("w%i",fAOD->GetRunNumber()));
    else fWeights = (AliGFWWeights*)fWeightList->FindObject(Form("w%i_SystFlag%i",fAOD->GetRunNumber(),fCurrSystFlag));
    if(!fWeights)
    {
        printf("Weights could not be found in list!\n");
        return kFALSE;
    }
    fWeights->CreateNUA();
    return kTRUE;
    
}
void AliAnalysisDecorrTask::FillWeights()
{
    int NParts(fAOD->GetNumberOfTracks());
    if(NParts < 1) { return; }
    double Vz = fAOD->GetPrimaryVertex()->GetZ();

    if(!fFillQA && !fFillWeights) { return; }
    int NumTracksB = 0;
    int NumTracksA = 0;
    int NumTracksHighPt = 0;
    for(int index(0); index < NParts; ++index)
    {   
        Float_t dcaxy = 0.0;
        Float_t dcaz = 0.0;
        Float_t tpcchi2 = 0.0;
        Float_t tpcchi2percls = 0.0;
        int ntpccls = 0;
        int nitscls = 0;
        Float_t phi = 0.0;
        Float_t eta = 0.0;
        Float_t pt = 0.0;
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(index));
        if(!track) { continue; }
        bool pass = IsTrackSelected(track);
        pt = track->Pt(); 
        phi = track->Phi();
        eta = track->Eta(); 

        if(fFillQA)
        {
            track->GetImpactParameters(dcaxy,dcaz);
            tpcchi2=track->GetTPCchi2();
            ntpccls=track->GetTPCNcls();
            tpcchi2percls = (ntpccls==0)?0.0:tpcchi2/ntpccls;
            nitscls=track->GetNcls(0);
            hDCAxyB->Fill(pt,dcaxy); if(pass) hDCAxyA->Fill(pt,dcaxy);
            hTPCchi2B->Fill(tpcchi2percls); if(pass) hTPCchi2A->Fill(tpcchi2percls);
            hTPCclsB->Fill(ntpccls); if(pass) hTPCclsA->Fill(ntpccls);
            hITSclsB->Fill(nitscls); if(pass) hITSclsA->Fill(nitscls);
            hPtPhiEtaB->Fill(pt,phi,eta); if(pass) hPtPhiEtaA->Fill(pt,phi,eta);
            if(pass) NumTracksB++;
            if(fRedTracks && pass && index < fTrackprevent*NParts) 
            { 
                NumTracksA++;
            }
            if(fRequireHighPtTracks && pass) { 
                if(pt >= fHighPtCut) ++NumTracksHighPt; 
                if(NumTracksHighPt >= fNHighPtTracks) hNumHighPtTracksA->Fill(NumTracksHighPt);
            }
        }
    }
    hNumTracksB->Fill(NumTracksB);
    hNumTracksA->Fill(NumTracksA);
}
void AliAnalysisDecorrTask::UserExec(Option_t *)
{
    if(!fInitTask) { AliFatal("Something went wrong! Task not initialized"); return; }
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliVEventHandler* inputHandler = (AliVEventHandler*) (man->GetInputEventHandler());
    if(fIsMC) {
        fMCEvent = inputHandler->MCEvent();
        if (!fMCEvent) return;
    }

    if(fOnTheFly) { fMCEvent = getMCEvent(); }
    else
    {   
        fEvent = dynamic_cast<AliVEvent*>(InputEvent());
        fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
        if(!fAOD) { return; }
        if(!fEvent) { AliFatal("fEvent not found!"); return; }
    }
    double lCent = (fOverrideCentrality>0)?fOverrideCentrality:getCentrality();
    if(lCent<0) return;
    if(!fOnTheFly)
    {
        if(!IsEventSelected(lCent)) { return; }
    }
    
    if(!fOnTheFly && fFillQA) FillWeights();

    if(!fUseWeightsOne && !fOnTheFly && !fIsMC) {
        bool loaded = LoadWeights();
        if(!loaded) { AliFatal("\n\n\n\n\n\n\n\n Weights could not be loaded \n\n\n\n\n\n\n\n"); return; }
    }
    fIndexSampling = GetSamplingIndex();
    
    Int_t iNumTask = fVecCorrTask.size();
    fTaskCounter = 0;
    int iPtRefLS = 0;
    int iPtRefOS = 0;
    for(Int_t iTask(0); iTask < iNumTask; ++iTask)
    {
        const AliDecorrFlowCorrTask* const task = fVecCorrTask.at(iTask);
        if(!task) { AliError("AliDecorrFlowCorrTask does not exist"); return; }
        fCorrName = task->fsName;
        fCorrOrder = task->fiNumHarm;
        lSampleList = (TList*)fFlowList->FindObject(fCorrName);
        InitTypes(task);
        if(bRef) refData.push_back(StoreRef());
        if(bDiff||bPtA||bPtRef) ptaData.push_back(StorePtA());
        if(bPtB) ptbData.push_back(StorePtB());
        if(fCorrOrder==4&&bPtRef) { if(task->fiHarm[1] > 0) iPtRefLS = fTaskCounter; else iPtRefOS = fTaskCounter; }
        FillRPvectors(task);
        if(bRef) { CalculateCorrelations(task, lCent, -1.0, -1.0, bRef, kFALSE, kFALSE, kFALSE, kFALSE); } 
        if(!(bDiff || bPtA || bPtRef || bPtB)) continue;  
        double lPt,lPtLow,lPtHigh;
        for(int iPtA(1); iPtA <= NPtBin; ++iPtA)
        {
            lPt = fPtAxis->GetBinCenter(iPtA);
            lPtLow = fPtAxis->GetBinLowEdge(iPtA);
            lPtHigh = fPtAxis->GetBinUpEdge(iPtA);
            int TrackCounter = 0;
            TrackCounter = FillPOIvectors(task, lPtLow, lPtHigh);     
            if(fRequireTwoPart) {
                if(TrackCounter > 1) CalculateCorrelations(task, lCent, lPt, -1.0, kFALSE, bDiff, kFALSE, kFALSE, kFALSE);
                CalculateCorrelations(task, lCent, lPt, -1.0, kFALSE, kFALSE, bPtA, bPtRef, kFALSE); }
            else {
                CalculateCorrelations(task, lCent, lPt, -1.0, kFALSE, bDiff, bPtA, bPtRef, kFALSE); }                 
            if(!bPtB || lPt > 5.0 || lCent > fCentMax) continue;   //Save cpu by restricting double pt loops to central and semicentral centralities and low pt     
            for(int iPtB(1); iPtB <= NPtBin; ++iPtB)
            { 
                if(iPtB == iPtA) bEqualPt = kTRUE; else bEqualPt = kFALSE;
                double lPtB = fPtAxis->GetBinCenter(iPtB);
                double lPtBLow = fPtAxis->GetBinLowEdge(iPtB);
                double lPtBHigh = fPtAxis->GetBinUpEdge(iPtB);
                FillPtBvectors(task, lPtBLow, lPtBHigh);
                CalculateCorrelations(task, lCent, lPt, lPtB, kFALSE, kFALSE, kFALSE, kFALSE, bPtB); 
            } 
        }
        ++fTaskCounter;        
    }
    fV0MMulti->Fill(lCent);
    fMulti->Fill(lCent);
    PostData(1, fFlowList);
    PostData(2, fQA);
}
void AliAnalysisDecorrTask::CalculateCorrelations(const AliDecorrFlowCorrTask* const task, const double &lCent, const double &lpta, const double &lptb, Bool_t bRef, Bool_t bDiff, Bool_t bPtA, Bool_t bPtRef, Bool_t bPtB)
{
        double dGap = task->fdGaps[0]; 
        bool bHasGap = (dGap > -1.0)?kTRUE:kFALSE;
        TComplex c[2] = {TComplex(0.0,0.0,kFALSE),TComplex(0.0,0.0,kFALSE)};

        if(bRef) { 
            CalculateGFW(c,task,0,bHasGap); 
            FillProfiles(c,lCent,""); 
        }
        if(bDiff) { 
            CalculateGFW(c,task,1,bHasGap); 
            FillProfiles(c,lCent,lpta,"_diff"); 
        }
        if(bPtA) { 
            CalculateGFW(c,task,2,bHasGap); 
            FillProfiles(c,lCent,lpta,"_PtA"); 
        }
        if(bPtRef)
        {
            CalculateGFW(c,task,3,bHasGap);
            FillProfiles(c,lCent,lpta,"_PtRef");
            if(fCorrOrder==4 && task->fiHarm[1] > 0)
            {
                CalculateGFW(c,task,3,bHasGap,true);
                FillProfiles(c,lCent,lpta,"_PtRefP");
            }
        }        
        if(bPtB)
        {
            CalculateGFW(c,task,4,bHasGap);
            FillProfiles(c,lCent,lpta,lptb,"_PtAPtB");
            if(fCorrOrder==4 && task->fiHarm[1] > 0)
            {
                CalculateGFW(c,task,4,bHasGap,true);
                FillProfiles(c,lCent,lpta,lptb,"_PtAPtBP");
            }
        }
        
    return;
}
void AliAnalysisDecorrTask::CalculateGFW(TComplex c[2], const AliDecorrFlowCorrTask* const task, int type, bool gap, bool switchreg)
{ 
    if(gap)
    {
        switch(type)
        {
            case 0 : 
                if(fCorrOrder==2) { c[0] = TwoGap10(0,0); c[1] = TwoGap10(task->fiHarm[0],task->fiHarm[1]); }
                else { c[0] = FourGap10(0,0,0,0); c[1] = FourGap10(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                break;
            case 1 :
                if(fCorrOrder==2) { c[0] = TwoDiffGap10M(0,0); c[1] = TwoDiffGap10M(task->fiHarm[0],task->fiHarm[1]); }
                else { c[0] = FourDiffGap10M(0,0,0,0); c[1] = FourDiffGap10M(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                break;
            case 2 :
                if(fCorrOrder==2) { c[0] = TwoDiffGap10_Pt(0,0); c[1] = TwoDiffGap10_Pt(task->fiHarm[0],task->fiHarm[1]); }
                else { c[0] = FourDiffGap10_PtA_PtA(0,0,0,0); c[1] = FourDiffGap10_PtA_PtA(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                break;
            case 3 : 
                if(fCorrOrder==2) { AliWarning("CalculateGFW: PtRef for 2pc is just dn{2}!"); return; }
                else {
                    if(task->fiHarm[1] > 0) {
                        if(switchreg) {
                            c[0] = FourGapP_2Diff_2Ref(0,0,0,0);
                            c[1] = FourGapP_2Diff_2Ref(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                        else {
                            c[0] = FourGapM_2Diff_2Ref(0,0,0,0);
                            c[1] = FourGapM_2Diff_2Ref(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); } 
                    }
                    else {
                        c[0] = FourGap_2Diff_2Ref_OS(0,0,0,0);
                        c[1] = FourGap_2Diff_2Ref_OS(task->fiHarm[0],task->fiHarm[2],task->fiHarm[1],task->fiHarm[3]); 
                    }
                }
                break;
            case 4 :
                if(fCorrOrder==2) { c[0] = TwoDiffGap10_PtA_PtB(0,0); c[1] = TwoDiffGap10_PtA_PtB(task->fiHarm[0],task->fiHarm[1]); }
                else { 
                    if(task->fiHarm[1] > 0) { 
                        if(switchreg) { 
                            c[0] = FourDiffGap10P_PtA_PtB(0,0,0,0); c[1] = FourDiffGap10P_PtA_PtB(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); } 
                        else { 
                            c[0] = FourDiffGap10M_PtA_PtB(0,0,0,0); c[1] = FourDiffGap10M_PtA_PtB(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); } }
                    else { 
                        c[0] = FourDiffGap10_OS_PtA_PtB(0,0,0,0);
                        c[1] = FourDiffGap10_OS_PtA_PtB(task->fiHarm[0],task->fiHarm[2],task->fiHarm[1],task->fiHarm[3]); }
                }   
                break;
        }
    }
    else
    {
        switch(type)
        {
            case 0 :             
                if(fCorrOrder==2) { c[0] = Two(0,0); c[1] = Two(task->fiHarm[0],task->fiHarm[1]); }
                else if(fCorrOrder==4) { c[0] = Four(0,0,0,0); c[1] = Four(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                else if(fCorrOrder==6) { c[0] = Six(0,0,0,0,0,0); c[1] = Six(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]); }
                else if(fCorrOrder==8) { c[0] = Eight(0,0,0,0,0,0,0,0); c[1] = Eight(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]); }
                break;
            case 1 :
                if(fCorrOrder==2) { c[0] = TwoDiff(0,0); c[1] = TwoDiff(task->fiHarm[0],task->fiHarm[1]); }
                else { c[0] = FourDiff(0,0,0,0); c[1] = FourDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                break;
            case 2 :
                if(fCorrOrder==2) {c[0] = TwoDiff_PtA(0,0); c[1] = TwoDiff_PtA(task->fiHarm[0],task->fiHarm[1]); }
                else { c[0] = FourDiff_PtA_PtA(0,0,0,0); c[1] = FourDiff_PtA_PtA(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                break;
            case 3 :
                if(fCorrOrder==2) { AliWarning("CalculateGFW: PtRef for 2pc is just dn{2}!"); return; }
                else { c[0] = Four_2Diff_2Ref(0,0,0,0); c[1] = Four_2Diff_2Ref(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]); }
                break;
            case 4 :
                if(fCorrOrder==2) { c[0] = TwoDiff_PtA_PtB(0,0); c[1] = TwoDiff_PtA_PtB(task->fiHarm[0],task->fiHarm[1]); }
                else { c[0] = FourDiff_PtA_PtB(0,0,0,0); c[1] = FourDiff_PtA_PtB(task->fiHarm[0],task->fiHarm[2],task->fiHarm[1],task->fiHarm[3]);}
                break;   
        }
    }

    return;
}
void AliAnalysisDecorrTask::FillProfiles(TComplex c[2], const double &lCent, TString suffix)
{
    double d[2], val;
    d[0] = c[0].Re();
    d[1] = c[1].Re();
    val = 0.0;
    bool fill = kFALSE;

    if(d[0] > 0.0) {fill = kTRUE; val = d[1]/d[0]; }
    if(fill && TMath::Abs(val) > 1.0) { fill = kFALSE; }
    if(!fill) { return; }
    TProfile* prof = (TProfile*)lSampleList->FindObject(Form("%s%s",fCorrName.Data(),suffix.Data()));
    if(!prof) { AliError(Form("Profile %s%s not found",fCorrName.Data(),suffix.Data())); return; }
    prof->Fill(lCent,val,d[0]);
    //refData[fTaskCounter].val = val;
    if(fSampling) 
    {
        TProfile* prof_sample = (TProfile*)lSampleList->FindObject(Form("%s%s_sample%d",fCorrName.Data(),suffix.Data(),fIndexSampling));
        if(!prof_sample) { AliError(Form("Profile %s%s_sample%d not found",fCorrName.Data(),suffix.Data(),fIndexSampling)); return; }
        prof_sample->Fill(lCent,val,d[0]);
    }
    return;
}
void AliAnalysisDecorrTask::FillProfiles(TComplex c[2], const double &lCent, const double &lpta, TString suffix)
{
    int iPtA = fPtAxis->FindBin(lpta);
    if(iPtA<1 || iPtA > NPtBin) return;
    --iPtA;
    double d[2], val;
    d[0] = c[0].Re();
    d[1] = c[1].Re();
    val = 0.0;
    bool fill = kFALSE;
    if(d[0] > 0.0) {fill = kTRUE; val = d[1]/d[0]; } 
    if(fill && TMath::Abs(val) > 1.0) { fill = kFALSE; }
    if(!fill) { return; }
    TProfile2D* prof = (TProfile2D*)lSampleList->FindObject(Form("%s%s",fCorrName.Data(),suffix.Data()));
    if(!prof) { AliError(Form("Profile %s%s not found",fCorrName.Data(),suffix.Data())); return; }
    prof->Fill(lCent,lpta,val,d[0]);
    
    //if(suffix.Contains("_diff")) ptaData[fTaskCounter].val[0][iPtA] = val;
    //if(suffix.Contains("_PtA")) ptaData[fTaskCounter].val[1][iPtA] = val;
    //if(suffix.Contains("_PtRef")) ptaData[fTaskCounter].val[2][iPtA] = val;
    if(fSampling) 
    {
        TProfile2D* prof_sample = (TProfile2D*)lSampleList->FindObject(Form("%s%s_sample%d",fCorrName.Data(),suffix.Data(),fIndexSampling));
        if(!prof_sample) { AliError(Form("Profile %s%s_sample%d not found",fCorrName.Data(),suffix.Data(),fIndexSampling)); return; }
        prof_sample->Fill(lCent,lpta,val,d[0]);
    }
    return;
}
void AliAnalysisDecorrTask::FillProfiles(TComplex c[2], const double &lCent, const double &lpta, const double lptb, TString suffix)
{
    int iPtA = fPtAxis->FindBin(lpta);
    int iPtB = fPtAxis->FindBin(lptb);
    if(iPtA <1 || iPtA > NPtBin || iPtB <1 || iPtB > NPtBin) return;
    --iPtA;
    --iPtB;
    double d[2], val;
    d[0] = c[0].Re();
    d[1] = c[1].Re();
    val = 0.0;
    bool fill = kFALSE;

    if(d[0] > 0.0) {fill = kTRUE; val = d[1]/d[0]; }
    if(fill && TMath::Abs(val) > 1.0) { fill = kFALSE; }
    if(!fill) { return; }
    TProfile3D* prof = (TProfile3D*)lSampleList->FindObject(Form("%s%s",fCorrName.Data(),suffix.Data()));
    if(!prof) { AliError(Form("Profile %s%s not found",fCorrName.Data(),suffix.Data())); return; }
    prof->Fill(lCent,lpta,lptb,val,d[0]);
    //ptbData[fTaskCounter].val[iPtA][iPtB] = val;
    if(fSampling) 
    {
        TProfile3D* prof_sample = (TProfile3D*)lSampleList->FindObject(Form("%s%s_sample%d",fCorrName.Data(),suffix.Data(),fIndexSampling));
        if(!prof_sample) { AliError(Form("Profile %s%s_sample%d not found",fCorrName.Data(),suffix.Data(),fIndexSampling)); return; }
        prof_sample->Fill(lCent,lpta,lptb,val,d[0]);
    }
    return;
}
void AliAnalysisDecorrTask::InitTypes(const AliDecorrFlowCorrTask* const task)
{
    bRef = task->fbDoRef;
    bDiff = task->fbDoDiff;
    bPtA = task->fbDoPtA;
    bPtRef = task->fbDoPtRef; 
    bPtB = task->fbDoPtB; 
    return;
}
void AliAnalysisDecorrTask::FillRPvectors(const AliDecorrFlowCorrTask* const task)
{
    ResetFlowVector(Qvector);
    ResetFlowVector(Qvector10P);
    ResetFlowVector(Qvector10M);

    bool bIsRP; 
    bool bHasGap;
    double dGap = task->fdGaps[0]; 
    if(dGap > -1.0) { bHasGap = kTRUE; } else { bHasGap = kFALSE; }
    double dEtaLimit = 0.5*dGap;
    double dVz = (fOnTheFly)?0.0:fAOD->GetPrimaryVertex()->GetZ();
    
    if(fOnTheFly)
    {
        Int_t iNumTracks = fMCEvent->GetNumberOfPrimaries();
        if(iNumTracks < 1) { return; }
        for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++) 
        {
            AliMCParticle* track = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(iTrack));

            if(!track) { continue; }    
            bIsRP = IsWithinRP(track);
            if (!bIsRP) { continue; }

            if(bOnlyPrimariesAndCh && (!track->IsPhysicalPrimary() || track->Charge())) continue;
            if(bUseLikeSign && (iSign <0)?(track->Charge() >= 0):(track->Charge() <= 0)) continue;
            hCharge->Fill(track->Charge());

            double dPhi = track->Phi();
            double dEta = track->Eta();

            //Calculating weights    
            double dWeight = 1.0;
            
            //Filling Q-vectors for RPs
            if(!bHasGap)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++) 
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        Qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    } //End for iPower
                }  //End for iHarm
            }
            // RFP in positive and negative eta acceptance
            if(dEta > dEtaLimit && bHasGap)
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
            else if(dEta < -dEtaLimit && bHasGap)
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
    }
    else if(fIsMC)
    {
        TClonesArray* tca = (TClonesArray*)fAOD->FindListObject("mcparticles");
        Int_t nPrim = tca->GetEntries();
        AliAODMCParticle* lpart;
        if(nPrim < 1 ) { return; }
        for(Int_t iPart(0); iPart < nPrim; iPart++) 
        {
            lpart = (AliAODMCParticle*)tca->At(iPart);
            if(!lpart->IsPhysicalPrimary()) continue;
            if(lpart->Charge()==0) continue;
            double dPhi = lpart->Phi();
            double dEta = lpart->Eta();
            double dPt = lpart->Pt();
            if (TMath::Abs(dEta) > 0.8) continue;
            
            if (dPt < fRFPsPtMin || dPt > fRFPsPtMax) continue;
            if(fRedTracks && iPart > fTrackprevent*nPrim) continue;
            if(fAbsEtaMax > 0.0 && Abs(dEta) > fAbsEtaMax) continue;
            //Calculating weights
            double dWeight = 1.0;
            
            //Filling Q-vectors for RPs
            if(!bHasGap)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++) 
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        Qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    } //End for iPower
                }  //End for iHarm
            }
            // RFP in positive and negative eta acceptance
            if(dEta > dEtaLimit && bHasGap)
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
            else if(dEta < -dEtaLimit && bHasGap)
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
        }
    }
    else 
    {
        int NTracks(fAOD->GetNumberOfTracks());
        if(NTracks < 1 ) { return; }
        for(Int_t iTrack(0); iTrack < NTracks; iTrack++) 
        {
            AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
            if(!track || !IsTrackSelected(track)) { continue; }
            if(fRedTracks && iTrack > fTrackprevent*NTracks) { continue; }
            bIsRP = IsWithinRP(track);
            if (!bIsRP) { continue; }

            double dPhi = track->Phi();
            double dEta = track->Eta();
            //double dPt = track->Pt();

            //Calculating weights    
            double dWeight = fUseWeightsOne?1.0:GetWeights(dPhi, dEta, dVz);
            if(dWeight <= 0.0) { dWeight = 1.0; }
            
            //Filling Q-vectors for RPs
            if(!bHasGap)
            {
                for(Int_t iHarm(0); iHarm < fNumHarms; iHarm++) 
                {
                    for(Int_t iPower(0); iPower < fNumPowers; iPower++)
                    {
                        Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
                        Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
                        Qvector[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    } //End for iPower
                }  //End for iHarm
            }
            // RFP in positive and negative eta acceptance
            if(dEta > dEtaLimit && bHasGap)
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
            else if(dEta < -dEtaLimit && bHasGap)
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
    }
    
    return;
}
int AliAnalysisDecorrTask::FillPOIvectors(const AliDecorrFlowCorrTask* const task, const double dPtLow, const double dPtHigh)
{

    ResetFlowVector(pvector);
    ResetFlowVector(pvector10M);
    ResetFlowVector(pvector10P);
    ResetFlowVector(qvector);
 
    bool bHasGap;
    double dGap = task->fdGaps[0]; 
    if(dGap > -1.0) { bHasGap = kTRUE; } else { bHasGap = kFALSE; }
    double dEtaLimit = 0.5*dGap;

    double dVz = (fOnTheFly)?0.0:fAOD->GetPrimaryVertex()->GetZ();
    Int_t TrackCounter = 0;
    if(fOnTheFly)
    {
        Int_t iNumTracks = fMCEvent->GetNumberOfPrimaries();
        if(iNumTracks < 1) { return 0; }
        for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++) 
        {
            AliMCParticle* track = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(iTrack));
            if(!track) { continue; }
            if(bOnlyPrimariesAndCh && (!track->IsPhysicalPrimary() || track->Charge())) continue;
            if(bUseLikeSign && (iSign <0)?(track->Charge() >= 0):(track->Charge() <= 0)) continue;

            double dPt = track->Pt();
            double dPhi = track->Phi();
            double dEta = track->Eta();

            //Check for overlap with RPs
            bool bIsWithinRP = IsWithinRP(track);

            //Check that particle is in POI range
            bool bIsWithinPOI = IsWithinPOI(track);

            if(!bIsWithinPOI) { continue; }

            //Load weights
            double dWeight = 1.0;

            //POI with no eta gap
            if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
            {
                TrackCounter++;
                if(!bHasGap)
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
                }
                //POI with eta gap
                if(dEta > dEtaLimit && bHasGap)
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
                if (dEta < -dEtaLimit && bHasGap)
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
    }
    else if(fIsMC)
    {
        TClonesArray* tca = (TClonesArray*)fAOD->FindListObject("mcparticles");
        Int_t nPrim = tca->GetEntries();
        AliAODMCParticle* lpart;
        if(nPrim < 1 ) { return 0; }

        for(Int_t iPart(0); iPart < nPrim; iPart++) 
        {
            lpart = (AliAODMCParticle*)tca->At(iPart);
            if(!lpart->IsPhysicalPrimary()) continue;
            if(lpart->Charge()==0) continue;
            double dPt = lpart->Pt();
            double dPhi = lpart->Phi();
            double dEta = lpart->Eta();
        
            if(fRedTracks && iPart > fTrackprevent*nPrim) continue;
            if(fAbsEtaMax > 0.0 && Abs(dEta) > fAbsEtaMax) continue;

            //Check for overlap with RPs
            bool bIsWithinRP = kFALSE;
            if (dPt > fRFPsPtMin || dPt < fRFPsPtMax) bIsWithinRP = kTRUE;
            

            //Check that particle is in POI range
            bool bIsWithinPOI = kFALSE; 
            if(dPt > fPOIsPtmin || dPt < fPOIsPtmax) bIsWithinPOI = kTRUE;

            if(!bIsWithinPOI) { continue; }

            //Load weights
            double dWeight = 1.0;

            //POI with no eta gap
            if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
            {
                TrackCounter++;
                if(!bHasGap)
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
                }
                //POI with eta gap
                if(dEta > dEtaLimit && bHasGap)
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
                if (dEta < -dEtaLimit && bHasGap)
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
        
    }
    else
    {
        int NTracks(fAOD->GetNumberOfTracks());
        if(NTracks < 1) { return 0; }
        for(int iTrack(0); iTrack < NTracks; ++iTrack)
        {   
            AliAODTrack* POItrack = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
            if(!POItrack || !IsTrackSelected(POItrack)) { continue; }
            if(fRedTracks && iTrack > fTrackprevent*NTracks) { continue; }
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
            double dWeight = fUseWeightsOne?1.0:GetWeights(dPhi, dEta, dVz);
            if(dWeight <= 0.0) { dWeight = 1.0; }

            //POI with no eta gap
            if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
            {
                TrackCounter++;
                if(!bHasGap)
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
                }
                //POI with eta gap
                if(dEta > dEtaLimit && bHasGap)
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
                if (dEta < -dEtaLimit && bHasGap)
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
    }
    
    return TrackCounter;
}
void AliAnalysisDecorrTask::FillPtBvectors(const AliDecorrFlowCorrTask* const task, const double dPtLow, const double dPtHigh)
{

    ResetFlowVector(pvectorPtB);
    ResetFlowVector(pvectorPtB10M);
    ResetFlowVector(pvectorPtB10P);
    ResetFlowVector(qvectorPtB);


    bool bHasGap;
    double dGap = task->fdGaps[0]; 
    if(dGap > -1.0) { bHasGap = kTRUE; } else { bHasGap = kFALSE; }
    double dEtaLimit = 0.5*dGap;

    double dVz = (fOnTheFly)?0.0:fAOD->GetPrimaryVertex()->GetZ();
    if(fOnTheFly)
    {
        Int_t iNumTracks = fMCEvent->GetNumberOfPrimaries();
        if(iNumTracks < 1) { return; }
        for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++) 
        {
            AliMCParticle* track = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(iTrack));
            if(!track) { continue; }
            if(bOnlyPrimariesAndCh && (!track->IsPhysicalPrimary() || track->Charge())) continue;
            if(bUseLikeSign && (iSign <0)?(track->Charge() >= 0):(track->Charge() <= 0)) continue;

            double dPt = track->Pt();
            double dPhi = track->Phi();
            double dEta = track->Eta();

            //if (Abs(dEta) < dEtaLimit) { continue; } Why the hell is this here?

            //Check for overlap with RPs
            bool bIsWithinRP = IsWithinRP(track);

            //Check that particle is in POI range
            bool bIsWithinPOI = IsWithinPOI(track);

            if(!bIsWithinPOI) { continue; }

            //Load weights
            double dWeight = 1.0;

            //POI with no eta gap
            if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
            {
                if(!bHasGap)
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
                }
                //POI with eta gap
                if(dEta > dEtaLimit && bHasGap)
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
                if (dEta < -dEtaLimit && bHasGap)
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
    }
    else if(fIsMC)
    {
        TClonesArray* tca = (TClonesArray*)fAOD->FindListObject("mcparticles");
        Int_t nPrim = tca->GetEntries();
        AliAODMCParticle* lpart;
        if(nPrim < 1 ) { return; }
        for(Int_t iPart(0); iPart < nPrim; iPart++) 
        {
            lpart = (AliAODMCParticle*)tca->At(iPart);
            if(!lpart->IsPhysicalPrimary()) continue;
            if(lpart->Charge()==0) continue;

            if(fRedTracks && iPart > fTrackprevent*nPrim) { continue; }
            double dPt = lpart->Pt();
            double dPhi = lpart->Phi();
            double dEta = lpart->Eta();

            //if (Abs(dEta) < dEtaLimit) { continue; } Why the hell is this here?

            //Check for overlap with RPs
            bool bIsWithinRP = kFALSE;
            if (dPt > fRFPsPtMin || dPt < fRFPsPtMax) bIsWithinRP = kTRUE;
            

            //Check that particle is in POI range
            bool bIsWithinPOI = kFALSE; 
            if(dPt > fPOIsPtmin || dPt < fPOIsPtmax) bIsWithinPOI = kTRUE;

            if(!bIsWithinPOI) { continue; }

            //Load weights
            double dWeight = 1.0;

            //POI with no eta gap
            if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
            {
                if(!bHasGap)
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
                }
                //POI with eta gap
                if(dEta > dEtaLimit && bHasGap)
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
                if (dEta < -dEtaLimit && bHasGap)
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
    }
    else
    {
        int NTracks(fAOD->GetNumberOfTracks());
        if(NTracks < 1) { return; }
        for(int iTrack(0); iTrack < NTracks; ++iTrack)
        {   
            AliAODTrack* POItrack = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
            if(!POItrack || !IsTrackSelected(POItrack)) { continue; }
            if(fRedTracks && iTrack > fTrackprevent*NTracks) { continue; }
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
            double dWeight = fUseWeightsOne?1.0:GetWeights(dPhi, dEta, dVz);
            if(dWeight <= 0.0) { dWeight = 1.0; }

            //POI with no eta gap
            if(dPt > dPtLow && dPt <= dPtHigh)      //Added = to <= 
            {
                if(!bHasGap)
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
                }
                //POI with eta gap
                if(dEta > dEtaLimit && bHasGap)
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
                if (dEta < -dEtaLimit && bHasGap)
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
    }
    
 

    return;
}
bool AliAnalysisDecorrTask::IsWithinRP(const AliVParticle* track) const
{
    if(fAbsEtaMax > 0.0 && Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fRFPsPtMax > 0.0 && track->Pt() < fRFPsPtMin) { return kFALSE; }
    if(fRFPsPtMax > 0.0 && track->Pt() > fRFPsPtMax) { return kFALSE; }

    return kTRUE;
}
bool AliAnalysisDecorrTask::IsWithinPOI(const AliVParticle* track) const
{
    if(fAbsEtaMax > 0.0 && Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(fPOIsPtmax > 0.0 && track->Pt() < fPOIsPtmin) { return kFALSE; }
    if(fPOIsPtmin > 0.0 && track->Pt() > fPOIsPtmax) { return kFALSE; }

    return kTRUE;
}
int AliAnalysisDecorrTask::GetSamplingIndex() const
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
double AliAnalysisDecorrTask::getCentrality()
{
  if(fOnTheFly)
  {
    vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
    vector<double>::iterator it = upper_bound(b.begin(),b.end(),fImpactParameterMC);
    double l_cent = (fImpactParameterMC<0)?-1.0:(centralitymap[b[it-b.begin()]]+centralitymap[b[it-b.begin()-1]])/2.0;
    return l_cent;
  }
  else 
  {
    AliMultSelection *l_MultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    if(!l_MultSel) { printf("MultSelection not found\n"); return -1.0; }
    double l_cent = l_MultSel->GetMultiplicityPercentile(fCentEstimator.Data());
    return l_cent;
  }
  return -1.0;
}
AliMCEvent *AliAnalysisDecorrTask::getMCEvent() {
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
  for(auto&& listObject: *ltgen){
    genName = Form("%s",listObject->GetName());
    if (genName.Contains("Hijing")) {
      headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
      break;
      }
    }
  }
  else headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      fImpactParameterMC = headerH->ImpactParameter();
  }
  return ev;
}
bool AliAnalysisDecorrTask::IsEventSelected(Double_t lCent)
{
  fhEventSel->Fill(0);
  if(!fEventCuts.AcceptEvent(fAOD)) return kFALSE;
  fhEventSel->Fill(1);
  if(!CheckTrigger(lCent)) return kFALSE;
  fhEventSel->Fill(2);
  if(lCent < 0) return kFALSE;
  //if(fCentMin > 0 && lCent < fCentMin) return kFALSE; 
  //if(fCentMax > 0 && lCent > fCentMax) return kFALSE; 
  /*
  if(!bIs2018Data)
  {
      if(!(fSelectMask & fTrigger)) { return kFALSE; }
  }
  else
  {
    if((dPercentile<10) || (dPercentile>30 && dPercentile<50)){
      if(!(fSelectMask & (AliVEvent::kCentral|AliVEvent::kSemiCentral|fTrigger))) { return kFALSE; }
    }
    else{
      if(!(fSelectMask & fTrigger)) { return kFALSE;}
    }
  }
  */
  fhEventSel->Fill(3);
  if(lCent > 100 || lCent < 0) { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1; }
  //if(fEventRejectAddPileUp && lCent > 0 && lCent < 10 && IsEventRejectedAddPileUp(fCentralPileupCut)) return kFALSE;
  //else if(fEventRejectAddPileUp && lCent > 10 && IsEventRejectedAddPileUp(fDefaultPileupCut)) return kFALSE;
  fhEventSel->Fill(4);
  if(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ()) > fPVtxCutZ) return kFALSE;
  fhEventSel->Fill(5);
  if(fRequireHighPtTracks)
  {
    int NTracks = fAOD->GetNumberOfTracks();
    int tmp = 0;
    for(int iTrack(0); iTrack < NTracks; ++iTrack)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        if(!track || !IsTrackSelected(track)) continue;
        if(track->Pt() > fHighPtCut) ++tmp;
    }
    if(tmp < fNHighPtTracks) return kFALSE;
  }
  fhEventSel->Fill(6);

  return kTRUE;
}
bool AliAnalysisDecorrTask::CheckTrigger(double lCent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fTrigger&fSelMask)) { return kFALSE; }; 
  if(fSelMask&(fTrigger&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; 
  if((fSelMask&fTrigger&AliVEvent::kCentral) && lCent>10) {return kFALSE; }; 
  if((fSelMask&fTrigger&AliVEvent::kSemiCentral) && (lCent<30 || lCent>50)) {return kFALSE; }; 
  return kTRUE;
}
bool AliAnalysisDecorrTask::IsTrackSelected(const AliAODTrack* track) const
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
      if(fCutDCAxyMax > 0.0 && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0]+dDCAXYZ[1]*dDCAXYZ[1]) > (0.0105+0.0350/pow(track->Pt(),1.1))*fCutDCAxyMax) { return kFALSE; }
  }
  if(bUseLikeSign)
  {
      if(!(track->Charge() == iSign)) { return kFALSE; }
  }
  if(fChi2Cut > 0.0 && track->GetTPCchi2()/track->GetTPCNcls() > fChi2Cut) { return kFALSE; }

  return kTRUE;
}
bool AliAnalysisDecorrTask::IsEventRejectedAddPileUp(const int fPileupCut) const
{
  // Check for additional pile-up rejection in Run 2 Pb-Pb collisions (15o, 17n)
  // based on multiplicity correlations
  // ***************************************************************************

  Bool_t bIs17n = kFALSE;
  Bool_t bIs15o = kFALSE;
  Bool_t bIs18qr = kFALSE;

  Int_t iRunNumber = fAOD->GetRunNumber();
  if(iRunNumber >= 244824 && iRunNumber <= 246994) { bIs15o = kTRUE; }
  else if(iRunNumber == 280235 || iRunNumber == 20234) { bIs17n = kTRUE; }
  else if(iRunNumber >= 295585 && iRunNumber <= 297595 ) { bIs18qr = kTRUE; }
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

  if(bIs18qr)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > fPileupCut) { return kTRUE; }
      
  }

  // QA Plots
  if(fFillQA) {
    fhQAEventsfMult32vsCentr->Fill(v0Centr, multTrk);
    fhQAEventsMult128vsCentr->Fill(v0Centr, multTPC128);
    fhQAEventsfMultTPCvsTOF->Fill(multTPC32, multTOF);
    fhQAEventsfMultTPCvsESD->Fill(multTPC128, multESD);
  }

  return kFALSE;
}
double AliAnalysisDecorrTask::GetWeights(double dPhi, double dEta, double dVz)
{
    double dWeight = 1.0;
    dWeight = fWeights->GetNUA(dPhi, dEta, dVz);
    return dWeight;  
}
double *AliAnalysisDecorrTask::GetBinsFromAxis(TAxis *inax) {
  Int_t lBins = inax->GetNbins();
  Double_t *retBins = new Double_t[lBins+1];
  for(Int_t i=0;i<lBins;i++)
    retBins[i] = inax->GetBinLowEdge(i+1);
  retBins[lBins] = inax->GetBinUpEdge(lBins);
  return retBins;
}

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
TComplex AliAnalysisDecorrTask::QGap10P(int n, int p)
{

	if(n>=0) return Qvector10P[n][p];
  else return TComplex::Conjugate(Qvector10P[-n][p]);

}
TComplex AliAnalysisDecorrTask::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[-n][p]);

}
TComplex AliAnalysisDecorrTask::pGap10M(int n, int p)
{

	if(n>=0) return pvector10M[n][p];
	else return TComplex::Conjugate(pvector10M[-n][p]);

}
TComplex AliAnalysisDecorrTask::pGap10P(int n, int p)
{

	if(n>=0) return pvector10P[n][p];
	else return TComplex::Conjugate(pvector10P[-n][p]);

}
TComplex AliAnalysisDecorrTask::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[-n][p]);

}
TComplex AliAnalysisDecorrTask::qGap10M(int n, int p)
{

    if(n>=0) return pvector10M[n][p];
    else return TComplex::Conjugate(pvector10M[-n][p]);

}
TComplex AliAnalysisDecorrTask::qGap10P(int n, int p)
{

    if(n>=0) return pvector10P[n][p];
    else return TComplex::Conjugate(pvector10P[-n][p]);

}
TComplex AliAnalysisDecorrTask::pPtA(int n, int p)
{

    if(n>=0) return pvector[n][p];
    else return TComplex::Conjugate(pvector[-n][p]);

}
TComplex AliAnalysisDecorrTask::pPtB(int n, int p)
{

    if(n>=0) return pvectorPtB[n][p];
    else return TComplex::Conjugate(pvectorPtB[-n][p]);

}
TComplex AliAnalysisDecorrTask::qPtA(int n, int p)
{

    if(n>=0) return qvector[n][p];
    else return TComplex::Conjugate(qvector[-n][p]);

}
TComplex AliAnalysisDecorrTask::qPtB(int n, int p)
{

    if(n>=0) return qvectorPtB[n][p];
    else return TComplex::Conjugate(qvector[-n][p]);

}
TComplex AliAnalysisDecorrTask::pPtBGap10M(int n, int p)
{

    if(n>=0) return pvectorPtB10M[n][p];
    else return TComplex::Conjugate(pvectorPtB10M[-n][p]);

}
TComplex AliAnalysisDecorrTask::pPtBGap10P(int n, int p)
{

    if(n>=0) return pvectorPtB10P[n][p];
    else return TComplex::Conjugate(pvectorPtB10P[-n][p]);

}
void AliAnalysisDecorrTask::ResetFlowVector(TComplex (&array)[fNumHarms][fNumPowers])
{
  for(Int_t iHarm(0); iHarm < fNumHarms; ++iHarm) {
    for(Int_t iPower(0); iPower < fNumPowers; ++iPower) {
      array[iHarm][iPower](0.0,0.0);
    }
  }
  return;

}
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
//cn{2} in subevents
TComplex AliAnalysisDecorrTask::Two_SubP(int n1, int n2)
{
	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1) - QGap10P(n1+n2,2);
    return formula;
}
TComplex AliAnalysisDecorrTask::Two_SubM(int n1, int n2)
{
	TComplex formula = QGap10M(n1,1)*QGap10M(n2,1) - QGap10M(n1+n2,2);
    return formula;
}

//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff(int n1, int n2)
{
    TComplex formula = p(n1,1)*Q(n2,1) - q(n1+n2,2);
    return formula;
}
//dn{2} in subevents
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubP(int n1, int n2)
{
    TComplex formula = pGap10P(n1,1)*QGap10P(n2,1) - qGap10P(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubM(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*QGap10M(n2,1) - qGap10M(n1+n2,2);
    return formula;
}
//dn{2} with gap between POI and RP
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
//Two particles from same pt in subevents
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubM_PtA(int n1, int n2)
{
    TComplex formula = pGap10M(n1,1)*pGap10M(n2,1) - pGap10M(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubP_PtA(int n1, int n2)
{
    TComplex formula = pGap10P(n1,1)*pGap10P(n2,1) - pGap10P(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubM_PtB(int n1, int n2)
{
    TComplex formula = pPtBGap10M(n1,1)*pPtBGap10M(n2,1) - pPtBGap10M(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubP_PtB(int n1, int n2)
{
    TComplex formula = pPtBGap10P(n1,1)*pPtBGap10P(n2,1) - pPtBGap10P(n1+n2,2);
    return formula;

}
//2 particles from different pt in sub event
TComplex AliAnalysisDecorrTask::TwoDiff_SubM_PtA_PtB(int n1, int n2)
{
    TComplex formula;
    if(bEqualPt) formula = pGap10M(n1,1)*pPtBGap10M(n2,1) - pGap10M(n1+n2,2);
    else formula = pGap10M(n1,1)*pPtBGap10M(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::TwoDiff_SubP_PtA_PtB(int n1, int n2)
{
    TComplex formula;
    if(bEqualPt) formula = pGap10P(n1,1)*pPtBGap10P(n2,1) - pGap10P(n1+n2,2);
    else formula = pGap10P(n1,1)*pPtBGap10P(n2,1);
    return formula;

}
//____________________________________________________________________
//____________________________________________________________________
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
//_____________________________________________________________________
TComplex AliAnalysisDecorrTask::Four_2Diff_2Ref(int n1, int n2, int n3, int n4)
{

    TComplex formula = p(n1,1)*Q(n2,1)*p(n3,1)*Q(n4,1)-q(n1+n2,2)*p(n3,1)*Q(n4,1)-Q(n2,1)*q(n1+n3,2)*Q(n4,1)
    - p(n1,1)*q(n2+n3,2)*Q(n4,1)+2.*q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*p(n3,1)*q(n1+n4,2)
    + q(n2+n3,2)*q(n1+n4,2)-p(n1,1)*p(n3,1)*Q(n2+n4,2)+q(n1+n3,2)*Q(n2+n4,2)
    + 2.*p(n3,1)*q(n1+n2+n4,3)-p(n1,1)*Q(n2,1)*q(n3+n4,2)+q(n1+n2,2)*q(n3+n4,2)
    + 2.*Q(n2,1)*q(n1+n3+n4,3)+2.*p(n1,1)*q(n2+n3+n4,3)-6.*q(n1+n2+n3+n4,4);
    //TComplex *out = (TComplex*) &formula;
    //return out;
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourGapM_2Diff_2Ref(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_SubM_PtA(n1,n2)*Two_SubP(n3,n4);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourGapP_2Diff_2Ref(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_SubP_PtA(n1,n2)*Two_SubM(n3,n4);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisDecorrTask::FourGap_2Diff_2Ref_OS(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_SubM(n1,n2)*TwoDiff_SubP(n3,n4);
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
TComplex AliAnalysisDecorrTask::FourDiffGap10_PtA_PtA(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiffGap10_Pt(n1,n2)*TwoDiffGap10_Pt(n3,n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10M_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_SubM_PtA(n1, n2)*TwoDiff_SubP_PtB(n3, n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10P_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_SubP_PtA(n1, n2)*TwoDiff_SubM_PtB(n3, n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10_OS_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiff_SubP_PtA_PtB(n1, n2)*TwoDiff_SubM_PtA_PtB(n3, n4);
    return formula;
}
/*
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiffGap10_PtA_PtB(n1, n2)*TwoDiffGap10_PtA_PtB(n3, n4);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisDecorrTask::FourDiffGap10_OS_PtA_PtB(int n1, int n2, int n3, int n4)
{
    TComplex formula = TwoDiffGap10_Pt(n1, n2)*TwoDiffGap10_PtB(n3, n4);
    return formula;
}
*/
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

