#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "AliAnalysisTaskAntiProtons_vs_RT_pPb.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif


//____________________________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_RT_pPb *AddTaskAntiProtons_vs_RT_pPb (
                                                Bool_t    isITSrecalib = kTRUE,
                                                Bool_t    isMC  = kFALSE,
                                                Bool_t    isPion = kFALSE,
                                                UInt_t    triggerType = (AliVEvent::kINT7),//(AliVEvent::kINT7|AliVEvent::kHighMultV0)
                                                Double_t  multMin = 0.0,
                                                Double_t  multMax = 100.0,
                                                Double_t  pt_min_leading = 5.0,
                                                Double_t  average_Nch_Transv = 13.73)  {
    
    
    //Get nsigmaITS recalibration map
    TFile *inputFile_data = TFile::Open("alien:///alice/cern.ch/user/m/mrasa/ITS_recalib_map_protons/ITSrecalibrationMaps_Data.root");
    TH2F *hITSnsigma_Mean_data = (TH2F*) inputFile_data->Get("hITSnsigma_Mean");
    TH2F *hITSnsigma_Width_data = (TH2F*) inputFile_data->Get("hITSnsigma_Width");

    TFile *inputFile_mc = TFile::Open("alien:///alice/cern.ch/user/m/mrasa/ITS_recalib_map_protons/ITSrecalibrationMaps_MC.root");
    TH2F *hITSnsigma_Mean_mc = (TH2F*) inputFile_mc->Get("hITSnsigma_Mean");
    TH2F *hITSnsigma_Width_mc = (TH2F*) inputFile_mc->Get("hITSnsigma_Width");

    TH2F *hITSnsigma_Mean;
    TH2F *hITSnsigma_Width;
    
    if(!isMC)
    {
        hITSnsigma_Mean = hITSnsigma_Mean_data;
        hITSnsigma_Width = hITSnsigma_Width_data;
    }
    else if(isMC)
    {
        hITSnsigma_Mean = hITSnsigma_Mean_mc;
        hITSnsigma_Width = hITSnsigma_Width_mc;
    }

    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    if ( isMC && !isPion) fileName += ":antiprotons_RT_MC";
    if (!isMC && !isPion) fileName += ":antiprotons_RT_Data";
    if ( isMC && isPion) fileName += ":pions_RT_MC";
    if (!isMC && isPion) fileName += ":pions_RT_Data";

    
    //Container Name
    Int_t index_ptLeadingTrackMin = (Int_t)(10.0*pt_min_leading);
    TString results_cont = Form("Results_ptLeadingTrackMin_%d",index_ptLeadingTrackMin);
    TString qaplots_cont = Form("QAplots_ptLeadingTrackMin_%d",index_ptLeadingTrackMin);

    
    //Analysis Task
    AliAnalysisTaskAntiProtons_vs_RT_pPb *task = new AliAnalysisTaskAntiProtons_vs_RT_pPb ("task_antiprotons_vs_RT");
    task -> SelectCollisionCandidates (triggerType);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetRunningMode            (isITSrecalib);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetTriggerType            (triggerType);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetMultiplicityInterval   (multMin,multMax);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetAverageTransverseMult  (average_Nch_Transv);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetInputData              (isMC);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetParticleType           (isPion);         
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetMinPtLeadingTrack      (pt_min_leading);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetITSRecalibrationMaps   (hITSnsigma_Mean, hITSnsigma_Width);
    mgr  -> AddTask(task);
    mgr  -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr  -> ConnectOutput(task,1,mgr->CreateContainer(results_cont, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr  -> ConnectOutput(task,2,mgr->CreateContainer(qaplots_cont, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
//____________________________________________________________________________________________________________________________________________

