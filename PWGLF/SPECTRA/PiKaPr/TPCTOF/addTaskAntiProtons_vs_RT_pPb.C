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
AliAnalysisTaskAntiProtons_vs_RT_pPb *addTaskAntiProtons_vs_RT_pPb (
                                                Bool_t    isITSrecalib = kTRUE,
                                                Bool_t    isMC  = kFALSE,
                                                UInt_t    triggerType = (AliVEvent::kINT7),//(AliVEvent::kINT7|AliVEvent::kHighMultV0)
                                                Double_t  multMin = 0.0,
                                                Double_t  multMax = 100.0,
                                                Double_t  pt_min_leading = 5.0,
                                                Double_t  average_Nch_Transv = 13.73)  {
    
    
    //Get nsigmaITS recalibration map
    TFile *inputFile;
    TH2F *hITSnsigma_Mean;
    TH2F *hITSnsigma_Width;

    //temporary files, to be changed as soon as the corrected recalibrated maps are availables
    if(!isMC)
    {
        inputFile = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/nSigmaITS_Recalib/ITSrecalibrationMaps_data.root");
        hITSnsigma_Mean = (TH2F*) inputFile->Get("hITSnsigma_Mean");
        hITSnsigma_Width = (TH2F*) inputFile->Get("hITSnsigma_Width");
    }
    else if(isMC)
    {
        inputFile = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/nSigmaITS_Recalib/ITSrecalibrationMaps_mc.root");
        hITSnsigma_Mean = (TH2F*) inputFile->Get("hITSnsigma_Mean");
        hITSnsigma_Width = (TH2F*) inputFile->Get("hITSnsigma_Width");
    }
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    if ( isMC) fileName += ":antiprotons_RT_MC";
    if (!isMC) fileName += ":antiprotons_RT_Data";
    
    
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
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetMinPtLeadingTrack      (pt_min_leading);
    task -> AliAnalysisTaskAntiProtons_vs_RT_pPb::SetITSRecalibrationMaps   (hITSnsigma_Mean, hITSnsigma_Width);
    mgr  -> AddTask(task);
    mgr  -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr  -> ConnectOutput(task,1,mgr->CreateContainer(results_cont, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr  -> ConnectOutput(task,2,mgr->CreateContainer(qaplots_cont, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
//____________________________________________________________________________________________________________________________________________

