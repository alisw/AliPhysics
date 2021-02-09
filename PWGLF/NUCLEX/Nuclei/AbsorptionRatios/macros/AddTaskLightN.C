//
//  AddTaskLightN.C
//

#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskLightN.h"
#include "AliLightNTrackCuts.h"
#include "AliLightNEventCuts.h"
#include <TString.h>
#include <TList.h>
#include "TROOT.h"
#include "TSystem.h"
#endif


AliAnalysisTaskLightN* AddTaskLightN(Bool_t isMC, const char* suffix = "" ) {
    
    // the manager is static, so get the existing manager via the static method
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
    if (!mgr){
        printf("No analysis manager to connect to!\n");
        return nullptr;
    }
    // just to see if all went well, check if the input event handler has been
    // connected
    if (!mgr->GetInputEventHandler()) {
        printf("This task requires an input event handler!\n");
        return nullptr;
    }

    AliLightNEventCuts *evtCutsParticle = AliLightNEventCuts::StandardCutsRun1();
    bool DCAPlots=true;
    bool CPAPlots=false;
    bool CombSigma=false;
    bool ContributionSplitting=false;
    bool ContributionSplittingDaug=false;
    if(isMC){
        ContributionSplitting=false;
        ContributionSplittingDaug=false;
    }else{
        ContributionSplitting=false;
        ContributionSplittingDaug=false;
    }
    
    //Do the standard analysis using the suffix: "default" and "0"
    //The other parameters, 1,2,... etc. are for systemactic variations
    
    //Track Cuts Proton
    AliLightNTrackCuts *ProtonTrackCuts = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    ProtonTrackCuts->SetCutCharge(1);
    if(strcmp (suffix,"0") == 0)ProtonTrackCuts->SetPID(AliPID::kProton, 0.7,3.,1e30);
    //Systematic track cuts
    if(strcmp (suffix,"1") == 0)ProtonTrackCuts->SetNClsTPC(60);
    if(strcmp (suffix,"2") == 0)ProtonTrackCuts->SetNClsITS(1);
    if(strcmp (suffix,"3") == 0)ProtonTrackCuts->SetTPCRatioCut(0.5);
    if(strcmp (suffix,"4") == 0)ProtonTrackCuts->SetChi2perNDFCut(6);
    if(strcmp (suffix,"5") == 0)ProtonTrackCuts->SetTPCCrossedRowsCut(50);
    if(strcmp (suffix,"6") == 0)ProtonTrackCuts->SetNClsTPC(80);
    if(strcmp (suffix,"7") == 0)ProtonTrackCuts->SetNClsITS(3);
    if(strcmp (suffix,"8") == 0)ProtonTrackCuts->SetTPCRatioCut(0.9);
    if(strcmp (suffix,"9") == 0)ProtonTrackCuts->SetChi2perNDFCut(3);
    if(strcmp (suffix,"10") == 0)ProtonTrackCuts->SetTPCCrossedRowsCut(110);
    if(strcmp (suffix,"11") == 0){
        ProtonTrackCuts->SetNClsTPC(60);
        ProtonTrackCuts->SetNClsITS(1);
    }
    if(strcmp (suffix,"12") == 0){
        ProtonTrackCuts->SetTPCRatioCut(0.5);
        ProtonTrackCuts->SetChi2perNDFCut(6);
    }
    if(strcmp (suffix,"13") == 0){
        ProtonTrackCuts->SetTPCCrossedRowsCut(50);
        ProtonTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"14") == 0){
        ProtonTrackCuts->SetNClsITS(3);
        ProtonTrackCuts->SetTPCRatioCut(0.9);
    }
    if(strcmp (suffix,"15") == 0){
        ProtonTrackCuts->SetChi2perNDFCut(3);
        ProtonTrackCuts->SetTPCCrossedRowsCut(110);
    }
    if(strcmp (suffix,"16") == 0){
        ProtonTrackCuts->SetNClsTPC(60);
        ProtonTrackCuts->SetNClsITS(1);
        ProtonTrackCuts->SetTPCRatioCut(0.5);
    }
    if(strcmp (suffix,"17") == 0){
        ProtonTrackCuts->SetChi2perNDFCut(6);
        ProtonTrackCuts->SetTPCCrossedRowsCut(50);
        ProtonTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"18") == 0){
        ProtonTrackCuts->SetNClsITS(3);
        ProtonTrackCuts->SetTPCRatioCut(0.9);
        ProtonTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"19") == 0){
        ProtonTrackCuts->SetNClsTPC(60);
        ProtonTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"20") == 0){
        ProtonTrackCuts->SetNClsITS(1);
        ProtonTrackCuts->SetTPCCrossedRowsCut(110);
    }
    
    //Systematic Secondary correction
    if(strcmp (suffix,"21") == 0){
        ProtonTrackCuts->SetDCAVtxXY(0.5);
        ProtonTrackCuts->SetDCAVtxZ(1.0);
    }
    if(strcmp (suffix,"22") == 0){
        ProtonTrackCuts->SetDCAVtxXY(0.05);
        ProtonTrackCuts->SetDCAVtxZ(.1);
    }
    
    //Systematic PID cuts
    if(strcmp (suffix,"23") == 0)ProtonTrackCuts->SetPID(AliPID::kProton, 0.7,2.,2.);
    if(strcmp (suffix,"24") == 0)ProtonTrackCuts->SetPID(AliPID::kProton, 0.7,4.,4.);
    
    
    AliLightNTrackCuts *AntiProtonTrackCuts = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiProtonTrackCuts->SetCutCharge(-1);
    if(strcmp (suffix,"0") == 0)AntiProtonTrackCuts->SetPID(AliPID::kProton, 0.7,3.,1e30);
    //Systematic track cuts
    if(strcmp (suffix,"1") == 0)AntiProtonTrackCuts->SetNClsTPC(60);
    if(strcmp (suffix,"2") == 0)AntiProtonTrackCuts->SetNClsITS(1);
    if(strcmp (suffix,"3") == 0)AntiProtonTrackCuts->SetTPCRatioCut(0.5);
    if(strcmp (suffix,"4") == 0)AntiProtonTrackCuts->SetChi2perNDFCut(6);
    if(strcmp (suffix,"5") == 0)AntiProtonTrackCuts->SetTPCCrossedRowsCut(50);
    if(strcmp (suffix,"6") == 0)AntiProtonTrackCuts->SetNClsTPC(80);
    if(strcmp (suffix,"7") == 0)AntiProtonTrackCuts->SetNClsITS(3);
    if(strcmp (suffix,"8") == 0)AntiProtonTrackCuts->SetTPCRatioCut(0.9);
    if(strcmp (suffix,"9") == 0)AntiProtonTrackCuts->SetChi2perNDFCut(3);
    if(strcmp (suffix,"10") == 0)AntiProtonTrackCuts->SetTPCCrossedRowsCut(110);
    if(strcmp (suffix,"11") == 0){
        AntiProtonTrackCuts->SetNClsTPC(60);
        AntiProtonTrackCuts->SetNClsITS(1);
    }
    if(strcmp (suffix,"12") == 0){
        AntiProtonTrackCuts->SetTPCRatioCut(0.5);
        AntiProtonTrackCuts->SetChi2perNDFCut(6);
    }
    if(strcmp (suffix,"13") == 0){
        AntiProtonTrackCuts->SetTPCCrossedRowsCut(50);
        AntiProtonTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"14") == 0){
        AntiProtonTrackCuts->SetNClsITS(3);
        AntiProtonTrackCuts->SetTPCRatioCut(0.9);
    }
    if(strcmp (suffix,"15") == 0){
        AntiProtonTrackCuts->SetChi2perNDFCut(3);
        AntiProtonTrackCuts->SetTPCCrossedRowsCut(110);
    }
    if(strcmp (suffix,"16") == 0){
        AntiProtonTrackCuts->SetNClsTPC(60);
        AntiProtonTrackCuts->SetNClsITS(1);
        AntiProtonTrackCuts->SetTPCRatioCut(0.5);
    }
    if(strcmp (suffix,"17") == 0){
        AntiProtonTrackCuts->SetChi2perNDFCut(6);
        AntiProtonTrackCuts->SetTPCCrossedRowsCut(50);
        AntiProtonTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"18") == 0){
        AntiProtonTrackCuts->SetNClsITS(3);
        AntiProtonTrackCuts->SetTPCRatioCut(0.9);
        AntiProtonTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"19") == 0){
        AntiProtonTrackCuts->SetNClsTPC(60);
        AntiProtonTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"20") == 0){
        AntiProtonTrackCuts->SetNClsITS(1);
        AntiProtonTrackCuts->SetTPCCrossedRowsCut(110);
    }
    
    //Systematic Secondary correction
    if(strcmp (suffix,"21") == 0){
        AntiProtonTrackCuts->SetDCAVtxXY(0.5);
        AntiProtonTrackCuts->SetDCAVtxZ(1.0);
    }
    if(strcmp (suffix,"22") == 0){
        AntiProtonTrackCuts->SetDCAVtxXY(0.05);
        AntiProtonTrackCuts->SetDCAVtxZ(.1);
    }
    
    //Systematic PID cuts
    if(strcmp (suffix,"23") == 0)AntiProtonTrackCuts->SetPID(AliPID::kProton, 0.7,2.,2.);
    if(strcmp (suffix,"24") == 0)ProtonTrackCuts->SetPID(AliPID::kProton, 0.7,4.,4.);
    
    
    //Track Cuts Deuteron
    AliLightNTrackCuts *DeuteronTrackCuts = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    DeuteronTrackCuts->SetCutCharge(1);
    //Systematic track cuts
    if(strcmp (suffix,"1") == 0)DeuteronTrackCuts->SetNClsTPC(60);
    if(strcmp (suffix,"2") == 0)DeuteronTrackCuts->SetNClsITS(1);
    if(strcmp (suffix,"3") == 0)DeuteronTrackCuts->SetTPCRatioCut(0.5);
    if(strcmp (suffix,"4") == 0)DeuteronTrackCuts->SetChi2perNDFCut(6);
    if(strcmp (suffix,"5") == 0)DeuteronTrackCuts->SetTPCCrossedRowsCut(50);
    if(strcmp (suffix,"6") == 0)DeuteronTrackCuts->SetNClsTPC(80);
    if(strcmp (suffix,"7") == 0)DeuteronTrackCuts->SetNClsITS(3);
    if(strcmp (suffix,"8") == 0)DeuteronTrackCuts->SetTPCRatioCut(0.9);
    if(strcmp (suffix,"9") == 0)DeuteronTrackCuts->SetChi2perNDFCut(3);
    if(strcmp (suffix,"10") == 0)DeuteronTrackCuts->SetTPCCrossedRowsCut(110);
    if(strcmp (suffix,"11") == 0){
        DeuteronTrackCuts->SetNClsTPC(60);
        DeuteronTrackCuts->SetNClsITS(1);
    }
    if(strcmp (suffix,"12") == 0){
        DeuteronTrackCuts->SetTPCRatioCut(0.5);
        DeuteronTrackCuts->SetChi2perNDFCut(6);
    }
    if(strcmp (suffix,"13") == 0){
        DeuteronTrackCuts->SetTPCCrossedRowsCut(50);
        DeuteronTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"14") == 0){
        DeuteronTrackCuts->SetNClsITS(3);
        DeuteronTrackCuts->SetTPCRatioCut(0.9);
    }
    if(strcmp (suffix,"15") == 0){
        DeuteronTrackCuts->SetChi2perNDFCut(3);
        DeuteronTrackCuts->SetTPCCrossedRowsCut(110);
    }
    if(strcmp (suffix,"16") == 0){
        DeuteronTrackCuts->SetNClsTPC(60);
        DeuteronTrackCuts->SetNClsITS(1);
        DeuteronTrackCuts->SetTPCRatioCut(0.5);
    }
    if(strcmp (suffix,"17") == 0){
        DeuteronTrackCuts->SetChi2perNDFCut(6);
        DeuteronTrackCuts->SetTPCCrossedRowsCut(50);
        DeuteronTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"18") == 0){
        DeuteronTrackCuts->SetNClsITS(3);
        DeuteronTrackCuts->SetTPCRatioCut(0.9);
        DeuteronTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"19") == 0){
        DeuteronTrackCuts->SetNClsTPC(60);
        DeuteronTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"20") == 0){
        DeuteronTrackCuts->SetNClsITS(1);
        DeuteronTrackCuts->SetTPCCrossedRowsCut(110);
    }
    
    //Systematic Secondary correction
    if(strcmp (suffix,"21") == 0){
        DeuteronTrackCuts->SetDCAVtxXY(0.5);
        DeuteronTrackCuts->SetDCAVtxZ(1.0);
    }
    if(strcmp (suffix,"22") == 0){
        DeuteronTrackCuts->SetDCAVtxXY(0.05);
        DeuteronTrackCuts->SetDCAVtxZ(.1);
    }

    //Systematic PID cuts
    if(strcmp (suffix,"23") == 0){
        DeuteronTrackCuts->SetPID(AliPID::kDeuteron, 1.4,2.,1e30);
        DeuteronTrackCuts->SetCutITSPID(-1.,1e30,true);
    }
    if(strcmp (suffix,"24") == 0){
        DeuteronTrackCuts->SetPID(AliPID::kDeuteron, 1.4,4.,1e30);
        DeuteronTrackCuts->SetCutITSPID(-3.,1e30,true);
    }
    
    //No ITS PID
    if(strcmp (suffix,"25") == 0)DeuteronTrackCuts->SetCutITSPID(-2.,1e30,false);
    
    //ITS PID investigation
    if(strcmp (suffix,"26") == 0){
        DeuteronTrackCuts->SetPtRange(0.1,1.4);
        DeuteronTrackCuts->SetCutITSPID(-2.,1e30,true);
    }
    if(strcmp (suffix,"27") == 0){
        DeuteronTrackCuts->SetPtRange(0.1,1.4);
        DeuteronTrackCuts->SetCutITSPID(-15.,-2.,true);
    }
    
    
    AliLightNTrackCuts *AntiDeuteronTrackCuts = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiDeuteronTrackCuts->SetCutCharge(-1);
    //Systematic track cuts
    if(strcmp (suffix,"1") == 0)AntiDeuteronTrackCuts->SetNClsTPC(60);
    if(strcmp (suffix,"2") == 0)AntiDeuteronTrackCuts->SetNClsITS(1);
    if(strcmp (suffix,"3") == 0)AntiDeuteronTrackCuts->SetTPCRatioCut(0.5);
    if(strcmp (suffix,"4") == 0)AntiDeuteronTrackCuts->SetChi2perNDFCut(6);
    if(strcmp (suffix,"5") == 0)AntiDeuteronTrackCuts->SetTPCCrossedRowsCut(50);
    if(strcmp (suffix,"6") == 0)AntiDeuteronTrackCuts->SetNClsTPC(80);
    if(strcmp (suffix,"7") == 0)AntiDeuteronTrackCuts->SetNClsITS(3);
    if(strcmp (suffix,"8") == 0)AntiDeuteronTrackCuts->SetTPCRatioCut(0.9);
    if(strcmp (suffix,"9") == 0)AntiDeuteronTrackCuts->SetChi2perNDFCut(3);
    if(strcmp (suffix,"10") == 0)AntiDeuteronTrackCuts->SetTPCCrossedRowsCut(110);
    if(strcmp (suffix,"11") == 0){
        AntiDeuteronTrackCuts->SetNClsTPC(60);
        AntiDeuteronTrackCuts->SetNClsITS(1);
    }
    if(strcmp (suffix,"12") == 0){
        AntiDeuteronTrackCuts->SetTPCRatioCut(0.5);
        AntiDeuteronTrackCuts->SetChi2perNDFCut(6);
    }
    if(strcmp (suffix,"13") == 0){
        AntiDeuteronTrackCuts->SetTPCCrossedRowsCut(50);
        AntiDeuteronTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"14") == 0){
        AntiDeuteronTrackCuts->SetNClsITS(3);
        AntiDeuteronTrackCuts->SetTPCRatioCut(0.9);
    }
    if(strcmp (suffix,"15") == 0){
        AntiDeuteronTrackCuts->SetChi2perNDFCut(3);
        AntiDeuteronTrackCuts->SetTPCCrossedRowsCut(110);
    }
    if(strcmp (suffix,"16") == 0){
        AntiDeuteronTrackCuts->SetNClsTPC(60);
        AntiDeuteronTrackCuts->SetNClsITS(1);
        AntiDeuteronTrackCuts->SetTPCRatioCut(0.5);
    }
    if(strcmp (suffix,"17") == 0){
        AntiDeuteronTrackCuts->SetChi2perNDFCut(6);
        AntiDeuteronTrackCuts->SetTPCCrossedRowsCut(50);
        AntiDeuteronTrackCuts->SetNClsTPC(80);
    }
    if(strcmp (suffix,"18") == 0){
        AntiDeuteronTrackCuts->SetNClsITS(3);
        AntiDeuteronTrackCuts->SetTPCRatioCut(0.9);
        AntiDeuteronTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"19") == 0){
        AntiDeuteronTrackCuts->SetNClsTPC(60);
        AntiDeuteronTrackCuts->SetChi2perNDFCut(3);
    }
    if(strcmp (suffix,"20") == 0){
        AntiDeuteronTrackCuts->SetNClsITS(1);
        AntiDeuteronTrackCuts->SetTPCCrossedRowsCut(110);
    }
    
    //Systematic Secondary correction
    if(strcmp (suffix,"21") == 0){
        AntiDeuteronTrackCuts->SetDCAVtxXY(0.5);
        AntiDeuteronTrackCuts->SetDCAVtxZ(1.0);
    }
    if(strcmp (suffix,"22") == 0){
        AntiDeuteronTrackCuts->SetDCAVtxXY(0.05);
        AntiDeuteronTrackCuts->SetDCAVtxZ(.1);
    }
    
    //Systematic PID cuts
    if(strcmp (suffix,"23") == 0){
        AntiDeuteronTrackCuts->SetPID(AliPID::kDeuteron, 1.4,2.,1e30);
        AntiDeuteronTrackCuts->SetCutITSPID(-1.,1e30,true);
    }
    if(strcmp (suffix,"24") == 0){
        AntiDeuteronTrackCuts->SetPID(AliPID::kDeuteron, 1.4,4.,1e30);
        AntiDeuteronTrackCuts->SetCutITSPID(-3.,1e30,true);
    }
    
    //No ITS PID
    if(strcmp (suffix,"25") == 0)AntiDeuteronTrackCuts->SetCutITSPID(-2.,1e30,false);
    
    
    TString TaskName = Form("LightN_%s",suffix);
    AliAnalysisTaskLightN *task=new AliAnalysisTaskLightN(TaskName.Data(),isMC);
    
    /*if(CentEst == "kInt7"){
        task->SelectCollisionCandidates(AliVEvent::kINT7);
        task->SetMVPileUp(kTRUE);
        }
    }else if(CentEst == "kMB"){
        task->SelectCollisionCandidates(AliVEvent::kMB);
        task->SetMVPileUp(kFALSE);
    }else if(CentEst == "kHighMultV0"){
        task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
        task->SetMVPileUp(kFALSE);
    }*/

    task->SetDebugLevel(0);
    task->SetEvtCutQA(true);
    task->SetTrackBufferSize(2500);
    task->SetEventCutsParticle(evtCutsParticle);
    task->SetTrackCutsProton(ProtonTrackCuts);
    task->SetAntiTrackCutsProton(AntiProtonTrackCuts);
    task->SetTrackCutsDeuteron(DeuteronTrackCuts);
    task->SetAntiTrackCutsDeuteron(AntiDeuteronTrackCuts);

    mgr->AddTask(task);
    TString file = AliAnalysisManager::GetCommonFileName();
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(task, 0, cinput);

    AliAnalysisDataContainer *coutputQA;
    TString QAName = Form("QA_%s",suffix);
    coutputQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                     QAName.Data(), TList::Class(),
                                     AliAnalysisManager::kOutputContainer,
                                     Form("%s:%s", file.Data(), QAName.Data()));
    mgr->ConnectOutput(task, 1, coutputQA);
    
    AliAnalysisDataContainer *coutputEvtCuts;
    TString EvtCutsName = Form("EvtCuts_%s",suffix);
    coutputEvtCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                          EvtCutsName.Data(), TList::Class(),
                                          AliAnalysisManager::kOutputContainer,
                                          Form("%s:%s", file.Data(), EvtCutsName.Data()));
    mgr->ConnectOutput(task, 2, coutputEvtCuts);
    
    AliAnalysisDataContainer *couputTrkCutsProton;
    TString TrackCutsNameProton = Form("ProtonTrackCuts_%s",suffix);
    couputTrkCutsProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                               TrackCutsNameProton.Data(), TList::Class(),
                                               AliAnalysisManager::kOutputContainer,
                                               Form("%s:%s", file.Data(), TrackCutsNameProton.Data()));
    mgr->ConnectOutput(task, 3, couputTrkCutsProton);
    
    AliAnalysisDataContainer *coutputAntiTrkCutsProton;
    TString AntiTrackCutsNameProton = Form("AntiProtonTrackCuts_%s",suffix);
    coutputAntiTrkCutsProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                    AntiTrackCutsNameProton.Data(), TList::Class(),
                                                    AliAnalysisManager::kOutputContainer,
                                                    Form("%s:%s", file.Data(), AntiTrackCutsNameProton.Data()));
    mgr->ConnectOutput(task, 4, coutputAntiTrkCutsProton);
    
    AliAnalysisDataContainer *couputTrkCutsDeuteron;
    TString TrackCutsNameDeuteron = Form("DeuteronTrackCuts_%s",suffix);
    couputTrkCutsDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                 TrackCutsNameDeuteron.Data(), TList::Class(),
                                                 AliAnalysisManager::kOutputContainer,
                                                 Form("%s:%s", file.Data(), TrackCutsNameDeuteron.Data()));
    mgr->ConnectOutput(task, 5, couputTrkCutsDeuteron);
    
    AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron;
    TString AntiTrackCutsNameDeuteron = Form("AntiDeuteronTrackCuts_%s",suffix);
    coutputAntiTrkCutsDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                      AntiTrackCutsNameDeuteron.Data(), TList::Class(),
                                                      AliAnalysisManager::kOutputContainer,
                                                      Form("%s:%s", file.Data(), AntiTrackCutsNameDeuteron.Data()));
    mgr->ConnectOutput(task, 6, coutputAntiTrkCutsDeuteron);
    

    
    if (isMC) {
        AliAnalysisDataContainer *coutputTrkCutsMCProton;
        TString TrkCutsMCNameProton = Form("ProtonTrkCutsMC_%s",suffix);
        coutputTrkCutsMCProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                      TrkCutsMCNameProton.Data(), TList::Class(),
                                                      AliAnalysisManager::kOutputContainer,
                                                      Form("%s:%s", file.Data(), TrkCutsMCNameProton.Data()));
        mgr->ConnectOutput(task, 7, coutputTrkCutsMCProton);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsMCProton;
        TString AntiTrkCutsMCNameProton = Form("AntiProtonTrkCutsMC_%s",suffix);
        coutputAntiTrkCutsMCProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                          AntiTrkCutsMCNameProton.Data(), TList::Class(),
                                                          AliAnalysisManager::kOutputContainer,
                                                          Form("%s:%s", file.Data(), AntiTrkCutsMCNameProton.Data()));
        mgr->ConnectOutput(task, 8, coutputAntiTrkCutsMCProton);
        
        AliAnalysisDataContainer *coutputTrkCutsMCDeuteron;
        TString TrkCutsMCNameDeuteron = Form("DeuteronTrkCutsMC_%s",suffix);
        coutputTrkCutsMCDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                        TrkCutsMCNameDeuteron.Data(), TList::Class(),
                                                        AliAnalysisManager::kOutputContainer,
                                                        Form("%s:%s", file.Data(), TrkCutsMCNameDeuteron.Data()));
        mgr->ConnectOutput(task, 9, coutputTrkCutsMCDeuteron);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsMCDeuteron;
        TString AntiTrkCutsMCNameDeuteron = Form("AntiDeuteronTrkCutsMC_%s",suffix);
        coutputAntiTrkCutsMCDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                            AntiTrkCutsMCNameDeuteron.Data(), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s:%s", file.Data(), AntiTrkCutsMCNameDeuteron.Data()));
        mgr->ConnectOutput(task, 10, coutputAntiTrkCutsMCDeuteron);
        
        
    }
    return task;
}
//#endif
