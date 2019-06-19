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
    
    //Track Cuts Proton
    AliLightNTrackCuts *ProtonTrackCuts = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    ProtonTrackCuts->SetCutCharge(1);
    AliLightNTrackCuts *AntiProtonTrackCuts = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiProtonTrackCuts->SetCutCharge(-1);
    
    //Track Cuts Deuteron
    AliLightNTrackCuts *DeuteronTrackCuts = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    DeuteronTrackCuts->SetCutCharge(1);
    AliLightNTrackCuts *AntiDeuteronTrackCuts = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiDeuteronTrackCuts->SetCutCharge(-1);
    
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
