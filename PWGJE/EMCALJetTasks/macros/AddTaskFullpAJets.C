AliAnalysisTaskFullpAJets *AddTaskFullpAJets(const Double_t jetRadius=0.4)
{
    const char *usedTracks = "PicoTracks";
    const char *usedClusters = "caloClusters";
    const char *outClusName = "caloClustersCorr";
    const char *centEst = "V0A";
    const Double_t minTrackPt=0.15;
    const Double_t minClusterPt=0.30;
    Double_t scaleFactor=1.45; // Obtained from previous runs...
    
    // Some constants for the jet finders
    const Int_t cKT=0;
    const Int_t cANTIKT=1;
    const Int_t cFULLJETS=0;
    const Int_t cCHARGEDJETS=1;
    const Int_t cNEUTRALJETS=2;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        Error("AddTaskJetCommon","No analysis manager found.");
        return 0;
    }

    // Determine the int of the jet radius for naming purposes
    Int_t drjet=Int_t(100*jetRadius);
    if (drjet%10 == 0)
    {
        drjet/=10;
    }
    
    TString taskName = Form("AnalysisFullpAJetsR%d",drjet);
    TString listName = Form("ListR%d",drjet);
    TString fileName = Form("%s:FullpAJets", AliAnalysisManager::GetCommonFileName());
    
    // Jet finders (RECONSTRUCTED DATA)
    TString tmpTaskName("");
    AliEmcalJetTask* jetFinderTask = NULL;

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    
    // ########## CHARGED JETS ##########
    jetFinderTask = AddTaskEmcalJet(usedTracks,"",cANTIKT,jetRadius,cCHARGEDJETS,minTrackPt,minClusterPt);

    jetFinderTask = AddTaskEmcalJet(usedTracks,"",cKT,jetRadius,cCHARGEDJETS,minTrackPt,minClusterPt);

    // ########## FULL JETS ##########
    // last two settings are for min pt tracks/clusters
    // anti-kT
    jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cANTIKT,jetRadius,cFULLJETS,minTrackPt,minClusterPt);

    // kT
    jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cKT,jetRadius,cFULLJETS,minTrackPt,minClusterPt);

    // Add User Task
    AliAnalysisTaskFullpAJets *task = new AliAnalysisTaskFullpAJets(taskName);
    task->SetRjet(drjet);
    task->SetCentralityTag(centEst);
    task->SetScaleFactor(scaleFactor);
    task->SetTrackPtCut(minTrackPt);
    task->SetClusterPtCut(minClusterPt);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    
    mgr->AddTask(task);

    AliAnalysisDataContainer *coutput = mgr->CreateContainer(listName,TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput);

    return task;

}
