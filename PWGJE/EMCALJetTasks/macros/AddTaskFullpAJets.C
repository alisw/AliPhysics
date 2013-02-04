// AddTaskFullpAJets.C 2013-02-03 16:58:39Z yaldo

void AddTaskFullpAJets()
{
    const char *usedTracks="PicoTracks";
    const char *usedClusters="CaloClusters";
    const char *outClusName="CaloClustersCorr";
    const Double_t hadcorr=2.0;
    const Double_t minTrackPt=0.15;
    const Double_t minClusterPt=0.30;
    const Double_t minChargedJetPt=0.15;
    const Double_t minFullJetPt=0.15;
    const Double_t Eexcl=0.00;
    const Double_t phiMatch=0.03;
    const Double_t etaMatch=0.015;
    
    // Some constants for the jet finders
    const Int_t cKT                 = 0;
    const Int_t cANTIKT             = 1;
    const Int_t cFULLJETS           = 0;
    const Int_t cCHARGEDJETS        = 1;
    const Int_t cNEUTRALJETS        = 2;
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        Error("AddTaskJetCommon","No analysis manager found.");
        return 0;
    }

    // Jet finders (RECONSTRUCTED DATA)
    TString tmpTaskName("");
    AliEmcalJetTask* jetFinderTask = NULL;
    
    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    // ########## CHARGED JETS ##########
    // R=0.2
    jetFinderTask = AddTaskEmcalJet(usedTracks,"",cANTIKT,0.2,1,minTrackPt,minClusterPt);
    RequestMemory(jetFinderTask,250*1024);//more memory

    // R=0.4
    jetFinderTask = AddTaskEmcalJet(usedTracks,"",cANTIKT,0.4,1,minTrackPt,minClusterPt);
    RequestMemory(jetFinderTask,250*1024);//more memory

    // ########## FULL JETS ##########
    // last two settings are for min pt tracks/clusters
    // R=0.2, anti-kT
    jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cANTIKT,0.2,cFULLJETS,minTrackPt,minClusterPt);
    RequestMemory(jetFinderTask,250*1024);//more memory

    // R=0.2 kT
    jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cKT,0.2,cFULLJETS,minTrackPt,minClusterPt);
    RequestMemory(jetFinderTask,250*1024);//more memory

    // R=0.4
    jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cANTIKT,0.4,cFULLJETS,minTrackPt,minClusterPt);
    RequestMemory(jetFinderTask,250*1024);//more memory
    
    // R=0.4 kT
    jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cKT,0.4,cFULLJETS,minTrackPt,minClusterPt);
    RequestMemory(jetFinderTask,250*1024);//more memory

    // Add User Tasks'
    // Run with R=0.2
    AliAnalysisTaskFullpAJets *task1 = new AliAnalysisTaskFullpAJets("FileR2");
    mgr->AddTask(task1);
    task1->SetR_JET(2);
    AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("R2List",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "FullpAJetsR2.root");
    mgr->ConnectInput(task1,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task1,1,coutput1);
    RequestMemory(task1,250*1024);//more memory

    // Run with R=0.4
    AliAnalysisTaskFullpAJets *task2 = new AliAnalysisTaskFullpAJets("FileR4");
    mgr->AddTask(task2);
    task2->SetR_JET(4);
    AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("R4List",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "FullpAJetsR4.root");
    mgr->ConnectInput(task2,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task2,1,coutput2);
    RequestMemory(task2,250*1024);//more memory
}
