AliAnalysisTaskFullpAJets *AddTaskFullpAJets(const char* proj_name, const Double_t jetRadius=0.4, Bool_t IsMC=kFALSE, const char* track_name="PicoTracks", const char* clus_name="caloClusters", const char* corrclus_name="caloClustersCorr", const char* mcpart_name="MCParticles", const char* Centrality_name="V0A", Double_t scaleFactor = 1.28, Double_t nefJetCut = 1.0, Bool_t doNEF=kFALSE, Bool_t signalTrackBias=kFALSE, Bool_t doTrackQA=kFALSE, Bool_t doClusterQA=kFALSE, Int_t calcRhoJet=0, Bool_t doNEFSignalOnly=kTRUE, Bool_t doVertexRCut=kTRUE, Bool_t isMCParticleLevel=kFALSE)
{
    char *usedTracks = track_name;
    char *usedClusters = clus_name;
    char *outClusName = corrclus_name;
    char *usedMCParticles = mcpart_name;
    char *centEst = Centrality_name;
    char *projName = proj_name;
    const Double_t minTrackPt=0.15;
    const Double_t minClusterPt=0.30;
    const Double_t minMCPartPt=0.00;
    Double_t scaleFac = scaleFactor; // Obtained from previous runs...
    Double_t NEFSignalJetCut = nefJetCut; // Require signal jet to not exceed a Neutral Energy Fraction of this setting...
    
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
    TString listName = Form("List%sR%d",projName,drjet);
    TString fileName = Form("%s:FullpAJets", AliAnalysisManager::GetCommonFileName());
    
    // Jet finders (RECONSTRUCTED DATA)
    TString tmpTaskName("");
    AliEmcalJetTask* jetFinderTask = NULL;

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    
    // Add User Task
    AliAnalysisTaskFullpAJets *task = new AliAnalysisTaskFullpAJets(taskName);

    // Used for physics selection
    task->SetUseAliAnaUtils(kTRUE);
    task->DoVertexRCut(doVertexRCut);

    if (IsMC == kTRUE)
    {
        task->SetTrackName(usedMCParticles);
        task->SetClusterName("");
        task->SetTrackPtCut(minMCPartPt);
        task->SetClusterPtCut(minMCPartPt);
        task->SetMCParticleLevel(isMCParticleLevel);
        
        // ########## CHARGED JETS ##########
        jetFinderTask = AddTaskEmcalJet(usedMCParticles,"",cKT,jetRadius,cCHARGEDJETS,minMCPartPt,minMCPartPt,0.01,1,"Jet");
        task->SetkTChargedJetName(jetFinderTask->GetName());
        
        jetFinderTask = AddTaskEmcalJet(usedMCParticles,"",cANTIKT,jetRadius,cCHARGEDJETS,minMCPartPt,minMCPartPt,0.01,1,"Jet");
        task->SetAkTChargedJetName(jetFinderTask->GetName());
        
        // ########## FULL JETS ##########
        // No Full jets or clusters are used if run over MCParticles!
        task->SetkTFullJetName("");
        task->SetAkTFullJetName("");
    }
    else
    {
        task->SetTrackName(usedTracks);
        task->SetClusterName(outClusName);
        task->SetTrackPtCut(minTrackPt);
        task->SetClusterPtCut(minClusterPt);

        // ########## CHARGED JETS ##########
        jetFinderTask = AddTaskEmcalJet(usedTracks,"",cKT,jetRadius,cCHARGEDJETS,minTrackPt,minClusterPt,0.01,1,"Jet");
        task->SetkTChargedJetName(jetFinderTask->GetName());
        
        jetFinderTask = AddTaskEmcalJet(usedTracks,"",cANTIKT,jetRadius,cCHARGEDJETS,minTrackPt,minClusterPt,0.01,1,"Jet");
        task->SetAkTChargedJetName(jetFinderTask->GetName());
        
        // ########## FULL JETS ##########
        jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cKT,jetRadius,cFULLJETS,minTrackPt,minClusterPt,0.01,1,"Jet");
        task->SetkTFullJetName(jetFinderTask->GetName());

        jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cANTIKT,jetRadius,cFULLJETS,minTrackPt,minClusterPt,0.01,1,"Jet");
        task->SetAkTFullJetName(jetFinderTask->GetName());
    }

    task->SetRjet(drjet);
    task->SetCentralityTag(centEst);
    task->SetScaleFactor(scaleFac);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetNColl(7);
    task->SetNEFSignalJetCut(NEFSignalJetCut);
    task->DoNEFCalibration(doNEF);
    task->DoNEFSignalOnly(doNEFSignalOnly);
    task->SetJetChargeBias(signalTrackBias);
    task->DoTrackQA(doTrackQA);
    task->DoClusterQA(doClusterQA);
    task->CalculateRhoJet(calcRhoJet);
    
    mgr->AddTask(task);

    AliAnalysisDataContainer *coutput = mgr->CreateContainer(listName,TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput);

    return task;
}
