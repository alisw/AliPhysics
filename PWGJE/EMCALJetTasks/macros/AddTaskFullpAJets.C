AliAnalysisTaskFullpAJets *AddTaskFullpAJets(const char* proj_name, const Double_t jetRadius = 0.4, Bool_t IsMC = kFALSE, const char* track_name = "PicoTracks", const char* clus_name = "caloClusters", const char* corrclus_name = "caloClustersCorr", const char* mcpart_name = "MCParticles", TString Centrality_name = "V0A", Double_t scaleFactor = 1.28, Double_t nefJetCut = 1.0, Bool_t doNEF = kFALSE, Bool_t signalTrackBias = kFALSE, Bool_t doTrackQA = kFALSE, Bool_t doClusterQA = kFALSE, Int_t calcRhoJet = 0, Bool_t doNEFSignalOnly = kTRUE, Bool_t doVertexRCut = kTRUE, Bool_t isMCParticleLevel = kFALSE, Double_t jetRAccept = 0.4, Bool_t doTHnSparse = kFALSE, Bool_t doJetRhoDensity = kFALSE, Int_t setTriggerClass = 0, Bool_t do3DHistos = kFALSE)
{
    char *usedTracks = track_name;
    char *usedClusters = clus_name;
    char *outClusName = corrclus_name;
    char *usedMCParticles = mcpart_name;
    char *projName = proj_name;
    const Double_t minTrackPt=0.15;
    const Double_t minClusterPt=0.30;
    const Double_t minMCPartPt=0.00;
    
    TString centEst = Centrality_name;
    
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
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    
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
        jetFinderTask = AddTaskEmcalJet(usedMCParticles,"",cKT,jetRadius,cCHARGEDJETS,minMCPartPt,minMCPartPt,0.005,1,"Jet");
        task->SetkTChargedJetName(jetFinderTask->GetName());
        
        jetFinderTask = AddTaskEmcalJet(usedMCParticles,"",cANTIKT,jetRadius,cCHARGEDJETS,minMCPartPt,minMCPartPt,0.005,1,"Jet");
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
        jetFinderTask = AddTaskEmcalJet(usedTracks,"",cKT,jetRadius,cCHARGEDJETS,minTrackPt,minClusterPt,0.005,1,"Jet");
        task->SetkTChargedJetName(jetFinderTask->GetName());
        
        jetFinderTask = AddTaskEmcalJet(usedTracks,"",cANTIKT,jetRadius,cCHARGEDJETS,minTrackPt,minClusterPt,0.005,1,"Jet");
        task->SetAkTChargedJetName(jetFinderTask->GetName());
        
        // ########## FULL JETS ##########
        jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cKT,jetRadius,cFULLJETS,minTrackPt,minClusterPt,0.005,1,"Jet");
        task->SetkTFullJetName(jetFinderTask->GetName());
        
        jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cANTIKT,jetRadius,cFULLJETS,minTrackPt,minClusterPt,0.005,1,"Jet");
        task->SetAkTFullJetName(jetFinderTask->GetName());
    }
    
    // Set Trigger
    // setTriggerClass == 0 -> MB or kINT7 trigger doesn't need to be set
    if (setTriggerClass == 1) // i.e. EMCal Jet Trigger kJ1, pT,jet > 20 GeV
    {
        task->SetTrigClass("J1");
    }
    else if (setTriggerClass == 2) // i.e. EMCal Jet Trigger kJ2, pT,jet > 10 GeV
    {
        task->SetTrigClass("J2");
    }
    else if (setTriggerClass == 3) // i.e. EMCal Gamma Trigger kG1
    {
        task->SetTrigClass("G1");
    }
    else if (setTriggerClass == 4) // i.e. EMCal Gamma Trigger kG2
    {
        task->SetTrigClass("G2");
    }
    else if (setTriggerClass == 5) // i.e. EMCal Trigger kL0
    {
        task->SetTrigClass("L0");
    }
    else if (setTriggerClass == -1) // i.e. EMCal Jet Trigger not defined
    {
        task->SetTrigClass("ND");
    }
    
    task->SetRjet(drjet);
    task->SetJetRAcceptance(jetRAccept);
    task->SetCentralityTag(centEst.Data());
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
    task->DoTHnSparse(doTHnSparse);
    task->DoJetRhoDensity(doJetRhoDensity);
    task->Do3DPlotting(do3DHistos);
    
    mgr->AddTask(task);
    
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(listName,TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput);
    
    return task;
}
