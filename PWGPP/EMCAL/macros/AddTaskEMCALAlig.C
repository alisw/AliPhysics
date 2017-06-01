/// \file AddTaskEmcalAlig
/// \ingroup EMCALPerformanceMacros
/// \brief Configuration analysis task on alignment checks.
///
/// The configuration of alignment checks analysis
///
/// \author Henrique Zanoli <Henrique.Zanoli@cern.ch>, University of Sao Paulo and Utrecht University
///

AliAnalysisTaskEMCALAlig* AddTaskEmcalAlig(
                                           const char *ntracks            = "usedefault",
                                           const char *nclusters          = "usedefault",
                                           const char* ncells             = "usedefault",
                                           const char *suffix             = "",
                                           Int_t TriggerMode              = 0,
                                           const char *sRunPeriod         = "lhc16k"
                                           )
{
    /* TriggerMode 
     default = INT7
     1 = GA
     2 = JE
     3 = INT7 | GA | JE
     */
    
   UInt_t Selection;
    
    switch (TriggerMode) {
        case 0:
            Selection = AliVEvent::kINT7;
            break;
        case 1:
            Selection = AliVEvent::kEMCEGA;
            break;
        case 2:
            Selection = AliVEvent::kEMCEJE;
            break;
        case 3:
            Selection = AliVEvent::kINT7 | AliVEvent::kEMCEJE | AliVEvent::kEMCEGA;
            break;
            
    }
    
    //AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);

    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalAlig", "No analysis manager to connect to.");
        return 0;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    AliVEventHandler* handler = mgr->GetInputEventHandler();
    if (!handler)
    {
        ::Error("AddTaskEmcalAlig", "This task requires an input event handler");
        return 0;
    }
    
    enum EDataType_t {
        kUnknown,
        kESD,
        kAOD
    };
    
    EDataType_t dataType = kUnknown;
    
    
    if (handler->InheritsFrom("AliESDInputHandler")) {
        dataType = kESD;
    }
    else if (handler->InheritsFrom("AliAODInputHandler")) {
        dataType = kAOD;
    }
    
    //-------------------------------------------------------
    // Init the task and do settings
    //-------------------------------------------------------
    
    TString trackName(ntracks);
    TString clusName(nclusters);
    TString cellName(ncells);
    
    if (trackName == "usedefault") {
        if (dataType == kESD) {
            trackName = "Tracks";
        }
        else if (dataType == kAOD) {
            trackName = "tracks";
        }
        else {
            trackName = "";
        }
    }
    
    if (clusName == "usedefault") {
        if (dataType == kESD) {
            clusName = "CaloClusters";
        }
        else if (dataType == kAOD) {
            clusName = "caloClusters";
        }
        else {
            clusName = "";
        }
    }
    
    if (cellName == "usedefault") {
        if (dataType == kESD) {
            cellName = "EMCALCells";
        }
        else if (dataType == kAOD) {
            cellName = "emcalCells";
        }
        else {
            cellName = "";
        }
    }
    
    TString name("AliAnalysisTaskEMCALAlig");
    if (!trackName.IsNull()) {
        name += "_";
        name += trackName;
    }
    if (!clusName.IsNull()) {
        name += "_";
        name += clusName;
    }
    if (!cellName.IsNull()) {
        name += "_";
        name += cellName;
    }
    if (strcmp(suffix,"") != 0) {
        name += "_";
        name += suffix;
    }
    
    AliAnalysisTaskEMCALAlig* sampleTask = new AliAnalysisTaskEMCALAlig(name);
    sampleTask->SetCaloCellsName(cellName);
    sampleTask->SetVzRange(-10,10);
    
    if (trackName == "mcparticles") {
        AliMCParticleContainer* mcpartCont = sampleTask->AddMCParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
        AliTrackContainer* trackCont = sampleTask->AddTrackContainer(trackName);
    }
    else if (!trackName.IsNull()) {
        sampleTask->AddParticleContainer(trackName);
    }
    sampleTask->AddClusterContainer(clusName);
    
    //Cluster cuts
    sampleTask->GetClusterContainer(0)->SetClusECut(0.);
    sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
    sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
    sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kNonLinCorr);
    
    //Track cuts
    sampleTask->GetParticleContainer(0)->SetParticlePtCut(1.2);
    sampleTask->GetTrackContainer(0)->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
    sampleTask->GetTrackContainer(0)->SetAODFilterBits(AliAODTrack::kTrkGlobalNoDCA);
    sampleTask->GetTrackContainer(0)->SetEtaLimits(-0.7, 0.7);
    
    sampleTask->SetHistoBins(600, 0, 300);
    sampleTask->SelectCollisionCandidates(Selection);
    
    
    //-------------------------------------------------------
    // Final settings, pass to manager and set the containers
    //-------------------------------------------------------
    
    mgr->AddTask(sampleTask);
    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                              TList::Class(),AliAnalysisManager::kOutputContainer,
                                                              Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput  (sampleTask, 0,  cinput1 );
    mgr->ConnectOutput (sampleTask, 1, coutput1 );
    
    return sampleTask;
}
