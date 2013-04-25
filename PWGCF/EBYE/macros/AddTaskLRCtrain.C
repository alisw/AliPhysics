AliAnalysisTaskLRC *AddTaskLRCtrain(
        Bool_t RunKine = kFALSE
        , Bool_t isIonsAnalysis = kFALSE
        , TString strPrefixTaskName = "testTask"
        , TString strRunMode = "default"
        , int nEtaWindows = 8
        , int nPhiWindows = 8
        , double etaWinWidth = 0.2
        , double etaWindowStep = 0.1
        , double phiWinWidth = TMath::Pi()/4
        , int nEventPhiRotations = 8
        , double ptMin = 0.3
        , double ptMax = 1.5
        , int multMin = 0
        , int multMax = 100000
        , double vertexZmax = 10.0
        )
{
    gROOT->LoadMacro("$ALICE_ROOT/PWGCF/EBYE/macros/configLRCAnalysis.C");
    //gROOT->LoadMacro("AliAnalysisTaskIA.cxx+g");
    
    // A. Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskLRC", "No analysis manager to connect to.");
        return NULL;
    }
    
    // B. Check the analysis type using the event handlers connected to the analysis
    //    manager. The availability of MC handler cann also be checked here.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        Error("AddTaskLRC", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType();
    cout << " # TaskLRC - input :" << type << "\n";
    
    // ########## Task LRC
    AliAnalysisTaskLRC* taskLRC = 0x0;
    taskLRC = createLRCtaskSkeleton("Task_LRC_" + strPrefixTaskName /*+ "_MB_Global2"*/,RunKine);
    //addAliLRCProcessors(taskLRC,windowsConfigId, nPhiSectors,gapAsPartOfPi);//, whichParticleToFill, whichPidCondition);
    
    tuneEtaPhiWindows( taskLRC
                       , nEtaWindows
                       , nPhiWindows
                       , etaWinWidth
                       , etaWindowStep
                       , 0
                       , phiWinWidth
                       );
    
    taskLRC->SetNumberOfPhiSectors( nEventPhiRotations );
    taskLRC->SetAnalysisLevel( type );
    
    if ( type == "AOD" )
        taskLRC->SetAODtrackCutBit( 128 );//128 );//128 );
    else if ( type == "ESD" )
        taskLRC->SetTrackCuts(createAliLRCcuts( "StandardTPCOnlyTrackCuts" )); // "Global2" ));

    taskLRC->SetVtxDiamond(0.4,0.4,vertexZmax);
    taskLRC->SetMaxPtLimit(ptMax);
    taskLRC->SetMinPtLimit(ptMin);
    //taskLRC->SetMinNumberOfSPDtracklets( systConfig->minNtracklets );
    taskLRC->SetNchCuts( multMin, multMax ); // !!! for mult studies
    
    int maxMultInHistograms = ( !isIonsAnalysis ? 15 : 70 );
    setHistMultRange( taskLRC,  0, 0, maxMultInHistograms );
    setHistPtRange( taskLRC, ptMin, ptMax, 0.01, 5 );
    
    //if PbPb or pPb analysis
    taskLRC->SetIonsAnalysis(isIonsAnalysis);
    // to be set in the train wagon config: taskLRC->SetCentralityClass(min, max);
    
    taskLRC->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented
    mgr->AddTask(taskLRC);
    
    configureLRCtaskOutput(taskLRC, /*strOutputRootFolder,*/ strPrefixTaskName, strRunMode);
    
    return taskLRC;
}



