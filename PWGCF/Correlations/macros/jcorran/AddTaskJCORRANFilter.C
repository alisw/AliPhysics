const TString kInputData = "AOD";
const TString kJCORRANInputFormat = "AOD"; // ESD, AODout, AODin

//_____________________________________________________________________
AliAnalysisTask *AddTaskJCORRANFilter(Bool_t IsMC = kFALSE, Int_t beamtype = 1){
    // Load Custom Configuration and parameters
    // override values with parameters
	cout <<"AddTaskJCORRANFilter:: IsMC = "<< IsMC <<"\t beamtype = "<< beamtype <<endl;
    UInt_t triggerSelection;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    // AOD handler
    //AliAODInputHandler *aodHandler = new AliAODInputHandler();
    //mgr->SetInputEventHandler(aodHandler);

    //================================================
    // TASKS
    //================================================

    // set the trigger selection
    triggerSelection =  AliVEvent::kAny;
    //										| AliVEvent::kHighMult 
    //										| AliVEvent::kEMCEGA;
    //										| AliVEvent::kEMCEJE
    //										| AliVEvent::kEMC1
    //										| AliVEvent::kEMC7
    //                    | AliVEvent::kCentral
    //                    | AliVEvent::kSemiCentral; 
    //============================
    //   JCORRANTask
    //===========================

    //==== JCORRAN TASK
    AliJCORRANTask *jctask = new AliJCORRANTask("PWGCFJCORRANTask",kJCORRANInputFormat);
    jctask->SetDebugLevel(0);

    AliJRunHeader * hdr = new AliJRunHeader;
    hdr->SetIsMC( IsMC );
    hdr->SetBeamTypeI( beamtype ); // 0:pp 1:PbPb
    hdr->SetWithoutSDD(false);
    hdr->SetRunType("LHC13c");
    hdr->SetInputFormat( 1 ); // 0: ESD;
    hdr->SetRefitESDVertexTracks(kFALSE);
    hdr->SetStoreEventPlaneSource(kFALSE);
    hdr->SetStoreTPCTrackBitMask( 1<<7 ); // TODO : For what?
    hdr->SetStoreGCGTrackBitMask( 2^30-1 );
    hdr->SetStoreEMCalInfo( kFALSE ); // error with lego now.


    jctask->SetJRunHeader( hdr );

    AliJFilter *jfilter = jctask->GetFilter();
    jfilter->SetAliJRunHeader( hdr );
    jfilter->SetTrackThreshold( 0 );

    //==event selection
    jctask->SelectCollisionCandidates( triggerSelection );  //Apply offline trigger selection by AliPhysicsSelectionTask

    mgr->AddTask((AliAnalysisTask*) jctask);

    //==== Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

// Connect input/output
    mgr->ConnectInput(jctask, 0, cinput);

	return jctask;
}

