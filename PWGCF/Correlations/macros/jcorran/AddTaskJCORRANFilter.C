const Bool_t IsMC = kFALSE; //With real data kMC = kFALSE
const TString kInputData = "AOD";
const TString kJCORRANInputFormat = "AOD"; // ESD, AODout, AODin
const Bool_t kDoStoreJOD = kTRUE;

//_____________________________________________________________________
AliAnalysisTask *AddTaskJCORRANFilter(){
    // Load Custom Configuration and parameters
    // override values with parameters
    UInt_t triggerSelection;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    // AOD handler
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);

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
    jctask->SetOutputAODName("jcorran.root");
    jctask->SetDebugLevel(0);
    jctask->SetDoStoreJOD( kDoStoreJOD );
    cout<<"DEBUG 4"<< jctask->GetDoStoreJOD() <<"\t"<<kDoStoreJOD<<endl;

    AliJRunHeader * hdr = new AliJRunHeader;
    hdr->SetIsMC( IsMC );
    hdr->SetBeamTypeI( 0 ); // 0:pp 1:PbPb
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
    jfilter->SetClusterThreshold( 0 );
    jfilter->SetTrackThreshold( 0 );

    //==event selection
    jctask->SelectCollisionCandidates( triggerSelection );  //Apply offline trigger selection by AliPhysicsSelectionTask

    mgr->AddTask((AliAnalysisTask*) jctask);

    //==== Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);

// Connect input/output
	mgr->ConnectInput(jctask, 0, cinput);
	if( kDoStoreJOD ){ 
		AliAnalysisDataContainer *runinfoOutput = mgr->CreateContainer("RunInfo",  TList::Class(), AliAnalysisManager::kOutputContainer, "jcorran.root");
		AliAnalysisDataContainer *treeOutput = mgr->CreateContainer("JODTree",  TTree::Class(), AliAnalysisManager::kOutputContainer, "jcorran.root");
		mgr->ConnectOutput(jctask, 1, treeOutput );
		mgr->ConnectOutput(jctask, 2, runinfoOutput );
	}

	return jctask;
}

