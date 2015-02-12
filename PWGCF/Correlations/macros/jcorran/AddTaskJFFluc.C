//_____________________________________________________________________
AliAnalysisTask *AddTaskJFFluc(TString taskName,Bool_t IsMC = kFALSE, Int_t FilterBit = 768 , double eta_min, double eta_max, int debuglevel){
    // Load Custom Configuration and parameters
    // override values with parameters
	cout <<"AddTaskJFFluc:: IsMC = "<< IsMC <<endl;

	//==== Get the pointer to the Analyis mgr
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	/*
    //=== AOD handler
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
	*/

	//==== Check the analysis type using the event handlers connected to the analysis mgr	
	if (!mgr->GetInputEventHandler() ){
		::Error("AddTaskFFluc", "This task requires an input event handler" );
		return NULL;
	}


    //==== JCORRAN TASK
    int CollisionCandidates = AliVEvent::kCentral;
    AliJFFlucTask *FFtask = new AliJFFlucTask( taskName , CollisionCandidates, IsMC );

	FFtask->SetFFlucTaskName( taskName ) ;
	FFtask->SetIsMC( IsMC );
	FFtask->SetTestFilterBit( FilterBit ) ;
	FFtask->SetEtaRange( eta_min, eta_max);
	FFtask->SetDebugLevel( debuglevel ) ; 
	//FFtask->SelectCollisionCandidates( AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral ) ; 
	//FFTask->SelectCollisionCandidates( CollisionCandidates );
	//FFtask->SelectCollisionCandidates( AliVEvent::kCentral ) ; 

	//==== Add task
	mgr->AddTask((AliAnalysisTask*) FFtask);

	//==== Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	AliAnalysisDataContainer *FFhist = mgr->CreateContainer(Form("%scontainer",FFtask->GetName()),    TDirectory::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:%s",  AliAnalysisManager::GetCommonFileName(), FFtask->GetName()));

	//==== Connect input/output
	mgr->ConnectInput(FFtask, 0, cinput);
	mgr->ConnectOutput(FFtask, 1, FFhist);

	return FFtask;
}

