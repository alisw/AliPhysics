//_____________________________________________________________________
AliAnalysisTask *AddTaskJXt(TString taskName="JFFluc",Bool_t IsMC = kFALSE, Int_t FilterBit = 32 , double etaRange = 0.8,
		double minPtIsol = 2.0, int effMode = 0, int debuglevel = 0, char* suffix=""){
    // Load Custom Configuration and parameters
    // override values with parameters
    // surfix in last arguments are added for subwagons.
    TString combinedName = Form("%s%s", taskName.Data(), suffix);
    //
	cout <<"AddTaskJXt:: IsMC = "<< IsMC <<endl;
	cout <<"AddTaskJXt:: Efficiency Corr Mod? (0:no, 1:period, 2:run#, 3:auto) = " << effMode << endl;

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
    AliJXtTask *xTtask = new AliJXtTask( combinedName.Data(), FilterBit, IsMC );

	xTtask->SetXtTaskName( combinedName.Data() ) ;
	xTtask->SetIsMC( IsMC );
	xTtask->SetTestFilterBit( FilterBit ) ;
	xTtask->SetEtaRange( etaRange );
	xTtask->SetDebugLevel( debuglevel ) ; 

	//==== Add task
	mgr->AddTask((AliAnalysisTask*) xTtask);

	//==== Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *FFhist = mgr->CreateContainer(Form("%scontainer",xTtask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:%s",  AliAnalysisManager::GetCommonFileName(), xTtask->GetName()));

	//==== Connect input/output
	mgr->ConnectInput(xTtask, 0, cinput);
	mgr->ConnectOutput(xTtask, 1, FFhist);
	return xTtask;
}

