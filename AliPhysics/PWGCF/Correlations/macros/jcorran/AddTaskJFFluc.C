//_____________________________________________________________________
//AliAnalysisTask *AddTaskJFFluc(TString taskName = "JFFluc", Bool_t IsMC = kFALSE, Bool_t IsWeakExclude = kFALSE, Bool_t IsCentFlat = kFALSE, Int_t FilterBit = 768 , double eta_min = 0.4, double eta_max = 0.8, double pt_min = 0.2, double pt_max = 5.0, int effMode = 0, Bool_t IsPhiModule = kFALSE, TString InFileNameNUE  = "", int debuglevel = 0, TString suffix = ""){
AliAnalysisTask *AddTaskJFFluc(TString taskName = "JFFluc", UInt_t flags = 0, Int_t FilterBit = 768 , double eta_min = 0.4, double eta_max = 0.8, double pt_min = 0.2, double pt_max = 5.0, int effMode = 0, int debuglevel = 0, TString suffix = ""){
    TString combinedName = taskName;
	if(suffix.Length() > 0)
		combinedName += "_"+suffix;
    //
	cout<<"AddTaskJFFluc::flags = "<<flags<<endl;
	//cout <<"AddTaskJFFluc:: IsMC = "<< IsMC <<endl;
	//cout <<"AddTaskJFFluc:: IsWeakExclude = "<< IsWeakExclude << endl;
	//cout <<"AddTaskJFFluc:: Force to Cent flatting for LHC11h? = " << IsCentFlat << endl;
	//cout <<"AddTaskJFFluc:: Efficiency Corr Mod? (0:no, 1:period, 2:run#, 3:auto) = " << effMode << endl;
	
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
    AliJFFlucTask *FFtask = new AliJFFlucTask(combinedName.Data());
   
	FFtask->SetFFlucTaskName( combinedName.Data() ) ;
	FFtask->AddFlags(flags);
	FFtask->SetTestFilterBit( FilterBit ) ;
	FFtask->SetEtaRange( eta_min, eta_max);
	FFtask->SetDebugLevel( debuglevel ) ; 
	FFtask->SetPtRange( pt_min, pt_max);
	FFtask->SetEffConfig( effMode, FilterBit); 
	FFtask->SetParticleCharge( 0 );  
    //int CollisionCandidates = AliVEvent::kCentral;
	//FFtask->SelectCollisionCandidates( AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral ) ; 
	//FFTask->SelectCollisionCandidates( CollisionCandidates );
	//FFtask->SelectCollisionCandidates( AliVEvent::kCentral ) ; 

	//==== Add task
	mgr->AddTask((AliAnalysisTask*) FFtask);

	//==== Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *FFhist = mgr->CreateContainer(Form("%scontainer",FFtask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:%s",  AliAnalysisManager::GetCommonFileName(), FFtask->GetName()));

	//==== Connect input/output
	mgr->ConnectInput(FFtask, 0, cinput);
	mgr->ConnectOutput(FFtask, 1, FFhist);
	
	return FFtask;
}

