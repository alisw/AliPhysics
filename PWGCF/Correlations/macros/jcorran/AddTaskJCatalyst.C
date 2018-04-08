//_____________________________________________________________________
AliAnalysisTask *AddTaskJCatalyst(TString taskName = "JCatalyst", UInt_t flags = 0, Int_t FilterBit = 768 , double eta_min = -0.8, double eta_max = 0.8, double pt_min = 0.0, double pt_max = 100.0, int debuglevel = 0){
    cout<<"AddTaskJCatalyst::flags = "<<flags<<endl;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    //==== Check the analysis type using the event handlers connected to the analysis mgr	
    if (!mgr->GetInputEventHandler() ){
	    ::Error("AddTaskFFluc", "This task requires an input event handler" );
	    return NULL;
    }

    //==== JCORRAN TASK
    AliJCatalystTask *FFtask = new AliJCatalystTask(taskName.Data());

    //TODO: test flags for call()
    FFtask->SetJCatalystTaskName( taskName.Data() ) ;
    FFtask->AddFlags(flags);
    FFtask->SetTestFilterBit( FilterBit ) ;
    FFtask->SetEtaRange( eta_min, eta_max);
    FFtask->SetDebugLevel( debuglevel ) ; 
    FFtask->SetPtRange( pt_min, pt_max);
    FFtask->SetParticleCharge( 0 );  // 0 : all charged particles 1 : positive -1 :negative 
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

