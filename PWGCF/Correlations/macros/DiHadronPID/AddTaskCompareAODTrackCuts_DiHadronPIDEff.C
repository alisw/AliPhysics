AliAnalysisTaskCompareAODTrackCuts* AddTaskCompareAODTrackCuts_DiHadronPIDEff(
	Double_t MinCentrality = 5.,
	Double_t MaxCentrality = 0.,
	const char* CentralityEstimator = "V0M",
	Double_t maxVertexZ = 7.,
	Double_t maxEta = 0.8,
	Double_t minAssociatedPt = 0.2,
	Double_t maxAssociatedPt = 5.0,
	Double_t minTriggerPt = 5.,
	Double_t maxTriggerPt = 10.,
	Bool_t requestAllSingleTrackHistos = kTRUE,
	Int_t FilterMaskTrigger = 7,
	Int_t FilterMaskAssociated = 10,
	Bool_t isPbPb = kTRUE,
	Bool_t isMC = kTRUE,
	Int_t DebugLevel = 0,
	const char* outputFileName = 0,
	const char* containerName = "DiHadronPIDEff",
	const char* folderName = "PWGCF_DiHadronPID") 

{

	// Get a pointer to the analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        cout<<"AddTaskDiHadronPID.C -> No analysis manager found."<<endl;
        return 0x0;
    }

    // Creating the analysis task.
    AliAnalysisTaskCompareAODTrackCuts* EfficiencyTask = new AliAnalysisTaskCompareAODTrackCuts(containerName);
    EfficiencyTask->SetVerbose(kFALSE);
    EfficiencyTask->SetMC(isMC);
    EfficiencyTask->SetCalculateTOFMismatch(kTRUE);

        // Configure and add Event Cuts.
	AliAODEventCutsDiHadronPID* eventcuts = new AliAODEventCutsDiHadronPID("EventCuts");
	eventcuts->SetTrigger(AliVEvent::kMB);
	eventcuts->SetCentrality(MaxCentrality, MinCentrality);
	eventcuts->SetMaxVertexZ(maxVertexZ);
	eventcuts->SetCentralityEstimator(CentralityEstimator);
	eventcuts->SetIsPbPb(isPbPb);
	eventcuts->SetDebugLevel(DebugLevel);
	EfficiencyTask->SetEventCuts(eventcuts);

	// Configure and add track cuts for trigger.
	AliAODTrackCutsDiHadronPID* triggercuts = new AliAODTrackCutsDiHadronPID("TrackCutsTrigger");
	triggercuts->SetIsMC(isMC);
	triggercuts->SetFilterMask(FilterMaskTrigger);
	triggercuts->SetPtRange(minTriggerPt,maxTriggerPt);
	triggercuts->SetMaxEta(maxEta);
	triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllCharged);
	triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPositive);
	triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegative);
	if (requestAllSingleTrackHistos) {
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllPion);
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPosPion);
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegPion);					
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllKaon);
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPosKaon);
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegKaon);		
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllProton);
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPosProton);
		triggercuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegProton);
	}	
	triggercuts->SetDebugLevel(DebugLevel);
	EfficiencyTask->AddTrackCuts(triggercuts);

	// Configure and add track cuts for associateds.
	AliAODTrackCutsDiHadronPID* associatedscuts = new AliAODTrackCutsDiHadronPID("TrackCutsAssociated");
	associatedscuts->SetIsMC(isMC);
	associatedscuts->SetFilterMask(FilterMaskAssociated);
	associatedscuts->SetPtRange(minAssociatedPt,maxAssociatedPt);
	associatedscuts->SetMaxEta(maxEta);
	ULong_t associatedflags = (UInt_t)(AliAODTrack::kTOFout)|(UInt_t)(AliAODTrack::kTIME);	
	associatedscuts->SetDemandFlags(associatedflags);
	associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllCharged);
	associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPositive);
	associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegative);
	if (requestAllSingleTrackHistos) {
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllPion);
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPosPion);
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegPion);					
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllKaon);
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPosKaon);
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegKaon);		
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kAllProton);
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kPosProton);
		associatedscuts->RequestQAHistos(AliAODTrackCutsDiHadronPID::kNegProton);
	}
	associatedscuts->SetDebugLevel(DebugLevel);
	EfficiencyTask->AddTrackCuts(associatedscuts);

	// Add the task.
	mgr->AddTask(EfficiencyTask);
    
	// Data containers.
	AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
	mgr->ConnectInput(EfficiencyTask, 0, cinput); 
	
	if (!outputFileName) {outputFileName = AliAnalysisManager::GetCommonFileName();}
	
	AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(containerName, TList::Class(),
                         AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
	
	mgr->ConnectOutput(EfficiencyTask,1,coutput1);
	
	return EfficiencyTask;

}
