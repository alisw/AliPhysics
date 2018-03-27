AliBSDiJetTask * AddAliBSDiJetTask(TString taskname, bool isAA, UInt_t triggMask, Double_t leadingParticlePtMin, AliEmcalJetTask* jetFinderTask, AliEmcalJetTask* jetFinderTaskkt, bool kineon, TString option){
	AliBSDiJetTask *task = new AliBSDiJetTask(taskname, option.Data());
	if (!option.Contains("MC")) {
		task->SelectCollisionCandidates( triggMask ); // WARNING: default is AliVEvent::kEMCEJE. Double check it!
	}
	task->SetIsAA( isAA );
	if (option.Contains("MC")) task->SetIsMC(1);

	task->SetLeadingParticlePtMin( leadingParticlePtMin );

	AliJetContainer *jetCont = task->AddJetContainer(jetFinderTask->GetName(),taskname.Contains("FullJet") ? "EMCALfid" : "TPCfid");
	AliJetContainer *jetContKt = 0x0;
	AliParticleContainer *trackCont = 0x0;
	AliClusterContainer *clusterCont = 0x0;
	if (kineon) {
		AliMCParticleContainer *mcCont = task->AddMCParticleContainer("mcparticles");
	}
	else {
		TString trackcontname="tracks";
		//if (option.Contains("12a15e") || option.Contains("16j5")) trackcontname = "PicoTracksMer";
		trackCont = task->AddParticleContainer(trackcontname.Data());
		if (taskname.Contains("FullJet")) clusterCont = task->AddClusterContainer("caloClusters");
		jetContKt = task->AddJetContainer(jetFinderTaskkt->GetName(),taskname.Contains("FullJet") ? "EMCALfid" : "TPCfid");
	}

	if (kineon) jetCont->ConnectParticleContainer( mcCont );
	else {
		jetCont->ConnectParticleContainer( trackCont );
		if (taskname.Contains("FullJet")) jetCont->ConnectClusterContainer(clusterCont);
		jetContKt->ConnectParticleContainer ( trackCont );
		jetContKt->SetLeadingHadronType( 0 );

	}
	//jetCont->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
	cout<<"jetFinderTask->GetRadius() : "<<jetFinderTask->GetRadius()<<endl;
	if (jetFinderTask->GetRadius() >= 0.4)  jetCont->SetPercAreaCut( 0.6 );
	if (!taskname.Contains("FullJet")) jetCont->SetLeadingHadronType( 0 ); //leading hadron set to h+-
	//jetCont->SetJetPtCut( 5 );

	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(taskname, TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s", AliAnalysisManager::GetCommonFileName()));

	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	return task;
}
