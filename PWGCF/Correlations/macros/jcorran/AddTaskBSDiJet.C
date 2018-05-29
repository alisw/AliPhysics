#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C>
#endif

AliBSDiJetTask * AddTaskBSDiJet(TString taskname, bool isAA, Double_t leadingParticlePtMin, TString option ){

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskJJet", "No analysis manager to connect to.");
		return NULL;
	}
#if !defined(__CLING__)
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
#endif

		// Own class
	AliBSDiJetTask *task = new AliBSDiJetTask(taskname, option);
	task->SetIsAA( isAA );
	task->SetLeadingParticlePtMin( leadingParticlePtMin );

	AliEmcalJetTask *jetFinderTask;
	AliEmcalJetTask *jetFinderTaskkt;
	TString trackcontname="tracks";
	jetFinderTask = AddTaskEmcalJet(trackcontname.Data(),"", AliJetContainer::antikt_algorithm,0.4,AliJetContainer::kChargedJet,0.15, 0.300, 0.005, AliJetContainer::pt_scheme, "Jets", 0. );
	jetFinderTaskkt = AddTaskEmcalJet(trackcontname.Data(),"",AliJetContainer::kt_algorithm,0.2,AliJetContainer::kChargedJet,0.15, 0.300, 0.005, AliJetContainer::pt_scheme, "JetsKt", 0. );

	AliJetContainer *jetCont = task->AddJetContainer(jetFinderTask->GetName(),"TPCfid");
	AliJetContainer *jetContKt = task->AddJetContainer(jetFinderTaskkt->GetName(),"TPCfid");

	AliParticleContainer *trackCont = 0x0;
	trackCont = task->AddParticleContainer(trackcontname.Data());
	jetCont->ConnectParticleContainer( trackCont );
	jetContKt->ConnectParticleContainer ( trackCont );
	jetContKt->SetLeadingHadronType( 0 );

	//jetCont->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
	cout<<"jetFinderTask->GetRadius() : "<<jetFinderTask->GetRadius()<<endl;
	if (jetFinderTask->GetRadius() >= 0.4)  jetCont->SetPercAreaCut( 0.6 );
	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(taskname, TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s", AliAnalysisManager::GetCommonFileName()));

	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	return task;
}
