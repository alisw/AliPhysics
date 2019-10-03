#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C>
#endif
	
class AliAnalysisTaskEmcalEmbeddingHelper;


AliBSDiJetTask * AddTaskBSDiJet(TString taskname, bool isAA, Double_t leadingParticlePtMin, TString option ){

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskJJet", "No analysis manager to connect to.");
		return NULL;
	}

	if (0){
		gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
		AliEmcalPicoTrackMaker* picomaker = AddTaskEmcalPicoTrackMaker("PicoTracks","tracks");
		gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromAOD.C");
		AliJetEmbeddingFromAODTask* embTask = AddTaskJetEmbeddingFromAOD(
				"Tracks"
				, ""
				, ""
				, ""
				, "alien:///alice/cern.ch/user/k/kimb/filedata.txt"
				, "aodTree"
				, "tracks"
				, "" , ""
				, ""
				, "lhc10h"
				, kFALSE
				, 0
				, 10
				, AliVEvent::kMB
				, 0,0, 1234567890, "Embedding"
				);
		embTask->SetMarkMC(99999);
		embTask->SetZVertexCut(10);
		//embTask->SetMaxVertexDist(2);
		embTask->SetRandomAccess(true);
		//TString trackEmbDistr = embTask->GetOutTrackName();
		gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskMergeBranches.C");
		AliJetModelMergeBranches* brmergTask = AddTaskMergeBranches("PicoTracks","Tracks","Mer","");
		brmergTask->SetCopyArray(1);
	}

	if (option.Contains("Emb")){
		AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
		// Set the file pattern. This example uses ptHardBin 4 of LHC12a15e_fix.
		// The pT hard bin and anchor run can also be set by adding a printf() wild card to the string (ie %d)
		// See the documentation of AliAnalysisTaskEmcalEmbeddingHelper::GetFilenames()
		//embeddingHelper->SetInputFilename("alice/data/2015/LHC15o/");
		//embeddingHelper->SetFilePattern("alien:///alice/data/2015/LHC15o/000246844/pass1/AOD194/");
		//embeddingHelper->SetFilePattern("alien:///alice/sim/2016/LHC16j5/15/246488/AOD200/");
		//if (option.Contains("LHC15o")) gSystem->Exec(Form("alien_find /alice/sim/2016/LHC16j5/ AliAOD.root | grep AOD200 | perl -nle'print \"alien://\".$_' | sort -R | head -300 > embfile.txt"));
		if (option.Contains("LHC15o")) gSystem->Exec(Form("alien_find /alice/sim/2016/LHC16j5/ AliAOD.root | grep AOD200 | perl -nle'print \"alien://\".$_' > embfile.txt"));
		embeddingHelper->SetFileListFilename("./embfile.txt");

		//embeddingHelper->SetFilePattern("alien:///alice/sim/2016/LHC16j5/%d/%d/AOD200/");
		//embeddingHelper->SetFileListFilename("alien:///alice/cern.ch/user/k/kimb/Emb_LHC15o.txt");
		// If the embedded file is an ESD, then set:
		embeddingHelper->SetAOD();
		// Add additional configure as desired.
		// For some examples...
		// ... randomly select which file to start from:
		embeddingHelper->SetRandomFileAccess(kTRUE);
		// ... Start from a random event within each file
		//embeddingHelper->SetRandomEventNumberAccess(kTRUE);
		embeddingHelper->SetTriggerMask(AliVEvent::kINT7);
		embeddingHelper->SetZVertexCut(10);
		//embeddingHelper->SetMaxVertexDistance(2);
		//embeddingHelper->SetCentralityRange(centmin,centmax);
		//embeddingHelper->SetUseManualInternalEventCuts(true) ;
		// ... Set pt hard bin properties
		//embeddingHelper->SetPtHardBin(5);
		//embeddingHelper->SetNPtHardBins(20);
		// etc..
		// As your last step, always initialize the helper!
		embeddingHelper->Initialize();
		AliAnalysisTaskBSEmbedding * preprocess = AliAnalysisTaskBSEmbedding::AddTaskBSEmbedding();
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
	AliEmcalJetTask *jetFinderTaskkine;
	//TString trackcontname = option.Contains("Emb") ? "PicoTracksMer" : "tracks";
	TString trackcontname  = "tracks";
	TString mccontname="mcparticles";
	jetFinderTask = AddTaskEmcalJet("usedefault","", AliJetContainer::antikt_algorithm,0.4,AliJetContainer::kChargedJet,0.15, 0.300, 0.005, AliJetContainer::pt_scheme, "Jets", 0. );
	jetFinderTaskkt = AddTaskEmcalJet("usedefault","",AliJetContainer::kt_algorithm,0.2,AliJetContainer::kChargedJet,0.15, 0.300, 0.005, AliJetContainer::pt_scheme, "Jets", 0. );

	if (option.Contains("Emb")) {
		jetFinderTaskkine = AddTaskEmcalJet(mccontname.Data(),"", AliJetContainer::antikt_algorithm,0.4,AliJetContainer::kChargedJet,0.15, 0.300, 0.005, AliJetContainer::pt_scheme, "Jets", 0. );
	}


	AliJetContainer *jetCont = task->AddJetContainer(jetFinderTask->GetName(),"TPC");
	AliJetContainer *jetContKt = task->AddJetContainer(jetFinderTaskkt->GetName(),"TPC");

	AliParticleContainer *trackCont = 0x0;
	trackCont = task->AddParticleContainer(trackcontname.Data());
	if (option.Contains("Emb")) trackCont -> SetIsEmbedding(true);
	jetFinderTask->AdoptParticleContainer( trackCont );
	jetFinderTaskkt->AdoptParticleContainer ( trackCont );
	//jetContKt->SetLeadingHadronType( 0 );

	AliMCParticleContainer *mcCont = 0x0;
	if (option.Contains("Emb")) {
		mcCont= task->AddMCParticleContainer (mccontname.Data());
		AliJetContainer *jetContMC = task->AddJetContainer(jetFinderTaskkine->GetName(),"TPC");
		jetFinderTaskkine->AdoptMCParticleContainer( mcCont );
		
		//jetFinderTaskkine -> GetParticleContainer(0)->SetIsEmbedding(true);
	}

	//jetCont->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
	cout<<"jetFinderTask->GetRadius() : "<<jetFinderTask->GetRadius()<<endl;
	if (jetFinderTask->GetRadius() >= 0.4)  jetCont->SetPercAreaCut( 0.6 );
	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(taskname, TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s", AliAnalysisManager::GetCommonFileName()));
	if (option.Contains("Emb")) task -> SetIsMC(true);

	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	return task;
}
