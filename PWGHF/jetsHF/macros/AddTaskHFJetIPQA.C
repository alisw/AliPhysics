AliAnalysisTaskHFJetIPQA* AddTaskHFJetIPQA(
		Bool_t IsESD = kFALSE,
		const char *ntracks            = "Tracks",
		const char *nclusters           = "",
		const char *njets              = "Jets",
		const char *nrho               = "			",
		Double_t jetradius =0.4,
		Bool_t isMC = kFALSE,
		const char * type = "TPC",
		const char *taskname           = "AliAnalysisTaskEmcalJetBJetTaggingIP",
		const char *njetsMC              = "Jets",
		const char *nrhoMC               = "RhoMC",
		TString PathToWeights = 	"alien:///alice/cern.ch/user/l/lfeldkam/Weights.root",
		const char* suffix = ""
)
{


	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
		return NULL;
	}


	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler())
	{
		::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
		return NULL;
	}




	TFile* filecorrectionfactors;
	Bool_t useweighting = kTRUE;

	if(isMC){

		if( PathToWeights.EqualTo("") ) {
			useweighting = kFALSE;
		} else {
			filecorrectionfactors=TFile::Open(PathToWeights.Data());
			if(!filecorrectionfactors ||(filecorrectionfactors&& !filecorrectionfactors->IsOpen())){
				AliFatal("Input weight file not found");
				return 0x0;
			}
		}
	}


	TString name(taskname);

	TString combinedName;
	combinedName.Form("%s%s", name.Data(),suffix);

	AliAnalysisTaskHFJetIPQA* jetTask = new AliAnalysisTaskHFJetIPQA(combinedName);





	if(isMC && filecorrectionfactors){
		TH1F * h[20] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};

	    const char * nampart[20] = {"pi0","eta","etap","rho","phi","omega","k0s","lambda","pi","kaon","proton","D0","Dp","Dsp","Ds","lambdac","bplus","b0","lambdab","bsp"};
	    for(int i = 0;i<20;++i){
			h[i]=(TH1F*)filecorrectionfactors->Get(nampart[i]);
	    }
	    Printf("Done Reading");
		for (int i = 0 ; i<20;++i) if(h[i]==0) return 0x0;
		jetTask->SetUseMonteCarloWeighingLinus(h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10],h[11],h[12],h[13],h[14],h[15],h[16],h[17],h[18],h[19]);
	    Printf("Weights written");

	}

	AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
	AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

	TString strType(type);

	AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);

	if(jetCont) {
		jetCont->SetRhoName(nrho);

		jetCont->ConnectParticleContainer(trackCont);

		jetCont->ConnectClusterContainer(clusterCont);
	}

	if(isMC)
	{
		AliJetContainer *jetContMC = jetTask->AddJetContainer(njetsMC,strType,jetradius);

		if(jetContMC) {

			jetContMC->SetRhoName(nrhoMC);

			jetContMC->SetIsParticleLevel(kTRUE);

			jetContMC->SetMaxTrackPt(1000);
		}
	}
	//-------------------------------------------------------
	//  Configure analysis task
	//-------------------------------------------------------

	jetTask->SetIsPythia(isMC);

	DefineCutsTaskpp(jetTask,-1.,100);
	if(IsESD) {
		jetTask->SetRunESD();
		  AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
			   esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
			   esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
			   esdTrackCutsH->SetDCAToVertex2D(kTRUE);
		 jetTask->SetESDCuts(new AliESDtrackCuts(*esdTrackCutsH));
	}

	//-------------------------------------------------------
	// Final settings, pass to manager and set the containers
	//-------------------------------------------------------

	mgr->AddTask(jetTask);


	// Create containers for input/output
	AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

	TString contname(combinedName);
	contname += "_histos";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
			TList::Class(),AliAnalysisManager::kOutputContainer,
			Form("%s", AliAnalysisManager::GetCommonFileName()));

	mgr->ConnectInput  (jetTask, 0,  cinput1 );

	mgr->ConnectOutput (jetTask, 1, coutput1 );

	return jetTask;
}

Bool_t DefineCutsTaskpp(AliAnalysisTaskHFJetIPQA *task, Float_t minC, Float_t maxC)
{
	// define cuts for task
	AliRDHFJetsCuts *cuts=task->GetJetCutsHF();
	// jets
	cuts->SetJetRadius(0.4); // this cut does nothing
	cuts->SetMaxEtaJet(0.5);//0.9-R
	cuts->SetMinPtJet(5.);
	cuts->SetMaxPtJet(250.);
	// Set centrality
	cuts->SetMinCentrality(minC);
	cuts->SetMaxCentrality(maxC);
	cuts->SetUsePhysicsSelection(kFALSE);
	cuts->SetOptPileup(1);
	cuts->ConfigurePileupCuts(5,0.8);
	cuts->SetTriggerClass("");
	cuts->SetTriggerMask(AliVEvent::kMB);
	cuts->PrintAll();
	// pPb minbias only
	return kTRUE;
}


 //
