
AliAnalysisTaskHFJetIPQA* AddTaskHFJetIPQA(
		const char *ntracks            = "PicoTracks",
		const char *nclusters           = "",
		const char *njets              = "Jets",
		const char *nrho               = "",
		Double_t jetradius =0.4,
		Bool_t isMC = kFALSE,
		const char * type = "TPC",
		const char *taskname           = "AliAnalysisTaskEmcalJetBJetTaggingIP",
		const char *njetsMC              = "Jets",
		const char *nrhoMC               = "RhoMC",
		TString PathToWeights = 	" alien:///alice/cern.ch/user/l/lfeldkam/nonHFEcorrect.root",
		const char *NamePionHist = 	 "hRatio_10f6_pion",
		const char *NameEtaHist = 	 "hRatio_10f6_eta",
		const char *NameOmegaHist =  "hRatio_10f6_omega",
		const char *NamePhiHist = 	 "hRatio_10f6_phi",
		const char *NameEtapHist = 	 "hRatio_10f6_etap",
		const char *NameRhoHist = 	 "hRatio_10f6_rho",
		const char *NameKaonHist = 	 "hRatio_10f6_kaon",
		const char *NameK0sHist = 	 "hRatio_10f6_k0s",
		const char *NameLambdaHist = "hRatio_10f6_lambda",
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
	if (strcmp(njets,"")) {
		name += "_";
		name += njets;
	}
	if (strcmp(nrho,"")) {
		name += "_";
		name += nrho;
	}
	if (strcmp(type,"")) {
		name += "_";
		name += type;
	}
	TString combinedName;
	combinedName.Form("%s%s", name.Data(),suffix);

	AliAnalysisTaskHFJetIPQA* jetTask = new AliAnalysisTaskHFJetIPQA(combinedName);



	const Double_t binLimit[45] =  {
			0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,
			0.295573,0.333397,0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,
			1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,
			4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20
	};

	if(isMC && filecorrectionfactors){
		TH1D * h[9] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
		h[0]=(TH1D*)filecorrectionfactors->Get(NamePionHist);
		h[1]=(TH1D*)filecorrectionfactors->Get(NameEtaHist);
		h[2]=(TH1D*)filecorrectionfactors->Get(NameOmegaHist);
		h[3]=(TH1D*)filecorrectionfactors->Get(NamePhiHist);
		h[4]=(TH1D*)filecorrectionfactors->Get(NameEtapHist);
		h[5]=(TH1D*)filecorrectionfactors->Get(NameRhoHist);
		h[6]=(TH1D*)filecorrectionfactors->Get(NameKaonHist);
		h[7]=(TH1D*)filecorrectionfactors->Get(NameK0sHist);
		h[8]=(TH1D*)filecorrectionfactors->Get(NameLambdaHist);


		for (int i = 0 ; i<9;++i) if(h[i]==0) return 0x0;
		jetTask->SetUseMonteCarloWeighing(h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8]);

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


