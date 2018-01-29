void AddAnalysisTaskPP(Bool_t isMC = kFALSE, UInt_t offlineTriggerMask, TString description, TString suff = "", TString badmap = "")
{
	AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();

	// Setup Selections
	TList * selections = new TList();

	AliPP13ClusterCuts cuts_pi0 = AliPP13ClusterCuts::GetClusterCuts();
	AliPP13ClusterCuts cuts_eta = AliPP13ClusterCuts::GetClusterCuts();
	cuts_eta.fAsymmetryCut = 0.7;


	if (!isMC)
	{
		// TODO: Add plain selections
		AliPP13SelectionWeights & data_weights = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kData);
		AliPP13SelectionWeights & data_weights_plain = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kPlain);

		selections->Add(new AliPP13PhysPhotonSelection("Phys", "Physics Selection", cuts_pi0, &data_weights));
		selections->Add(new AliPP13PhysPhotonSelection("PhysPlain", "Physics Selection no TOF cut efficiency", cuts_pi0, &data_weights_plain));
		selections->Add(new AliPP13PhotonTimecutStudySelection("Time", "Testing Timing Selection", cuts_pi0, &data_weights));

		selections->Add(new AliPP13PhysPhotonSelection("Eta", "Physics Selection for eta meson", cuts_eta, &data_weights));
		selections->Add(new AliPP13PhysPhotonSelection("EtaPlain", "Physics Selection for eta meson no TOF cut efficiency", cuts_eta, &data_weights_plain));
		selections->Add(new AliPP13PhotonTimecutStudySelection("EtaTime", "Testing Timing Selection for eta meson", cuts_eta, &data_weights));

		selections->Add(new AliPP13TagAndProbeSelection("TagAndProbleTOF", "Cluster p_{T} Selection", cuts_pi0, &data_weights));
		selections->Add(new AliPP13NonlinearitySelection("PhysNonlinEst", "Physics efficiency for neutral particles", cuts_pi0, &data_weights));
		selections->Add(new AliPP13NonlinearitySelection("PhysNonlinEstPlain", "Physics efficiency for neutral particles", cuts_pi0, &data_weights_plain));

		selections->Add(new AliPP13QualityPhotonSelection("Qual", "Cluster quality Selection", cuts_pi0, &data_weights));
		selections->Add(new AliPP13PhotonSpectrumSelection("PhotonsTime", "Cluster p_{T} Selection with timing cut", cuts_pi0, &data_weights, 10., 3.));
		
		selections->Add(new AliPP13PhotonSpectrumSelection("Photons", "Cluster p_{T} Selection", cuts_pi0, &data_weights));
		selections->Add(new AliPP13PhotonSpectrumSelection("PhotonsPlain", "Cluster p_{T} Selection", cuts_pi0, &data_weights_plain));

	}	
	else
	{
		AliPP13SelectionWeightsMC & mc_weights = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kMC);
		AliPP13SelectionWeightsMC & mc_weights_only = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kMC);

		// Nonlinearity for zs 20 Run2Default (Daiki's approximation)
		// The pi^0 peak is misplaced in this fit: A * 1.03274e+00 (global energy scale)
		// Calculated for the updated version for the corrected Data
		mc_weights.fNonA = -0.020025549129372242;
		mc_weights.fNonSigma = 1.1154536660217529;
	    mc_weights.fNonGlobal = 1.0493128193171741;		

		mc_weights_only.fNonGlobal = 1.0;
		mc_weights_only.fNonA = 0.0;

		selections->Add(new AliPP13EfficiencySelectionMC("PhysEff", "Physics efficiency for neutral particles fully corrected", cuts_pi0, &mc_weights));
		selections->Add(new AliPP13EfficiencySelectionMC("PhysEffPlain", "Physics efficiency for neutral particles, no nonlinearity", cuts_pi0, &mc_weights_only));

		selections->Add(new AliPP13NonlinearitySelection("PhysNonlin", "Physics nonlinearity for neutral particles", cuts_pi0, &mc_weights, kTRUE));
		selections->Add(new AliPP13NonlinearitySelection("PhysNonlinRaw", "Raw nonlinearity for neutral particles", cuts_pi0, &mc_weights_only, kTRUE));

		selections->Add(new AliPP13NonlinearityScanSelection("PhysNonlinScan", "Physics efficiency for neutral particles", cuts_pi0, &mc_weights));
		selections->Add(new AliPP13MesonSelectionMC("MCStudy", "MC Selection with timing cut", cuts_pi0, &mc_weights));	
	}


	if (isMC)
	{
		selections->Add(new AliPP13PhysPhotonSelection("Phys", "Physics Selection", cuts_pi0));
		selections->Add(new AliPP13PhotonTimecutStudySelection("Time", "Testing Timing Selection", cuts_pi0));
		selections->Add(new AliPP13TagAndProbeSelection("TagAndProbleTOF", "Cluster P_{t} Selection", cuts_pi0));

		selections->Add(new AliPP13PhysPhotonSelection("Eta", "Physics Selection for eta meson", cuts_eta));
		selections->Add(new AliPP13PhotonTimecutStudySelection("EtaTime", "Testing Timing Selection for eta meson", cuts_eta));

		selections->Add(new AliPP13QualityPhotonSelection("Qual", "Cluster quality Selection", cuts_pi0));
		selections->Add(new AliPP13PhotonSpectrumSelection("Photons", "Cluster P_{t} Selection", cuts_pi0));
		selections->Add(new AliPP13PhotonSpectrumSelection("PhotonsTime", "Cluster P_{t} Selection with timing cut", cuts_pi0, 10., 3.));
	}

	// Nonlinearity for zs 20 Run2Default (Daiki's approximation)
	// The pi^0 peak is misplaced in this fit: A * 1.03274e+00 (global energy scale)
	// Calculated for the updated version for the corrected Data
	Float_t nonlin_a = -0.020025549129372242;
	Float_t nonlin_b = 1.1154536660217529;
	Float_t ge_scale = 1.0493128193171741;

	Float_t weigh_a = -1.063;
	Float_t weigh_b = 0.855;


	if (isMC)
	{
		selections->Add(new AliPP13PhysPhotonSelectionMC("PhysNonlin", "Corrected for nonlinearity Physics Selection", cuts_pi0, nonlin_a, nonlin_b, ge_scale));
		selections->Add(new AliPP13PhysPhotonSelectionMC("PhysRaw", "Raw Physics Selection", cuts_pi0));

		selections->Add(new AliPP13MesonSelectionMC("MCStudy", "MC Selection with timing cut", cuts_pi0,
		                nonlin_a, nonlin_b, ge_scale,
		                weigh_a, weigh_b));

		selections->Add(new AliPP13QualityPhotonSelection("Qual", "Cluster quality Selection", cuts_pi0));

		if (suff.Contains("Only") && IsJetJetMC(description, isMC))
			selections->Add(new AliPP13PythiaInfoSelection("PythiaInfo", "Cross section and ntrials for a pthard bin."));
	}

	// Setup task
	AliAnalysisTaskPP13 * task = new AliAnalysisTaskPP13("PhosProtons", selections);

	if ( !badmap.IsNull() )
		task->SetBadMap(badmap);

	task->SelectCollisionCandidates(offlineTriggerMask);
	mgr->AddTask(task);


	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	AliAnalysisDataContainer * coutput = 0;
	for (Int_t i = 0; i < task->GetSelections()->GetEntries(); ++ i)
	{
		AliPP13PhotonSelection * fSel = dynamic_cast<AliPP13PhotonSelection *> (task->GetSelections()->At(i));
		fSel->SetTitle(description);
		cout << fSel->GetTitle() << endl;

		coutput = mgr->CreateContainer(fSel->GetName() + suff,
		                               TList::Class(),
		                               AliAnalysisManager::kOutputContainer,
		                               AliAnalysisManager::GetCommonFileName());
		mgr->ConnectOutput(task, i + 1, coutput);
	}
}

Bool_t IsJetJetMC(TString description, Bool_t isMC)
{
	if (!isMC)
		return kFALSE;

	if (description.Contains("Jet-Jet"))
		return kTRUE;


	cout << "Not Using PythiaInfo!!! " << endl;

	// Don't include Jet-Jet counters by default
	return kFALSE;
}