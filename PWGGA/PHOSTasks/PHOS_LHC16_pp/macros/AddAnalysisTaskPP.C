AliAnalysisTaskPP13 * AddAnalysisTaskPP(
	Bool_t isMC = kFALSE,
	TString description = "",
	TString suff = ""
)
{
	// Restore the analysis manager
	//
	AliAnalysisManager * manager = AliAnalysisManager::GetAnalysisManager();

	// Setup Analysis Selections
	//
	TList * selections = new TList();

	AliPP13ClusterCuts cuts_pi0 = AliPP13ClusterCuts::GetClusterCuts();
	AliPP13ClusterCuts cuts_eta = AliPP13ClusterCuts::GetClusterCuts();
	cuts_eta.fAsymmetryCut = 0.7;

	AliPP13ClusterCuts cuts_no_tof = AliPP13ClusterCuts::GetClusterCuts();
	cuts_no_tof.fTimingCut = 999;

	AliPP13ClusterCuts cuts_wide_tof = AliPP13ClusterCuts::GetClusterCuts();
	cuts_wide_tof.fTimingCut = 30.0e-9;


	if (!isMC)
	{
		// TODO: Add plain selections
		AliPP13SelectionWeights & data_weights = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kData);
		AliPP13SelectionWeights & data_weights_plain = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kPlain);

		selections->Add(new AliPP13SpectrumSelection("Phys", "Physics Selection", cuts_pi0, &data_weights));
		selections->Add(new AliPP13SpectrumSelection("PhysWide", "Physics Selection, 30 ns cut", cuts_wide_tof, &data_weights));
		// selections->Add(new AliPP13SpectrumSelection("PhysPlain", "Physics Selection no TOF cut efficiency", cuts_pi0, &data_weights_plain));
		// selections->Add(new AliPP13PhotonTimecutStudySelection("Time", "Testing Timing Selection", cuts_pi0, &data_weights));

		selections->Add(new AliPP13SpectrumSelection("Eta", "Physics Selection for eta meson", cuts_eta, &data_weights));
		// selections->Add(new AliPP13SpectrumSelection("EtaPlain", "Physics Selection for eta meson no TOF cut efficiency", cuts_eta, &data_weights_plain));
		selections->Add(new AliPP13TagAndProbeSelection("TagAndProbleTOF", "Timing cut 12.5 ns", cuts_pi0, &data_weights_plain));
		selections->Add(new AliPP13TagAndProbeSelection("TagAndProbleTOFWide", "Timing cut 30 ns", cuts_wide_tof, &data_weights_plain));
		selections->Add(new AliPP13QualityPhotonSelection("Qual", "Cluster quality Selection", cuts_no_tof, &data_weights));
		selections->Add(new AliPP13EpRatioSelection("EpRatio", "E/p ratio selection for electrons", cuts_pi0, &data_weights));

		delete &data_weights;
		delete &data_weights_plain;
	}


	if (isMC)
	{
		AliPP13SelectionWeightsMC & mc_weights = dynamic_cast<AliPP13SelectionWeightsMC &>(AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kMC));

		selections->Add(new AliPP13EfficiencySelectionMC("PhysEff", "Physics efficiency for neutral particles fully corrected", cuts_pi0, &mc_weights));
		selections->Add(new AliPP13EfficiencySelectionMC("PhysEffEta", "Physics efficiency for neutral particles with asymmetry fully corrected", cuts_eta, &mc_weights));
		selections->Add(new AliPP13QualityPhotonSelection("Qual", "Cluster quality Selection", cuts_pi0, &mc_weights));
		selections->Add(new AliPP13MesonSelectionMC("MCStudy", "MC Selection with timing cut", cuts_pi0, &mc_weights));
		selections->Add(new AliPP13EpRatioSelection("EpRatio", "E/p ratio selection for electrons", cuts_pi0, &mc_weights));

		AliPP13SelectionWeightsMC & mc_weights_feeddown = dynamic_cast<AliPP13SelectionWeightsMC &>(AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kFeeddown));
		selections->Add(new AliPP13FeeddownSelection("FeeddownSelection", "FeeddownSelection", cuts_pi0, &mc_weights_feeddown));

		delete & mc_weights;
		delete & mc_weights_feeddown;
	}

	// Setup the main task
	//
	AliAnalysisTaskPP13 * task = new AliAnalysisTaskPP13("PhosProtons", selections);
	manager->AddTask(task);
	manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
	AliAnalysisDataContainer * coutput = 0;
	for (Int_t i = 0; i < task->GetSelections()->GetEntries(); ++ i)
	{
		AliPP13PhysicsSelection * fSel = dynamic_cast<AliPP13PhysicsSelection *> (task->GetSelections()->At(i));
		fSel->SetTitle(description);
		cout << fSel->GetTitle() << endl;

		coutput = manager->CreateContainer(
			fSel->GetName() + suff,
			TList::Class(),
			AliAnalysisManager::kOutputContainer,
			AliAnalysisManager::GetCommonFileName()
		);

		manager->ConnectOutput(task, i + 1, coutput);
	}
	return task;
}
