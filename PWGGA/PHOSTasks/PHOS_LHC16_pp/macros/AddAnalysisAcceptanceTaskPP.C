AliAnalysisTaskPP13 * AddAnalysisAcceptanceTaskPP(
	Bool_t isMC = kFALSE,
	TString description = "",
	TString suff = "",
	Int_t minDistanceMaximum = 5,
	Float_t scale = 1 // Scale parameter in cm
)
{
	AliAnalysisManager * manager = AliAnalysisManager::GetAnalysisManager();

	// Setup Analysis Selections
	//
	TList * selections = new TList();

	AliPP13SelectionWeights & data_weights = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kData);
	for (Int_t i = 0; i < minDistanceMaximum; ++i)
	{
		AliPP13ClusterCuts cuts_pi0 = AliPP13ClusterCuts::GetClusterCuts();
		cuts_pi0.fMinimalDistance = i * scale;

		AliPP13ClusterCuts cuts_eta = AliPP13ClusterCuts::GetClusterCuts();
		cuts_eta.fAsymmetryCut = 0.7;
		cuts_eta.fMinimalDistance = i * scale;

		// TODO: Add plain selections
		selections->Add(new AliPP13SpectrumSelectionSimple(Form("Phys%d", i), "Physics Selection", cuts_pi0, &data_weights));
		selections->Add(new AliPP13SpectrumSelectionSimple(Form("Eta%d", i), "Physics Selection for eta meson", cuts_eta, &data_weights));
	}
	delete &data_weights;

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
