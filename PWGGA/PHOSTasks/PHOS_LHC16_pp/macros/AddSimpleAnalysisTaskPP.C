AliAnalysisTaskPP13 * AddSimpleAnalysisTaskPP(
	Bool_t isMC = kFALSE,
	TString description = "",
	TString suff = ""
)
{
	// Setup Analysis Selections
	//
	TList * selections = new TList();

	AliPP13ClusterCuts cuts_pi0 = AliPP13ClusterCuts::GetClusterCuts();
	AliPP13ClusterCuts cuts_eta = AliPP13ClusterCuts::GetClusterCuts();
	cuts_eta.fAsymmetryCut = 0.7;


	AliPP13SelectionWeights & data_weights = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kData);

	selections->Add(new AliPP13SpectrumSelectionSimple("SimplePhys", "Physics Selection", cuts_pi0, &data_weights));
	selections->Add(new AliPP13SpectrumSelectionSimple("SimpleEta", "Physics Selection for eta meson", cuts_eta, &data_weights));

	delete &data_weights;

	// Setup the main task
	//
	AliAnalysisTaskPP13 * task = new AliAnalysisTaskPP13("PhosProtons", selections);

	AliAnalysisManager * manager = AliAnalysisManager::GetAnalysisManager();
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
