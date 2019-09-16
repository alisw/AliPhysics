AliAnalysisTaskPP13 * AddTOFAnalysisTask(
	Bool_t isMC = kFALSE,
	TString description = "",
	TString suff = ""
)
{
	// Setup Analysis Selections
	//

	AliPP13ClusterCuts cuts_pi0 = AliPP13ClusterCuts::GetClusterCuts();
	AliPP13SelectionWeights & data_weights_plain = AliPP13SelectionWeights::Init(AliPP13SelectionWeights::kPlain);

	TList * selections = new TList();
	selections->Add(new AliPP13TagAndProbeSelection("TagAndProbleTOF", "Cluster p_{T} Selection", cuts_pi0, &data_weights_plain));
	delete &data_weights_plain;

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
