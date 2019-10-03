void AddTask_jklein_gtusim()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskTRDgtuSim *task = 0x0;

  // Int_t deltaAlpha[] = { 11, 16, 21 };
  Int_t deltaAlpha[] = { 21 };
  Int_t nDeltaAlpha = sizeof(deltaAlpha) / sizeof(deltaAlpha[0]);
  // Int_t deltaY[] = { 9, 12, 15, 21, 27 };
  Int_t deltaY[] = { 18 };
  Int_t nDeltaY = sizeof(deltaY) / sizeof(deltaY[0]);

  for (Int_t iDeltaY = 0; iDeltaY < nDeltaY; ++iDeltaY)
    for (Int_t iDeltaAlpha = 0; iDeltaAlpha < nDeltaAlpha; ++iDeltaAlpha) {
      Int_t label = - deltaY[iDeltaY] * 100 - deltaAlpha[iDeltaAlpha];

      TString name = TString::Format("gtusim_%i_%i", deltaY[iDeltaY], deltaAlpha[iDeltaAlpha]);

      AliAnalysisTaskTRDgtuSim *task = new AliAnalysisTaskTRDgtuSim(name.Data());
      task->SetDeltaY(deltaY[iDeltaY]);
      task->SetDeltaAlpha(deltaAlpha[iDeltaAlpha]);
      task->SetLabel(label);

      mgr->AddTask(task);

      AliAnalysisDataContainer *coutput =
	mgr->CreateContainer(Form("hist_%s", name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,
			     Form("%s:TRD_gtusim", AliAnalysisManager::GetCommonFileName()));

      mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(task, 1, coutput);
    }

  return;
}
