AliAnalysisTaskSEpPbCorrelationsJetV2Kine *AddTaskSEpPbCorrelationsJetV2Kine(TString sMode = "TPCTPCFMDA")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskKineAmpt.C::AddTaskKineAmpt", "No analysis manager to connect to");
    return nullptr;
  }
//=============================================================================

  AliAnalysisTaskSEpPbCorrelationsJetV2Kine *task = new AliAnalysisTaskSEpPbCorrelationsJetV2Kine("AliAnalysisTaskSEpPbCorrelationsJetV2Kine");

  Double_t trigPtLimits[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0};
  const Int_t nBinTrigPt = sizeof(trigPtLimits) / sizeof(Double_t) - 1;
  task->SetTrigPtBinning(nBinTrigPt, trigPtLimits);

  Double_t assocPtLimits[] = {0.5, 1., 1.5, 50.};
  const Int_t nBinAssocPt = sizeof(assocPtLimits) / sizeof(Double_t) - 1;
  task->SetAssocPtBinning(nBinAssocPt, assocPtLimits);

  task->SetAnaMode(sMode);
  task->SetAssoCut(1.0);
 
  task->SetPtMin(0.5);
  task->SetPtMax(50.);

  task->SetCen1(0);
  task->SetCen2(10);


  // Input File
  TGrid::Connect("alien://");
  TFile *file = TFile::Open("alien:///alice/cern.ch/user/s/sitang/AMPT_Centrality_Calibration/Centrality.root");

  //TFile *file = TFile::Open("./Centrality/Centrality.root");
  TH1D *h_Charge  = (TH1D*)file->Get("hChargeV0A"); h_Charge->SetDirectory(0);
  file->Close();

  // create input container
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("h_centrality",
                                    TH1D::Class(),
                                    AliAnalysisManager::kInputContainer);
  cinput1->SetData(h_Charge);

  mgr->AddTask(task);
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput(task, 1, cinput1);
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("listKineAmpt_%s",sMode.Data()),
                                                   TList::Class(),
                                                   AliAnalysisManager::kOutputContainer,
                                                   AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
