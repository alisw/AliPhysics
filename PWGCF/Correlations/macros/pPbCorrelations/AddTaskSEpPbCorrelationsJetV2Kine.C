AliAnalysisTaskSEpPbCorrelationsJetV2Kine *AddTaskSEpPbCorrelationsJetV2Kine(TString sMode = "TPCTPCFMDA", TString sNameList = "0_10", TString sEst = "V0A")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskKineAmpt.C::AddTaskKineAmpt", "No analysis manager to connect to");
    return nullptr;
  }
//=============================================================================

  AliAnalysisTaskSEpPbCorrelationsJetV2Kine *task = new AliAnalysisTaskSEpPbCorrelationsJetV2Kine("AliAnalysisTaskSEpPbCorrelationsJetV2Kine");

  Double_t trigPtLimits[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0};
  //Double_t trigPtLimits[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
  //Double_t trigPtLimits[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0};
  const Int_t nBinTrigPt = sizeof(trigPtLimits) / sizeof(Double_t) - 1;
  task->SetTrigPtBinning(nBinTrigPt, trigPtLimits);

  Double_t assocPtLimits[] = {0.5, 1., 1.5, 50.};
  const Int_t nBinAssocPt = sizeof(assocPtLimits) / sizeof(Double_t) - 1;
  task->SetAssocPtBinning(nBinAssocPt, assocPtLimits);

  task->SetMuonDecay("B");
  task->SetTPCEtaLimits(1.5); // Less than 2

  task->SetAnaMode(sMode);
  task->SetAssoCut(1.0);
 
  task->SetPtMin(0.5);
  task->SetPtMax(20.);

  task->SetCen1(0);
  task->SetCen2(10);

  task->SetIsPbp(kTRUE);
  task->SetEst(sEst);

  task->SetCorrBeam(kTRUE);

  // Input File
  TGrid::Connect("alien://");
  TFile *file = TFile::Open("alien:///alice/cern.ch/user/s/sitang/AMPT_Centrality_Calibration/Centrality.root");

  //TFile *file = TFile::Open("./Centrality_502_Muonv2/Centrality.root"); // Old AMPT
  //TFile *file = TFile::Open("./Centrality/Centrality.root");  // w/o Hadronic scattering
  //TFile *file = TFile::Open("./Centrality_HS/Centrality.root"); // w/ Hardronic scattering
 // TFile *file = TFile::Open("./Centrality_816_MuonV2/Centrality.root"); // 8.16 TeV w/ Hardronic scattering
  //TFile *file = TFile::Open("./Centrality_816_MuonV2_Pbp/Centrality.root"); // 8.16 TeV w/ Hardronic scattering
  TH1D *h_Charge  = (TH1D*)file->Get(Form("hCharge%s",sEst.Data())); h_Charge->SetName(Form("hCharge%s_%s_%s",sEst.Data(), sMode.Data(),sNameList.Data())); h_Charge->SetDirectory(0);
  file->Close();

  // create input container
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer(Form("h_centrality_%s_%s",sMode.Data(), sNameList.Data()),
                                    TH1D::Class(),
                                    AliAnalysisManager::kInputContainer);
  cinput1->SetData(h_Charge);

  mgr->AddTask(task);
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput(task, 1, cinput1);
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("listKineAmpt_%s_%s",sMode.Data(),sNameList.Data()),
                                                   TList::Class(),
                                                   AliAnalysisManager::kOutputContainer,
                                                   AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
