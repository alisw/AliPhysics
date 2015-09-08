AlidNdPtAnalysisPbPbAOD *AddTask_dNdPt_PbPbAOD( UInt_t uTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral , 
						Double_t dNCrossedRowsTPC = 100,
						Int_t iFilterBit = AliAODTrack::kTrkGlobal,
						char *contName = "dNdPtPbPbAOD",
						Double_t dNClustersTPC = 0,
						Bool_t bDoCutTPCLength = kTRUE,
						Double_t dPrefactorLengthInTPC = 0.85,
						char *centEstimator = "V0M",
						char *suffix = ""
											  )
{
  // Creates, configures and attaches to the train a cascades check task.
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_dNdPt_PbPbAOD", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_dNdPt_PbPbAOD", "This task requires an input event handler");
    return NULL;
  }   
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Create and configure the task
  TString combinedName;
  combinedName.Form("%s%s",contName, suffix);
  ::Info("AddTask_dNdPt_PbPbAOD",Form("Name of Task: %s", combinedName.Data()));
  AlidNdPtAnalysisPbPbAOD *task = new AlidNdPtAnalysisPbPbAOD(combinedName.Data());
  //   UInt_t triggerMask = AliVEvent::kMB;
  //   triggerMask |= AliVEvent::kCentral;
  //   triggerMask |= AliVEvent::kSemiCentral;
  
  task->SelectCollisionCandidates(uTriggerMask);
  
  task->SetCentralityEstimator(centEstimator);
  
//   char *cDisableTrigger = "",
//   
//   task->DisableOnlineTriggerStrings(cDisableTrigger);
  
  task->SetEventplaneSelector("Q");
  
  task->SetCutMinNCrossedRowsTPC(dNCrossedRowsTPC);
  task->SetCutMinNClustersTPC(dNClustersTPC);
  task->SetCutLengthInTPCPtDependent(bDoCutTPCLength);
  task->SetPrefactorLengthInTPCPtDependent(dPrefactorLengthInTPC);
  //task->SetCutLengthInTPCPtDependent();
  //
//   task->EnableRelativeCuts();
//   task->SetCutPercMinNClustersTPC(1.0);
//   task->SetCutPercMinNCrossedRowsTPC(1.0);
  
  Double_t binsPtCheck[] = {0., 0.15, 1.0, 5.0, 10.0, 50.0, 200.0}; 
  Int_t nBinPtCheck = sizeof(binsPtCheck)/sizeof(Double_t);
  task->SetBinsPtCheck(nBinPtCheck, binsPtCheck);
  
  Double_t binsPtCorr[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 200.0};
  Int_t nBinPtCorr = sizeof(binsPtCorr)/sizeof(Double_t);
  task->SetBinsPtCorr(nBinPtCorr, binsPtCorr);
  
//   Double_t binsPhi[] = {0 ,0.10472 ,0.20944 ,0.314159 ,0.418879 ,0.523599 ,0.628319 ,0.733038 ,0.837758 ,0.942478 ,1.0472 ,1.15192 ,1.25664 ,1.36136 ,1.46608 ,1.5708 ,1.67552 ,1.78024 ,1.88496 ,1.98968 ,2.0944 ,2.19911 ,2.30383 ,2.40855 ,2.51327 ,2.61799 ,2.72271 ,2.82743 ,2.93215 ,3.03687 ,3.14159 ,3.24631 ,3.35103 ,3.45575 ,3.56047 ,3.66519 ,3.76991 ,3.87463 ,3.97935 ,4.08407 ,4.18879 ,4.29351 ,4.39823 ,4.50295 ,4.60767 ,4.71239 ,4.81711 ,4.92183 ,5.02655 ,5.13127 ,5.23599 ,5.34071 ,5.44543 ,5.55015 ,5.65487 ,5.75959 ,5.86431 ,5.96903 ,6.07375 ,6.17847 ,6.28319};
  
//   Double_t binsPhi[] = {-1.*TMath::Pi(), -2.97625, -2.8109, -2.64555, -2.4802, -2.31486, -2.14951, -1.98416, -1.81882, -1.65347, -1.48812, -1.32278, -1.15743, -0.992082, -0.826735, -0.661388, -0.496041, -0.330694, -0.165347, 0, 0.165347, 0.330694, 0.496041, 0.661388, 0.826735, 0.992082, 1.15743, 1.32278, 1.48812, 1.65347, 1.81882, 1.98416, 2.14951, 2.31486, 2.4802, 2.64555, 2.8109, 2.97625, TMath::Pi()};
  
//   Double_t binsPhi[] = {-1.*TMath::Pi(), -0.75*TMath::Pi(), -0.5*TMath::Pi(), -0.25*TMath::Pi(), 0, 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0.75*TMath::Pi(), TMath::Pi()};
//   Double_t binsPhi[] = {
// 	-2.*TMath::Pi(), -15./8.*TMath::Pi(), -14./8.*TMath::Pi(), -13./8.*TMath::Pi(), -12./8.*TMath::Pi(), -11./8.*TMath::Pi(), -10./8.*TMath::Pi(), -9./8.*TMath::Pi(),
// 	-1.*TMath::Pi(), -7./8.*TMath::Pi(), -6./8.*TMath::Pi(), -5./8.*TMath::Pi(),  -4./8.*TMath::Pi(), -3./8.*TMath::Pi(), -2./8.*TMath::Pi(), -1./8.*TMath::Pi(),  0, 
// 	1./8.*TMath::Pi(), 2./8.*TMath::Pi(), 3./8.*TMath::Pi(), 4./8.*TMath::Pi(), 5./8.*TMath::Pi(), 6./8.*TMath::Pi(), 7./8.*TMath::Pi(), 1.*TMath::Pi(),
// 	10./8.*TMath::Pi(), 11./8.*TMath::Pi(), 12./8.*TMath::Pi(), 13./8.*TMath::Pi(), 14./8.*TMath::Pi(), 15./8.*TMath::Pi(), 2.*TMath::Pi()
//   };

//   Double_t binsPhi[] = {
// 	-2.*TMath::Pi(), -1.75*TMath::Pi(), -1.5*TMath::Pi(), -1.25*TMath::Pi(), 
// 	-1.*TMath::Pi(), -0.75*TMath::Pi(), -0.5*TMath::Pi(), -0.25*TMath::Pi(), 
// 	0, 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0.75*TMath::Pi(), TMath::Pi(),
// 	1.25*TMath::Pi(), 1.5*TMath::Pi(), 1.75*TMath::Pi(), 2.*TMath::Pi()
//   };
// Double_t binsPhi[] = {
//   0, 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0.75*TMath::Pi(), TMath::Pi()
// };
Double_t binsPhi[] = {
  0, 1./16.*TMath::Pi(), 2./16.*TMath::Pi(), 3./16.*TMath::Pi(), 4./16.*TMath::Pi(), 5./16.*TMath::Pi(), 6./16.*TMath::Pi(), 7./16.*TMath::Pi(), 8./16.*TMath::Pi()
};

// Double_t binsPhi[] = {
//   0, 1./8.*TMath::Pi(), 2./8.*TMath::Pi(), 3./8.*TMath::Pi(), 4./8.*TMath::Pi()
// };

  Int_t nBinPhi = sizeof(binsPhi)/sizeof(Double_t);
  task->SetBinsDeltaphi(nBinPhi, binsPhi);
  
//   Double_t binsZv[] = {-30.,-10.,-7.5, -5., -2.5 ,0.,2.5, 5., 7.5, 10.,30.};
  Double_t binsZv[] = {-30.,-10.,-8., -6., -4., -2., 0., 2., 4., 6., 8., 10., 30.};
  Int_t nBinsZv = sizeof(binsZv)/sizeof(Double_t);
  task->SetBinsZv(nBinsZv, binsZv);
    
  task->SetFilterBit(iFilterBit);
//   task->SetRequireHybridTracking();
//   task->SetNContributorsVertex(1);
  
  ::Info("AddTask_dNdPt_PbPbAOD",Form("CrossedRowCut set to %.0f", task->GetCutMinNCrossedRowsTPC()));
  ::Info("AddTask_dNdPt_PbPbAOD",Form("FilterBit set to %d", task->GetFilterBit()));
  
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", combinedName.Data()), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
  
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}   
