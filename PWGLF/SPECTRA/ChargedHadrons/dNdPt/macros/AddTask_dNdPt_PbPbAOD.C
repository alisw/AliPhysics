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
  
  // Get Analysis Type, can be "ESD" or "AOD"
  //==============================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType();
  
  //==============================================================================
  // Create and configure the task
  //==============================================================================
  
  // create combined Name from basic name of container and suffix
  TString combinedName;
  combinedName.Form("%s%s",contName, suffix);

  // print combined name of task
  ::Info("AddTask_dNdPt_PbPbAOD",Form("Name of Task: %s", combinedName.Data()));
  
  // create Task itself
  AlidNdPtAnalysisPbPbAOD *task = new AlidNdPtAnalysisPbPbAOD(combinedName.Data());
  
  // select a given trigger mask
  task->SelectCollisionCandidates(uTriggerMask);
  
  // select centrality estimator
  task->SetCentralityEstimator(centEstimator);
  
  // enable possibility do deselect certain trigger strings
  //   char *cDisableTrigger = "";
  //   task->DisableOnlineTriggerStrings(cDisableTrigger);
  
  // select eventplane determinator
  task->SetEventplaneSelector("Q");

  // select filter bit for analysis
  task->SetFilterBit(iFilterBit);
  
  // enable use of hybrid tracks
  // task->SetRequireHybridTracking();
  
  // cut on minimum number of contributors to the vertex
  // task->SetNContributorsVertex(1);
  
  // cut on minimum number of crossed rows in the TPC
  task->SetCutMinNCrossedRowsTPC(dNCrossedRowsTPC);
  
  // cut on minimum number of clusters in the TPC
  task->SetCutMinNClustersTPC(dNClustersTPC);
  
  // enable cut on length in TPC, dependent on pT
  task->SetCutLengthInTPCPtDependent(bDoCutTPCLength);
  
  // prefactor for the minimum length cut in the TPC
  task->SetPrefactorLengthInTPCPtDependent(dPrefactorLengthInTPC);
  
  // enable relative cuts
  // task->EnableRelativeCuts();
  
  // use certain fraction of Nclusters in the TPC
  // task->SetCutPercMinNClustersTPC(1.0);
  
  // use certain fraction of NCrossedRows in the TPC
  // task->SetCutPercMinNCrossedRowsTPC(1.0);
  
  // define pT bins for cross checks
  Double_t binsPtCheck[] = {0., 0.15, 1.0, 5.0, 10.0, 50.0, 200.0};
  Int_t nBinPtCheck = sizeof(binsPtCheck)/sizeof(Double_t);
  task->SetBinsPtCheck(nBinPtCheck, binsPtCheck);
  
  // define pT bins for MC corrections
  Double_t binsPtCorr[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 200.0};
  Int_t nBinPtCorr = sizeof(binsPtCorr)/sizeof(Double_t);
  task->SetBinsPtCorr(nBinPtCorr, binsPtCorr);
  
  // define phi bins
  Double_t binsPhi[] = {
    0, 1./16.*TMath::Pi(), 2./16.*TMath::Pi(), 3./16.*TMath::Pi(), 4./16.*TMath::Pi(), 5./16.*TMath::Pi(), 6./16.*TMath::Pi(), 7./16.*TMath::Pi(), 8./16.*TMath::Pi()
  };
  Int_t nBinPhi = sizeof(binsPhi)/sizeof(Double_t);
  task->SetBinsDeltaphi(nBinPhi, binsPhi);
  
  // define Zv bins
  Double_t binsZv[] = {-30.,-10.,-8., -6., -4., -2., 0., 2., 4., 6., 8., 10., 30.};
  Int_t nBinsZv = sizeof(binsZv)/sizeof(Double_t);
  task->SetBinsZv(nBinsZv, binsZv);
  
  // print some information
  ::Info("AddTask_dNdPt_PbPbAOD",Form("CrossedRowCut set to %.0f", task->GetCutMinNCrossedRowsTPC()));
  ::Info("AddTask_dNdPt_PbPbAOD",Form("FilterBit set to %d", task->GetFilterBit()));
  
  // add task to manager
  mgr->AddTask(task);
  
  // create output container
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", combinedName.Data()), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
  
  // connect common input container to task
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  
  // connect created output container
  mgr->ConnectOutput(task, 1, coutput);
  
  // return task in case something needs to be changed from macro
  return task;
}   
