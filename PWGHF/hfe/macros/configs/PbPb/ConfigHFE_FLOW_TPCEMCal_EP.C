AliAnalysisTaskFlowTPCEMCalEP* ConfigHFE_FLOW_TPCEMCal_EP(Bool_t useMC){
  //
  // HFE standard task configuration
  //

  Bool_t kAnalyseTaggedTracks = kTRUE;
  
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsEMCAL","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(100);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetMinNClustersITS(3);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetVertexRange(10.);
  hfecuts->SetTOFPIDStep(kFALSE);
  hfecuts->SetPtRange(2, 50);
  hfecuts->SetMaxImpactParam(1,2);
  
  AliAnalysisTaskFlowTPCEMCalEP *task = new AliAnalysisTaskFlowTPCEMCalEP("HFE v2");
  printf("task ------------------------ %p\n ", task);
  task->SetHFECuts(hfecuts);
  task->SetInvariantMassCut(0.05);

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TPC", 0);
  pid->AddDetector("EMCAL", 1);
  // change E/p cuts
  AliHFEpidEMCAL *emcpid = pid->AliHFEpid::GetDetPID(AliHFEpid::kEMCALpid);
  emcpid->SetEoPMax(1.2);
  emcpid->SetEoPMim(0.9);

  Double_t params[4];
  char *cutmodel;
  if(useMC){
	  // Monte-Carlo needs modelling of the falling mean with momentum at low momentum
	  // for high momentum it is consistent with a flat -0.94
	  cutmodel = "[0]*TMath::Exp([1]*x) + [2] + [3]*x";
	  Double_t paramsMC[4] = {0.7174, -1.588, -0.9395, 0.0246};
	  for(int ipar = 0; ipar < 4; ipar++) params[ipar] = paramsMC[ipar];
  } else {
	  // Data is consistent with a flat 0.12
	  cutmodel = "pol0";
	  //params[0] = -0.0015;
	  //params[0] = -3.0;
	  //params[0] = -0.05; //sigma min
	  params[0] = -1.0; //sigma min
  }
  pid->ConfigureTPCdefaultCut(cutmodel, params,3.0); 

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
//  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
