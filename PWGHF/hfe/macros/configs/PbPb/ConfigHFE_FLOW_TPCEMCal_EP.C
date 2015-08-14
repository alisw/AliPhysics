AliAnalysisTaskFlowTPCEMCalEP* ConfigHFE_FLOW_TPCEMCal_EP(Bool_t useMC, 
 Double_t fP2_lowPtEta0010 = 1., Double_t fP3_lowPtEta0010 = 1., Double_t fP4_lowPtEta0010 = 1., 
 Double_t fP2_highPtEta0010 = 1., Double_t fP3_highPtEta0010 = 1., Double_t fP4_highPtEta0010 = 1., 
 Double_t fP2_lowPtPi00010 = 1., Double_t fP3_lowPtPi00010 = 1., Double_t fP4_lowPtPi00010 = 1., 
 Double_t fP2_highPtPi00010 = 1., Double_t fP3_highPtPi00010 = 1., Double_t fP4_highPtPi00010 = 1., 
 Double_t fP2_lowPtEta1020 = 1., Double_t fP3_lowPtEta1020 = 1., Double_t fP4_lowPtEta1020 = 1., 
 Double_t fP2_highPtEta1020 = 1., Double_t fP3_highPtEta1020 = 1., Double_t fP4_highPtEta1020 = 1., 
 Double_t fP2_lowPtPi01020 = 1., Double_t fP3_lowPtPi01020 = 1., Double_t fP4_lowPtPi01020 = 1., 
 Double_t fP2_highPtPi01020 = 1., Double_t fP3_highPtPi01020 = 1., Double_t fP4_highPtPi01020 = 1., 
 Double_t fP2_lowPtEta2040 = 1., Double_t fP3_lowPtEta2040 = 1., Double_t fP4_lowPtEta2040 = 1., 
 Double_t fP2_highPtEta2040 = 1., Double_t fP3_highPtEta2040 = 1., Double_t fP4_highPtEta2040 = 1., 
 Double_t fP2_lowPtPi02040 = 1., Double_t fP3_lowPtPi02040 = 1., Double_t fP4_lowPtPi02040 = 1., 
 Double_t fP2_highPtPi02040 = 1., Double_t fP3_highPtPi02040 = 1., Double_t fP4_highPtPi02040 = 1.)
{
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
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetVertexRange(10.);
  hfecuts->SetTOFPIDStep(kFALSE);
  hfecuts->SetPtRange(1.5, 50);
  hfecuts->SetMaxImpactParam(1,2);
  
  AliAnalysisTaskFlowTPCEMCalEP *task = new AliAnalysisTaskFlowTPCEMCalEP("HFE v2");
  printf("task ------------------------ %p\n ", task);
  task->SetHFECuts(hfecuts);
  
  task->SetP2_lowPtEta0010(fP2_lowPtEta0010); 
  task->SetP3_lowPtEta0010(fP3_lowPtEta0010); 
  task->SetP4_lowPtEta0010(fP4_lowPtEta0010); 

  task->SetP2_highPtEta0010(fP2_highPtEta0010); 
  task->SetP3_highPtEta0010(fP3_highPtEta0010); 
  task->SetP4_highPtEta0010(fP4_highPtEta0010); 

  task->SetP2_lowPtPi00010(fP2_lowPtPi00010); 
  task->SetP3_lowPtPi00010(fP3_lowPtPi00010); 
  task->SetP4_lowPtPi00010(fP4_lowPtPi00010); 

  task->SetP2_highPtPi00010(fP2_highPtPi00010); 
  task->SetP3_highPtPi00010(fP3_highPtPi00010); 
  task->SetP4_highPtPi00010(fP4_highPtPi00010); 

  task->SetP2_lowPtEta1020(fP2_lowPtEta1020); 
  task->SetP3_lowPtEta1020(fP3_lowPtEta1020); 
  task->SetP4_lowPtEta1020(fP4_lowPtEta1020); 

  task->SetP2_highPtEta1020(fP2_highPtEta1020); 
  task->SetP3_highPtEta1020(fP3_highPtEta1020); 
  task->SetP4_highPtEta1020(fP4_highPtEta1020); 

  task->SetP2_lowPtPi01020(fP2_lowPtPi01020); 
  task->SetP3_lowPtPi01020(fP3_lowPtPi01020); 
  task->SetP4_lowPtPi01020(fP4_lowPtPi01020); 

  task->SetP2_highPtPi01020(fP2_highPtPi01020); 
  task->SetP3_highPtPi01020(fP3_highPtPi01020); 
  task->SetP4_highPtPi01020(fP4_highPtPi01020); 

  task->SetP2_lowPtEta2040(fP2_lowPtEta2040); 
  task->SetP3_lowPtEta2040(fP3_lowPtEta2040); 
  task->SetP4_lowPtEta2040(fP4_lowPtEta2040); 

  task->SetP2_highPtEta2040(fP2_highPtEta2040); 
  task->SetP3_highPtEta2040(fP3_highPtEta2040); 
  task->SetP4_highPtEta2040(fP4_highPtEta2040); 

  task->SetP2_lowPtPi02040(fP2_lowPtPi02040); 
  task->SetP3_lowPtPi02040(fP3_lowPtPi02040); 
  task->SetP4_lowPtPi02040(fP4_lowPtPi02040); 

  task->SetP2_highPtPi02040(fP2_highPtPi02040); 
  task->SetP3_highPtPi02040(fP3_highPtPi02040); 
  task->SetP4_highPtPi02040(fP4_highPtPi02040); 

  

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
//   pid->AddDetector("TPC", 0);

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
