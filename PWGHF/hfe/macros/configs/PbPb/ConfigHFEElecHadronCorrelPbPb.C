AliAnalysisTaskElecHadronCorrel* ConfigHFEElecHadronCorrelPbPb(Bool_t useMC,
                                                               Bool_t EventTrigSelMB=kTRUE,
                                                               Bool_t TrigSelCen = kTRUE,
                                                               Double_t CentMin = 0,
                                                               Double_t CentMax = 7,
                                                               Double_t TPCNsigMinE = -2,
                                                               Double_t TPCNsigMaxE = 2,
                                                               Double_t TPCNsigMinH = -10,
                                                               Double_t TPCNsigMaxH = -3.5,
                                                               Double_t SSM02Min = 0.03,  
                                                               Double_t SSM02Max = 0.5,
                                                               Double_t SSM20Min = 0.03,  
                                                               Double_t SSM20Max = 0.3,
                                                               Double_t Disp = 1,
                                                               Double_t EovPMin = 0.8,    
                                                               Double_t EovPMax = 1.2,    
                                                               Double_t InvM = 0.1,
                                                               TString ContNameExt = "Central",
                                                               TString TaskName="hfeCorrl"){

//
  // HFE standard task configuration
  //

  Bool_t kAnalyseTaggedTracks = kTRUE;
  
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsEMCAL","HFE Standard Cuts");
//  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(100);
  hfecuts->SetMinNClustersITS(3);
  hfecuts->SetMinNTrackletsTRD(0);
  hfecuts->SetMinRatioTPCclusters(0.6);
   // hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetRequireITSPixel();
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny); 
  hfecuts->SetMaxChi2perClusterITS(-1);
  hfecuts->SetMaxChi2perClusterTPC(3.5);
  hfecuts->SetCheckITSLayerStatus(kFALSE); // shud be put back
//  hfecuts->UnsetVertexRequirement();
  hfecuts->SetVertexRange(10.);
  hfecuts->SetRequireSigmaToVertex();
  //hfecuts->SetSigmaToVertex(10);
  hfecuts->SetTOFPIDStep(kFALSE);
//  hfecuts->SetQAOn();
 hfecuts->SetPtRange(0, 30);

 TString taskName = TaskName;

  AliAnalysisTaskElecHadronCorrel *task = new AliAnalysisTaskElecHadronCorrel(taskName);
  printf("task ------------------------ %p\n ", task);
  task->SetHFECuts(hfecuts);
//  task->SetRemovePileUp(kTRUE);
//  task->SetInvariantMassCut(0.1);

  task->SetEventTriggerSelectionMB(EventTrigSelMB);
  task->SetTriggerSelection(TrigSelCen);
  task->SetCentralityParameters(CentMin, CentMax, "V0M");
  task->SetInvariantMassCut(InvM);
  task->SetTPCnsigmaCutsElecSelection(TPCNsigMinE,TPCNsigMaxE);
  task->SetTPCnsigmaCutsHadSelection(TPCNsigMinH,TPCNsigMaxH);
  task->SetShowerShapeCutsM02(SSM02Min,SSM02Max);
  task->SetShowerShapeCutsM20(SSM20Min,SSM20Max);
  task->SetShowerShapeCutsDisp(0,Disp);
  task->SetEovPCuts(EovPMin,EovPMax);
  task->SetRejectKinkMother(kTRUE);

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TPC", 0);
  pid->AddDetector("EMCAL", 1);
  /*
  // change E/p cuts
  AliHFEpidEMCAL *emcpid = pid->AliHFEpid::GetDetPID(AliHFEpid::kEMCALpid);
  emcpid->SetEoPMax(1.2);
  emcpid->SetEoPMim(0.8);

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
	  params[0] = -1.0; //sigma min
  }
//  pid->ConfigureTPCdefaultCut(cutmodel, params,3.0); 
  pid->ConfigureTPCasymmetric(0,30,-1,3.0); 
*/
  printf("*************************************\n");
  printf("Configuring standard Task:\n");
//  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
