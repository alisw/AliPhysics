AliAnalysisTaskFlowTPCEMCalEP* ConfigHFE_FLOW_TPCEMCal_EP(Bool_t useMC, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut,Double_t AssPtCut, Int_t AssTPCnCut, Int_t  ITSncut, Bool_t AssITSrefitCut, Int_t TPCnCut, Bool_t UseNewEP, Bool_t UseTender,Int_t period, TString passV0, TString passTPC,Bool_t TimeCut,Bool_t WeightSyst,Bool_t SystTOFcut, Double_t CutM02, Double_t CutM20, Bool_t SScut, Bool_t EnablePileupRejVZEROTPCout){
  //
  // HFE standard task configuration
  //

  Bool_t kAnalyseTaggedTracks = kTRUE;
  
  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsEMCAL","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCnCut);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetMaxChi2perClusterTPC(3.5);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  //hfecuts->SetMinNClustersITS(ITSncut);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetVertexRange(10.);
  hfecuts->SetTOFPIDStep(kFALSE);
  hfecuts->SetPtRange(1.5, 50);
  hfecuts->SetMaxImpactParam(2.4,3.2); // radial, z
  
  AliAnalysisTaskFlowTPCEMCalEP *task = new AliAnalysisTaskFlowTPCEMCalEP("HFE v2");
  printf("task ------------------------ %p\n ", task);
  task->SetHFECuts(hfecuts);
  task->SetAssPtCut(AssPtCut);
  task->SetAssTPCnCut(AssTPCnCut);
  task->SetTPCnCut(TPCnCut);
  task->SetITSnCut(ITSncut);
  task->SetAssITSrefitCut(AssITSrefitCut);
  task->SetPeriod(period);
  task->SetEP(UseNewEP);
  task->SetTender(UseTender);
  task->SetCorrectionPassV0(passV0);
  task->SetCorrectionPassTPC(passTPC);
  task->SetSigmaITScut(SigmaITScut);
  task->SetSigmaTOFcut(SigmaTOFcut);
  task->SetSigmaTPCcut(SigmaTPCcut);
  task->SetTimeCut(TimeCut);
  task->SetWeightSyst(WeightSyst);
  task->SetSystTOFcut(SystTOFcut);
  task->SetCutM02(CutM02);
  task->SetCutM20(CutM20);
  task->SetSScut(SScut);
  task->SetPileUpCut(EnablePileupRejVZEROTPCout);
    
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
