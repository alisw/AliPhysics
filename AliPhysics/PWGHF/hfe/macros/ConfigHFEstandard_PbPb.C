AliAnalysisTaskHFE* ConfigHFEstandard_PbPb(Bool_t useMC){
  //
  // HFE standard task configuration
  //

  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(120);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);

  hfecuts->SetMinNClustersITS(4);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);

  //hfecuts->UnsetVertexRequirement();
  hfecuts->SetMaxImpactParam(3.,3.);

  hfecuts->SetVertexRange(10.);
  //hfecuts->SetMaxChi2perClusterITS(36);
  //hfecuts->SetSigmaToVertex(10);
  hfecuts->SetTOFPIDStep(kTRUE);
  //hfecuts->SetTOFMISMATCHStep(kTRUE);
  //hfecuts->SetTPCPIDCleanUpStep(kTRUE);
  hfecuts->SetQAOn();
  //hfecuts->SetMinNTrackletsTRD(5);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("HFEanalysisStandard");
  task->SetHFECuts(hfecuts);
  task->SetPbPbAnalysis(kTRUE);
  //task->SetRemovePileUp(kTRUE);
  task->GetPIDQAManager()->SetHighResolutionHistos();

  // Define Variables
  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt");
  vm->AddVariable("eta");
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");
  vm->AddVariable("centrality");

  if(!useMC){

      for(Int_t a=0;a<12;a++)
      {
	  TF1 *hBackground = new TF1("hadronicBackgroundFunction","TMath::Exp([0]/x + [1])", 0., 20.);
	  hBackground->SetParameter(0, -43.87);
	  hBackground->SetParameter(1, 2.85);
	  task->SetBackGroundFactorsFunction(hBackground,a);
      }


  }

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);
  pid->AddDetector("TOF", 0);
  pid->AddDetector("TPC", 1);
//  pid->ConfigureTPCrejection();

  if(!useMC){

      Double_t params_centr_0_5[1];
      Double_t params_centr_5_10[1];
      Double_t params_centr_10_20[1];
      Double_t params_centr_20_30[1];
      Double_t params_centr_per[1];
      params_centr_0_5[0]=0.16;  // cut tuned for 0-10%
      params_centr_5_10[0]=0.16; // cut tuned for 0-10%
      params_centr_10_20[0]=0.29;
      params_centr_20_30[0]=0.38;
      params_centr_per[0]=0.44;
      char *cutmodel;
      cutmodel="pol0";


      for(Int_t a=0;a<11;a++)
      {
	  if(a>3)  pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_per,3.0);
	  if(a==0) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_0_5,3.0);    //  0-5%
	  if(a==1) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_5_10,3.0);    //  5-10%
	  if(a==2) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_10_20,3.0);    //  10-20%
	  if(a==3) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_20_30,3.0);    //  20-30%
      }


  }


  // V0 tagged tracks
  AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
  v0trackCuts->CreateStandardCuts();
  v0trackCuts->SetMinNClustersTPC(120);
  v0trackCuts->SetMinRatioTPCclusters(0.6);
  v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  v0trackCuts->SetMinNClustersITS(1);
  v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  v0trackCuts->SetCheckITSLayerStatus(kFALSE);
  v0trackCuts->UnsetVertexRequirement();
  //v0trackCuts->SetMaxChi2perClusterITS(36);
  //hfecuts->SetSigmaToVertex(10);
  v0trackCuts->SetTOFPIDStep(kTRUE);
//  v0trackCuts->SetTOFMISMATCHStep(kTRUE);
  //v0trackCuts->SetTPCPIDCleanUpStep(kTRUE);
  v0trackCuts->SetQAOn();
  
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
  task->SetTaggedTrackCuts(v0trackCuts);
  task->SetCleanTaggedTrack(kFALSE);
  
  // QA
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
//  task->SetFillSignalOnly(kFALSE);    // for DE pluging for MC
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);



  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
