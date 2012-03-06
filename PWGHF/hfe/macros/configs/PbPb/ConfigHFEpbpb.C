AliAnalysisTaskHFE* ConfigHFEpbpb(Bool_t useMC=kFALSE, Bool_t beauty=kFALSE, UChar_t Sample=10,
				  UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
				  Double_t TPCclRatio = 0.6, Double_t TPCclshared = 1.1,
				  UChar_t ITScl=3,  Double_t ITSchi2perclusters=99999999.,
                                  Double_t dcaxy=1000.0, Double_t dcaz=2000.0,
				  Double_t TPCs=0., Double_t TPCu=3.09, 
				  Double_t TOFs=3.,Double_t IpSig=3., TString appendix){
  //
  // HFE standard task configuration
  //
    Bool_t kAnalyseTaggedTracks = kTRUE;

  printf("\n String settings: %s \n",appendix);

  AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts pbpb TOF TPC");
  hfecuts->CreateStandardCuts();

  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinRatioTPCclusters(TPCclRatio);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetFractionOfSharedTPCClusters(TPCclshared);

  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMaxChi2perClusterITS(ITSchi2perclusters);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);

  hfecuts->SetIPcutParam(0,0,0,IpSig,kTRUE);
  if(useMC && beauty) hfecuts->SetProductionVertex(0,100,0,100);

  hfecuts->SetMaxImpactParam(dcaxy,dcaz);

  // event cuts
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);

  // others
  hfecuts->SetTOFPIDStep(kTRUE);
  //hfecuts->SetMaxChi2perClusterITS(36);
  //hfecuts->SetTOFMISMATCHStep(kTRUE);
  //hfecuts->SetTPCPIDCleanUpStep(kTRUE);
  hfecuts->SetQAOn();

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(appendix);
  task->SetHFECuts(hfecuts);
  task->SetPbPbAnalysis(kTRUE);
  task->SetRemovePileUp(kTRUE);
  task->GetPIDQAManager()->SetHighResolutionHistos();
  
  // Define Variables
  Double_t ptbinning[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
  //Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
  //      		    1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5.,
  //      		    5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 
			     0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  AliHFEvarManager *vm = task->GetVarManager();
  //vm->AddVariable("pt");
  //vm->AddVariable("eta");
  // vm->AddVariable("pt", 35, ptbinning);
  vm->AddVariable("pt", 18, ptbinning);
  vm->AddVariable("eta", 16, etabinning);
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
  //pid->ConfigureTPCrejection();
  
  if(!useMC){
    
    Double_t params_centr_0_5[1];
    Double_t params_centr_5_10[1];
    Double_t params_centr_10_20[1];
    Double_t params_centr_20_30[1];
    Double_t params_centr_per[1];
    params_centr_0_5[0]=  0.16; // cut tuned for 0-10% // TPCs
    params_centr_5_10[0]= 0.16; // cut tuned for 0-10%
    params_centr_10_20[0]= 0.29;
    params_centr_20_30[0]= 0.38;
    params_centr_per[0]=   0.44;
    char *cutmodel;
    cutmodel="pol0";
    
    
    for(Int_t a=0;a<11;a++)
      {
	if(a>3)  pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_per,TPCu);
	if(a==0) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_0_5,TPCu);    //  0-5%
	if(a==1) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_5_10,TPCu);    //  5-10%
	if(a==2) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_10_20,TPCu);    //  10-20%
	if(a==3) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_20_30,TPCu);    //  20-30%
      }
  }
  pid->ConfigureTOF(TOFs);

  if(kAnalyseTaggedTracks)
  {
      // V0 tagged tracks
      AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
      v0trackCuts->CreateStandardCuts();

      v0trackCuts->SetMinNClustersTPC(TPCcl);
      v0trackCuts->SetMinNClustersTPCPID(TPCclPID);
      v0trackCuts->SetFractionOfSharedTPCClusters(TPCclshared);
      v0trackCuts->SetMinRatioTPCclusters(TPCclRatio);
      v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
      v0trackCuts->SetMinNClustersITS(1);
      v0trackCuts->SetMaxChi2perClusterITS(ITSchi2perclusters);
      v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
      v0trackCuts->SetCheckITSLayerStatus(kFALSE);
      v0trackCuts->UnsetVertexRequirement();
      //v0trackCuts->SetMaxChi2perClusterITS(36);
      //hfecuts->SetSigmaToVertex(10);
      v0trackCuts->SetTOFPIDStep(kTRUE);
      //v0trackCuts->SetTOFMISMATCHStep(kTRUE);
      //v0trackCuts->SetTPCPIDCleanUpStep(kTRUE);
      v0trackCuts->SetQAOn();

      task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
      task->SetTaggedTrackCuts(v0trackCuts);
      task->SetCleanTaggedTrack(kFALSE);
  }
  
  // QA
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  //task->SetFillSignalOnly(kFALSE);    // for DE pluging for MC
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
  if(useMC) task->SetDebugStreaming();

  printf("*************************************\n");
  printf("Configuring task PbPb \n"); 
  //if(isLHC10) printf("Configuring TPC1 Task 2010 :\n");
  //if(isLHC11) printf("Configuring TPC1 Task 2011 :\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
