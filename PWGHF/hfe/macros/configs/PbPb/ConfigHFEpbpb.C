AliAnalysisTaskHFE* ConfigHFEpbpb(Bool_t useMC=kFALSE, Bool_t beauty=kFALSE,
				  UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
				  Double_t TPCclRatio = 0.6, Double_t TPCclshared = 1.1,
				  UChar_t ITScl=3,  Double_t ITSchi2perclusters=99999999.,
                                  Double_t dcaxy=1.0, Double_t dcaz=2.0,
				  Double_t TOFs=3.,Double_t IpSig=3., Int_t itspixelcut=AliHFEextraCuts::kFirst, TString appendix,
				  Float_t prodlow=0., Float_t prodhigh=100., Int_t addflag=0., Int_t ptbin=0,
				  Int_t nondefaultcentr=0, Float_t* arraycentr=NULL,
				  Double_t* tpcdEdxcut=NULL,Double_t tpcu=3.0){
  //
  // HFE standard task configuration
  //

  AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts pbpb TOF TPC");
  hfecuts->CreateStandardCuts();

  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinRatioTPCclusters(TPCclRatio);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetFractionOfSharedTPCClusters(TPCclshared);

  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMaxChi2perClusterITS(ITSchi2perclusters);
  hfecuts->SetCutITSpixel(itspixelcut);
  hfecuts->SetCheckITSLayerStatus(kFALSE);

  hfecuts->SetIPcutParam(0,0,0,IpSig,kTRUE,kTRUE);
  // if(useMC && beauty) hfecuts->SetProductionVertex(prodlow,prodhigh,prodlow,prodhigh);
  if(useMC) hfecuts->SetProductionVertex(prodlow,prodhigh,prodlow,prodhigh);

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
  task->SetRemovePileUp(kFALSE);
  task->GetPIDQAManager()->SetHighResolutionHistos();
 
  
  if((nondefaultcentr!=0) && arraycentr) {
    for(Int_t i=0;i<12;i++)
      {
	task->SetPbPbUserCentralityLimit(kTRUE);
	task->SetPbPbUserCentralityArray(i,arraycentr[i]);
	
      }
  }
  
  // Define Variables

  if(ptbin==0)
  {
      Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
      1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5.,
      5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  }

  if(ptbin==1)
  {
      Double_t ptbinning[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
  }


  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 
			     0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  AliHFEvarManager *vm = task->GetVarManager();
  Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
  Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;

  //  printf("ptbinning: %i \n",sizept);

  //vm->AddVariable("pt");
  //vm->AddVariable("eta");
  vm->AddVariable("pt", sizept, ptbinning);
  // vm->AddVariable("pt", 18, ptbinning);
  //  vm->AddVariable("eta", 16, etabinning);
  vm->AddVariable("eta", sizeeta, etabinning);
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
    
    // 0-10% 0.16
    // 10-20% 0.29
    // 20-30% 0.38
    // 30-40% 0.44
    Double_t paramsTPCdEdxcut[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    if((nondefaultcentr!=0) && tpcdEdxcut) memcpy(paramsTPCdEdxcut,tpcdEdxcut,sizeof(paramsTPCdEdxcut));
    
    char *cutmodel;
    cutmodel="pol0";
    
    for(Int_t a=0;a<11;a++)
    {
     //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
	Double_t tpcparam[1]={paramsTPCdEdxcut[a]};

	pid->ConfigureTPCcentralityCut(a,cutmodel,tpcparam,tpcu);
    }

  }
  pid->ConfigureTOF(TOFs);

  AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
  if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs);

  // QA
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  //task->SetFillSignalOnly(kFALSE);    // for DE pluging for MC
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
  if(useMC && addflag==1) task->SetDebugStreaming();

  printf("*************************************\n");
  printf("Configuring task PbPb \n"); 
  //if(isLHC10) printf("Configuring TPC1 Task 2010 :\n");
  //if(isLHC11) printf("Configuring TPC1 Task 2011 :\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
