TF1* GetEtaCorrection(TString listname){
  
  TString etaMap="$ALICE_PHYSICS/PWGDQ/dielectron/files/EtaCorrMaps.root";
  
  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigHFEpbpb","Eta map not found: %s",etaMap.Data());
    return 0;
  }

  TFile f(etaMap.Data());
  if (!f.IsOpen()) return 0;
  gROOT->cd();
  TList *keys=f.GetListOfKeys();

  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(listname)){
      printf("Using Eta Correction Function: %s\n",kName.Data());
      return (TF1*)f.Get(kName.Data());
    }
  }
  return 0;
}


AliAnalysisTaskHFE* ConfigHFEpbpb(Bool_t isaod,
				  Bool_t useMC,
				  TString appendix,
				  Int_t aodfilter=-1,
				  Int_t clusterdef=AliHFEextraCuts::kFound, Int_t clusterrdef=AliHFEextraCuts::kFoundOverFindable, 
				  UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
				  Double_t TPCclRatio = 0.6, Double_t TPCclshared = 1.1,
				  Bool_t rejectkinkmother = kFALSE,
				  UChar_t ITScl=3,  Double_t ITSchi2perclusters=99999999.,
				  Int_t itspixelcut=AliHFEextraCuts::kBoth,
                                  Double_t dcaxy=1.0, Double_t dcaz=2.0,
				  Bool_t usetof=kFALSE,
				  Double_t TOFs=3.,
				  Bool_t etacor=kFALSE,TString listname="",
				  Double_t* tpcdEdxcutlow=NULL,Double_t* tpcdEdxcuthigh=NULL,
				  Float_t prodlow=0., Float_t prodhigh=100., 
				  Bool_t kNoPhotonic = kFALSE,
				  Int_t nondefaultcentr=0, Float_t* arraycentr=NULL,
				  Int_t ptbin=0){
  //
  // HFE standard task configuration
  //

  for(Int_t k=0; k < 12; k++) {
    printf("TPC dEdx cut low %f, high %f for %d\n",tpcdEdxcutlow[k],tpcdEdxcuthigh[k],k);
  }


  AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts pbpb TOF TPC");
  hfecuts->CreateStandardCuts();

  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinRatioTPCclusters(TPCclRatio);
  hfecuts->SetTPCmodes(clusterdef, clusterrdef);
  hfecuts->SetFractionOfSharedTPCClusters(TPCclshared);

  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMaxChi2perClusterITS(ITSchi2perclusters);
  hfecuts->SetCutITSpixel(itspixelcut);
  hfecuts->SetCheckITSLayerStatus(kFALSE);

 if(useMC) hfecuts->SetProductionVertex(prodlow,prodhigh,prodlow,prodhigh);
 
  hfecuts->SetMaxImpactParam(dcaxy,dcaz);

  // event cuts
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);

  // others
   if(usetof) hfecuts->SetTOFPIDStep(kTRUE);
  //hfecuts->SetMaxChi2perClusterITS(36);
  //hfecuts->SetTOFMISMATCHStep(kTRUE);
  //hfecuts->SetTPCPIDCleanUpStep(kTRUE);
  hfecuts->SetQAOn();

  // Background Subtraction
  AliHFENonPhotonicElectron *backe = new AliHFENonPhotonicElectron(appendix,"Background subtraction");
  AliESDtrackCuts *hfeBackgroundCuts = new AliESDtrackCuts();
  hfeBackgroundCuts->SetName("backgroundcuts");
  //hfeBackgroundCuts->SetAcceptKinkDaughters(kFALSE);
  hfeBackgroundCuts->SetRequireTPCRefit(kTRUE);
  hfeBackgroundCuts->SetRequireITSRefit(kTRUE);
  hfeBackgroundCuts->SetMinNClustersITS(2);
  hfeBackgroundCuts->SetEtaRange(-0.8,0.8);
  hfeBackgroundCuts->SetRequireSigmaToVertex(kTRUE);
  hfeBackgroundCuts->SetMaxChi2PerClusterTPC(3.5);
  hfeBackgroundCuts->SetMinNClustersTPC(100);
  hfeBackgroundCuts->SetPtRange(0.3,1e10);
  AliHFEpid *pidbackground = backe->GetPIDBackground();
  if(useMC) pidbackground->SetHasMCData(kTRUE);
  pidbackground->AddDetector("TPC", 1);
  pidbackground->ConfigureTPCasymmetric(0.0,9999.,-3.,3.);
  backe->SetHFEBackgroundCuts(hfeBackgroundCuts);

  // task
  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(appendix);
  task->SetHFECuts(hfecuts);
  task->SetHFEBackgroundSubtraction(backe);
  task->SetPbPbAnalysis(kTRUE);
  task->SetRemovePileUp(kFALSE);
  task->SetRejectKinkMother(rejectkinkmother);
  task->GetPIDQAManager()->SetHighResolutionHistos();
  if(useMC) task->SetHasMCData(kTRUE); // necessary for AOD
  printf("AOD filter %d On/OFF?\n",aodfilter);
  if(aodfilter > 0) {
    printf("ON AOD filter %d\n",aodfilter);
    task->SetUseFilterAOD(kTRUE);
    task->SetFilter(aodfilter);
  }
  else {
    task->SetUseFilterAOD(kFALSE);
  }
 
  
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
  if(usetof){
    pid->AddDetector("TOF", 0);
    pid->AddDetector("TPC", 1);
    //pid->ConfigureTPCrejection();
  }
  else {
    pid->AddDetector("TPC", 0);
  }

  if(usetof) pid->ConfigureTOF(TOFs);

  if(!useMC){
    Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    if(tpcdEdxcutlow) memcpy(paramsTPCdEdxcutlow,tpcdEdxcutlow,sizeof(paramsTPCdEdxcutlow));
    
    Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
    if(tpcdEdxcuthigh) memcpy(paramsTPCdEdxcuthigh,tpcdEdxcuthigh,sizeof(paramsTPCdEdxcuthigh));

    char *cutmodel;
    cutmodel="pol0";
    
    for(Int_t a=0;a<11;a++)
      {
	//   cout << a << " " << paramsTPCdEdxcut[a] << endl;
	Double_t tpcparamlow[1]={paramsTPCdEdxcutlow[a]};
	Float_t tpcparamhigh=paramsTPCdEdxcuthigh[a];
	pid->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlow,tpcparamhigh);
      }
    
    if(etacor)
      { 
	// Apply eta correction
	AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
	TF1 *etacorrection = GetEtaCorrection(listname);
	if(etacorrection) tpcpid->SetEtaCorrection(etacorrection);
      }
  }
  
  // QA
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  //task->SetFillSignalOnly(kFALSE);    // for DE pluging for MC
  if(kNoPhotonic) task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectron);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
  //if(useMC && addflag==1) task->SetDebugStreaming();

  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);

  task->SetFillNoCuts(kTRUE);

  if(isaod) {
    task->SetAODAnalysis();
    task->SetApplyCutAOD(kTRUE);
  }
 
  printf("*************************************\n");
  printf("Configuring task PbPb \n"); 
  //if(isLHC10) printf("Configuring TPC1 Task 2010 :\n");
  //if(isLHC11) printf("Configuring TPC1 Task 2011 :\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
