Bool_t ReadContaminationFunctions(TString filename, TF1 **functions, double sigma){
  TFile *in = TFile::Open(Form("$ALICE_ROOT/PWGHF/hfe/configs/pPb/%s", filename.Data()));
  gROOT->cd();
  int isig = static_cast<int>(sigma * 100.);
  printf("Getting hadron background for the sigma cut: %d\n", isig);
  bool status = kTRUE;
  for(int icent = 0; icent < 12; icent++){
    functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_%d_%d", isig, icent)));
    if(functions[icent]) printf("Config for centrality class %d found\n", icent);
    else{
      printf("Config for the centrality class %d not found\n", icent);
      status = kFALSE;
    }
  }
  delete in;
  return status;
}

AliAnalysisTaskHFE* ConfigHFEpPb(Bool_t useMC, Bool_t isAOD, TString appendix, 
				 UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
				 UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000.,
				 Double_t* tpcdEdxcutlow=NULL,Double_t* tpcdEdxcuthigh=NULL,
				 Double_t TOFs=3., Int_t TOFmis=0,
				 Int_t itshitpixel = 0, Int_t icent = 1,
         Double_t etami=-0.8, Double_t etama=0.8){
  
  Bool_t kAnalyseTaggedTracks = isAOD ? kFALSE : kTRUE;
  
  //***************************************//
  //        Setting up the HFE cuts        //
  //***************************************//
  
  AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts for pPb");

  //hfecuts->SetQAOn();
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(itshitpixel);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetEtaRange(etami,etama);
  if(isAOD) hfecuts->SetAODFilterBit(4);
  
  //if((iPixelAny==AliHFEextraCuts::kAny) || (iPixelAny==AliHFEextraCuts::kSecond))     
  //hfecuts->SetProductionVertex(0,7,0,7);
 
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);
  // New pPb cuts (February 2013)
  hfecuts->SetUseCorrelationVertex();
  hfecuts->SetSPDVtxResolutionCut();

  // TOF settings:
  Int_t usetof=0;
  Bool_t kTOFmis=kFALSE;
  if (TOFs>0.){
    usetof = 1;
    printf("CONFIGURATION FILE: TOF is used \n");
    hfecuts->SetTOFPIDStep(kTRUE);
    printf("CONFIGURATION FILE: TOF PID step is requested !!!! \n");
    if (TOFmis>0){
      kTOFmis = kTRUE;
      printf("CONFIGURATION FILE: TOF mismatch rejection is set ON \n");
    }
  }

  //***************************************//
  //        Setting up the task            //
  //***************************************//

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(appendix);
  printf("task %p\n", task);
  task->SetpPbAnalysis();
  if(!isAOD) task->SetRemoveFirstEventInChunk();
  task->SetRemovePileUp(kFALSE);
  task->SetHFECuts(hfecuts);
  task->GetPIDQAManager()->SetHighResolutionHistos();

  // Determine the centrality estimator
  task->SetCentralityEstimator("V0A");
  if (icent == 2) task->SetCentralityEstimator("V0M");
  else if (icent == 3) task->SetCentralityEstimator("CL1");
  else if (icent == 4) task->SetCentralityEstimator("ZNA");

  //***************************************//
  //          Variable manager             //
  //***************************************//

  // Define Variables
  Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};

  //Double_t etabinning[9] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};
  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
  Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;

  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt", sizept, ptbinning);
  vm->AddVariable("eta", sizeeta, -0.8,0.8);
  vm->AddVariable("phi",21, -0, 2*TMath::Pi());
  vm->AddVariable("charge");
  vm->AddVariable("source");
  vm->AddVariable("centrality");

  // For the moment, remove the part dedicated to the background subtraction.
  // It will be implemented in a different way, reading it from a root file.

  //***************************************//
  //          Configure the PID            //
  //***************************************//

  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);

  if (usetof){
    pid->AddDetector("TOF", 0);
    pid->AddDetector("TPC", 1);
  } else {
    pid->AddDetector("TPC", 0);
  }
  
  // Configure TPC PID
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
  }

  // Configure TOF PID
  if (usetof){
    pid->ConfigureTOF(TOFs);
    AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
    if (kTOFmis){
      tofpid->SetRejectTOFmismatch();
    }
  }
  // To make different upper TOF cut to see contamination effect
  // The below two lines should be removed after this check
  //AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
  //if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs); //only to check the assymmetric tof cut

  // Load hadron background
  if(!useMC){
    Bool_t status = kTRUE;
    TF1 *hBackground[12];
    status = ReadContaminationFunctions("hadronContamination_pPbTPCTOF_forwardEta.root", hBackground, tpcdEdxcutlow[0]);
    for(Int_t a=0;a<12;a++) {
      //  printf("back %f \n",p0[a]);
      if(status) task->SetBackGroundFactorsFunction(hBackground[a],a);
      else printf("not all background functions found\n");
    }
  }

  //***************************************//
  //          V0 tagged tracks             //
  //***************************************//
  if(kAnalyseTaggedTracks){
    AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    v0trackCuts->CreateStandardCuts();
    v0trackCuts->SetMinNClustersTPC(TPCcl);  
    v0trackCuts->SetMinNClustersTPCPID(TPCclPID);
    v0trackCuts->SetMinRatioTPCclusters(0.6);
    v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    v0trackCuts->SetMinNClustersITS(1);
    v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
    v0trackCuts->SetCheckITSLayerStatus(kFALSE);
    v0trackCuts->UnsetVertexRequirement();
    //hfecuts->SetSigmaToVertex(10);
    v0trackCuts->SetTOFPIDStep(kTRUE);
    v0trackCuts->SetQAOn();

    task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
    task->SetTaggedTrackCuts(v0trackCuts);
    task->SetCleanTaggedTrack(kTRUE);
  }

  // QA
  printf("task %p\n", task);
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);    
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
