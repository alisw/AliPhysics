Bool_t ReadContaminationFunctions(TString filename, TF1 **functions, double sigma){
    TFile *in = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
  gROOT->cd();
  int isig = static_cast<int>(sigma * 100.);
  printf("Getting hadron background for the sigma cut: %d\n", isig);
  bool status = kTRUE;
  for(int icent = 0; icent < 12; icent++){
    functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_%d_%d", isig, icent)));
    if(functions[icent])
    {
//	new TCanvas;
//        functions[0]->Draw();
	printf("Config for centrality class %d found\n", icent);
    }
    else{
      printf("Config for the centrality class %d not found\n", icent);
      status = kFALSE;
    }
  }
  delete in;
  return status;
}

AliAnalysisTaskHFE* ConfigHFEpPbTRD(Bool_t useMC, Bool_t isAOD, TString appendix,
				 UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
				 UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000.,
				 Double_t* tpcdEdxcutlow=NULL,Double_t* tpcdEdxcuthigh=NULL,
				 Double_t TOFs=3., Int_t TOFmis=0,
				 Int_t itshitpixel = 0, Int_t icent=1,
				 Double_t etami=-0.8, Double_t etama=0.8,
				 Int_t TRDtrigger=1, Int_t trdpidmethod=1, Int_t TRDtl=6, Int_t TRDeff=4,
				 TString detector){
  
    Bool_t kAnalyseTaggedTracks = kTRUE;

    // TRD settings
    Float_t eeff[6] = {0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
    Int_t eeffint[6] = {70, 75, 80, 85, 90, 95};
    //  if(TRDeff >= 6 || TRDtl < 4 || TRDtl > 6) return NULL;
    if(TRDeff >= 6 || TRDtl > 6) return NULL;
    printf("TRD settings: %i %f \n",TRDtl, eeff[TRDeff]);
  
  //***************************************//
  //        Setting up the HFE cuts        //
  //***************************************//
  
  AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts for pPb with TRD");

 // hfecuts->SetQAOn();
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  hfecuts->SetCutITSpixel(itshitpixel);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
//  hfecuts->SetEtaRange(etami,etama);
  if(TRDtl>0) hfecuts->SetMinNTrackletsTRD(TRDtl, kTRUE);   // number of trd tracklets
  if(isAOD) hfecuts->SetAODFilterBit(4);
  
  //if((iPixelAny==AliHFEextraCuts::kAny) || (iPixelAny==AliHFEextraCuts::kSecond))     
  //hfecuts->SetProductionVertex(0,7,0,7);
 
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);
  // New pPb cuts (February 2013)
  hfecuts->SetUseCorrelationVertex();
  hfecuts->SetSPDVtxResolutionCut();
  hfecuts->SetpApileupCut();

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

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEtask%s",appendix.Data()));
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

  Bool_t activateTRDTrigger=kFALSE;
  if(TRDtrigger>0) activateTRDTrigger=kTRUE;
  task->SetTRDTrigger(activateTRDTrigger,TRDtrigger);

  //***************************************//
  //          Variable manager             //
  //***************************************//
  // Define Variables
  Double_t ptbinning[41] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30.};
  //Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  //Double_t ptbinning[53] = {0., 0.1, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};

  //Double_t etabinning[9] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};
  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  //Double_t etabinning[33] = {-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05, 
  //0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8};

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
  pid->SetDetectorsForAnalysis(detector);
  
  // Configure TPC PID
  // do the identical thing in data and MC
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
  pid->ConfigureTPCdefaultCut(cutmodel,paramsTPCdEdxcutlow,paramsTPCdEdxcuthigh[0]);



  // Apply eta correction
  AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
  // Setting to cut on absolute dEdx (not sigma)
  //tpcpid->UsedEdx();
  //TF1 *etacorrection = GetEtaCorrection();
  //if(etacorrection) tpcpid->SetEtaCorrection(etacorrection);

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

  if(TRDtl>0){
      AliHFEpidTRD *trdpid = pid->GetDetPID(AliHFEpid::kTRDpid);
      if(trdpidmethod==2) trdpid->SetTRD2DPID();
      trdpid->SetElectronEfficiency(eeff[TRDeff]);   // efficiency
      trdpid->SetNTracklets(TRDtl);      // ntracklets threshold
      trdpid->SetCutNTracklets(TRDtl, kTRUE);
  }


  // Load hadron background
  if(!useMC){
    Bool_t status = kTRUE;
    TF1 *hBackground[12];
    //status = ReadContaminationFunctions("hadronContamination_pPbTPCTOF_forwardEta.root", hBackground, tpcdEdxcutlow[0]);
    //status = ReadContaminationFunctions("hadronContamination_pPbTPCTOF_eta66.root", hBackground, tpcdEdxcutlow[0]);
    if(TRDtl==5){
	status = ReadContaminationFunctions("hadroncontamination_TRDTPC_pPb_eta06_5tracklets.root", hBackground, tpcdEdxcutlow[0]);
        printf("hadroncont5tracklets \n");
    }
    if(TRDtl==6)
    {
	status = ReadContaminationFunctions("hadroncontamination_TRDTPC_pPb_eta06_6tracklets.root", hBackground, tpcdEdxcutlow[0]);
	printf("hadroncont6tracklets \n");
    }
    for(Int_t a=0;a<12;a++) {
      //printf("back %f \n",hBackground[a]);
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
    if(usetof) v0trackCuts->SetTOFPIDStep(kTRUE);
    if(TRDtl>0)v0trackCuts->SetMinNTrackletsTRD(TRDtl,kTRUE); // condition for TRD tracklets
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
