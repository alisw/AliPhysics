Bool_t ReadContaminationFunctionsBeauty(TString filename, TF1 **functions, double sigma){
  TFile *in = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
  gROOT->cd();
  int isig = static_cast<int>(sigma * 100.);
  if (isig == -44) isig = -42;
  if (isig == 6) isig = 9;
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

Bool_t ReadContaminationFunctions(TString filename, TF1 **functions, double sigma){
  TFile *in = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/%s", filename.Data()));
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


AliAnalysisTaskHFE* ConfigHFEnpepPb(Bool_t useMC, Bool_t isAOD, Bool_t isBeauty, TString appendix,
                UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
                UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000., 
                Double_t* tpcdEdxcutlow=NULL, Double_t* tpcdEdxcuthigh=NULL, 
                Double_t TOFs=3., Int_t TOFmis=0, 
                Int_t itshitpixel = 0, 
		Int_t spdcheck = 0,
		Int_t icent, 
		Double_t etami=-0.8, Double_t etama=0.8,
		Double_t phimi=-1., Double_t phima=-1.,
		Double_t assETAm=-0.8, Double_t assETAp=0.8,
		Double_t assMinPt=0.2, Int_t assITS=2, 
                Int_t assTPCcl=100, Int_t assTPCPIDcl=80, 
                Double_t assDCAr=1.0, Double_t assDCAz=2.0, 
                Double_t *assTPCSminus=NULL, Double_t *assTPCSplus=NULL, 
                Double_t assITSpid=-3., 
		Double_t assTOFs=3.,
                Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE, Int_t weightlevelback = -1, 
		Double_t etadalwei = 0.,
	        Bool_t nonPhotonicElectronBeauty = kFALSE, Bool_t ipCharge = kFALSE, Bool_t ipOpp = kFALSE, Double_t *ipPar=NULL)
{
  Bool_t kAnalyseTaggedTracks = kTRUE;
  Bool_t kApplyPreselection = kFALSE;

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
  //if (spdcheck>0){
  //hfecuts->SetCheckITSLayerStatus(kTRUE);
  //printf("\n\nWe check the live status of the pixels\n\n");
  //}
  hfecuts->SetEtaRange(etami,etama);
  if(phimi >= 0. && phima >= 0) hfecuts->SetPhiRange(phimi,phima);
  hfecuts->SetRejectKinkDaughters();
  hfecuts->SetAcceptKinkMothers();
  if(isAOD) hfecuts->SetAODFilterBit(4);
  
  //if((iPixelAny==AliHFEextraCuts::kAny) || (iPixelAny==AliHFEextraCuts::kSecond))     
 
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);
  // New pPb cuts (February 2013)
  hfecuts->SetUseCorrelationVertex();
  hfecuts->SetSPDVtxResolutionCut();
  hfecuts->SetpApileupCut();

  Bool_t ipSig = kFALSE;

  if(!ipPar) hfecuts->SetIPcutParam(0.0054,0.078,-0.56,0,ipSig,ipCharge,ipOpp); // reference
  else hfecuts->SetIPcutParam(ipPar[0],ipPar[1],ipPar[2],0,ipSig,ipCharge,ipOpp);
  /*if(ipsys==1) hfecuts->SetIPcutParam(0.0054,0.057,-0.66,0,ipSig,ipCharge,ipOpp);
  if(ipsys==2) hfecuts->SetIPcutParam(0.012,0.088,-0.65,0,ipSig,ipCharge,ipOpp);*/


  if(isBeauty) hfecuts->SetProductionVertex(0,100,0,100);

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
  task->SetRejectKinkMother(kFALSE);

  // Determine the centrality estimator
  task->SetCentralityEstimator("V0A");
  if (icent == 2) task->SetCentralityEstimator("V0M");
  else if (icent == 3) task->SetCentralityEstimator("CL1");
  else if (icent == 4) task->SetCentralityEstimator("ZNA");

  //***************************************//
  //        Prepare preselection           //
  // This mimics the ESD->AOD filter in    //
  // case of the ESD analysis and selects  //
  // only tracks which will be selected in //
  // the AOD analysis with the given filter//
  // bit. Not to be applied for AODS.      //
  // For pPb the cuts used are (bit 4)     //
  // esdTrackCutsHG0 from file $ALICE_ROOT///
  // ANALYSIS/macros/AddTaskESDFilter.C    //
  //***************************************//

  if(kApplyPreselection){    
    AliESDtrackCuts* esdfilter = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    esdfilter->SetMaxDCAToVertexXY(2.4);
    esdfilter->SetMaxDCAToVertexZ(3.2);
    esdfilter->SetDCAToVertex2D(kTRUE);

    task->SetHFECutsPreselect(esdfilter);
    printf("Put a preselection cut\n");
    task->SetFillNoCuts(kTRUE);
  }

  //***************************************//
  //          Variable manager             //
  //***************************************//
  // Define Variables
  Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
  Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;

  AliHFEvarManager *vm = task->GetVarManager();
  vm->AddVariable("pt", sizept, ptbinning);
  vm->AddVariable("eta", sizeeta, -0.8,0.8);
  //vm->AddVariable("phi",18, -0, 2*TMath::Pi());
  vm->AddVariable("phi",3, -0, 2*TMath::Pi());
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
  // do the identical thing in data and MC
  Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  if(tpcdEdxcutlow) memcpy(paramsTPCdEdxcutlow,tpcdEdxcutlow,sizeof(paramsTPCdEdxcutlow));

  Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  if(tpcdEdxcuthigh) memcpy(paramsTPCdEdxcuthigh,tpcdEdxcuthigh,sizeof(paramsTPCdEdxcuthigh));

  char *cutmodel;

  if(useMC){ // constant (default) cut for MC
      cutmodel="pol0(0)";
      Double_t params[1];
      params[0]=paramsTPCdEdxcutlow[0];
      pid->ConfigureTPCdefaultCut(cutmodel, params,tpcdEdxcuthigh[0]);
  } else { // correct for mean shift in data
      cutmodel="min(pol1(0),pol0(2))";
      Double_t params[3];
      //params[0]=-0.12; params[1]=0.14; params[2]=0.09;
      params[0]=-0.21 + paramsTPCdEdxcutlow[0];
      params[1]=0.14;
      params[2]=paramsTPCdEdxcutlow[0];
      pid->ConfigureTPCdefaultCut(cutmodel, params,tpcdEdxcuthigh[0]);
  }
  /*
  char *cutmodel;
  cutmodel="pol0";

  for(Int_t a=0;a<11;a++){
    // Not necessary anymore, since the pPb case is handled similarly to the pp case
    //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
    Double_t tpcparamlow[1]={paramsTPCdEdxcutlow[a]};
    Float_t tpcparamhigh=paramsTPCdEdxcuthigh[a];
    pid->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlow,tpcparamhigh);
  }
  pid->ConfigureTPCdefaultCut(cutmodel,paramsTPCdEdxcutlow,paramsTPCdEdxcuthigh[0]); // After introducing the pPb flag, pPb is merged with pp and this line defines the cut
  */

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
    if(isAOD==1) {
	if (usetof)  status = ReadContaminationFunctions("hadroncontamination_AOD139_TOFPID_pPb_eta06.root", hBackground, tpcdEdxcutlow[0]);
	else { 
	  if (spdcheck == 0){
	    status = ReadContaminationFunctions("hadroncontamination_AOD139_noTOFPID_pPb_eta06.root", hBackground, tpcdEdxcutlow[0]);
	  } else if (spdcheck == 1){ 
	    status = ReadContaminationFunctions("hadroncontamination_noTOFPID_pPb_AOD_eta06_TPCcut0_envelope_minsys.root", hBackground, tpcdEdxcutlow[0]);	    
	  } else if (spdcheck == 2){ 
	    status = ReadContaminationFunctions("hadroncontamination_noTOFPID_pPb_AOD_eta06_TPCcut0_envelope_maxsys.root", hBackground, tpcdEdxcutlow[0]);	    
	  } else if (spdcheck == 3){ 
	    status = ReadContaminationFunctions("hadroncontamination_AOD139_noTOFPID_pPb_eta06_polyFit.root", hBackground, tpcdEdxcutlow[0]);	    
	  }
	}
    }
    else if (isBeauty==1) {
	if (spdcheck == 0){
	    printf("hadron cont standard beauty");

	    if(TOFs==0){
		printf("hadron cont no TOF\n");
		status = ReadContaminationFunctions("hadroncontamination_noTOF_pPb_Beauty_ESD_eta06.root", hBackground, tpcdEdxcutlow[0]);
	    }
	    else{
                printf("hadron cont with TOF\n");
		status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2.root", hBackground, tpcdEdxcutlow[0]);
	    }
	} else if (spdcheck == 1){
	    printf("hadron cont min beauty");
            if(TOFs==0){
		printf("hadron cont no TOF\n");
		status = ReadContaminationFunctions("hadroncontamination_noTOF_pPb_Beauty_ESD_eta06_envelopemin.root", hBackground, tpcdEdxcutlow[0]);
	    }
            status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2_envelope_minsys.root", hBackground, tpcdEdxcutlow[0]);
	} else if (spdcheck == 2){
	    printf("hadron cont max beauty");
	    if(TOFs==0){
		printf("hadron cont no TOF\n");
		status = ReadContaminationFunctions("hadroncontamination_noTOF_pPb_Beauty_ESD_eta06_envelopemax.root", hBackground, tpcdEdxcutlow[0]);
	    }
            status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2_envelope_maxsys.root", hBackground, tpcdEdxcutlow[0]);
	} else if (spdcheck == 3){
	    printf("hadron cont pol3 beauty");
            status = ReadContaminationFunctionsBeauty("hadroncontamination_ESD_Beauty_TOFPID_pPb_eta06_iter2_pol3.root", hBackground, tpcdEdxcutlow[0]);
	}
    }
  else  status = ReadContaminationFunctions("hadroncontamination_TOFTPC_pPb_eta06_newsplines_try3.root", hBackground, tpcdEdxcutlow[0]);
    for(Int_t a=0;a<12;a++) {
      //printf("back %f \n",hBackground[a]);
      if(status) task->SetBackGroundFactorsFunction(hBackground[a],a);
      else printf("not all background functions found\n");
    }
  }

  //***************************************//
  //       Configure NPE plugin            //
  //***************************************//

  AliHFENonPhotonicElectron *backe = new AliHFENonPhotonicElectron(Form("HFEBackGroundSubtractionPID2%s",appendix.Data()),"Background subtraction");  //appendix
    //Setting the Cuts for the Associated electron-pool
  AliHFEcuts *hfeBackgroundCuts = new AliHFEcuts(Form("HFEBackSub%s",appendix.Data()),"Background sub Cuts");
  //  hfeBackgroundCuts->SetEtaRange(assETA);
  hfeBackgroundCuts->SetEtaRange(assETAm,assETAp);
  hfeBackgroundCuts->SetPtRange(assMinPt,20.);

  hfeBackgroundCuts->SetMaxChi2perClusterTPC(4);
  hfeBackgroundCuts->SetMinNClustersITS(assITS);
  hfeBackgroundCuts->SetMinNClustersTPC(assTPCcl);
  hfeBackgroundCuts->SetMinNClustersTPCPID(assTPCPIDcl);
  hfeBackgroundCuts->SetMaxImpactParam(assDCAr,assDCAz);
  if(isAOD) hfeBackgroundCuts->SetAODFilterBit(0); // Standard TPC only tracks
  hfeBackgroundCuts->SetQAOn();			        // QA

  AliHFEpid *pidbackground = backe->GetPIDBackground();
  if(useMC) pidbackground->SetHasMCData(kTRUE);

  if (etadalwei>0){
    printf("\n\n\n\n\n WEIGHTS FOR ETA DALITZ !!!!!!!!!!!!!!!!!! \n\n\n\n");
    backe->SetEtaDalitzWeightFactor(etadalwei);  
  }   

  if (assTOFs>0.){
    pidbackground->AddDetector("TOF", 0);
    pidbackground->AddDetector("TPC", 1);
  } else {
    pidbackground->AddDetector("TPC", 0);
  }

  Double_t paramsTPCdEdxcutlowAssoc[12] ={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
  if(assTPCSminus) memcpy(paramsTPCdEdxcutlowAssoc,assTPCSminus,sizeof(paramsTPCdEdxcutlowAssoc));

  Double_t paramsTPCdEdxcuthighAssoc[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  if(assTPCSplus) memcpy(paramsTPCdEdxcuthighAssoc,assTPCSplus,sizeof(paramsTPCdEdxcuthighAssoc));
    
  char *cutmodelAssoc;
  cutmodelAssoc="pol0";
  for(Int_t a=0;a<11;a++){
    // Not necessary anymore, since the pPb case is handled similarly to the pp case
    //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
    Double_t tpcparamlow[1]={paramsTPCdEdxcutlowAssoc[a]};
    Float_t tpcparamhigh=paramsTPCdEdxcuthighAssoc[a];
    pidbackground->ConfigureTPCcentralityCut(a,cutmodelAssoc,tpcparamlow,tpcparamhigh);
  }
  pidbackground->ConfigureTPCdefaultCut(cutmodelAssoc,paramsTPCdEdxcutlowAssoc,paramsTPCdEdxcuthighAssoc[0]); // After introducing the pPb flag, pPb is merged with pp and this line defines the cut
  //backe->GetPIDBackgroundQAManager()->SetHighResolutionHistos();

  if (assTOFs>0.){
    pidbackground->ConfigureTOF(TOFs);
  }

  backe->SetHFEBackgroundCuts(hfeBackgroundCuts);

  // Selection of associated tracks for the pool
  if(useCat1Tracks) backe->SelectCategory1Tracks(kTRUE);
  if(useCat2Tracks){
    backe->SelectCategory2Tracks(kTRUE);
    backe->SetITSMeanShift(-0.5);
    backe->SetITSnSigmaHigh(assITSpid);
    Double_t assITSminus = -1.0 * assITSpid;
    backe->SetITSnSigmaLow(assITSminus);
    //backe->SetminPt(assMinPt);
  }

  // apply opening angle cut to reduce file size
  backe->SetMaxInvMass(0.3);
  backe->SetPtBinning(sizept, ptbinning);
  backe->SetEtaBinning(sizeeta, etabinning);
  //backe->SetAnaPairGen(kTRUE,2);
  //backe->SetDisplayMCStack();
  // MC weight
  if(useMC) {
    //printf("test put weight %d\n",weightlevelback);
    if((weightlevelback >=0) && (weightlevelback < 3)) backe->SetWithWeights(weightlevelback);
  }
  task->SetHFEBackgroundSubtraction(backe);

  //task->SetWeightHist(); 
  //if(useMC) task->SetDebugStreaming();
  //task->SetCalcContamBeauty(kTRUE);

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
  if(isBeauty){
  if(nonPhotonicElectronBeauty) task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectronBeauty);
  }
  else task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectron);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n");
  return task;
}
