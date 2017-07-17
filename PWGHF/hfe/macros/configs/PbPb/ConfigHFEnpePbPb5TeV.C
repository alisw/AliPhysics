TF1* GetEtaCorrection(){
  TString list=gSystem->Getenv("LIST");

  TString etaMap="$TRAIN_ROOT/hfe_HFE/EtaCorrMapsTest.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");

  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigPbPb2010_Cent","Eta map not found: %s",etaMap.Data());
    return 0;
  }

  TFile f(etaMap.Data());
  if (!f.IsOpen()) return 0;
  gROOT->cd();
  TList *keys=f.GetListOfKeys();

  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(list)){
      printf("Using Eta Correction Function: %s\n",kName.Data());
      return (TF1*)f.Get(kName.Data());
    }
  }
  return 0;
}

  // ***** Background selection for PbPb 5TeV *****                                // TOF sigma and ITS sigma added
Bool_t ReadContaminationFunctions(TString filename, TF1 **functions, double sigma, double TOFs, double ITSs){
  //TFile *in = TFile::Open(Form("$TRAIN_ROOT/util/hfe/%s", filename.Data()));   // GSI version 
  TFile *in = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/%s", filename.Data()));   // GRID version 
  gROOT->cd();
  //int isig = static_cast<int>(sigma * 100.);  // original
  int isig      = static_cast<int>(sigma * 1000.);   
  int nTOFsigma = static_cast<int>(TOFs*10);
  int nITSsigma = static_cast<int>(ITSs*10);

  printf("Getting hadron background for the sigma cut: %d\n", isig);
  printf("Getting hadron background for TOF sigma (INTEGER*10): %d\n", nTOFsigma);
  printf("Getting hadron background for ITS sigma (INTEGER*10): %d\n", nITSsigma);
  bool status = kTRUE;

  for(int icent = 0; icent < 12; icent++){
    //functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_%d_%d", isig, icent)));        // original
    if(isig<0)  // --- case of negative low TPC cut ---
    {
        int isigSignSwitched = 0-isig;  // sign switched 
        //cout << " *** isig<0 *** " << endl;
        functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_ITS%d_TOF%d_m%d_%d", nITSsigma, nTOFsigma, isigSignSwitched, icent)));
        //printf("function[%d] name = hback_ITS%d_TOF%d_m%d_%d\n",icent, nITSsigma, nTOFsigma, isigSignSwitched, icent);
    }
    else       functions[icent] = dynamic_cast<TF1 *>(in->Get(Form("hback_ITS%d_TOF%d_%d_%d", nITSsigma, nTOFsigma, isig, icent))); 
    if(functions[icent]) printf("Config for centrality class %d found - function name: %s\n", icent, functions[icent]->GetName());
    else{
      printf("Config for the centrality class %d not found\n", icent);
      status = kFALSE;
    }

  }
  delete in;
  return status;
}
  // **********************************************

AliAnalysisTaskHFE* ConfigHFEnpePbPb5TeV(Bool_t useMC, Bool_t isAOD, TString appendix,
				     UChar_t TPCcl=70, UChar_t TPCclPID = 80, 
				     UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000., 
				     Double_t* tpcdEdxcutlow=NULL, Double_t* tpcdEdxcuthigh=NULL, 
				     Double_t TOFs=3., Int_t TOFmis=0, 
				     Double_t ITSs=0.,
				     Int_t itshitpixel = 0, Double_t itsChi2PerClusters, Double_t tpcClShared,
				     Bool_t etacor = kFALSE, Bool_t multicor = kFALSE, Bool_t toflast = kFALSE,
				     Double_t etami=-0.8, Double_t etama=0.8,
				     Double_t assETAm=-0.8, Double_t assETAp=0.8,
				     Int_t assITS=2, 
				     Int_t assTPCcl=100, Int_t assTPCPIDcl=80, 
				     Double_t assDCAr=1.0, Double_t assDCAz=2.0, 
				     Double_t *assTPCSminus=NULL, Double_t *assTPCSplus=NULL, 
				     Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE, 
                                     Int_t weightlevelback = -1, 
                                     Double_t assMinpT = 0.1,  // associated particle minimum pT syst. (mfaggin, 14th July 2017)
                                     Bool_t releasemcvx = kFALSE,
				     Bool_t nondefaultcentr = kFALSE,Bool_t ipCharge = kFALSE, Bool_t ipOpp = kFALSE,
				     Bool_t usekfparticle = kFALSE
                                     // ----- Asymmetric ITS cut (mfaggin, June 26th 2017) -----

                                     // --------------------------------------------------------
                                     )
{
  Bool_t kAnalyseTaggedTracks = kFALSE;
  Bool_t kApplyPreselection = kFALSE;

  Bool_t isBeauty = kFALSE;

  //***************************************//
  //        Setting up the HFE cuts        //
  //***************************************//

  AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts for PbPb");
  //hfecuts->SetQAOn();
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(TPCcl);
  hfecuts->SetMinNClustersTPCPID(TPCclPID);
  hfecuts->SetMinNClustersITS(ITScl);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFoundAll, AliHFEextraCuts::kFoundAllOverFindable);
  hfecuts->SetCutITSpixel(itshitpixel);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMaxChi2perClusterITS(itsChi2PerClusters);
  hfecuts->SetEtaRange(etami,etama);
  hfecuts->SetFractionOfSharedTPCClusters(tpcClShared);
  hfecuts->SetAcceptKinkMothers();
  if(isAOD) hfecuts->SetAODFilterBit(2);  

  if((itshitpixel==AliHFEextraCuts::kAny) || (itshitpixel==AliHFEextraCuts::kSecond))     
  hfecuts->SetProductionVertex(0,7,0,7);
 
  hfecuts->SetMaxImpactParam(DCAxy,DCAz);
  hfecuts->SetUseMixedVertex(kTRUE);
  hfecuts->SetVertexRange(10.);

  Bool_t ipSig = kFALSE;

  hfecuts->SetIPcutParam(0.0054,0.078,-0.56,0,ipSig,ipCharge,ipOpp);
  if(isBeauty || releasemcvx) hfecuts->SetProductionVertex(0,100,0,100);

  // TOF settings:
  Int_t usetof=0;
  Bool_t kTOFmis=kFALSE;
  if (TOFs>0.){
    usetof = 1;
    printf("CONFIGURATION FILE: TOF is used \n");
    hfecuts->SetTOFPIDStep(kTRUE);
    if(useMC && (!isAOD)) hfecuts->SetMatchTOFLabel(kTRUE);
    printf("CONFIGURATION FILE: TOF PID step is requested !!!! \n");
    if (TOFmis>0){
      kTOFmis = kTRUE;
      printf("CONFIGURATION FILE: TOF mismatch rejection is set ON \n");
    }
  }

  
  // ITS settings:
  Int_t useits=0;
  if (ITSs>0.){
    useits = 1;
    printf("CONFIGURATION FILE: ITS is used \n");
  }

  //***************************************//
  //        Setting up the task            //
  //***************************************//

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEtask%s",appendix.Data()));
  printf("task %p\n", task);
  task->SetPbPbAnalysis();
  task->SetRemovePileUp(kFALSE);
  task->SetHFECuts(hfecuts);
  task->GetPIDQAManager()->SetHighResolutionHistos();
  task->SetRejectKinkMother(kFALSE);
  
  // Determine the centrality estimator
  task->SetCentralityEstimator("V0M");

  // Get weights
  task->SetWeightHist();

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
    AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdTrackCutsH->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    task->SetHFECutsPreselect(esdTrackCutsH);
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

  if (usetof && (useits<1)){
    if(toflast){
      pid->AddDetector("TPC", 0);
      pid->AddDetector("TOF", 1);
    } else {
      pid->AddDetector("TOF", 0);
      pid->AddDetector("TPC", 1);
    }
  } else if(usetof && (useits>0)){
    if(toflast){
      pid->AddDetector("ITS", 0);
      pid->AddDetector("TPC", 1);
      pid->AddDetector("TOF", 2);
     
    } else {
      pid->AddDetector("TOF", 0);
      pid->AddDetector("ITS", 1);
      pid->AddDetector("TPC", 2);
    }
  }
  else {
    pid->AddDetector("TPC", 0);
  }
  
  // Configure TPC PID
  // do the identical thing in data and MC

  Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  if(tpcdEdxcutlow) memcpy(paramsTPCdEdxcutlow,tpcdEdxcutlow,sizeof(paramsTPCdEdxcutlow));

  Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  if(tpcdEdxcuthigh) memcpy(paramsTPCdEdxcuthigh,tpcdEdxcuthigh,sizeof(paramsTPCdEdxcuthigh));

  char *cutmodel;
  cutmodel="pol0";

  for(Int_t a=0;a<11;a++){
    // Not necessary anymore, since the PbPb case is handled similarly to the pp case
    //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
    Double_t tpcparamlow[1]={paramsTPCdEdxcutlow[a]};
    Float_t tpcparamhigh=paramsTPCdEdxcuthigh[a];
    pid->ConfigureTPCcentralityCut(a,cutmodel,tpcparamlow,tpcparamhigh);
  }

  if(!useMC){
    AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTOFpid);
    if(etacor){ 
	  // Apply eta correction
	  TF1 *etacorrection = GetEtaCorrection();
	  if(etacorrection) tpcpid->SetEtaCorrection(etacorrection);
    }
    if(multicor){
	  TF1 *centralityCorrection = new TF1("centralityCorrection", "pol1", 0., 10000.);
	  centralityCorrection->SetParameter(0, 1.0);
	  centralityCorrection->SetParameter(1, -0.00002);
	  tpcpid->SetCentralityCorrection(centralityCorrection);
    }
  }

  // Configure TOF PID
  if (usetof){
    pid->ConfigureTOF(TOFs);
    AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
    if (kTOFmis){
      tofpid->SetRejectTOFmismatch();
    }
    if(toflast) tofpid->SetGenerateTOFmismatch(); // Makes only sense if TOF is the last detector
  }

  // Configure ITS PID
  if (useits>0){
    AliHFEpidITS *itspid = pid->GetDetPID(AliHFEpid::kITSpid);
    //itspid->SetITSnSigma(1.);
    itspid->SetITSnSigma(ITSs); // ***** modified 11/06/2017 (mfaggin)
  }

  // To make different upper TOF cut to see contamination effect
  // The below two lines should be removed after this check
  //AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
  //if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs); //only to check the assymmetric tof cut


  // ***** Background selection for PbPb 5TeV *****
  // Load hadron background
  if(!useMC){
    Bool_t status = kTRUE;
    TF1 *hBackground[12];                                                                                // TOF sigma and ITS sigma added
    status = ReadContaminationFunctions("hadronContamination_PbPb5TeV.root", hBackground, tpcdEdxcutlow[0], TOFs, ITSs);
    for(Int_t a=0;a<12;a++) {
      //printf("back %f \n",hBackground[a]);
      if(status) task->SetBackGroundFactorsFunction(hBackground[a],a);
      else printf("not all background functions found\n");
    }
  }
  // **********************************************


  // Load user defined centrality bin for beauty analysis
  if(nondefaultcentr) {
      Float_t arraycentr[12] = {0.,20., 25., 30., 35., 40., 80.,85., 90., 95.,100.,105.};
      for(Int_t i=0;i<12;i++)
      {
          task->SetPbPbUserCentralityLimit(kTRUE);
          task->SetPbPbUserCentralityArray(i,arraycentr[i]);

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

  //hfeBackgroundCuts->SetPtRange(0.1,1e10);    // old hardcoded minimum pT (mfaggin, 14th July 2017)
  hfeBackgroundCuts->SetPtRange(assMinpT,1e10); // associated particle minimum pT syst. (mfaggin, 14th July 2017)

  hfeBackgroundCuts->SetMaxChi2perClusterTPC(4);
  hfeBackgroundCuts->SetMinNClustersITS(assITS);
  hfeBackgroundCuts->SetMinNClustersTPC(assTPCcl);
  hfeBackgroundCuts->SetMinNClustersTPCPID(assTPCPIDcl);
  hfeBackgroundCuts->SetMaxImpactParam(assDCAr,assDCAz);
  if(isAOD) hfeBackgroundCuts->SetAODFilterBit(0);
  if(usekfparticle) backe->SetAlgorithmMA(kTRUE);
  //hfeBackgroundCuts->SetQAOn();			        // QA

  AliHFEpid *pidbackground = backe->GetPIDBackground();
  if(useMC) pidbackground->SetHasMCData(kTRUE);
  pidbackground->AddDetector("TPC", 0);
  Double_t paramsTPCdEdxcutlowAssoc[12] ={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
  if(assTPCSminus) memcpy(paramsTPCdEdxcutlowAssoc,assTPCSminus,sizeof(paramsTPCdEdxcutlowAssoc));

  Double_t paramsTPCdEdxcuthighAssoc[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  if(assTPCSplus) memcpy(paramsTPCdEdxcuthighAssoc,assTPCSplus,sizeof(paramsTPCdEdxcuthighAssoc));
    
  char *cutmodelAssoc;
  cutmodelAssoc="pol0";
  for(Int_t a=0;a<11;a++){
    //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
    Double_t tpcparamlow[1]={paramsTPCdEdxcutlowAssoc[a]};
    Float_t tpcparamhigh=paramsTPCdEdxcuthighAssoc[a];
    pidbackground->ConfigureTPCcentralityCut(a,cutmodelAssoc,tpcparamlow,tpcparamhigh);
  }
  //backe->GetPIDBackgroundQAManager()->SetHighResolutionHistos();
  backe->SetHFEBackgroundCuts(hfeBackgroundCuts);

  // Selection of associated tracks for the pool
  if(useCat1Tracks) backe->SelectCategory1Tracks(kTRUE);
  if(useCat2Tracks){
    backe->SelectCategory2Tracks(kTRUE);
    backe->SetITSMeanShift(-0.5);
  }

  // apply opening angle cut to reduce file size
  backe->SetMaxInvMass(0.3);
  backe->SetStudyRadius(kTRUE);
  backe->SetPtBinning(sizept, ptbinning);
  backe->SetEtaBinning(sizeeta, etabinning);
  // MC weight
  if(useMC) {
    //printf("test put weight %d\n",weightlevelback);
    if((weightlevelback >=0) && (weightlevelback < 3)) backe->SetWithWeights(weightlevelback);
  }

  task->SetHFEBackgroundSubtraction(backe);

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
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectron);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->PrintStatus();
  pid->PrintStatus();
  printf("*************************************\n");
  return task;
}
