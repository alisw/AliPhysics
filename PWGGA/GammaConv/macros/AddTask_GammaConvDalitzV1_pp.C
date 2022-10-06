
void AddTask_GammaConvDalitzV1_pp(  Int_t trainConfig = 1,  //change different set of cuts
                                    Int_t isMC   = 0, //0 Exp data, 1 MC data, 2 JJ MC for improves.
                                    TString photonCutNumberV0Reader       = "",
                                    TString periodNameV0Reader            = "",
                                    Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                    Bool_t    enableElecdEdxPostCalibration = kFALSE,     //PostCalibration
                                    Bool_t    ForceRPostCalibration = kFALSE,   //Works only with Postcalibration.
                                    Bool_t    lightVersion = kFALSE,                      //lightVersion Systematic.
                                    TString   fileNameExternalInputs        = "",
                                    Int_t   enableMatBudWeightsPi0          =  0,              // 1 = three radial bins, 2 = 10 radial bins
                                    TString   additionalTrainConfig         = "0"
         ) {





  Int_t isHeavyIon = 0;
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode




  AliCutHandlerPCM cuts;
  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");

  TString addTaskName                 = "AddTask_GammaConvDalitzV1_pp";
  TString sAdditionalTrainConfig  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "","", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO:"<< addTaskName.Data()<< " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi()   << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF","",addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;
  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights,addTaskName);
  TString strTrackMatcherRunningMode         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM","",addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i",addTaskName.Data(), trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

   //========= Check whether PID Reponse is there ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    Error(Form("%s_%i",addTaskName.Data(), trainConfig), "No PID response has been initialized aborting.");
    return;
  }


  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }


  if( !(AliDalitzElectronSelector*)mgr->GetTask("ElectronSelector") ){
    AliDalitzElectronSelector *fElectronSelector = new AliDalitzElectronSelector("ElectronSelector");
    //ConfigV0ReaderV1(fV0ReaderV1,ConvCutnumber,IsHeavyIon);
    // Set AnalysisCut Number
    AliDalitzElectronCuts *fElecCuts=0;
    TString ElecCuts = "30105400000003300000";

    if( ElecCuts!=""){
      fElecCuts= new AliDalitzElectronCuts(ElecCuts.Data(),ElecCuts.Data());
      if(fElecCuts->InitializeCutsFromCutString(ElecCuts.Data())){
        fElectronSelector->SetDalitzElectronCuts(fElecCuts);
        fElecCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    fElectronSelector->Init();
    mgr->AddTask(fElectronSelector);
    //connect input fElectronSelector

    mgr->ConnectInput (fElectronSelector,0,cinput);
  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //            find input container
  AliAnalysisTaskGammaConvDalitzV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvDalitzV1(Form("GammaConvDalitzV1_%i",trainConfig));
  task->SetIsHeavyIon(0);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
    if (enableQAMesonTask && lightVersion){
    cout << "\n\n**************************************************************" << endl;
    cout << "******* QA and Light version are incompatible, changing lightVersion to kFALSE........" << endl;
    cout << "**************************************************************\n\n" << endl;
    lightVersion=kFALSE;
  }
  task->SetDoLightVersion(lightVersion);



  if(trainConfig == 1){
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253202221710", "0263103100900000"); // standard cut number for pp7TeV
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253202221710", "0163103100900000");
    // train configs 2 to 4 for estimation of systematics of standard cut number pp7TeV
  } else if (trainConfig == 2) {
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "10475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "30475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253201221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253203221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253202121710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253302221710", "0263103100900000");
  } else if (trainConfig == 3) {
    cuts.AddCutPCMDalitz("00000113", "00200009360300007900004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007200004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007100004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300001800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300002800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200049360300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200019360300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20475400253202222710", "0263103100900000");
  } else if (trainConfig == 4) {
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20575400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20775400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009260300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009660300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20425400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009360300007800004000", "20407200253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009320300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00000113", "00200009305100007800004000", "20475400253202221710", "0263103100900000");
  } else if (trainConfig == 102) {
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "10475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "30475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "20475400253201221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "20475400253203221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "20475400253202121710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009360300007800004000", "20475400253202321710", "0263103100900000");
  } else if (trainConfig == 103) {
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "20475400253202301710", "0263103100900000");
  } else if (trainConfig == 104) {  // to be used with MBW
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "20475400253202301710", "0263103100900000");
  } else if (trainConfig == 105) {  // nominal B
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "10477400233202321710", "0263103100900000"); //kpiMinMomdedxSigmaTPCCut = 0.3,
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "10478400233202321710", "0263103100900000"); //kpiMinMomdedxSigmaTPCCut = 0.25,
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "10477400233202321510", "0263103100900000"); // massCut < 0.35 GeV/c^2
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "10477400233202321910", "0263103100900000"); // if pT < 1 GeV  massCut < 0.25 GeV/c^2 if pT > 1 GeV massCut < 0.35 GeV/c^2
  } else if (trainConfig == 106) {  // low B
    cuts.AddCutPCMDalitz("00010113", "00200089266300008854404000", "10477400233202301710", "0263103100900000"); //prim electron Pt > 0.075 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCutPCMDalitz("00010113", "00200089266300008854404000", "10477400233202361710", "0263103100900000"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if (trainConfig == 107) {   // R  scan
    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00a00009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00b00009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00c00009266300008854404000", "20475400253202321710", "0263103100900000");
  } else if (trainConfig == 108) {  // nominal B
    cuts.AddCutPCMDalitz("00010113", "00200009266300008850404000", "10477400233202321710", "0263103100900000"); //kpiMinMomdedxSigmaTPCCut = 0.3,
    cuts.AddCutPCMDalitz("00010113", "00200009266300008850404000", "10478400233202321710", "0263103100900000"); //kpiMinMomdedxSigmaTPCCut = 0.25,
    cuts.AddCutPCMDalitz("00010113", "00200009266300008850404000", "10477400233202321510", "0263103100900000"); // massCut < 0.35 GeV/c^2
    cuts.AddCutPCMDalitz("00010113", "00200009266300008850404000", "10477400233202321910", "0263103100900000"); // if pT < 1 GeV  massCut < 0.25 GeV/c^2 if pT > 1 GeV massCut < 0.35 GeV/c^2
  } else if (trainConfig == 109) {  // low B
    cuts.AddCutPCMDalitz("00010113", "00200089266300008850404000", "10477400233202301710", "0263103100900000"); //prim electron Pt > 0.075 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCutPCMDalitz("00010113", "00200089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if ( trainConfig ==110 ) { // low B

    cuts.AddCutPCMDalitz("00010113", "00200089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCutPCMDalitz("00010113", "00a00089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt >
    cuts.AddCutPCMDalitz("00010113", "00b00089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt
    cuts.AddCutPCMDalitz("00010113", "00c00089266300008850404000", "10477400233202361710", "0263103100900000");

  }  else if (trainConfig == 117) {   // as iConfig R scan to be used with MBW

    cuts.AddCutPCMDalitz("00010113", "00200009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00a00009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00b00009266300008854404000", "20475400253202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "00c00009266300008854404000", "20475400253202321710", "0263103100900000");

  } else if ( trainConfig ==120 ) { // low B

    cuts.AddCutPCMDalitz("00010113", "00200089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCutPCMDalitz("00010113", "00a00089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt >
    cuts.AddCutPCMDalitz("00010113", "00b00089266300008850404000", "10477400233202361710", "0263103100900000"); //prim electron Pt
    cuts.AddCutPCMDalitz("00010113", "00c00089266300008850404000", "10477400233202361710", "0263103100900000");

 } else if (trainConfig == 202) {
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "20475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "10475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "30475400253202221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "20475400253201221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "20475400253203221710", "0263103100900000");
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "20475400253202121710", "0263103100900000");
    cuts.AddCutPCMDalitz("00074113", "00200009360300007800004000", "20475400253202321710", "0263103100900000");
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// 5TeV 2017 and 13TeV low Magnetic Field////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //NOTE for low Magnetic Field we will run the same trains that we use for normal Magnetic Field just with one difference on the
    // pT of electrons and positrons, So a pretty simple approximation we will change B_{Normal}/B_{Low}=0.5/0.2=2.5,
    //with this factor we will recalculate the pT for the primary and secondary like 0.125/2.5=0.05.
 }  else if (trainConfig == 300) {//Primary cut (No Standard yet!), there are the standard cut from Lucia at 5TeV plus Pedro cuts for 7 Tev or 5.02 TeV pPb on electrons.
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263710", "0152103500000000");

    ///////////////////////////////Study on electrons cuts////////////////////////////////////
//0->8, 0.05 to 0.02, position 7
 }  else if (trainConfig == 301) {//Variation on the dE/dx standar -4 to 5
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "205c6400263202263710", "0152103500000000");//Primary, change de/dx TPCsigma -3,5
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "206c6400263202263710", "0152103500000000");//Primary, change de/dx TPCsigma -4,4
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "208c6400263202263710", "0152103500000000");//Primary, change de/dx TPCsigma -2,3.5
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "209c6400263202263710", "0152103500000000");//Primary, change de/dx TPCsigma -3,4
  }  else if (trainConfig == 302) {  //Variation on Mass_{e^+e^-}.0pT & 0.015 mass, >1.0p, Primary, change <1T & 0.035 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263a10", "0152103500000000");//Primary, change <1.0pT & 0.02 mass, >1.0pT & 0.03 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263b10", "0152103500000000");//Primary, change <1.0pT & 0.027 mass, >1.0pT & 0.057 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263c10", "0152103500000000");//Primary, change < 0.02 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263d10", "0152103500000000");//Primary, change < 0.04 mass
  }  else if (trainConfig == 303) {  //Variation on Mass_{e^+e^-}.0pT & 0.015 mass, >1.0p, Primary, change <1T & 0.035 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263e10", "0152103500000000");//Primary, change < 0.06 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263f10", "0152103500000000");//Primary, change < 0.08 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263g10", "0152103500000000");//Primary, change < 0.15 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263h10", "0152103500000000");//Primary, change < 0.2 mass
  //////////////// systematic variations for 5 TeV 2017 (Lucia)/////////////////////////////

  } else if (trainConfig == 310){
    cuts.AddCutPCMDalitz("00010113", "00200089227300008250404000", "204c6400263202263710", "0152103500000000"); // eta < 0.9
    cuts.AddCutPCMDalitz("00010113", "0c200089227300008250404000", "204c6400263202263710", "0152103500000000"); // eta < 0.85
    cuts.AddCutPCMDalitz("00010113", "0d200089247300008250404000", "204c6400263202263710", "0152103500000000"); // pidEdx 3, 1
    cuts.AddCutPCMDalitz("00010113", "0d200089227100008250404000", "204c6400263202263710", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 311){
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008860404000", "204c6400263202263710", "0152103500000000"); // variation chi2 20 psi pair 0.2 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089227300002252404000", "204c6400263202263710", "0152103500000000"); // variation qT max 0.06 2D, asym var 1
    cuts.AddCutPCMDalitz("00010113", "0d200089227300002254404000", "204c6400263202263710", "0152103500000000"); // variation qT max 0.06 2D, asym vat 2 pt dep
    cuts.AddCutPCMDalitz("00010113", "0d200089227300002256404000", "204c6400263202263710", "0152103500000000"); // variation qT max 0.06 2D, asym var 3
  } else if (trainConfig == 312){
    cuts.AddCutPCMDalitz("00010113", "0d200089a27300008250904120", "204c6400263202263710", "0152103500000000"); //cosPA, 0.99 eta 0.9
    cuts.AddCutPCMDalitz("00010113", "0d200079227300008250404000", "204c6400263202263710", "0152103500000000"); // min pT no cut
    cuts.AddCutPCMDalitz("00010113", "0d200049227300008250404000", "204c6400263202263710", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCMDalitz("00010113", "0d200019227300008250404000", "204c6400263202263710", "0152103500000000"); // min pT 100 MeV
    cuts.AddCutPCMDalitz("00010113", "0d200088227300008250404000", "204c6400263202263710", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCMDalitz("00010113", "0d200086227300008250404000", "204c6400263202263710", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 313){
    cuts.AddCutPCMDalitz("00010113", "0d200089327300008250404000", "204c6400263202263710", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCMDalitz("00010113", "0d200089627300008250404000", "204c6400263202263710", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCMDalitz("00010113", "0d200089257300008250404000", "204c6400263202263710", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCMDalitz("00010113", "0d200089217300008250404000", "204c6400263202263710", "0152103500000000"); // pidEdx 0,-10
    cuts.AddCutPCMDalitz("00010113", "0d200089226300008250404000", "204c6400263202263710", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMDalitz("00010113", "0d200089227600008250404000", "204c6400263202263710", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
  } else if (trainConfig == 314){
    cuts.AddCutPCMDalitz("00010113", "0d200089227300002250404000", "204c6400263202263710", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089227300009250404000", "204c6400263202263710", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008260404000", "204c6400263202263710", "0152103500000000"); // Psi pair 0.05 2D, chi2 30.
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008280404000", "204c6400263202263710", "0152103500000000"); // Psi pair 0.2  2D, chi2 30.
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008850404000", "204c6400263202263710", "0152103500000000"); // chi2 20. psi pair 0.1 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008150404000", "204c6400263202263710", "0152103500000000"); // chi2 50. psi pair 0.1 2D
  } else if (trainConfig == 315){
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008254404000", "204c6400263202263710", "0152103500000000"); // Photon Asymmetry Cut
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250604000", "204c6400263202263710", "0152103500000000"); // CosPA 0.9
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250004000", "204c6400263202263710", "0152103500000000"); // no CosPA
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250400000", "204c6400263202263710", "0152103500000000"); // no double counting
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263710", "0152101500000000"); // meson alpha pt dep
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263202263710", "0152107500000000"); // meson alpha < 0.85
  }  else if (trainConfig == 325) {//Primary cut, +with no mass c psi pair cut cuts.ut,no psi cut, +with
    cuts.AddCutPCMDalitz("00010113", "0d200089227300008250404000", "204c6400263002263010", "0152103500000000");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// 5TeV 2017 and 13TeV Normal Magnetic Field////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  }  else if (trainConfig == 400) {////Primary cut, Standard 2021 for 5 TeV and 13TeV
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    ///////////////////////////////Study on electrons cuts////////////////////////////////////
  }  else if (trainConfig == 401) {//Variation on the dE/dx standar -4 to 5
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "205c6400863f02223710", "0152103500000000");//Primary, change de/dx TPCsigma -3,5
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "206c6400863f02223710", "0152103500000000");//Primary, change de/dx TPCsigma -4,4
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c64008a3f02223710", "0152103500000000");//Number of TPC cluster standar > 80.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c64008b3f02223710", "0152103500000000");//Number of TPC cluster standar > 60.
  }  else if (trainConfig == 402) {  //Variation on Mass_{e^+e^-}.0pT & 0.015 mass, >1.0p, Primary, change <1T & 0.035 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223a10", "0152103500000000");//Primary, change <1.0pT & 0.02 mass, >1.0pT & 0.03 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223b10", "0152103500000000");//Primary, change <1.0pT & 0.027 mass, >1.0pT & 0.057 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223c10", "0152103500000000");//Primary, change < 0.02 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223d10", "0152103500000000");//Primary, change < 0.04 mass
  }  else if (trainConfig == 403) {  //Variation on Mass_{e^+e^-}.0pT & 0.015 mass, >1.0p, Primary, change <1T & 0.035 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223e10", "0152103500000000");//Primary, change < 0.06 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223f10", "0152103500000000");//Primary, change < 0.08 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223g10", "0152103500000000");//Primary, change < 0.15 mass
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223h10", "0152103500000000");//Primary, change < 0.2 mass
  }  else if (trainConfig == 404) {//No Primary cut, we change a 0 for 8 on second cut, and study the psi pair cut.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863002273010", "0152103500000000");//No Psi pair cut
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863202273010", "0152103500000000");//0.60 Psi pair cut
  }  else if (trainConfig == 405) {//Primary cut, with DCAxy(Distance of Closest Approach on plane xy) standar pt dependance, single pT standar >0.125GeV
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02229710", "0152103500000000");//DCAxy change to 0.8cm
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f0222a710", "0152103500000000");//DCAxy change to 1.2cm
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02213710", "0152103500000000");//Single pT change to > 1.0GeV
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02233710", "0152103500000000");//Single pT change to > 1.5GeV
  } else if (trainConfig == 406) {//Primary cut, with (Number of TPC cluster found)/(Number of TPC cluster findable) standar > 0.6, and Number of TPC cluster standar > 70.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400853f02223710", "0152103500000000");//(Number of TPC cluster found)/(Number of TPC cluster findable) > 0.35
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400813f02223710", "0152103500000000");//(Number of TPC cluster found)/(Number of TPC cluster findable) > 0.7
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "208c6400863f02223710", "0152103500000000");//Primary, change de/dx TPCsigma -2,3.5
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "209c6400863f02223710", "0152103500000000");//Primary, change de/dx TPCsigma -3,4
  } else if (trainConfig == 407) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "20466400863f02223710", "0152103500000000");//dE/dx Sigma [2.5,-10] pions
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "20476400863f02223710", "0152103500000000");//dE/dx Sigma [2.0,-10] pions
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c5400863f02223710", "0152103500000000");//min pT 0.5 GeV/c pions
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c7400863f02223710", "0152103500000000");//min pT 0.3 GeV/c pions
  } else if (trainConfig == 408) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6200863f02223710", "0152103500000000");//max pT 5 GeV/c pions
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6600863f02223710", "0152103500000000");//max pT 2 GeV/c pions
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863102223710", "0152103500000000");//PsiPair 0.45
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863602223710", "0152103500000000");//PsiPair 0.65
  } else if (trainConfig == 409) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//New Standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204a6200863f02223710", "0152103500000000");//
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204b6600863f02223710", "0152103500000000");//
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863102223710", "0152103500000000");//
  //////////////// systematic variations for 5 TeV 2017 (Lucia)/////////////////////////////

  } else if (trainConfig == 410){
    cuts.AddCutPCMDalitz("00010113", "00m00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // eta < 0.9
    cuts.AddCutPCMDalitz("00010113", "0cm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // eta < 0.85
    cuts.AddCutPCMDalitz("00010113", "0dm00009f4730000dge0404000", "204c6400863f02223710", "0152103500000000"); // pidEdx 3, 1
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9710000dge0404000", "204c6400863f02223710", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 411){
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000d860404000", "204c6400863f02223710", "0152103500000000"); // variation chi2 20 psi pair 0.2 2D
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300002ge2404000", "204c6400863f02223710", "0152103500000000"); // variation qT max 0.06 2D, asym var 1
    cuts.AddCutPCMDalitz("00010113", "0dm00009a9730000dge0904120", "204c6400863f02223710", "0152103500000000"); //cosPA, 0.99 eta 0.9
    cuts.AddCutPCMDalitz("00010113", "0dm00079f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // min pT no cut
  } else if (trainConfig == 412){
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300002ge4404000", "204c6400863f02223710", "0152103500000000"); // variation qT max 0.06 2D, asym vat 2 pt dep
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300002ge6404000", "204c6400863f02223710", "0152103500000000"); // variation qT max 0.06 2D, asym var 3
    cuts.AddCutPCMDalitz("00010113", "0dm00049f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCMDalitz("00010113", "0dm00019f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // min pT 100 MeV
  } else if (trainConfig == 413){
    cuts.AddCutPCMDalitz("00010113", "0dm0000939730000dge0404000", "204c6400863f02223710", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCMDalitz("00010113", "0dm0000969730000dge0404000", "204c6400863f02223710", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCMDalitz("00010113", "0dm00009f5730000dge0404000", "204c6400863f02223710", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCMDalitz("00010113", "0dm00009f1730000dge0404000", "204c6400863f02223710", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 414){
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300002250404000", "204c6400863f02223710", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300009250404000", "204c6400863f02223710", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300008260404000", "204c6400863f02223710", "0152103500000000"); // Psi pair 0.05 2D, chi2 30.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300008280404000", "204c6400863f02223710", "0152103500000000"); // Psi pair 0.2  2D, chi2 30.
  } else if (trainConfig == 415){
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge4404000", "204c6400863f02223710", "0152103500000000"); // Photon Asymmetry Cut
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0604000", "204c6400863f02223710", "0152103500000000"); // CosPA 0.9
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0004000", "204c6400863f02223710", "0152103500000000"); // no CosPA
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0400000", "204c6400863f02223710", "0152103500000000"); // no double counting
  } else if (trainConfig == 416) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Standard with kBoth on electrons
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863002223710", "0152103500000000");//No PsiPair cut with kBoth on electrons
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263f02223710", "0152103500000000");//kAny with new Psipair
//     cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263502223710", "0152103500000000");//Standard a 0.06
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263002223710", "0152103500000000");//No PsiPair on electrons
  } else if (trainConfig == 417) {//Systematic for TPC cluster and pion energy
    cuts.AddCutPCMDalitz("00010113", "0dm00008f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCMDalitz("00010113", "0dm00006f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9630000dge0404000", "204c6400863f02223710", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9760000dge0404000", "204c6400863f02223710", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
  } else if (trainConfig == 418) {//Systematic for chi2 psi pair on Photons
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000d8e0404000", "204c6400863f02223710", "0152103500000000"); // chi2 20. psi pair 0.1 2D
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000d1e0404000", "204c6400863f02223710", "0152103500000000"); // chi2 50. psi pair 0.1 2D
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152101500000000"); // meson alpha pt dep
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152107500000000"); // meson alpha < 0.85
  } else if (trainConfig == 419) {//Cross chech of efficiency show for range on TPC.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");// Standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");// exclusion of range o transition between chambers of TPC 52 to 72 cm.
  } else if (trainConfig == 420) {//Test with mass cut removed on primary selection
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02273010", "0152103500000000");
    //New Standard+0.9 single pt cut
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "20476400863f02283010", "0152103500000000");
    //New Standard+1.1 single pt cut + pion rejection improve
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "20476450863f02223010", "0152103500000000");
    //New Standard+Low Rejection for kaons protons and pions, 2,2,2.5 sigmas + pion rejection improve
    cuts.AddCutPCMDalitz("00010113", "0dm00009f97300007ge0404000", "204c6400863f02223710", "0152103500000000");
    //New Standard with qt < 0.15 flat
  } else if (trainConfig == 421) {//Test with mass cut removed on primary selection + Event selection
    cuts.AddCutPCMDalitz("0008d113", "0dm00009f9730000dge0404000", "204c6400863f02223010", "0152103500000000");
    //New Standard
    cuts.AddCutPCMDalitz("0008d113", "0dm00009f9730000dge0404000", "204c6400863f02273010", "0152103500000000");
    //New Standard+0.9 single pt cut
    cuts.AddCutPCMDalitz("0008d113", "0dm00009f9730000dge0404000", "204c6400863f02243010", "0152103500000000");
    //New Standard+0.5 single pt cut
    cuts.AddCutPCMDalitz("0008d113", "0dm00009f9730000dge0404000", "204c6420863f02223010", "0152103500000000");
    //New Standard+Low Rejection for pions kaons and protons
   } else if (trainConfig == 430) {//Psipair Optimization.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863202223710", "0152103500000000");//Standard with kBoth on electrons
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263202223710", "0152103500000000");//Standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263702223710", "0152103500000000");//kAny 5% GG lineal
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263802223710", "0152103500000000");//kAny 8% GG lineal
    } else if (trainConfig == 431) {//Psipair Optimization.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863e02223710", "0152103500000000");//kBoth 5 % GG
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//kBoth 8 % GG
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863h02223710", "0152103500000000");//kBoth pt dependance
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263g02223710", "0152103500000000");//kAny pt dependance
    } else if (trainConfig == 432) {//Psipair Optimization.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400863j02223710", "0152103500000000");// Nominal B kBoth 10%
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263l02223710", "0152103500000000");// Nominal B kAny 10%
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263m02223710", "0152103500000000");// Nominal B kAny 8%
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263n02223710", "0152103500000000");// Nominal B kAny 5%
    } else if (trainConfig == 433) {//kANY and kBOTH with shared layer
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400a63f02223710", "0152103500000000");//kAny 8 % GG shared 0, mass standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400a63f02223c10", "0152103500000000");//kAny 8 % GG shared 0, mass 0.02
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400b63f02223710", "0152103500000000");//kAny 8 % GG shared 1, mass standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400b63f02223c10", "0152103500000000");//kAny 8 % GG shared 1, mas 0.02
    } else if (trainConfig == 434) {//kANY and kBOTH with shared layer
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400c63f02223710", "0152103500000000");//kBoth 8 % GG shared 0, mass standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400c63f02223c10", "0152103500000000");//kBoth 8 % GG shared 0, mass 0.02
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400d63f02223710", "0152103500000000");//kBoth 8 % GG shared 1, mass standard
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400d63f02223c10", "0152103500000000");//kBoth 8 % GG shared 1, mas 0.02

    //V0M mult slices
  } else if (trainConfig == 440) {//10% mult slices part I Update
    cuts.AddCutPCMDalitz("n0115113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-10%
    cuts.AddCutPCMDalitz("n1215113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//10-20%
    cuts.AddCutPCMDalitz("n2315113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//20-30%
    cuts.AddCutPCMDalitz("n3515113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//30-50%
  } else if (trainConfig == 441) {//10% mult slices part II Update
    cuts.AddCutPCMDalitz("n5715113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//50-70%
    cuts.AddCutPCMDalitz("n7a15113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//70-100%
  } else if (trainConfig == 442) {//1-5-10% mult slices Update
    cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-1%
    cuts.AddCutPCMDalitz("m1515113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//1-5%
    cuts.AddCutPCMDalitz("m5a15113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//5-10%
  } else if (trainConfig == 443) {//0.01% mult slices
    cuts.AddCutPCMDalitz("r0115113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-0.01%
    cuts.AddCutPCMDalitz("r1515113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.01-0.05%
    cuts.AddCutPCMDalitz("r5a15113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.05-0.1%
    cuts.AddCutPCMDalitz("r0a15113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-1%
  } else if (trainConfig == 444) {//0.01% mult slices, V0M high mult triggered Update
    cuts.AddCutPCMDalitz("r0178113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-0.01%      V0M high mult trigger
    cuts.AddCutPCMDalitz("r1578113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.01-0.05%   V0M high mult trigger
    cuts.AddCutPCMDalitz("r5a78113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.05-0.1%    V0M high mult trigger
    cuts.AddCutPCMDalitz("r0a78113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-1%         V0M high mult trigger
  } else if (trainConfig == 445) {//10% mult slices 5 TeV Mike
    cuts.AddCutPCMDalitz("n0210113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-20%
    cuts.AddCutPCMDalitz("n2410113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//20-40%
    cuts.AddCutPCMDalitz("n4610113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//40-60%
    cuts.AddCutPCMDalitz("n6a10113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//60-100%

    //SPD mult slices
    } else if (trainConfig == 450) {//10% mult slices part I
    cuts.AddCutPCMDalitz("p0110113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-10%
    cuts.AddCutPCMDalitz("p1210113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//10-20%
    cuts.AddCutPCMDalitz("p2310113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//20-30%
    cuts.AddCutPCMDalitz("p3510113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//30-50%
    } else if (trainConfig == 451) {//10% mult slices part II
    cuts.AddCutPCMDalitz("p5710113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//50-70%
    cuts.AddCutPCMDalitz("p7a10113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//70-100%
    } else if (trainConfig == 452) {//1-5-10% mult slices
    cuts.AddCutPCMDalitz("o0110113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-1%
    cuts.AddCutPCMDalitz("o1510113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//1-5%
    cuts.AddCutPCMDalitz("o5a10113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//5-10%
    } else if (trainConfig == 453) {//0.1% mult slices
    cuts.AddCutPCMDalitz("u0110113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-0.1%
    cuts.AddCutPCMDalitz("u1510113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.1-0.5%
    cuts.AddCutPCMDalitz("u5a10113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.5-1%
    cuts.AddCutPCMDalitz("u0a10113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-1%
  } else if (trainConfig == 454) {//0.1% mult slices, SPD High Mult trigger
    cuts.AddCutPCMDalitz("u0175113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-0.1%   SPD High Mult trigger
    cuts.AddCutPCMDalitz("u1575113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.1-0.5% SPD High Mult trigger
    cuts.AddCutPCMDalitz("u5a75113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0.5-1%   SPD High Mult trigger
    cuts.AddCutPCMDalitz("u0a75113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-1%     SPD High Mult trigger

  } else if (trainConfig == 455) {//10% mult slices 5 TeV Mike
    cuts.AddCutPCMDalitz("p0210113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//0-20%
    cuts.AddCutPCMDalitz("p2410113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//20-40%
    cuts.AddCutPCMDalitz("p4610113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//40-60%
    cuts.AddCutPCMDalitz("p6a10113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//60-100%

  } else if (trainConfig == 459) {//Selection of multiplicties for efficiency studies
    cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); //0-1%
    cuts.AddCutPCMDalitz("m5a15113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); //5-10%
    cuts.AddCutPCMDalitz("n2315113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); //20-30%
    cuts.AddCutPCMDalitz("n7a15113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000"); //70-100%

  } else if (trainConfig == 460) {//Exclusive MC train with no Cuts
    cuts.AddCutPCMDalitz("00115113", "0dm00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");// MC train with no cut in the event

///////////////////////////////// Multiplicity selection, ITS + TRD photons 0-100% ///////////////////////////
} else if (trainConfig == 461) {// Systematics
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000dge0474000", "0152103500000000"); // Multiplicity selection, ITS + TRD photons
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000dge0674000", "0152103500000000"); // cosPointing angle 0.9
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000dge0374000", "0152103500000000"); // cosPointing angle 0.75
  cuts.AddCutPCMDalitz("00015113", "0dm0000969730000dge0474000", "0152103500000000"); // photons nsig electron -4,5

} else if (trainConfig == 462) {// Systematics
  cuts.AddCutPCMDalitz("00015113", "0dm0000939730000dge0474000", "0152103500000000"); // photons nsig electron -2.5,4
  cuts.AddCutPCMDalitz("00015113", "0dm00009f5730000dge0474000", "0152103500000000"); // photons nsig pion 2,-10
  cuts.AddCutPCMDalitz("00015113", "0dm00009f1730000dge0474000", "0152103500000000"); // photons nsig pion 0,-10
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000age0474000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.99

} else if (trainConfig == 463) {// Systematics
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000ege0474000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.99
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000dfe0474000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2)
  cuts.AddCutPCMDalitz("00015113", "0dm00009f9730000dgg0474000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2)

///////////////////////////////// Multiplicity selection, ITS + TRD photons 50-100% ///////////////////////////
} else if (trainConfig == 461) {// Systematics
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000dge0474000", "0152103500000000"); // Multiplicity selection, ITS + TRD photons
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000dge0674000", "0152103500000000"); // cosPointing angle 0.9
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000dge0374000", "0152103500000000"); // cosPointing angle 0.75
  cuts.AddCutPCMDalitz("n5a15113", "0dm0000969730000dge0474000", "0152103500000000"); // photons nsig electron -4,5

} else if (trainConfig == 462) {// Systematics
  cuts.AddCutPCMDalitz("n5a15113", "0dm0000939730000dge0474000", "0152103500000000"); // photons nsig electron -2.5,4
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f5730000dge0474000", "0152103500000000"); // photons nsig pion 2,-10
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f1730000dge0474000", "0152103500000000"); // photons nsig pion 0,-10
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000age0474000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.99

} else if (trainConfig == 463) {// Systematics
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000ege0474000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.99
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000dfe0474000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2)
  cuts.AddCutPCMDalitz("n5a15113", "0dm00009f9730000dgg0474000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2)

///////////////////////////////// Multiplicity selection, ITS + TRD photons 0-1% ///////////////////////////
} else if (trainConfig == 461) {// Systematics
  cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dge0474000", "0152103500000000"); // Multiplicity selection, ITS + TRD photons
  cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dge0674000", "0152103500000000"); // cosPointing angle 0.9
  cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dge0374000", "0152103500000000"); // cosPointing angle 0.75
  cuts.AddCutPCMDalitz("m0115113", "0dm0000969730000dge0474000", "0152103500000000"); // photons nsig electron -4,5

} else if (trainConfig == 462) {// Systematics
  cuts.AddCutPCMDalitz("m015113", "0dm0000939730000dge0474000", "0152103500000000"); // photons nsig electron -2.5,4
  cuts.AddCutPCMDalitz("m015113", "0dm00009f5730000dge0474000", "0152103500000000"); // photons nsig pion 2,-10
  cuts.AddCutPCMDalitz("m015113", "0dm00009f1730000dge0474000", "0152103500000000"); // photons nsig pion 0,-10
  cuts.AddCutPCMDalitz("m015113", "0dm00009f9730000age0474000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.99

} else if (trainConfig == 463) {// Systematics
  cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000ege0474000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.99
  cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dfe0474000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2)
  cuts.AddCutPCMDalitz("m0115113", "0dm00009f9730000dgg0474000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  6XX for lowB,    65X  lowB and MBW ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  } else if (trainConfig == 606) {  // low B   eta<0.8
    cuts.AddCutPCMDalitz("00010113", "0d200089266300008850404000", "10477400234202361710", "0263103100900000"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if ( trainConfig ==607 ) { // low B
    cuts.AddCutPCMDalitz("00010113", "0da00089266300008850404000", "10477400234202361710", "0263103100900000"); //prim electron Pt >
    cuts.AddCutPCMDalitz("00010113", "0db00089266300008850404000", "10477400234202361710", "0263103100900000"); //prim electron Pt
    cuts.AddCutPCMDalitz("00010113", "0dc00089266300008850404000", "10477400234202361710", "0263103100900000");
    // to be used with weights
  } else if (trainConfig == 656) {  // low B   eta<0.8
    cuts.AddCutPCMDalitz("00010113", "0d200089266300008850404000", "10477400234202361710", "0263103100900000"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if ( trainConfig ==657 ) { // low B
    cuts.AddCutPCMDalitz("00010113", "0da00089266300008850404000", "10477400234202361710", "0263103100900000"); //prim electron Pt >
    cuts.AddCutPCMDalitz("00010113", "0db00089266300008850404000", "10477400234202361710", "0263103100900000"); //prim electron Pt
    cuts.AddCutPCMDalitz("00010113", "0dc00089266300008850404000", "10477400234202361710", "0263103100900000");


    //  7XX for NomB,    75X  nomB and MBW
  } else if (trainConfig == 701) {  //Cross check on new cuts, modification had been made on standard.
    cuts.AddCutPCMDalitz("00010113", "0d200009f9730000dge0404000", "204c6400263202223710", "0152103500000000");
    //New Standard cut, Nico Implemantation(Gamma)
    cuts.AddCutPCMDalitz("00010113", "0d200009227300008250404000", "204c4640263202223710", "0152103500000000");
    //Standard used untill 21 Octuber 2019
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263202223710", "0152103500000000");
    //New Standard cut, Nico Implemantation(Gamma) + Ana TPC range remove 55-72
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dg70404000", "204c6400263202223710", "0152103500000000");
    //New Standard cut, Nico Implemantation(Gamma) + Ana TPC range remove 55-72 + 0.15 line on qT(Armenteros)
    cuts.AddCutPCMDalitz("00010113", "0d200009227300008250404000", "204c6400263202223710", "0152103500000000");
    //Standard cut + Correction on Low Rejection and momentum for pions = New Standard
  } else if (trainConfig == 706) {   // Nominal eta<0.8 for primactronsry and secondary ele   //  gamma asymmetry cut removed on 15.04.2019
    cuts.AddCutPCMDalitz("00010113", "0d200009266300008850404000", "20475400254202321710", "0263103100900000");
  } else if (trainConfig == 707) {    // R scan
    cuts.AddCutPCMDalitz("00010113", "0da00009266300008850404000", "20475400254202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "0db00009266300008850404000", "20475400254202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "0dc00009266300008850404000", "20475400254202321710", "0263103100900000");
  } else if (trainConfig == 714) {    // R scan
    cuts.AddCutPCMDalitz("00010113", "0dh00009266300008850404000", "20475400254202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "0di00009266300008850404000", "20475400254202321710", "0263103100900000");
  } else if (trainConfig == 719) {  //Removing range 55-72 of TPC range.
    cuts.AddCutPCMDalitz("00010113", "0dm00009f9730000dge0404000", "204c6400263202223710", "0152103500000000");
  } else if (trainConfig == 720) {  //New standar?
    cuts.AddCutPCMDalitz("00010113", "0d200009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
  } else if (trainConfig == 721) {  //Study on range of TPC
    cuts.AddCutPCMDalitz("00010113", "0da00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 5-33.5
    cuts.AddCutPCMDalitz("00010113", "0db00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 33.5-72
    cuts.AddCutPCMDalitz("00010113", "0dc00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 72-180
  } else if (trainConfig == 722) {  //still Study on range of TPC
    cuts.AddCutPCMDalitz("00010113", "0dh00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 5-13
    cuts.AddCutPCMDalitz("00010113", "0di00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 13-33.5
    cuts.AddCutPCMDalitz("00010113", "0dj00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 33.5-55
  } else if (trainConfig == 723) {  //still still Study on range of TPC
    cuts.AddCutPCMDalitz("00010113", "0dk00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 55-72
    cuts.AddCutPCMDalitz("00010113", "0dl00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 72-95
    cuts.AddCutPCMDalitz("00010113", "0dg00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
    //Range TPC 95-180
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// to be used with weights equal to 70X +50 /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  } else if (trainConfig == 756) {   // Nominal eta<0.8 for primactronsry and secondary ele   //  gamma asymmetry cut removed on 15.04.2019
    cuts.AddCutPCMDalitz("00010113", "0d200009266300008850404000", "20475400254202321710", "0263103100900000");
  } else if (trainConfig == 757) {    // R scan
    cuts.AddCutPCMDalitz("00010113", "0da00009266300008850404000", "20475400254202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "0db00009266300008850404000", "20475400254202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "0dc00009266300008850404000", "20475400254202321710", "0263103100900000");

  } else if (trainConfig == 764) {    // R scan
    cuts.AddCutPCMDalitz("00010113", "0dh00009266300008850404000", "20475400254202321710", "0263103100900000");
    cuts.AddCutPCMDalitz("00010113", "0di00009266300008850404000", "20475400254202321710", "0263103100900000");
  } else if (trainConfig == 770) {  //New standar?
    cuts.AddCutPCMDalitz("00010113", "0d200009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");
  } else if (trainConfig == 771) {  //Study on range of TPC
    cuts.AddCutPCMDalitz("00010113", "0da00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 5-33.5i
    cuts.AddCutPCMDalitz("00010113", "0db00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 33.5-72
    cuts.AddCutPCMDalitz("00010113", "0dc00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 72-180
  } else if (trainConfig == 772) {  //still Study on range of TPC
    cuts.AddCutPCMDalitz("00010113", "0dh00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 5-13
    cuts.AddCutPCMDalitz("00010113", "0di00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 13-33.5
    cuts.AddCutPCMDalitz("00010113", "0dj00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 33.5-55
  } else if (trainConfig == 773) {  //still still Study on range of TPC
    cuts.AddCutPCMDalitz("00010113", "0dk00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 55-72
    cuts.AddCutPCMDalitz("00010113", "0dl00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 72-95
    cuts.AddCutPCMDalitz("00010113", "0dg00009f9730000dge0404000", "204c6400863f02223710", "0152103500000000");//Range TPC 95-180
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////Study region of R on LowB Field/////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//NOTE for low Magnetic Field we will run the same trains that we use for normal Magnetic Field just with one difference on the pT of electrons and positrons, So a pretty simple approximation we will change B_{Normal}/B_{Low}=0.5/0.2=2.5, with this factor we will recalculate the pT for the primary and secondary like 0.125/2.5=0.05.
///////////////////////////////////////////////////////////////////////////////////
}  else if (trainConfig == 900) {////Primary cut, Standard 2021 for 5 TeV and 13TeV
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");
    ///////////////////////////////Study on electrons cuts////////////////////////////////////
  }  else if (trainConfig == 901) {//Variation on the dE/dx standar -4 to 5
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "205c6400863d02263710", "0152103500000000");//Primary, change de/dx TPCsigma -3,5
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "206c6400863d02263710", "0152103500000000");//Primary, change de/dx TPCsigma -4,4
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c64008a3d02263710", "0152103500000000");//Number of TPC cluster standar > 80.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c64008b3d02263710", "0152103500000000");//Number of TPC cluster standar > 60.
  }  else if (trainConfig == 902) {  //Variation on Mass_{e^+e^-}.0pT & 0.015 mass, >1.0p, Primary, change <1T & 0.035 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263a10", "0152103500000000");//Primary, change <1.0pT & 0.02 mass, >1.0pT & 0.03 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263b10", "0152103500000000");//Primary, change <1.0pT & 0.027 mass, >1.0pT & 0.057 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263c10", "0152103500000000");//Primary, change < 0.02 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263d10", "0152103500000000");//Primary, change < 0.04 mass
  }  else if (trainConfig == 903) {  //Variation on Mass_{e^+e^-}.0pT & 0.015 mass, >1.0p, Primary, change <1T & 0.035 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263e10", "0152103500000000");//Primary, change < 0.06 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263f10", "0152103500000000");//Primary, change < 0.08 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263g10", "0152103500000000");//Primary, change < 0.15 mass
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263h10", "0152103500000000");//Primary, change < 0.2 mass
  }  else if (trainConfig == 904) {//No Primary cut, we change a 0 for 8 on second cut, and study the psi pair cut.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863002263010", "0152103500000000");//No Psi pair cut
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863202263010", "0152103500000000");//0.60 Psi pair cut
  }  else if (trainConfig == 905) {//Primary cut, with DCAxy(Distance of Closest Approach on plane xy) standar pt dependance, single pT standar >0.125GeV
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02269710", "0152103500000000");//DCAxy change to 0.8cm
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d0226a710", "0152103500000000");//DCAxy change to 1.2cm
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02213710", "0152103500000000");//Single pT change to > 1.0GeV
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02293710", "0152103500000000");//Single pT change to > 0.02 GeV
  } else if (trainConfig == 906) {//Primary cut, with (Number of TPC cluster found)/(Number of TPC cluster findable) standar > 0.6, and Number of TPC cluster standar > 70.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400853d02263710", "0152103500000000");//(Number of TPC cluster found)/(Number of TPC cluster findable) > 0.35
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400813d02263710", "0152103500000000");//(Number of TPC cluster found)/(Number of TPC cluster findable) > 0.7
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "208c6400863d02263710", "0152103500000000");//Primary, change de/dx TPCsigma -2,3.5
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "209c6400863d02263710", "0152103500000000");//Primary, change de/dx TPCsigma -3,4
  } else if (trainConfig == 907) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "20466400863d02263710", "0152103500000000");//dE/dx Sigma [2.5,-10] pions
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "20476400863d02263710", "0152103500000000");//dE/dx Sigma [2.0,-10] pions
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c5400863d02263710", "0152103500000000");//min pT 0.5 GeV/c pions
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c7400863d02263710", "0152103500000000");//min pT 0.3 GeV/c pions
  } else if (trainConfig == 908) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//Standard
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6200863d02263710", "0152103500000000");//max pT 5 GeV/c pions
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6600863d02263710", "0152103500000000");//max pT 2 GeV/c pions
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863102263710", "0152103500000000");//PsiPair 0.45
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863602263710", "0152103500000000");//PsiPair 0.65
  } else if (trainConfig == 909) {//Primary cut, dE/dx sigma, min pT, max pT for pions.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204a6200863d02263710", "0152103500000000");//
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204b6600863d02263710", "0152103500000000");//
  //////////////// systematic variations for Gammas low B /////////////////////////////
  } else if (trainConfig == 910){
    cuts.AddCutPCMDalitz("00010113", "00200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // eta < 0.9
    cuts.AddCutPCMDalitz("00010113", "0c200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // eta < 0.85
    cuts.AddCutPCMDalitz("00010113", "0d200089f4730000iih0404000", "204c6400863d02263710", "0152103500000000"); // pidEdx 3, 1
    cuts.AddCutPCMDalitz("00010113", "0d200089f9710000iih0404000", "204c6400863d02263710", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 911){
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000i860404000", "204c6400863d02263710", "0152103500000000"); // variation chi2 20 psi pair 0.2 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089f97300002ih2404000", "204c6400863d02263710", "0152103500000000"); // variation qT max 0.06 2D, asym var 1
    cuts.AddCutPCMDalitz("00010113", "0d200089a9730000iih0904120", "204c6400863d02263710", "0152103500000000"); //cosPA, 0.99 eta 0.9
    cuts.AddCutPCMDalitz("00010113", "0d200079f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // min pT no cut
  } else if (trainConfig == 912){
    cuts.AddCutPCMDalitz("00010113", "0d200089f97300002ih4404000", "204c6400863d02263710", "0152103500000000"); // variation qT max 0.06 2D, asym vat 2 pt dep
    cuts.AddCutPCMDalitz("00010113", "0d200089f97300002ih6404000", "204c6400863d02263710", "0152103500000000"); // variation qT max 0.06 2D, asym var 3
    cuts.AddCutPCMDalitz("00010113", "0d200049f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCMDalitz("00010113", "0d200019f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // min pT 100 MeV
  } else if (trainConfig == 913){
    cuts.AddCutPCMDalitz("00010113", "0d20008939730000iih0404000", "204c6400863d02263710", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCMDalitz("00010113", "0d20008969730000iih0404000", "204c6400863d02263710", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCMDalitz("00010113", "0d200089f5730000iih0404000", "204c6400863d02263710", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCMDalitz("00010113", "0d200089f1730000iih0404000", "204c6400863d02263710", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 914){
    cuts.AddCutPCMDalitz("00010113", "0d200089f97300002ih0404000", "204c6400863d02263710", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089f97300009ih0404000", "204c6400863d02263710", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000ii60404000", "204c6400863d02263710", "0152103500000000"); // Psi pair 0.05 2D, chi2 30.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000ii80404000", "204c6400863d02263710", "0152103500000000"); // Psi pair 0.2  2D, chi2 30.
  } else if (trainConfig == 915){
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih4404000", "204c6400863d02263710", "0152103500000000"); // Photon Asymmetry Cut
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0604000", "204c6400863d02263710", "0152103500000000"); // CosPA 0.9
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0004000", "204c6400863d02263710", "0152103500000000"); // no CosPA
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0400000", "204c6400863d02263710", "0152103500000000"); // no double counting
  } else if (trainConfig == 916) {// Study Low B Field kBoth for New PsiPair and kBoth (ITS hits)
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//Standard with kBoth on electrons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863002263710", "0152103500000000");//No PsiPair cut with kBoth on electrons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263d02263710", "0152103500000000");//kAny with new Psipair
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263002263710", "0152103500000000");//No PsiPair on electrons
  } else if (trainConfig == 917) {//Systematic for TPC cluster and pion energy
    cuts.AddCutPCMDalitz("00010113", "0d200008f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCMDalitz("00010113", "0d200006f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCMDalitz("00010113", "0d200089f9630000iih0404000", "204c6400863d02263710", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMDalitz("00010113", "0d200089f9760000iih0404000", "204c6400863d02263710", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
  } else if (trainConfig == 918) {//Systematic for chi2 psi pair on Photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000i8h0404000", "204c6400863d02263710", "0152103500000000"); // chi2 20. psi pair 0.1 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000i1h0404000", "204c6400863d02263710", "0152103500000000"); // chi2 50. psi pair 0.1 2D
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152101500000000"); // meson alpha pt dep
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152107500000000"); // meson alpha < 0.85
  } else if (trainConfig == 919) { // Study Low B Field
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCMDalitz("00010113", "0dd00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
//////////////////////////// TPC studies in Ranges///////////////////////////////////////////
  } else if (trainConfig == 920) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // all range 0 to 180
  } else if (trainConfig == 921) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0da00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 5-33.5
    cuts.AddCutPCMDalitz("00010113", "0db00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 33.5-72
    cuts.AddCutPCMDalitz("00010113", "0dc00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 72-180
  } else if (trainConfig == 922) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0dh00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 5-13
    cuts.AddCutPCMDalitz("00010113", "0di00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 13-33.5
    cuts.AddCutPCMDalitz("00010113", "0dj00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 33.5-55
  } else if (trainConfig == 923) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0dk00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 55-72
    cuts.AddCutPCMDalitz("00010113", "0dl00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 72-95
    cuts.AddCutPCMDalitz("00010113", "0dg00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 95-180
  } else if (trainConfig == 928) { // R 5-180  // Cat 1, cat 2+3   Meson Cat >=2
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0424000", "204c6400263202263710", "0152103500000000"); // Test Jpsi cat 1
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0454000", "204c6400263202263710", "0152103500000000"); // Test Jpsi cat 2
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0424000", "204c6400263202263710", "0152103500000000"); // Test Jpsi cat 3
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0454000", "204c6400263202263710", "0152103500000000"); // Test Jpsi cat 2+3
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0404000", "204c6400263202263710", "0152103520000000"); // Test Jpsi cat 1+2
  } else if (trainConfig == 930) { // Study Low B Field kBoth for New PsiPair and kBoth, Optimization at 5 and 8 % of contamiantion.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263a02263710", "0152103500000000"); //kAny 5 % of contamination from photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263b02263710", "0152103500000000"); //kAny 8 % of contamination from photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863c02263710", "0152103500000000"); //kBoth 5 % of contamination from photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); //kBoth 8 % of contamination from photons
  } else if (trainConfig == 931) { // Study Low Contamination.
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263i02263710", "0152103500000000"); //kAny 10 % of contamination from photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863k02263710", "0152103500000000"); //kBoth 10 % of contamination from photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263202263710", "0152103500000000"); //kAny Standard, contamination from photons
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863202263710", "0152103500000000"); //kBoth Standard, contamination from photons
  } else if (trainConfig == 932) {//kANY and kBOTH with shared layer204c6400863d02263710
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400a63d02263710", "0152103500000000");//kAny 8 % GG shared 0, mass standard
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400a63d02263c10", "0152103500000000");//kAny 8 % GG shared 0, mass 0.02
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400b63d02263710", "0152103500000000");//kAny 8 % GG shared 1, mass standard
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400b63d02263c10", "0152103500000000");//kAny 8 % GG shared 1, mas 0.02
  } else if (trainConfig == 933) {//kANY and kBOTH with shared layer
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400c63d02263710", "0152103500000000");//kBoth 8 % GG shared 0, mass standard
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400c63d02263c10", "0152103500000000");//kBoth 8 % GG shared 0, mass 0.02
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400d63d02263710", "0152103500000000");//kBoth 8 % GG shared 1, mass standard
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400d63d02263c10", "0152103500000000");//kBoth 8 % GG shared 1, mas 0.02

    //V0M mult slices
    } else if (trainConfig == 940) {//10% mult slices part I
    cuts.AddCutPCMDalitz("n0110113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-10%
    cuts.AddCutPCMDalitz("n1210113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//10-20%
    cuts.AddCutPCMDalitz("n2310113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//20-30%
    cuts.AddCutPCMDalitz("n3510113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//30-50%
    } else if (trainConfig == 941) {//10% mult slices part II
    cuts.AddCutPCMDalitz("n5710113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//50-70%
    cuts.AddCutPCMDalitz("n7a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//70-100%
    } else if (trainConfig == 942) {//1-5-10% mult slices
    cuts.AddCutPCMDalitz("m0110113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-1%
    cuts.AddCutPCMDalitz("m1510113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//1-5%
    cuts.AddCutPCMDalitz("m5a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//5-10%
    } else if (trainConfig == 943) {//0.01% mult slices
    cuts.AddCutPCMDalitz("r0110113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-0.01%
    cuts.AddCutPCMDalitz("r1510113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0.01-0.05%
    cuts.AddCutPCMDalitz("r5a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0.05-0.1%
    cuts.AddCutPCMDalitz("r0a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-1%
    } else if (trainConfig == 944) {//10% mult slices 5 TeV Mike
    cuts.AddCutPCMDalitz("n0210113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-20%
    cuts.AddCutPCMDalitz("n2410113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//20-40%
    cuts.AddCutPCMDalitz("n4610113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//40-60%
    cuts.AddCutPCMDalitz("n6a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//60-100%

    //SPD mult slices
    } else if (trainConfig == 950) {//10% mult slices part I
    cuts.AddCutPCMDalitz("p0110113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-10%
    cuts.AddCutPCMDalitz("p1210113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//10-20%
    cuts.AddCutPCMDalitz("p2310113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//20-30%
    cuts.AddCutPCMDalitz("p3510113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//30-50%
    } else if (trainConfig == 951) {//10% mult slices part II
    cuts.AddCutPCMDalitz("p5710113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//50-70%
    cuts.AddCutPCMDalitz("p7a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//70-100%
    } else if (trainConfig == 952) {//1-5-10% mult slices
    cuts.AddCutPCMDalitz("o0110113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-1%
    cuts.AddCutPCMDalitz("o1510113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//1-5%
    cuts.AddCutPCMDalitz("o5a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//5-10%
    } else if (trainConfig == 953) {//0.01% mult slices
    cuts.AddCutPCMDalitz("u0110113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-0.01%
    cuts.AddCutPCMDalitz("u1510113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0.01-0.05%
    cuts.AddCutPCMDalitz("u5a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0.05-0.1%
    cuts.AddCutPCMDalitz("u0a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-1%
    } else if (trainConfig == 954) {//10% mult slices 5 TeV Mike
    cuts.AddCutPCMDalitz("p0210113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//0-20%
    cuts.AddCutPCMDalitz("p2410113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//20-40%
    cuts.AddCutPCMDalitz("p4610113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//40-60%
    cuts.AddCutPCMDalitz("p6a10113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000");//60-100%

    //-----------------same as 6XX to be used with MBW extracted from 5TeV Nch
  } else if (trainConfig == 969) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400263202263710", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0404000", "204c6400263202263710", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCMDalitz("00010113", "0dd00089f9730000iih0404000", "204c6400263202263710", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 970) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // all range 0 to 180
  } else if (trainConfig == 971) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0da00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 5-33.5
    cuts.AddCutPCMDalitz("00010113", "0db00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 33.5-72
    cuts.AddCutPCMDalitz("00010113", "0dc00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 72-180
  } else if (trainConfig == 972) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0dh00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 5-13
    cuts.AddCutPCMDalitz("00010113", "0di00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 13-33.5
    cuts.AddCutPCMDalitz("00010113", "0dj00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 33.5-55
  } else if (trainConfig == 973) { // R 5-180
    cuts.AddCutPCMDalitz("00010113", "0dk00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 55-72
    cuts.AddCutPCMDalitz("00010113", "0dl00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 72-95
    cuts.AddCutPCMDalitz("00010113", "0dg00089f9730000iih0404000", "204c6400863d02263710", "0152103500000000"); // Range TPC 95-180
  } else if (trainConfig == 978) { // R 5-180  // Cat 1, cat 2+3, Meson Cat >= 2
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0424000", "204c6400263202263710", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCMDalitz("00010113", "0d200089f9730000iih0454000", "204c6400263202263710", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0424000", "204c6400263202263710", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0454000", "204c6400263202263710", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCMDalitz("00010113", "0dm00089f9730000iih0404000", "204c6400263202263710", "0152103520000000"); // eta < 0.8  // Test improved cuts
  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvDalitz! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList  *EventCutList = new TList();
  TList  *ConvCutList  = new TList();
  TList  *MesonCutList = new TList();
  TList  *ElecCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header2 = new TObjString("BOX");
  HeaderList->Add(Header2);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts   = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisConvPhotonsCuts   = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  ElecCutList->SetOwner(kTRUE);
  AliDalitzElectronCuts **analysisElecCuts   = new AliDalitzElectronCuts*[numberOfCuts];
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    TString cutName( Form("%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetElectronCut(i)).Data(),(cuts.GetMesonCut(i)).Data() ) );

    //Inizialtion of Events, Mesons, Conv Photons and Electrons cuts
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisConvPhotonsCuts[i] = new AliConversionPhotonCuts();
    analysisElecCuts[i] = new AliDalitzElectronCuts();

    ////////////////////////////////////
    //// LightOutputs for selection ////
    ////////////////////////////////////

    if (lightVersion){
        analysisEventCuts[i]->SetLightOutput(2);
        analysisConvPhotonsCuts[i]->SetLightOutput(2);
        analysisMesonCuts[i]->SetLightOutput(2);
        analysisElecCuts[i]->SetDoLightVersion(kTRUE);
    }

    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    //////////////////////////////////////
    //// Material Budget Weights flag ////
    //////////////////////////////////////

    if (enableMatBudWeightsPi0 > 0){
        if (isMC > 0){
            if (analysisConvPhotonsCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;}
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

    analysisConvPhotonsCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisConvPhotonsCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisConvPhotonsCuts[i]);
    analysisConvPhotonsCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");

    ////////////////////////////////////
    //// Elec dE/dx PostCalibration ////
    ////////////////////////////////////

    if (enableElecdEdxPostCalibration){
        if (isMC == 0){
            analysisConvPhotonsCuts[i]->SetDoElecDeDxPostCalibration(kTRUE);
            analysisElecCuts[i]->SetDoElecDeDxPostCalibrationPrimaryPair(kTRUE);
            if(ForceRPostCalibration){
              //Dependance in range 0-33 on R, and not TPC dependace
            analysisConvPhotonsCuts[i]->ForceTPCRecalibrationAsFunctionOfConvR();
            analysisElecCuts[i]->ForceTPCRecalibrationAsFunctionOfConvRPrimaryPair();
            }
        } else{
        cout << "ERROR enableElecdEdxPostCalibration set to True even if it is MC file. Automatically reset to 0"<< endl;
        analysisConvPhotonsCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
        analysisElecCuts[i]->SetDoElecDeDxPostCalibrationPrimaryPair(kFALSE);
      }
    }

    analysisElecCuts[i]->InitializeCutsFromCutString((cuts.GetElectronCut(i)).Data());
    ElecCutList->Add(analysisElecCuts[i]);
    analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,cutName);
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetElectronCutList(ElecCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }
  if (trainConfig == 325 || trainConfig == 425) {
    task->SetDoChicAnalysis(kTRUE);
  }
    ////////////////////////////////////
    /////////// QA Meson task //////////
    ////////////////////////////////////
  if (enableQAMesonTask) task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  //if (enableQAMesonTask) task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvDalitzV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvDalitzV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
