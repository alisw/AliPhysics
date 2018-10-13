
//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************

class CutHandlerConvDalitz{
  public:
    CutHandlerConvDalitz(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; elecCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; elecCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString elecCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvDalitz: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || elecCut.Length()!=20 ) {
        cout << "ERROR in CutHandlerConvDalitz: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; elecCutArray[nCuts]=elecCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvDalitz: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvDalitz: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConvDalitz: GetMesonCut wrong index i" << endl;return "";}}
    TString GetElecCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return elecCutArray[i]; else {cout << "ERROR in CutHandlerConvDalitz: GetElecCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* elecCutArray;
};


void AddTask_GammaConvDalitzV1_pp(  Int_t trainConfig = 1,  //change different set of cuts
                                    Bool_t isMC   = kFALSE, //run MC
                                    Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                    TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                    Int_t   enableMatBudWeightsPi0          = 0,              // 1 = three radial bins, 2 = 10 radial bins
                                    TString filenameMatBudWeights           = "MCInputFileMaterialBudgetWeights.root",
                                    TString periodNameV0Reader              = "",
                                    TString   additionalTrainConfig         = "0"
         ) {

  Int_t isHeavyIon = 0;
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_pp during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0){ rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    } else {
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout<< tempStr.Data()<<endl;

      if(tempStr.Contains("MaterialBudgetWeights") && enableMatBudWeightsPi0 > 0){
         if(tempStr.Contains("MaterialBudgetWeightsNONE")){
            enableMatBudWeightsPi0 = 0;
            cout << "INFO:  AddTask_GammaConvV1_pp materialBudgetWeights switched off signaled by additionalTrainConfigFlag" << endl;
         } else {
            TObjArray *fileNameMatBudWeightsArr = filenameMatBudWeights.Tokenize("/");
            if(fileNameMatBudWeightsArr->GetEntries()<1 ){
                cout<<"ERROR: AddTask_GammaConvV1_pp when reading material budget weights file name" << filenameMatBudWeights.Data()<< "'" << endl;
                return;
            }
            TObjString * oldMatObjStr = (TObjString*)fileNameMatBudWeightsArr->At( fileNameMatBudWeightsArr->GetEntries()-1);
            TString  oldfileName  = oldMatObjStr->GetString();
            TString  newFileName  = Form("MCInputFile%s.root",tempStr.Data());
            cout<<newFileName.Data()<<endl;
            if( oldfileName.EqualTo(newFileName.Data()) == 0 ){
              filenameMatBudWeights.ReplaceAll(oldfileName.Data(),newFileName.Data());
              cout << "INFO: AddTask_GammaConvV1_pp the material budget weights file has been change to " <<filenameMatBudWeights.Data()<<"'"<< endl;
          }
        }
      }
    }
  }

  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
  }
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;



  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvDalitzV1_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }



  //=========  Set Cutnumber for V0Reader ================================
  //TString cutnumber = "00000000000840010015000000";
  TString cutnumberPhoton     = "00000008400000000100000000";
  if (  periodNameV0Reader.CompareTo("LHC16f") == 0 || periodNameV0Reader.CompareTo("LHC17g")==0 || periodNameV0Reader.CompareTo("LHC18c")==0 ||
        periodNameV0Reader.CompareTo("LHC17d1") == 0  || periodNameV0Reader.CompareTo("LHC17d12")==0 ||
        periodNameV0Reader.CompareTo("LHC17h3")==0 || periodNameV0Reader.CompareTo("LHC17k1")==0 ||
        periodNameV0Reader.CompareTo("LHC17f8b") == 0 ||
        periodNameV0Reader.CompareTo("LHC16P1JJLowB") == 0 || periodNameV0Reader.CompareTo("LHC16P1Pyt8LowB") == 0 )
    cutnumberPhoton         = "00000088400000000100000000";

  TString cutnumberEvent  = "00000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){

    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);

    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(cutnumberPhoton!=""){
      fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

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


  CutHandlerConvDalitz cuts;

  if(trainConfig == 1){
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253202221710"); // standard cut number for pp7TeV
    cuts.AddCut("00000113", "00200009360300007800004000", "0163103100900000", "20475400253202221710");
    // train configs 2 to 4 for estimation of systematics of standard cut number pp7TeV
  } else if (trainConfig == 2) {
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "10475400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "30475400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253201221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253203221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253202121710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253302221710");
  } else if (trainConfig == 3) {
    cuts.AddCut("00000113", "00200009360300007900004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300007200004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300007100004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300001800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300002800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200049360300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200019360300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20475400253202222710");
  } else if (trainConfig == 4) {
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20575400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20775400253202221710");
    cuts.AddCut("00000113", "00200009260300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009660300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20425400253202221710");
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000", "20407200253202221710");
    cuts.AddCut("00000113", "00200009320300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00000113", "00200009305100007800004000", "0263103100900000", "20475400253202221710");
  } else if (trainConfig == 102) {
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "10475400253202221710");
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "30475400253202221710");
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "20475400253201221710");
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "20475400253203221710");
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "20475400253202121710");
    cuts.AddCut("00010113", "00200009360300007800004000", "0263103100900000", "20475400253202321710");
  } else if (trainConfig == 103) {
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "20475400253202301710");
  } else if (trainConfig == 104) {  // to be used with MBW
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "20475400253202301710");
  } else if (trainConfig == 105) {  // nominal B
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "10477400233202321710"); //kpiMinMomdedxSigmaTPCCut = 0.3,
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "10478400233202321710"); //kpiMinMomdedxSigmaTPCCut = 0.25,
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "10477400233202321510"); // massCut < 0.35 GeV/c^2
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "10477400233202321910"); // if pT < 1 GeV  massCut < 0.25 GeV/c^2 if pT > 1 GeV massCut < 0.35 GeV/c^2
  } else if (trainConfig == 106) {  // low B
    cuts.AddCut("00010113", "00200089266300008854404000", "0263103100900000", "10477400233202301710"); //prim electron Pt > 0.075 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCut("00010113", "00200089266300008854404000", "0263103100900000", "10477400233202361710"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if (trainConfig == 107) {   // R  scan
    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00a00009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00b00009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00c00009266300008854404000", "0263103100900000", "20475400253202321710");
  } else if (trainConfig == 108) {  // nominal B
    cuts.AddCut("00010113", "00200009266300008850404000", "0263103100900000", "10477400233202321710"); //kpiMinMomdedxSigmaTPCCut = 0.3,
    cuts.AddCut("00010113", "00200009266300008850404000", "0263103100900000", "10478400233202321710"); //kpiMinMomdedxSigmaTPCCut = 0.25,
    cuts.AddCut("00010113", "00200009266300008850404000", "0263103100900000", "10477400233202321510"); // massCut < 0.35 GeV/c^2
    cuts.AddCut("00010113", "00200009266300008850404000", "0263103100900000", "10477400233202321910"); // if pT < 1 GeV  massCut < 0.25 GeV/c^2 if pT > 1 GeV massCut < 0.35 GeV/c^2
  } else if (trainConfig == 109) {  // low B
    cuts.AddCut("00010113", "00200089266300008850404000", "0263103100900000", "10477400233202301710"); //prim electron Pt > 0.075 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCut("00010113", "00200089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if ( trainConfig ==110 ) { // low B

    cuts.AddCut("00010113", "00200089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCut("00010113", "00a00089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt >
    cuts.AddCut("00010113", "00b00089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt
    cuts.AddCut("00010113", "00c00089266300008850404000", "0263103100900000", "10477400233202361710");

  }  else if (trainConfig == 117) {   // as iConfig R scan to be used with MBW

    cuts.AddCut("00010113", "00200009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00a00009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00b00009266300008854404000", "0263103100900000", "20475400253202321710");
    cuts.AddCut("00010113", "00c00009266300008854404000", "0263103100900000", "20475400253202321710");

  } else if ( trainConfig ==120 ) { // low B

    cuts.AddCut("00010113", "00200089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
    cuts.AddCut("00010113", "00a00089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt >
    cuts.AddCut("00010113", "00b00089266300008850404000", "0263103100900000", "10477400233202361710"); //prim electron Pt
    cuts.AddCut("00010113", "00c00089266300008850404000", "0263103100900000", "10477400233202361710");

 } else if (trainConfig == 202) {
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253202221710");
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "10475400253202221710");
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "30475400253202221710");
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253201221710");
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253203221710");
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253202121710");
    cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253202321710");
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// 5TeV 2017 and 13TeV low Magnetic Field////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //NOTE for low Magnetic Field we will run the same trains that we use for normal Magnetic Field just with one difference on the
    // pT of electrons and positrons, So a pretty simple approximation we will change B_{Normal}/B_{Low}=0.5/0.2=2.5,
    //with this factor we will recalculate the pT for the primary and secondary like 0.125/2.5=0.05 (primary) and 0.05/2.5=0.02(se).
 }  else if (trainConfig == 300) {//Primary cut (No Standard yet!), there are the standard cut from Lucia at 5TeV plus Pedro cuts for 7 Tev or 5.02 TeV pPb on electrons.
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263710");

    ///////////////////////////////Study on electrons cuts////////////////////////////////////
//0->8, 0.05 to 0.02, position 7
 }  else if (trainConfig == 301) {//Variation on the dE/dx standar -4 to 5
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "205b6600263202263710");//Primary, change de/dx TPCsigma -3,5
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "206b6600263202263710");//Primary, change de/dx TPCsigma -4,4
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "208b6600263202263710");//Primary, change de/dx TPCsigma -2,3.5
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "209b6600263202263710");//Primary, change de/dx TPCsigma -3,4
  }  else if (trainConfig == 302) {  //Variation on Mass_{e^+e^-}, Primary, change <1.0pT & 0.015 mass, >1.0pT & 0.035 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263a10");//Primary, change <1.0pT & 0.02 mass, >1.0pT & 0.03 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263b10");//Primary, change <1.0pT & 0.027 mass, >1.0pT & 0.057 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263c10");//Primary, change < 0.02 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263d10");//Primary, change < 0.04 mass
  }  else if (trainConfig == 303) {  //Variation on Mass_{e^+e^-}, Primary, change <1.0pT & 0.015 mass, >1.0pT & 0.035 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263e10");//Primary, change < 0.06 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263f10");//Primary, change < 0.08 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263g10");//Primary, change < 0.15 mass
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263h10");//Primary, change < 0.2 mass
}  else if (trainConfig == 300) {//Primary cut (No Standard yet!), there are the standard cut from Lucia at 5TeV plus Pedro cuts for 7 Tev or 5.02 TeV pPb on electrons.
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202263710");
  //////////////// systematic variations for 5 TeV 2017 (Lucia)/////////////////////////////

  } else if (trainConfig == 310){
    cuts.AddCut("00010113", "00200089227300008250404000", "0152103500000000", "204b6600263202263710"); // eta < 0.9
    cuts.AddCut("00010113", "0c200089227300008250404000", "0152103500000000", "204b6600263202263710"); // eta < 0.85
    cuts.AddCut("00010113", "0d200089247300008250404000", "0152103500000000", "204b6600263202263710"); // pidEdx 3, 1
    cuts.AddCut("00010113", "0d200089227100008250404000", "0152103500000000", "204b6600263202263710"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 311){
    cuts.AddCut("00010113", "0d200089227300008860404000", "0152103500000000", "204b6600263202263710"); // variation chi2 20 psi pair 0.2 2D
    cuts.AddCut("00010113", "0d200089227300002252404000", "0152103500000000", "204b6600263202263710"); // variation qT max 0.06 2D, asym var 1
    cuts.AddCut("00010113", "0d200089227300002254404000", "0152103500000000", "204b6600263202263710"); // variation qT max 0.06 2D, asym vat 2 pt dep
    cuts.AddCut("00010113", "0d200089227300002256404000", "0152103500000000", "204b6600263202263710"); // variation qT max 0.06 2D, asym var 3
  } else if (trainConfig == 312){
    cuts.AddCut("00010113", "0d200089a27300008250904120", "0152103500000000", "204b6600263202263710"); //cosPA, 0.99 eta 0.9
    cuts.AddCut("00010113", "0d200079227300008250404000", "0152103500000000", "204b6600263202263710"); // min pT no cut
    cuts.AddCut("00010113", "0d200049227300008250404000", "0152103500000000", "204b6600263202263710"); // min pT 75 MeV
    cuts.AddCut("00010113", "0d200019227300008250404000", "0152103500000000", "204b6600263202263710"); // min pT 100 MeV
    cuts.AddCut("00010113", "0d200088227300008250404000", "0152103500000000", "204b6600263202263710"); // TPC cluster 35%
    cuts.AddCut("00010113", "0d200086227300008250404000", "0152103500000000", "204b6600263202263710"); // TPC cluster 70%
  } else if (trainConfig == 313){
    cuts.AddCut("00010113", "0d200089327300008250404000", "0152103500000000", "204b6600263202263710"); // edEdx -4,5
    cuts.AddCut("00010113", "0d200089627300008250404000", "0152103500000000", "204b6600263202263710"); // edEdx -2.5,4
    cuts.AddCut("00010113", "0d200089257300008250404000", "0152103500000000", "204b6600263202263710"); // pidEdx 2,-10
    cuts.AddCut("00010113", "0d200089217300008250404000", "0152103500000000", "204b6600263202263710"); // pidEdx 0,-10
    cuts.AddCut("00010113", "0d200089226300008250404000", "0152103500000000", "204b6600263202263710"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCut("00010113", "0d200089227600008250404000", "0152103500000000", "204b6600263202263710"); // pion nsig max mom 2.00 GeV/c
  } else if (trainConfig == 314){
    cuts.AddCut("00010113", "0d200089227300002250404000", "0152103500000000", "204b6600263202263710"); // qT max 0.06 2D
    cuts.AddCut("00010113", "0d200089227300009250404000", "0152103500000000", "204b6600263202263710"); // qT max 0.03 2D
    cuts.AddCut("00010113", "0d200089227300008260404000", "0152103500000000", "204b6600263202263710"); // Psi pair 0.05 2D, chi2 30.
    cuts.AddCut("00010113", "0d200089227300008280404000", "0152103500000000", "204b6600263202263710"); // Psi pair 0.2  2D, chi2 30.
    cuts.AddCut("00010113", "0d200089227300008850404000", "0152103500000000", "204b6600263202263710"); // chi2 20. psi pair 0.1 2D
    cuts.AddCut("00010113", "0d200089227300008150404000", "0152103500000000", "204b6600263202263710"); // chi2 50. psi pair 0.1 2D
  } else if (trainConfig == 315){
    cuts.AddCut("00010113", "0d200089227300008254404000", "0152103500000000", "204b6600263202263710"); // Photon Asymmetry Cut
    cuts.AddCut("00010113", "0d200089227300008250604000", "0152103500000000", "204b6600263202263710"); // CosPA 0.9
    cuts.AddCut("00010113", "0d200089227300008250004000", "0152103500000000", "204b6600263202263710"); // no CosPA
    cuts.AddCut("00010113", "0d200089227300008250400000", "0152103500000000", "204b6600263202263710"); // no double counting
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152101500000000", "204b6600263202263710"); // meson alpha pt dep
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152107500000000", "204b6600263202263710"); // meson alpha < 0.85
  }  else if (trainConfig == 325) {//Primary cut, +with no mass cut,no psi cut, +with psi pair cut cuts.
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263002263010");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// 5TeV 2017 and 13TeV Normal Magnetic Field////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }  else if (trainConfig == 400) {////Primary cut (No Standard yet!), there are the standard cut from Lucia at 5TeV plus Pedro cuts for 7 Tev or 5.02 TeV pPb on electrons.
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223710");

    ///////////////////////////////Study on electrons cuts////////////////////////////////////

  }  else if (trainConfig == 401) {//Variation on the dE/dx standar -4 to 5
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "205b6600263202223710");//Primary, change de/dx TPCsigma -3,5
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "206b6600263202223710");//Primary, change de/dx TPCsigma -4,4
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "208b6600263202223710");//Primary, change de/dx TPCsigma -2,3.5
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "209b6600263202223710");//Primary, change de/dx TPCsigma -3,4
  }  else if (trainConfig == 402) {  //Variation on Mass_{e^+e^-}, Primary, change <1.0pT & 0.015 mass, >1.0pT & 0.035 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223a10");//Primary, change <1.0pT & 0.02 mass, >1.0pT & 0.03 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223b10");//Primary, change <1.0pT & 0.027 mass, >1.0pT & 0.057 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223c10");//Primary, change < 0.02 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223d10");//Primary, change < 0.04 mass
  }  else if (trainConfig == 403) {  //Variation on Mass_{e^+e^-}, Primary, change <1.0pT & 0.015 mass, >1.0pT & 0.035 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223e10");//Primary, change < 0.06 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223f10");//Primary, change < 0.08 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223g10");//Primary, change < 0.15 mass
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000", "204b6600263202223h10");//Primary, change < 0.2 mass
  }  else if (trainConfig == 404) {//Primary cut, +with no mass cut,no psi cut, +with psi pair cut cuts, Pt e+e- Primary 0.9.
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263002273010");//No Psi pair cut
    cuts.AddCut("00010113", "0d200089227300008250404000", "0152103500000000", "204b6600263202273010");//0.60 Psi pair cut

  //////////////// systematic variations for 5 TeV 2017 (Lucia)/////////////////////////////

  } else if (trainConfig == 410){
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000", "204b6600263202223710"); // eta < 0.9
    cuts.AddCut("00010113", "0c200009227300008250404000", "0152103500000000", "204b6600263202223710"); // eta < 0.85
    cuts.AddCut("00010113", "0d200009247300008250404000", "0152103500000000", "204b6600263202223710"); // pidEdx 3, 1
    cuts.AddCut("00010113", "0d200009227100008250404000", "0152103500000000", "204b6600263202223710"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 411){
    cuts.AddCut("00010113", "0d200009227300008860404000", "0152103500000000", "204b6600263202223710"); // variation chi2 20 psi pair 0.2 2D
    cuts.AddCut("00010113", "0d200009227300002252404000", "0152103500000000", "204b6600263202223710"); // variation qT max 0.06 2D, asym var 1
    cuts.AddCut("00010113", "0d200009227300002254404000", "0152103500000000", "204b6600263202223710"); // variation qT max 0.06 2D, asym vat 2 pt dep
    cuts.AddCut("00010113", "0d200009227300002256404000", "0152103500000000", "204b6600263202223710"); // variation qT max 0.06 2D, asym var 3
  } else if (trainConfig == 412){
    cuts.AddCut("00010113", "0d200009a27300008250904120", "0152103500000000", "204b6600263202223710"); //cosPA, 0.99 eta 0.9
    cuts.AddCut("00010113", "0d200079227300008250404000", "0152103500000000", "204b6600263202223710"); // min pT no cut
    cuts.AddCut("00010113", "0d200049227300008250404000", "0152103500000000", "204b6600263202223710"); // min pT 75 MeV
    cuts.AddCut("00010113", "0d200019227300008250404000", "0152103500000000", "204b6600263202223710"); // min pT 100 MeV
    cuts.AddCut("00010113", "0d200008227300008250404000", "0152103500000000", "204b6600263202223710"); // TPC cluster 35%
    cuts.AddCut("00010113", "0d200006227300008250404000", "0152103500000000", "204b6600263202223710"); // TPC cluster 70%
  } else if (trainConfig == 413){
    cuts.AddCut("00010113", "0d200009327300008250404000", "0152103500000000", "204b6600263202223710"); // edEdx -4,5
    cuts.AddCut("00010113", "0d200009627300008250404000", "0152103500000000", "204b6600263202223710"); // edEdx -2.5,4
    cuts.AddCut("00010113", "0d200009257300008250404000", "0152103500000000", "204b6600263202223710"); // pidEdx 2,-10
    cuts.AddCut("00010113", "0d200009217300008250404000", "0152103500000000", "204b6600263202223710"); // pidEdx 0,-10
    cuts.AddCut("00010113", "0d200009226300008250404000", "0152103500000000", "204b6600263202223710"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCut("00010113", "0d200009227600008250404000", "0152103500000000", "204b6600263202223710"); // pion nsig max mom 2.00 GeV/c
  } else if (trainConfig == 414){
    cuts.AddCut("00010113", "0d200009227300002250404000", "0152103500000000", "204b6600263202223710"); // qT max 0.06 2D
    cuts.AddCut("00010113", "0d200009227300009250404000", "0152103500000000", "204b6600263202223710"); // qT max 0.03 2D
    cuts.AddCut("00010113", "0d200009227300008260404000", "0152103500000000", "204b6600263202223710"); // Psi pair 0.05 2D, chi2 30.
    cuts.AddCut("00010113", "0d200009227300008280404000", "0152103500000000", "204b6600263202223710"); // Psi pair 0.2  2D, chi2 30.
    cuts.AddCut("00010113", "0d200009227300008850404000", "0152103500000000", "204b6600263202223710"); // chi2 20. psi pair 0.1 2D
    cuts.AddCut("00010113", "0d200009227300008150404000", "0152103500000000", "204b6600263202223710"); // chi2 50. psi pair 0.1 2D
  } else if (trainConfig == 415){
    cuts.AddCut("00010113", "0d200009227300008254404000", "0152103500000000", "204b6600263202223710"); // Photon Asymmetry Cut
    cuts.AddCut("00010113", "0d200009227300008250604000", "0152103500000000", "204b6600263202223710"); // CosPA 0.9
    cuts.AddCut("00010113", "0d200009227300008250004000", "0152103500000000", "204b6600263202223710"); // no CosPA
    cuts.AddCut("00010113", "0d200009227300008250400000", "0152103500000000", "204b6600263202223710"); // no double counting
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152101500000000", "204b6600263202223710"); // meson alpha pt dep
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152107500000000", "204b6600263202223710"); // meson alpha < 0.85

    //  6XX for lowB,    65X  lowB and MBW
  } else if (trainConfig == 606) {  // low B   eta<0.8
    cuts.AddCut("00010113", "0d200089266300008850404000", "0263103100900000", "10477400234202361710"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if ( trainConfig ==607 ) { // low B
    cuts.AddCut("00010113", "0da00089266300008850404000", "0263103100900000", "10477400234202361710"); //prim electron Pt >
    cuts.AddCut("00010113", "0db00089266300008850404000", "0263103100900000", "10477400234202361710"); //prim electron Pt
    cuts.AddCut("00010113", "0dc00089266300008850404000", "0263103100900000", "10477400234202361710");


    // to be used with weights
  } else if (trainConfig == 656) {  // low B   eta<0.8
    cuts.AddCut("00010113", "0d200089266300008850404000", "0263103100900000", "10477400234202361710"); //prim electron Pt > 0.05 GeV & sec. electrons pt > 0.02 GeV
  } else if ( trainConfig ==657 ) { // low B
    cuts.AddCut("00010113", "0da00089266300008850404000", "0263103100900000", "10477400234202361710"); //prim electron Pt >
    cuts.AddCut("00010113", "0db00089266300008850404000", "0263103100900000", "10477400234202361710"); //prim electron Pt
    cuts.AddCut("00010113", "0dc00089266300008850404000", "0263103100900000", "10477400234202361710");


    //  7XX for NomB,    75X  nomB and MBW
  } else if (trainConfig == 706) {   // Nominal eta<0.8 for primary and secondary electrons
    cuts.AddCut("00010113", "0d200009266300008854404000", "0263103100900000", "20475400254202321710");
  } else if (trainConfig == 707) {    // R scan
    cuts.AddCut("00010113", "0da00009266300008854404000", "0263103100900000", "20475400254202321710");
    cuts.AddCut("00010113", "0db00009266300008854404000", "0263103100900000", "20475400254202321710");
    cuts.AddCut("00010113", "0dc00009266300008854404000", "0263103100900000", "20475400254202321710");
    // to be used with weights equal to 70X +50
  } else if (trainConfig == 756) {   // Nominal eta<0.8 for primary and secondary electrons
    cuts.AddCut("00010113", "0d200009266300008854404000", "0263103100900000", "20475400254202321710");
  } else if (trainConfig == 757) {    // R scan
    cuts.AddCut("00010113", "0da00009266300008854404000", "0263103100900000", "20475400254202321710");
    cuts.AddCut("00010113", "0db00009266300008854404000", "0263103100900000", "20475400254202321710");
    cuts.AddCut("00010113", "0dc00009266300008854404000", "0263103100900000", "20475400254202321710");

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
  AliConversionPhotonCuts **analysisCuts   = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  ElecCutList->SetOwner(kTRUE);
  AliDalitzElectronCuts **analysisElecCuts   = new AliDalitzElectronCuts*[numberOfCuts];
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    TString cutName( Form("%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetElecCut(i)).Data(),(cuts.GetMesonCut(i)).Data() ) );

    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    if (enableMatBudWeightsPi0 > 0){
        if (isMC > 0){
            if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,filenameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;}
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");

    analysisElecCuts[i] = new AliDalitzElectronCuts();
    analysisElecCuts[i]->InitializeCutsFromCutString((cuts.GetElecCut(i)).Data());
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
  //task->SetDoMesonAnalysis(kTRUE);
  if (enableQAMesonTask) task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
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
