
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
