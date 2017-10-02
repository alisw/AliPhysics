
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
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || elecCut.Length()!=20 ) {cout << "ERROR in CutHandlerConvDalitz: Incorrect length of cut string!" << endl; validCuts = false; return;}
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
                                    TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                    Int_t   enableMatBudWeightsPi0          = 0,              // 1 = three radial bins, 2 = 10 radial bins
                                    TString filenameMatBudWeights           = "MCInputFileMaterialBudgetWeights.root"
         ) {
  
  Int_t isHeavyIon = 0;
    
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
  TString cutnumberPhoton = "06000008400100007500000000";
  TString cutnumberEvent  = "00000103"; 
      
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
  
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

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
      fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
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
  } else if (trainConfig == 202) {
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253202221710");
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "10475400253202221710");
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "30475400253202221710");
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253201221710");
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253203221710");
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253202121710");
	cuts.AddCut("00074113", "00200009360300007800004000", "0263103100900000", "20475400253202321710");


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
  //if (enableQAMesonTask) task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
  //if (enableQAMesonTask) task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvDalitzV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
