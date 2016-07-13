void AddTask_GammaConvDalitzV1_pp(  Int_t trainConfig = 1,  //change different set of cuts
                                    Bool_t isMC   = kFALSE, //run MC 
                                    TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                    Int_t   enableMatBudWeightsPi0          = 0,              // 1 = three radial bins, 2 = 10 radial bins
                                    TString filenameMatBudWeights           = "MCInputFileMaterialBudgetWeights.root"
         ) {

  // ================= Load Librariers =================================
  gSystem->Load("libCore");  
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");  
  gSystem->Load("libPWGGAGammaConv");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  
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
  
  
  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 2;
  
  TString *eventCutArray   = new TString[numberOfCuts];
  TString *photonCutArray  = new TString[numberOfCuts];
  TString *MesonCutarray   = new TString[numberOfCuts];
  TString *ElecCutarray    = new TString[numberOfCuts];
  

  if(trainConfig == 1){
    //TOF PID                                           
    eventCutArray[0]="00000113"; photonCutArray[0] = "00200009366302007800000000"; MesonCutarray[0] = "0163103100900000";ElecCutarray[0] = "90478403253102621000";  //TOF[-3,5] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0
    eventCutArray[1]="00000113"; photonCutArray[1] = "00200009366302007800000000"; MesonCutarray[1] = "0163103100900000";ElecCutarray[1] = "90478404253102621000";  //TOF[-2,3] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0     
  } else if (trainConfig == 2) {
    //TOF PID
    eventCutArray[0]="00000113"; photonCutArray[0] = "00200009366302007800000000"; MesonCutarray[0] = "0163103100900000";ElecCutarray[0] = "90478403253102621000";  //TOF[-3,5] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0
    eventCutArray[1]="00000113"; photonCutArray[1] = "00200009366302007800000000"; MesonCutarray[1] = "0163103100900000";ElecCutarray[1] = "90478404253102621000";  //TOF[-2,3] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0     
  } else if (trainConfig == 3) {
    //TOF PID
    eventCutArray[0]="00000113"; photonCutArray[0] = "00200009366302007800000000"; MesonCutarray[0] = "0163103100900000";ElecCutarray[0] = "90478403253102621000";  //TOF[-3,5] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0
    eventCutArray[1]="00000113"; photonCutArray[1] = "00200009366302007800000000"; MesonCutarray[1] = "0163103100900000";ElecCutarray[1] = "90478404253102621000";  //TOF[-2,3] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0     
  } else if (trainConfig == 4) {
    //working electron cutnumber
    eventCutArray[0]="00000113"; photonCutArray[0] = "00200009360300007800004000"; MesonCutarray[0] = "0263103100900000";ElecCutarray[0] = "20475400253202221710";  //New Standard Only MB, standard pp7Tev cut dalitz
  } else {
    Error(Form("GammaConvDalitzV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

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
    TString cutName( Form("%s_%s_%s_%s",eventCutArray[i].Data(),photonCutArray[i].Data(),ElecCutarray[i].Data(),MesonCutarray[i].Data() ) );
    
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
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
    analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
  
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    
    analysisElecCuts[i] = new AliDalitzElectronCuts();
    analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data());
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
