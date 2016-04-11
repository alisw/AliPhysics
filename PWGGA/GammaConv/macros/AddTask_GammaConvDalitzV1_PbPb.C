void AddTask_GammaConvDalitzV1_PbPb(  Int_t   trainConfig               = 1,
                                      Bool_t  isMC                      = kFALSE, //run MC 
                                      Bool_t  enableQAMesonTask         = kFALSE, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                      Bool_t  enableDoMesonChic         = kFALSE, // enable additional Chic analysis
                                      TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                      Bool_t  doWeighting               = kFALSE,  //enable Weighting
                                      TString cutnumberAODBranch        = "0000000060084001001500000"
                                 ) {
  cout<<"Entro -1"<<endl;

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
  
  Int_t isHeavyIon = 1;
  cout<<"Entro 0"<<endl;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvDalitzV1_PbPb_%i",trainConfig), "No analysis manager found.");
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
          
  TString cutnumberPhoton = "00000008400100001500000000";
  TString cutnumberEvent  = "10000003";
  TString ElecCuts        = "90006200000002000000";            //Electron Cuts

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

  //================================================
  //========= Add Electron Selector ================


  if( !(AliDalitzElectronSelector*)mgr->GetTask("ElectronSelector") ){

    AliDalitzElectronSelector *fElectronSelector = new AliDalitzElectronSelector("ElectronSelector");

    // Set AnalysisCut Number
    AliDalitzElectronCuts *fElecCuts=0;

    if( ElecCuts!=""){
      fElecCuts= new AliDalitzElectronCuts(ElecCuts.Data(),ElecCuts.Data());
      if(fElecCuts->InitializeCutsFromCutString(ElecCuts.Data())){
      fElectronSelector->SetDalitzElectronCuts(fElecCuts);
      fElecCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    fElectronSelector->Init();
    mgr->AddTask(fElectronSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    //connect input V0Reader
    mgr->ConnectInput (fElectronSelector,0,cinput1);
  }

  cout<<"Entro"<<endl;
  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaConvDalitzV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvDalitzV1(Form("GammaConvDalitzV1_%i",trainConfig));
  task->SetIsHeavyIon(1);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);

  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 3;

  TString *eventCutArray     = new TString[numberOfCuts];
  TString *photonCutArray    = new TString[numberOfCuts];
  TString *ElecCutarray      = new TString[numberOfCuts];
  TString *MesonCutarray     = new TString[numberOfCuts];

  if ( trainConfig == 1 ) {
    eventCutArray[0]="12400013"; photonCutArray[0]  = "04200009297100007200000000"; MesonCutarray[0] = "0152204500900000"; ElecCutarray[0]  = "90516200255102521700"; //PbPb 20-40% kAny
    eventCutArray[1]="14600013"; photonCutArray[1]  = "04200009297100007200000000"; MesonCutarray[1] = "0152204500900000"; ElecCutarray[1]  = "90516200255102521700"; //PbPb 40-60% kAny
    eventCutArray[2]="16800013"; photonCutArray[2]  = "04200009297100007200000000"; MesonCutarray[2] = "0152204500900000"; ElecCutarray[2]  = "90516200255102521700"; //PbPb 60-80% kAny
  } else if ( trainConfig == 2 ) {
    eventCutArray[0]="52400013"; photonCutArray[0]  = "04200009297100003220000000"; MesonCutarray[0] = "0152208500900000"; ElecCutarray[0]  = "90516200255102521700"; //PbPb 20-40% kAny Alpha cut 0.6
    eventCutArray[1]="54600013"; photonCutArray[1]  = "04200009297100001200000000"; MesonCutarray[1] = "0152206500900000"; ElecCutarray[1]  = "90516200255102521700"; //PbPb 40-60% kAny Alpha cut 0.8      
    eventCutArray[2]="56800013"; photonCutArray[2]  = "04200009297100001200000000"; MesonCutarray[2] = "0152207500900000"; ElecCutarray[2]  = "90516200255102521700"; //PbPb 60-80% kAny Alpha cut 0.85
  } else if ( trainConfig == 3 ) {
    eventCutArray[0]="52400013"; photonCutArray[0]  = "04200009297100003220000000"; MesonCutarray[0] = "0152208500900000"; ElecCutarray[0]  = "90516200255102521710"; //PbPb 20-40% kAny Alpha cut 0.6
    eventCutArray[1]="54600013"; photonCutArray[1]  = "04200009297100001200000000"; MesonCutarray[1] = "0152206500900000"; ElecCutarray[1]  = "90516200255102521710"; //PbPb 40-60% kAny Alpha cut 0.8      
    eventCutArray[2]="56800013"; photonCutArray[2]  = "04200009297100001200000000"; MesonCutarray[2] = "0152207500900000"; ElecCutarray[2]  = "90516200255102521710"; //PbPb 60-80% kAny Alpha cut 0.85
  } else if ( trainConfig == 4 ) {
    eventCutArray[0]="52400023"; photonCutArray[0]  = "03200009297100003220000000"; MesonCutarray[0] = "0152301500900000"; ElecCutarray[0]  = "90516200255102521710"; //PbPb 20-40% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent ( 0.7, 1.2)
    eventCutArray[1]="54600023"; photonCutArray[1]  = "03200009297100001200000000"; MesonCutarray[1] = "0152301500900000"; ElecCutarray[1]  = "90516200255102521710"; //PbPb 40-60% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent ( 0.7, 1.2)
    eventCutArray[2]="56800023"; photonCutArray[2]  = "03200009297100001200000000"; MesonCutarray[2] = "0152302500900000"; ElecCutarray[2]  = "90516200255102521710"; //PbPb 60-80% kAny Gamma |Eta| < 0.80  only added signals alpha cut Pt dependent ( 0.80, 1.2)
  } else if ( trainConfig == 5 ) {
    eventCutArray[0]="52400013"; photonCutArray[0]  = "03200009297100003220000000"; MesonCutarray[0] = "0152301500900000"; ElecCutarray[0]  = "90516200255102521710"; //PbPb 20-40% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.7, 1.2)
    eventCutArray[1]="54600013"; photonCutArray[1]  = "03200009297100001200000000"; MesonCutarray[1] = "0152301500900000"; ElecCutarray[1]  = "90516200255102521710"; //PbPb 40-60% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.7, 1.2)  
    eventCutArray[2]="56800013"; photonCutArray[2]  = "03200009297100001200000000"; MesonCutarray[2] = "0152302500900000"; ElecCutarray[2]  = "90516200255102521710"; //PbPb 60-80% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.8, 1.2)
  } else if ( trainConfig == 6 ) {
    eventCutArray[0]="52400023"; photonCutArray[0]  = "03200009297100003220000000"; MesonCutarray[0] = "0152309500900000"; ElecCutarray[0]  = "90516200255102521710"; //PbPb 20-40% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent( 0.65, 1.2)
    eventCutArray[1]="54600023"; photonCutArray[1]  = "03200009297100001200000000"; MesonCutarray[1] = "0152309500900000"; ElecCutarray[1]  = "90516200255102521710"; //PbPb 40-60% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent( 0.65, 1.2)
    eventCutArray[2]="56800023"; photonCutArray[2]  = "03200009297100001200000000"; MesonCutarray[2] = "0152302500900000"; ElecCutarray[2]  = "90516200255102521710"; //PbPb 60-80% kAny Gamma |Eta| < 0.80  only added signals alpha cut Pt dependent( 0.80, 1.2)
  } else if ( trainConfig == 7 ) {
    eventCutArray[0]="52400013"; photonCutArray[0]  = "03200009297100003220000000"; MesonCutarray[0] = "0152309500900000"; ElecCutarray[0]  = "90516200255102521710"; //PbPb 20-40% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.65, 1.2)
    eventCutArray[1]="54600013"; photonCutArray[1]  = "03200009297100001200000000"; MesonCutarray[1] = "0152309500900000"; ElecCutarray[1]  = "90516200255102521710"; //PbPb 40-60% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.65, 1.2)
    eventCutArray[2]="56800013"; photonCutArray[2]  = "03200009297100001200000000"; MesonCutarray[2] = "0152302500900000"; ElecCutarray[2]  = "90516200255102521710"; //PbPb 60-80% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.80, 1.2)
  }


  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
  TList *MesonCutList = new TList();
  TList *ElecCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);

  //TObjString *Header3 = new TObjString("eta_2");
  //HeaderList->Add(Header3);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts          = new AliConvEventCuts*[numberOfCuts];
  
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts       = new AliConversionPhotonCuts*[numberOfCuts];

  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];

  ElecCutList->SetOwner(kTRUE);
  AliDalitzElectronCuts **analysisElecCuts     = new AliDalitzElectronCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    
    analysisEventCuts[i] = new AliConvEventCuts();    
    if( trainConfig == 1){
      if (i == 0 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
    } else if ( trainConfig == 2 || trainConfig == 3 || trainConfig == 5 || trainConfig == 7 ) {
      if (i == 0 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
    } else if ( trainConfig == 4 || trainConfig == 6 ) {
      if (i == 0 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
    }
    
    if( ! analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data()) ) {
      cout<<"ERROR: analysisEventCuts [" <<i<<"]"<<endl;
      return 0;
    }
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    if( ! analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data()) ) {
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    }      
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    }

    TString cutName( Form("%s_%s_%s_%s",eventCutArray[i].Data(),photonCutArray[i].Data(),ElecCutarray[i].Data(),MesonCutarray[i].Data() ) );
    analysisElecCuts[i] = new AliDalitzElectronCuts();
    if( !analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data())) {
      cout<< "ERROR:  analysisElecCuts [ " <<i<<" ] "<<endl;
      return 0;
    }  else { 
      ElecCutList->Add(analysisElecCuts[i]);
      analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetElectronCutList(ElecCutList);

  task->SetMoveParticleAccordingToVertex(kTRUE);


  if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);
  if(enableDoMesonChic) task->SetDoChicAnalysis(kTRUE);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvDalitzV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1Dalitz_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;
}
