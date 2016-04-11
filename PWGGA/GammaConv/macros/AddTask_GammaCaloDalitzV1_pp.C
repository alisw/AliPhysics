void AddTask_GammaCaloDalitzV1_pp(  Int_t trainConfig = 1,  //change different set of cuts
                  Bool_t isMC   = kFALSE, //run MC
                  Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                  Int_t enableQAPhotonTask = 0, // enable additional QA task
                  TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                  TString cutnumberAODBranch = "8000000060084000001500000", // cutnumber for AOD branch
                  Bool_t  enableExtendedMatching = kFALSE, 							//enable or disable extended matching histograms for conversion electrons <-> cluster
                  TString periodname = "LHC12f1x", 						// period name
                  Bool_t 	doWeighting = kFALSE,							// enables weighting
                  Bool_t 	enableV0findingEffi = kFALSE
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
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGGAGammaConv");

  Int_t isHeavyIon = 0;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }
  
  Printf("here \n");
  
  //=========  Set Cutnumber for V0Reader ================================
            //06000078400100001500000000
  TString cutnumberPhoton 	= "00000008400100001500000000";
  TString cutnumberEvent 		= "00000003";
  TString cutnumberElectron   = "90005400000002000000";            //Electron Cuts
  
  Bool_t doEtaShift = kFALSE;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
    if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    
        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

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
        fEventCuts->DoEtaShift(doEtaShift);
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

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

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

    if( cutnumberElectron!=""){
      fElecCuts= new AliDalitzElectronCuts(cutnumberElectron.Data(),cutnumberElectron.Data());
      if(fElecCuts->InitializeCutsFromCutString(cutnumberElectron.Data())){
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
  
  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaCaloDalitzV1 *task=NULL;
  task= new AliAnalysisTaskGammaCaloDalitzV1(Form("GammaCaloDalitz_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
    task->SetV0ReaderName(V0ReaderName);

  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 2;
  
  TString *eventCutArray   	= new TString[numberOfCuts];
  TString *photonCutArray  	= new TString[numberOfCuts];
  TString *clusterCutArray 	= new TString[numberOfCuts];
  TString *electronCutArray	= new TString[numberOfCuts];
  TString *mesonCutArray   	= new TString[numberOfCuts];
  
  //************************************************ EMCAL clusters **********************************************************
  if ( trainConfig == 1){ // min energy = 0.3 GeV/c																
      eventCutArray[0] = "00003113"; photonCutArray[0] = "00200009327002008250400000"; clusterCutArray[0] = "1111100053032230000"; electronCutArray[0] = "90475400233102621710"; mesonCutArray[0] = "0263103100000000"; //standart cut, kINT7
      eventCutArray[1] = "00051113"; photonCutArray[1] = "00200009327002008250400000"; clusterCutArray[1] = "1111100053032230000"; electronCutArray[1] = "90475400233102621710"; mesonCutArray[1] = "0263103100000000"; //standard cut, kEMC7
  } else if ( trainConfig == 2){     /****PHOS******/
      eventCutArray[0] = "00003113"; photonCutArray[0] = "00200009327002008250400000"; clusterCutArray[0] = "2444400048033200000"; electronCutArray[0] = "90475400233102621710"; mesonCutArray[0] = "0263103100000000"; //standart cut, kINT7
      eventCutArray[1] = "00061113"; photonCutArray[1] = "00200009327002008250400000"; clusterCutArray[1] = "2444400048033200000"; electronCutArray[1] = "90475400233102621710"; mesonCutArray[1] = "0263103100000000"; //standard cut, kPHI7
  } else if ( trainConfig == 3 ){  /*****EMCAL*****/
      eventCutArray[0] = "00003113"; photonCutArray[0] = "00200009327002008250400000"; clusterCutArray[0] = "1111200053032230000"; electronCutArray[0] = "90475400233102621710"; mesonCutArray[0] = "0263103100000000"; //standart cut, kINT7 with TRD
      eventCutArray[1] = "00003113"; photonCutArray[1] = "00200009327002008250400000"; clusterCutArray[1] = "1113100053032230000"; electronCutArray[1] = "90475400233102621710"; mesonCutArray[1] = "0263103100000000"; //standard cut, kINT7 No TRD
  }  else {
    Error(Form("GammaConvDalitzCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  TList *EventCutList     = new TList();
  TList *ConvCutList      = new TList();
  TList *ClusterCutList   = new TList();
  TList *ElectronCutList  = new TList();
  TList *MesonCutList     = new TList();

  TList *HeaderList = new TList();
  
  if (periodname.Contains("LHC12i3")){
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if ( periodname.CompareTo("LHC14e2b") == 0 ){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }	
  
  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (periodname.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (periodname.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  } 	
  if (periodname.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (periodname.Contains("LHC12f1a")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";
  } else if (periodname.Contains("LHC12f1b")){
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";
  } else if (periodname.Contains("LHC14e2a")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";			
  } else if (periodname.Contains("LHC14e2b")){	
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";				
  } else if (periodname.Contains("LHC14e2c")){		
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";					
  }	
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  
  ElectronCutList->SetOwner(kTRUE);
  AliDalitzElectronCuts   **analysisElectronCuts = new AliDalitzElectronCuts*[numberOfCuts];
    
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();  
    
        // definition of weighting input
    TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());
    Bool_t fAddedSignal = eventCutArray[i].EndsWith("2");
    
    TString mcInputNamePi0 = "";
    TString mcInputNameEta = "";
    
    if (fAddedSignal && (periodname.Contains("LHC12i3") || periodname.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }	
    
    if (doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());    
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
        
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
  
    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->InitializeCutsFromCutString(clusterCutArray[i].Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatching(enableExtendedMatching);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisElectronCuts[i] = new AliDalitzElectronCuts();
    if( !analysisElectronCuts[i]->InitializeCutsFromCutString(electronCutArray[i].Data())) {
      cout<< "ERROR:  analysisElectronCuts [ " <<i<<" ] "<<endl;
      return 0;
    }
    ElectronCutList->Add(analysisElectronCuts[i]);
    analysisElectronCuts[i]->SetFillCutHistograms("",kFALSE,electronCutArray[i].Data()); 
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetElectronCutList(numberOfCuts,ElectronCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
        //task->SetUseTHnSparse(isUsingTHnSparse);
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCaloDalitz_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloDalitz_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
