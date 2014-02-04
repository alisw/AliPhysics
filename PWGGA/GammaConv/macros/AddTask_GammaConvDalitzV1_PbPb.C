void AddTask_GammaConvDalitzV1_PbPb(   Int_t trainConfig = 1,
                                      Bool_t isMC   = kFALSE, //run MC 
                                      Bool_t enableQAMesonTask = kFALSE, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                      Bool_t enableDoMesonChic = kFALSE, // enable additional Chic analysis
                                      TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                      Bool_t doWeighting = kFALSE,  //enable Weighting
                                      TString cutnumberAODBranch = "0000000060084001001500000"
                                 ) {



  cout<<"Entro -1"<<endl;

 // ================= Load Librariers =================================
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");  
  gSystem->Load("libPWGGAGammaConv.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libTENDERSupplies.so");


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
			   
  TString ConvCutnumber = "1080000000084001001500000000";   //Online  V0 finder
  TString ElecCuts      = "9000620000000200000";            //Electron Cuts
                           //903162000550020210
                           //900054000000020000


  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
     AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");

     fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
     fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
     fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

     if (!mgr) {
        Error("AddTask_V0ReaderV1", "No analysis manager found.");
        return;
     }

     // Set AnalysisCut Number
     AliConversionCuts *fCuts=NULL;
     if( ConvCutnumber !=""){
        fCuts= new AliConversionCuts(ConvCutnumber.Data(),ConvCutnumber.Data());
        fCuts->SetPreSelectionCutFlag(kTRUE);
        if(fCuts->InitializeCutsFromCutString(ConvCutnumber.Data())){
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



  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 3;

  TString *ConvCutarray    = new TString[numberOfCuts];
  TString *ElecCutarray    = new TString[numberOfCuts];
  TString *MesonCutarray   = new TString[numberOfCuts];

   if ( trainConfig == 1 ) {

       ConvCutarray[0]  = "1240001042092971007200000000"; MesonCutarray[0] = "01522045009000"; ElecCutarray[0]  = "9051620025510252170"; //PbPb 20-40% kAny
       ConvCutarray[1]  = "1460001042092971007200000000"; MesonCutarray[1] = "01522045009000"; ElecCutarray[1]  = "9051620025510252170"; //PbPb 40-60% kAny
       ConvCutarray[2]  = "1680001042092971007200000000"; MesonCutarray[2] = "01522045009000"; ElecCutarray[2]  = "9051620025510252170"; //PbPb 60-80% kAny

    } else if ( trainConfig == 2 ) {

       ConvCutarray[0]  = "5240001042092971003220000000"; MesonCutarray[0] = "01522085009000"; ElecCutarray[0]  = "9051620025510252170"; //PbPb 20-40% kAny Alpha cut 0.6
       ConvCutarray[1]  = "5460001042092971001200000000"; MesonCutarray[1] = "01522065009000"; ElecCutarray[1]  = "9051620025510252170"; //PbPb 40-60% kAny Alpha cut 0.8      
       ConvCutarray[2]  = "5680001042092971001200000000"; MesonCutarray[2] = "01522075009000"; ElecCutarray[2]  = "9051620025510252170"; //PbPb 60-80% kAny Alpha cut 0.85
      
    } else if ( trainConfig == 3 ) {

       ConvCutarray[0]  = "5240001042092971003220000000"; MesonCutarray[0] = "01522085009000"; ElecCutarray[0]  = "9051620025510252171"; //PbPb 20-40% kAny Alpha cut 0.6
       ConvCutarray[1]  = "5460001042092971001200000000"; MesonCutarray[1] = "01522065009000"; ElecCutarray[1]  = "9051620025510252171"; //PbPb 40-60% kAny Alpha cut 0.8      
       ConvCutarray[2]  = "5680001042092971001200000000"; MesonCutarray[2] = "01522075009000"; ElecCutarray[2]  = "9051620025510252171"; //PbPb 60-80% kAny Alpha cut 0.85

    } else if ( trainConfig == 4 ) {

       ConvCutarray[0]  = "5240002032092971003220000000"; MesonCutarray[0] = "01523015009000"; ElecCutarray[0]  = "9051620025510252171"; //PbPb 20-40% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent ( 0.7, 1.2)
       ConvCutarray[1]  = "5460002032092971001200000000"; MesonCutarray[1] = "01523015009000"; ElecCutarray[1]  = "9051620025510252171"; //PbPb 40-60% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent ( 0.7, 1.2)
       ConvCutarray[2]  = "5680002032092971001200000000"; MesonCutarray[2] = "01523025009000"; ElecCutarray[2]  = "9051620025510252171"; //PbPb 60-80% kAny Gamma |Eta| < 0.80  only added signals alpha cut Pt dependent ( 0.80, 1.2)

    } else if ( trainConfig == 5 ) {

       ConvCutarray[0]  = "5240001032092971003220000000"; MesonCutarray[0] = "01523015009000"; ElecCutarray[0]  = "9051620025510252171"; //PbPb 20-40% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.7, 1.2)
       ConvCutarray[1]  = "5460001032092971001200000000"; MesonCutarray[1] = "01523015009000"; ElecCutarray[1]  = "9051620025510252171"; //PbPb 40-60% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.7, 1.2)  
       ConvCutarray[2]  = "5680001032092971001200000000"; MesonCutarray[2] = "01523025009000"; ElecCutarray[2]  = "9051620025510252171"; //PbPb 60-80% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.8, 1.2)
    } else if ( trainConfig == 6 ) {

       ConvCutarray[0]  = "5240002032092971003220000000"; MesonCutarray[0] = "01523095009000"; ElecCutarray[0]  = "9051620025510252171"; //PbPb 20-40% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent( 0.65, 1.2)
       ConvCutarray[1]  = "5460002032092971001200000000"; MesonCutarray[1] = "01523095009000"; ElecCutarray[1]  = "9051620025510252171"; //PbPb 40-60% kAny Gamma |Eta| < 0.65  only added signals alpha cut Pt dependent( 0.65, 1.2)
       ConvCutarray[2]  = "5680002032092971001200000000"; MesonCutarray[2] = "01523025009000"; ElecCutarray[2]  = "9051620025510252171"; //PbPb 60-80% kAny Gamma |Eta| < 0.80  only added signals alpha cut Pt dependent( 0.80, 1.2)

    } else if ( trainConfig == 7 ) {

       ConvCutarray[0]  = "5240001032092971003220000000"; MesonCutarray[0] = "01523095009000"; ElecCutarray[0]  = "9051620025510252171"; //PbPb 20-40% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.65, 1.2)
       ConvCutarray[1]  = "5460001032092971001200000000"; MesonCutarray[1] = "01523095009000"; ElecCutarray[1]  = "9051620025510252171"; //PbPb 40-60% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.65, 1.2)
       ConvCutarray[2]  = "5680001032092971001200000000"; MesonCutarray[2] = "01523025009000"; ElecCutarray[2]  = "9051620025510252171"; //PbPb 60-80% kAny Gamma |Eta| < 0.65 alpha cut Pt dependent ( 0.80, 1.2)
    }



  TList *ConvCutList  = new TList();
  TList *MesonCutList = new TList();
  TList *ElecCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);

 //TObjString *Header3 = new TObjString("eta_2");
 //HeaderList->Add(Header3);

  ConvCutList->SetOwner(kTRUE);
  AliConversionCuts **analysisCuts             = new AliConversionCuts*[numberOfCuts];

  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];

  ElecCutList->SetOwner(kTRUE);
  AliDalitzElectronCuts **analysisElecCuts     = new AliDalitzElectronCuts*[numberOfCuts];



  for(Int_t i = 0; i<numberOfCuts; i++){
     analysisCuts[i] = new AliConversionCuts();
     if( ! analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data()) ) {
           cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
           return 0;
     } else {
        ConvCutList->Add(analysisCuts[i]);
        analysisCuts[i]->SetFillCutHistograms("",kFALSE);
        if( trainConfig == 1){
            if (i == 0 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
            if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
            if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        } else if ( trainConfig == 2 || trainConfig == 3 || trainConfig == 5 || trainConfig == 7 ) {
            if (i == 0 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
            if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
            if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        } else if ( trainConfig == 4 || trainConfig == 6 ) {
            if (i == 0 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
            if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
            if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        }
     }

     analysisMesonCuts[i] = new AliConversionMesonCuts();
     if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
           cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
           return 0;
     } else {
       MesonCutList->Add(analysisMesonCuts[i]);
       analysisMesonCuts[i]->SetFillCutHistograms("");
     }

     TString cutName( Form("%s_%s_%s",ConvCutarray[i].Data(),ElecCutarray[i].Data(),MesonCutarray[i].Data() ) );
     analysisElecCuts[i] = new AliDalitzElectronCuts();
     if( !analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data())) {
           cout<< "ERROR:  analysisElecCuts [ " <<i<<" ] "<<endl;
           return 0;
     }  else { 
        ElecCutList->Add(analysisElecCuts[i]);
        analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
     }
     analysisCuts[i]->SetAcceptedHeader(HeaderList);

  }


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
