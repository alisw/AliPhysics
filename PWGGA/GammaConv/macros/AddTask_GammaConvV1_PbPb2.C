void AddTask_GammaConvV1_PbPb2(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Bool_t doWeighting = kFALSE,  //enable Weighting
                              TString cutnumberAODBranch = "1000000060084000001500000" 
                           ) {

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
   
   //=========  Set Cutnumber for V0Reader ================================
   TString cutnumber = "1000000000084001001500000000"; 
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
      if(cutnumber!=""){
         fCuts= new AliConversionCuts(cutnumber.Data(),cutnumber.Data());
         fCuts->SetPreSelectionCutFlag(kTRUE);
         if(fCuts->InitializeCutsFromCutString(cutnumber.Data())){
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
   //========= Add task to the ANALYSIS manager =====
   //================================================
   AliAnalysisTaskGammaConvV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
   task->SetIsHeavyIon(1);
   task->SetIsMC(isMC);
   // Cut Numbers to use in Analysis
   Int_t numberOfCuts = 1;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if (trainConfig == 1){ 
      cutarray[ 0] = "6010001042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; 
   } else if (trainConfig == 2) { 
      cutarray[ 0] = "6120001042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; 
   } else if (trainConfig == 3) { 
      cutarray[ 0] = "5010001042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; 
   } else if (trainConfig == 4) { 
      cutarray[ 0] = "5020001042092970023220000000"; mesonCutArray[ 0] = "01522045009000";    
   } else if (trainConfig == 5) { 
      cutarray[ 0] = "5120001042092970023220000000"; mesonCutArray[ 0] = "01522045009000";    
   } else if (trainConfig == 6) { 
      cutarray[ 0] = "5240001042092970023220000000"; mesonCutArray[ 0] = "01522045009000";       
   } else if (trainConfig == 7) {    
      cutarray[ 0] = "5460001042092970023220000000"; mesonCutArray[ 0] = "01522065009000"; 
   } else if (trainConfig == 8) {    
      cutarray[ 0] = "5480001042092970023220000000"; mesonCutArray[ 0] = "01522065009000";    
   } else if (trainConfig == 9) {    
      cutarray[ 0] = "5450001042092970023220000000"; mesonCutArray[ 0] = "01522065009000"; 
   } else if (trainConfig == 10) { 
      cutarray[ 0] = "5560001042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
   } else if (trainConfig == 11) { 
      cutarray[ 0] = "5680001042092970023220000000"; mesonCutArray[ 0] = "01522065009000";    
   } else if (trainConfig == 12) { 
      cutarray[ 0] = "5670001042092970023220000000"; mesonCutArray[ 0] = "01522065009000"; 
   } else if (trainConfig == 13) { 
      cutarray[ 0] = "5780001042092970023220000000"; mesonCutArray[ 0] = "01522065009000"; 
   } else if (trainConfig == 14) { 
      cutarray[ 0] = "4690001042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
   } else if (trainConfig == 15) { 
      cutarray[ 0] = "5890001042092970023220000000"; mesonCutArray[ 0] = "01522065009000";    
   } else  if (trainConfig == 16){ 
      cutarray[ 0] = "6010002032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; 
   } else if (trainConfig == 17) { 
      cutarray[ 0] = "6120001032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; 
   } else if (trainConfig == 18) { 
      cutarray[ 0] = "5010001032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; 
   } else if (trainConfig == 19) { 
      cutarray[ 0] = "5020001032092970023220000000"; mesonCutArray[ 0] = "01523045009000";    
   } else if (trainConfig == 20) { 
      cutarray[ 0] = "5120001032092970023220000000"; mesonCutArray[ 0] = "01523045009000";    
   } else if (trainConfig == 21) { 
      cutarray[ 0] = "5240001032092970023220000000"; mesonCutArray[ 0] = "01523045009000";       
   } else if (trainConfig == 22) {    
      cutarray[ 0] = "5460001032092970023220000000"; mesonCutArray[ 0] = "01523065009000"; 
   } else if (trainConfig == 23) {    
      cutarray[ 0] = "5480001032092970023220000000"; mesonCutArray[ 0] = "01523065009000";    
   } else if (trainConfig == 24) {    
      cutarray[ 0] = "5450001032092970023220000000"; mesonCutArray[ 0] = "01523065009000"; 
   } else if (trainConfig == 25) { 
      cutarray[ 0] = "5560001032092970023220000000"; mesonCutArray[ 0] = "01523065009000";
   } else if (trainConfig == 26) { 
      cutarray[ 0] = "5680001032092970023220000000"; mesonCutArray[ 0] = "01523065009000";    
   } else if (trainConfig == 27) { 
      cutarray[ 0] = "5670001032092970023220000000"; mesonCutArray[ 0] = "01523065009000"; 
   } else if (trainConfig == 28) { 
      cutarray[ 0] = "5780001032092970023220000000"; mesonCutArray[ 0] = "01523065009000"; 
   } else if (trainConfig == 29) { 
      cutarray[ 0] = "4690001032092970023220000000"; mesonCutArray[ 0] = "01523065009000";
   } else if (trainConfig == 30) { 
      cutarray[ 0] = "5890001032092970023220000000"; mesonCutArray[ 0] = "01523065009000";    
   } else if (trainConfig == 31) { 
      cutarray[ 0] = "5080001002092970023220000000"; mesonCutArray[ 0] = "01525065009000";    
   } else if (trainConfig == 32) { 
      cutarray[ 0] = "5080002002092970023220000000"; mesonCutArray[ 0] = "01525065009000";    
   } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
   }

   TList *ConvCutList = new TList();
   TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   TObjString *Header1 = new TObjString("BOX");
   HeaderList->Add(Header1);
//    TObjString *Header3 = new TObjString("eta_2");
//    HeaderList->Add(Header3);
   
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

   for(Int_t i = 0; i<numberOfCuts; i++){
      analysisCuts[i] = new AliConversionCuts();
      analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
      ConvCutList->Add(analysisCuts[i]);

      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
      analysisMesonCuts[i] = new AliConversionMesonCuts();
      analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
      analysisCuts[i]->SetAcceptedHeader(HeaderList);
   }

   task->SetConversionCutList(numberOfCuts,ConvCutList);
   task->SetMesonCutList(numberOfCuts,MesonCutList);
   task->SetMoveParticleAccordingToVertex(kTRUE);
   task->SetDoMesonAnalysis(kTRUE);
   task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
   task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return;

}
