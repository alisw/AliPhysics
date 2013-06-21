void AddTask_GammaConvV1_PbPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Bool_t enableQAMesonTask = kFALSE, //enable QA in AliAnalysisTaskGammaConvV1
                              Bool_t enableQAPhotonTask = kFALSE, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Bool_t doWeighting = kFALSE //enable Weighting
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
   TString cutnumber = "1000000000084001001500000"; 
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
   Int_t numberOfCuts = 4;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if(trainConfig == 1){ // Standard neutral pion cuts
      cutarray[ 0] = "3010001042092970023220000"; mesonCutArray[ 0] = "01522045009"; // 0-5%
      cutarray[ 1] = "3120001042092970023220000"; mesonCutArray[ 1] = "01522045009"; // 5-10%
      cutarray[ 2] = "1010001042092970023220000"; mesonCutArray[ 2] = "01522045009"; // 0-10%
      cutarray[ 3] = "1120001042092970023220000"; mesonCutArray[ 3] = "01522045009"; // 10-20%
   } else if (trainConfig == 2) { // Standard neutral pion cuts
      cutarray[ 0] = "1240001042092970023220000"; mesonCutArray[ 0] = "01522045009"; // 20-40%
      cutarray[ 1] = "1460001042092970023220000"; mesonCutArray[ 1] = "01522065009"; // 40-60%
      cutarray[ 2] = "1680001042092970023220000"; mesonCutArray[ 2] = "01522065009"; // 60-80%
      cutarray[ 3] = "1080001042092970023220000"; mesonCutArray[ 3] = "01522065009"; // 0-80%
   } else if (trainConfig == 3) { // Standard direct photon cuts
      cutarray[ 0] = "1020001042092370023700000"; mesonCutArray[ 0] = "01522045009"; // 0-20%
      cutarray[ 1] = "1040001042092370023700000"; mesonCutArray[ 1] = "01522045009"; // 0-40%
      cutarray[ 2] = "1480001042092370023700000"; mesonCutArray[ 2] = "01522065009"; // 40-80%
      cutarray[ 3] = "1080001042092370023700000"; mesonCutArray[ 3] = "01522045009"; // 0-80%
   } else if (trainConfig == 4) { // Standard direct photon cuts
      cutarray[ 0] = "1010001042092370023700000"; mesonCutArray[ 0] = "01522045009"; // 0-10%
      cutarray[ 1] = "1120001042092370023700000"; mesonCutArray[ 1] = "01522045009"; // 10-20%
      cutarray[ 2] = "1240001042092370023700000"; mesonCutArray[ 2] = "01522045009"; // 20-40%
      cutarray[ 3] = "1090001042092370023700000"; mesonCutArray[ 3] = "01522045009"; // 0-90%
   } else if (trainConfig == 5) { // Standard neutral pion cuts only added signals
      cutarray[ 0] = "3010002042092970023220000"; mesonCutArray[ 0] = "01522045009"; // 0-5%
      cutarray[ 1] = "3120002042092970023220000"; mesonCutArray[ 1] = "01522045009"; // 5-10%
      cutarray[ 2] = "1010002042092970023220000"; mesonCutArray[ 2] = "01522045009"; // 0-10%
      cutarray[ 3] = "1120002042092970023220000"; mesonCutArray[ 3] = "01522045009"; // 10-20%
   } else if (trainConfig == 6) { // Standard neutral pion cuts only added signals
      cutarray[ 0] = "1240002042092970023220000"; mesonCutArray[ 0] = "01522045009"; // 20-40%
      cutarray[ 1] = "1460002042092970023220000"; mesonCutArray[ 1] = "01522065009"; // 40-60%
      cutarray[ 2] = "1680002042092970023220000"; mesonCutArray[ 2] = "01522065009"; // 60-80%
      cutarray[ 3] = "1080002042092970023220000"; mesonCutArray[ 3] = "01522065009"; // 0-80%
   } else if (trainConfig == 7) { // Standard direct photon cuts
      cutarray[ 0] = "1020003042092370023700000"; mesonCutArray[ 0] = "01522045009"; // 0-20%
      cutarray[ 1] = "1040003042092370023700000"; mesonCutArray[ 1] = "01522045009"; // 0-40%
      cutarray[ 2] = "1480003042092370023700000"; mesonCutArray[ 2] = "01522065009"; // 40-80%
      cutarray[ 3] = "1080003042092370023700000"; mesonCutArray[ 3] = "01522045009"; // 0-80%
   } else if (trainConfig == 8) { // Standard direct photon cuts
      cutarray[ 0] = "1010003042092370023700000"; mesonCutArray[ 0] = "01522045009"; // 0-10%
      cutarray[ 1] = "1120003042092370023700000"; mesonCutArray[ 1] = "01522045009"; // 10-20%
      cutarray[ 2] = "1240003042092370023700000"; mesonCutArray[ 2] = "01522045009"; // 20-40%
      cutarray[ 3] = "1090003042092370023700000"; mesonCutArray[ 3] = "01522045009"; // 0-90%
   } else if(trainConfig == 9){ // Different trigger conditions
      cutarray[ 0] = "3010001042092970023220000"; mesonCutArray[ 0] = "01522045009"; // 0-5%
      cutarray[ 1] = "3014001042092970023220000"; mesonCutArray[ 1] = "01522045009"; // 0-5% kMB
      cutarray[ 2] = "3015001042092970023220000"; mesonCutArray[ 2] = "01522045009"; // 0-5% kSemiCentral
      cutarray[ 3] = "3016001042092970023220000"; mesonCutArray[ 3] = "01522045009"; // 0-5% kCentral
      
   } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
   }

   TList *ConvCutList = new TList();
   TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   TObjString *Header1 = new TObjString("pi0_1");
   HeaderList->Add(Header1);
   TObjString *Header3 = new TObjString("eta_2");
   HeaderList->Add(Header3);
   
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];


   for(Int_t i = 0; i<numberOfCuts; i++){
      if (trainConfig == 1 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0005", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0005");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0510", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0510"); 
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0010", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0010");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE,fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_1020", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_1020");
      } else if (trainConfig == 2){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_2040", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_2040");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_4060", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_4060");
        if (i==2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_6080", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_6080");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0080", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0080");
      } else if (trainConfig == 3){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0020", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0020");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0040", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0040");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_4080", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_4080");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0080", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0080");
      } else if (trainConfig == 4){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0010", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0010");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE,fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_1020", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_1020");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_2040", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_2040");
      } else if (trainConfig == 9){
         if (doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, fileNameInputForWeighting, "Pi0_Hijing_PbPb_2760GeV_0005", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0005");
         if (i == 1) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kMB, "AliVEvent::kMB" );
         if (i == 2) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kSemiCentral,"AliVEvent::kSemiCentral" );
         if (i == 3) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kCentral,"AliVEvent::kCentral" );
         
      } 
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
   if (enableQAMesonTask) task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
   if (enableQAMesonTask) task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return;

}
