void AddTask_GammaConvV1_PbPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Bool_t enableQAMesonTask = kFALSE, //enable QA in AliAnalysisTaskGammaConvV1
                              Bool_t enableQAPhotonTask = kFALSE, // enable additional QA task
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
   TString cutnumber = "100000000008400100150000000"; 
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
   Int_t numberOfCuts = 5;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if (trainConfig == 1){ // Standard cuts
      cutarray[ 0] = "601000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "502000104209297002322000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 2) { // Standard cuts
      cutarray[ 0] = "524000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "548000104209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
      cutarray[ 4] = "549000104209297002322000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 3) { // Standard cuts only added signals
      cutarray[ 0] = "601000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "502000204209297002322000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 4) { // Standard cuts only added signals
      cutarray[ 0] = "524000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "548000204209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
      cutarray[ 4] = "549000204209297002322000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 5){ // R-minCut 7.5 cm
      cutarray[ 0] = "601000104909297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104909297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104909297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104909297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "502000104909297002322000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 6) { // R-minCut 7.5 cm
      cutarray[ 0] = "524000104909297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104909297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104909297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "548000104909297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
      cutarray[ 4] = "549000104909297002322000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 7) {// R-minCut 7.5 cm
      cutarray[ 0] = "601000204909297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204909297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204909297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204909297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "502000204909297002322000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 8) { // R-minCut 7.5 cm
      cutarray[ 0] = "524000204009297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204909297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204909297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "548000204909297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
      cutarray[ 4] = "549000204909297002322000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%   
   } else if (trainConfig == 9){ // R-minCut 12.5 cm
      cutarray[ 0] = "601000104809297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104809297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104809297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104809297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "502000104809297002322000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 10) { // R-minCut 12.5 cm
      cutarray[ 0] = "524000104809297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104809297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104809297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "548000104809297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
      cutarray[ 4] = "549000104809297002322000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 11) {// R-minCut 12.5 cm
      cutarray[ 0] = "601000204809297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204809297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204809297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204809297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "502000204809297002322000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 12) { // R-minCut 12.5 cm
      cutarray[ 0] = "524000204009297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204809297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204809297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "548000204809297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
      cutarray[ 4] = "549000204809297002322000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%      
   } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
   }

   TList *ConvCutList = new TList();
   TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   TObjString *Header1 = new TObjString("pi0_1");
   HeaderList->Add(Header1);
//    TObjString *Header3 = new TObjString("eta_2");
//    HeaderList->Add(Header3);
   
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

   for(Int_t i = 0; i<numberOfCuts; i++){
      analysisCuts[i] = new AliConversionCuts();
      if (trainConfig == 1 ||trainConfig == 5 || trainConfig == 9 || trainConfig == 13 || trainConfig == 17 || trainConfig == 21 || trainConfig == 25 || trainConfig == 29 || trainConfig == 33 || trainConfig == 37 || trainConfig == 41 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
         if (i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
      } else if (trainConfig == 2 ||trainConfig == 6 || trainConfig == 10 || trainConfig == 14 || trainConfig == 18 || trainConfig == 22 || trainConfig == 26 || trainConfig == 30 || trainConfig == 34 || trainConfig == 38 || trainConfig == 42 ){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      } else if (trainConfig == 3 ||trainConfig == 7 || trainConfig == 11 || trainConfig == 15 || trainConfig == 19 || trainConfig == 23 || trainConfig == 27 || trainConfig == 31 || trainConfig == 35 || trainConfig == 39 || trainConfig == 43 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
          if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
      } else if (trainConfig == 4 ||trainConfig == 8 || trainConfig == 12 || trainConfig == 16 || trainConfig == 20 || trainConfig == 24 || trainConfig == 28 || trainConfig == 32 || trainConfig == 36 || trainConfig == 40 || trainConfig == 44 ){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      }
//      if (trainConfig == 45 ||trainConfig == 47 || trainConfig == 49 || trainConfig == 51 || trainConfig == 53 || trainConfig == 55  ){
//          if ((i == 0 || i == 1 || i == 2)&& doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
//          if ((i == 3 || i == 4 )&& doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
//      }
//      if (trainConfig == 46 ||trainConfig == 48 || trainConfig == 50 || trainConfig == 52 || trainConfig == 54 || trainConfig == 56  ){
//          if ((i == 0 || i == 1 || i == 2)&& doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
//          if ((i == 3 || i == 4 )&& doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
//      }
//      if (trainConfig == 57  ){
//          if (doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
//      }
//      if (trainConfig == 58){
//          if (doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
//      }
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
   if (enableQAPhotonTask) task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return;

}
