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
   Int_t numberOfCuts = 4;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if (trainConfig == 1){ // Standard neutral pion cuts
      cutarray[ 0] = "301000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "312000104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "101000104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "112000104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 2) { // Standard neutral pion cuts
      cutarray[ 0] = "124000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "146000104209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "168000104209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "108000104209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 0-80%
   } else if (trainConfig == 3) { // Standard direct photon cuts
      cutarray[ 0] = "102000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "104000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-40%
      cutarray[ 2] = "148000104209237002370000000"; mesonCutArray[ 2] = "01522065009000"; // 40-80%
      cutarray[ 3] = "108000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-80%
   } else if (trainConfig == 4) { // Standard direct photon cuts
      cutarray[ 0] = "101000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-10%
      cutarray[ 1] = "112000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 10-20%
      cutarray[ 2] = "124000104209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 20-40%
      cutarray[ 3] = "109000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-90%
   } else if (trainConfig == 5) { // Standard neutral pion cuts only added signals
      cutarray[ 0] = "301000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "312000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "101000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "112000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 6) { // Standard neutral pion cuts only added signals
      cutarray[ 0] = "124000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "146000204209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "168000204209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "108000204209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 0-80%
   } else if (trainConfig == 7) { // Standard direct photon cuts
      cutarray[ 0] = "102000304209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "104000304209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-40%
      cutarray[ 2] = "148000304209237002370000000"; mesonCutArray[ 2] = "01522065009000"; // 40-80%
      cutarray[ 3] = "108000304209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-80%
   } else if (trainConfig == 8) { // Standard direct photon cuts
      cutarray[ 0] = "101000304209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-10%
      cutarray[ 1] = "112000304209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 10-20%
      cutarray[ 2] = "124000304209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 20-40%
      cutarray[ 3] = "109000304209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-90%
   } else if(trainConfig == 9){ // Different trigger conditions
      cutarray[ 0] = "301000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "301400104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-5% kMB
      cutarray[ 2] = "301500104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-5% kSemiCentral
      cutarray[ 3] = "301600104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 0-5% kCentral
   } else if (trainConfig == 10) { // Standard direct photon cuts
      cutarray[ 0] = "102000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "104000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-40%
      cutarray[ 2] = "148000204209237002370000000"; mesonCutArray[ 2] = "01522065009000"; // 40-80%
      cutarray[ 3] = "108000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-80%
   } else if (trainConfig == 11) { // Standard direct photon cuts
      cutarray[ 0] = "101000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-10%
      cutarray[ 1] = "112000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 10-20%
      cutarray[ 2] = "124000204209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 20-40%
      cutarray[ 3] = "109000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-90%   
   } else if (trainConfig == 12){ // Standard neutral pion cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "601000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 13) { // Standard neutral pion cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "524000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "508000104209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 0-80%
   } else if (trainConfig == 14) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "504000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-40%
      cutarray[ 2] = "548000104209237002370000000"; mesonCutArray[ 2] = "01522065009000"; // 40-80%
      cutarray[ 3] = "508000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-80%
   } else if (trainConfig == 15) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "501000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-10%
      cutarray[ 1] = "512000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 10-20%
      cutarray[ 2] = "524000104209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 20-40%
      cutarray[ 3] = "509000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-90%
   } else if (trainConfig == 16) { // Standard neutral pion cuts only added signals, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "601000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 17) { // Standard neutral pion cuts only added signals, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "524000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "508000204209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 0-80%
   } else if (trainConfig == 18) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000304209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "504000304209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-40%
      cutarray[ 2] = "548000304209237002370000000"; mesonCutArray[ 2] = "01522065009000"; // 40-80%
      cutarray[ 3] = "508000304209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-80%
   } else if (trainConfig == 19) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "501000304209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-10%
      cutarray[ 1] = "512000304209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 10-20%
      cutarray[ 2] = "524000304209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 20-40%
      cutarray[ 3] = "509000304209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-90%
   } else if (trainConfig == 20) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "504000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-40%
      cutarray[ 2] = "548000204209237002370000000"; mesonCutArray[ 2] = "01522065009000"; // 40-80%
      cutarray[ 3] = "508000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-80%
   } else if (trainConfig == 21) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "501000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-10%
      cutarray[ 1] = "512000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 10-20%
      cutarray[ 2] = "524000204209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 20-40%
      cutarray[ 3] = "509000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 0-90%      
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
      analysisCuts[i] = new AliConversionCuts();
      if (trainConfig == 1 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M"); 
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 2){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i==2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");
      } else if (trainConfig == 3){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0040V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");
      } else if (trainConfig == 4){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 5 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M"); 
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 6){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i==2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      } else if (trainConfig == 9){
         if (doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kMB, "AliVEvent::kMB" );
         if (i == 2) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kSemiCentral,"AliVEvent::kSemiCentral" );
         if (i == 3) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kCentral,"AliVEvent::kCentral" );
         
      } else if (trainConfig == 10){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0040V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");
      } else if (trainConfig == 11){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 12 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 13){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i==2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");
      } else if (trainConfig == 14){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0040V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");
      } else if (trainConfig == 15){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 16 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 17){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i==2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");   
      } else if (trainConfig == 20){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0040V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0080V0M");
      } else if (trainConfig == 21){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      }
      
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
