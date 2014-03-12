void AddTask_GammaConvV1_PbPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Bool_t doWeighting = kFALSE,  //enable Weighting
                              TString cutnumberAODBranch = "1000000060084000001500000",
                              TString periodName = "LHC13d2"
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
   Int_t numberOfCuts = 5;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if (trainConfig == 1){ // Standard cuts
      cutarray[ 0] = "6010001042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "6120001042092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "5010001042092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "5120001042092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "5020001042092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 2) { // Standard cuts
      cutarray[ 0] = "5240001042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "5460001042092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "5680001042092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "5480001042092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
      cutarray[ 4] = "5490001042092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 3) { // Standard cuts only added signals
      cutarray[ 0] = "6010002042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "6120002042092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "5010002042092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "5120002042092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "5020002042092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 4) { // Standard cuts only added signals
      cutarray[ 0] = "5240002042092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "5460002042092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "5680002042092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "5480002042092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
      cutarray[ 4] = "5490002042092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 5){ // R-minCut 7.5 cm
      cutarray[ 0] = "6010001049092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "6120001049092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "5010001049092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "5120001049092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "5020001049092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 6) { // R-minCut 7.5 cm
      cutarray[ 0] = "5240001049092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "5460001049092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "5680001049092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "5480001049092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
      cutarray[ 4] = "5490001049092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 7) {// R-minCut 7.5 cm
      cutarray[ 0] = "6010002049092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "6120002049092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "5010002049092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "5120002049092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "5020002049092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 8) { // R-minCut 7.5 cm
      cutarray[ 0] = "5240002040092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "5460002049092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "5680002049092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "5480002049092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
      cutarray[ 4] = "5490002049092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%   
   } else if (trainConfig == 9){ // R-minCut 12.5 cm
      cutarray[ 0] = "6010001048092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "6120001048092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "5010001048092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "5120001048092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "5020001048092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 10) { // R-minCut 12.5 cm
      cutarray[ 0] = "5240001048092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "5460001048092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "5680001048092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "5480001048092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
      cutarray[ 4] = "5490001048092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%
   } else if (trainConfig == 11) {// R-minCut 12.5 cm
      cutarray[ 0] = "6010002048092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "6120002048092970023220000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "5010002048092970023220000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "5120002048092970023220000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
      cutarray[ 4] = "5020002048092970023220000000"; mesonCutArray[ 4] = "01522045009000"; // 0-20%
   } else if (trainConfig == 12) { // R-minCut 12.5 cm
      cutarray[ 0] = "5240002040092970023220000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "5460002048092970023220000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "5680002048092970023220000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "5480002048092970023220000000"; mesonCutArray[ 3] = "01522065009000"; // 20-40% 
      cutarray[ 4] = "5490002048092970023220000000"; mesonCutArray[ 4] = "01522065009000"; // 40-90%      
   } else  if (trainConfig == 13){ // eta 0.65 (new standard), y = 0.6 (new Standard)
      cutarray[ 0] = "6010001032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120001032092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010001032092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120001032092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020001032092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 14) {  // eta 0.65 (new standard), y = 0.6 (new Standard)
      cutarray[ 0] = "5240001032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460001032092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680001032092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480001032092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
      cutarray[ 4] = "5490001032092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else if (trainConfig == 15) { // eta 0.65 (new standard), y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "6010002032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120002032092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010002032092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120002032092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020002032092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 16) { // eta 0.65 (new standard), y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "5240002032092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460002032092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680002032092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480002032092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
      cutarray[ 4] = "5490002032092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else  if (trainConfig == 17){ // eta 0.6, y = 0.6 (new Standard)
      cutarray[ 0] = "6010001012092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120001012092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010001012092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120001012092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020001012092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 18) {  // eta 0.6, y = 0.6 (new Standard)
      cutarray[ 0] = "5240001012092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460001012092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680001012092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480001012092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
      cutarray[ 4] = "5490001012092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else if (trainConfig == 19) { // eta 0.6, y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "6010002012092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120002012092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010002012092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120002012092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020002012092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 20) { // eta 0.6, y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "5240002012092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460002012092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680002012092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480002012092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
      cutarray[ 4] = "5490002012092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else  if (trainConfig == 21){ // eta 0.7, y = 0.6 (new Standard)
      cutarray[ 0] = "6010001072092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120001072092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010001072092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120001072092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020001072092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 22) {  // eta 0.7, y = 0.6 (new Standard)
      cutarray[ 0] = "5240001072092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460001072092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680001072092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480001072092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
      cutarray[ 4] = "5490001072092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else if (trainConfig == 23) { // eta 0.7, y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "6010002072092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120002072092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010002072092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120002072092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020002072092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 24) { // eta 0.7, y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "5240002072092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460002072092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680002072092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480002072092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
      cutarray[ 4] = "5490002072092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else  if (trainConfig == 25){ // eta 0.5, y = 0.6 (new Standard)
      cutarray[ 0] = "6010001052092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120001052092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010001052092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120001052092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020001052092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 26) {  // eta 0.5, y = 0.6 (new Standard)
      cutarray[ 0] = "5240001052092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460001052092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680001052092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480001052092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 40-80%
      cutarray[ 4] = "5490001052092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else if (trainConfig == 27) { // eta 0.5, y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "6010002052092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 0-5%
      cutarray[ 1] = "6120002052092970023220000000"; mesonCutArray[ 1] = "01523045009000"; // 5-10%
      cutarray[ 2] = "5010002052092970023220000000"; mesonCutArray[ 2] = "01523045009000"; // 0-10%
      cutarray[ 3] = "5120002052092970023220000000"; mesonCutArray[ 3] = "01523045009000"; // 10-20%
      cutarray[ 4] = "5020002052092970023220000000"; mesonCutArray[ 4] = "01523045009000"; // 0-20%
   } else if (trainConfig == 28) { // eta 0.5, y = 0.6 (new Standard) cuts only added signals
      cutarray[ 0] = "5240002052092970023220000000"; mesonCutArray[ 0] = "01523045009000"; // 20-40%
      cutarray[ 1] = "5460002052092970023220000000"; mesonCutArray[ 1] = "01523065009000"; // 40-60%
      cutarray[ 2] = "5680002052092970023220000000"; mesonCutArray[ 2] = "01523065009000"; // 60-80%
      cutarray[ 3] = "5480002052092970023220000000"; mesonCutArray[ 3] = "01523065009000"; // 20-40% 
      cutarray[ 4] = "5490002052092970023220000000"; mesonCutArray[ 4] = "01523065009000"; // 40-90%
   } else  if (trainConfig == 29){ // eta 0.65 (new standard), y = 0.6 (new Standard) pt dependent alpha
      cutarray[ 0] = "6010001032092970023220000000"; mesonCutArray[ 0] = "01523015009000"; // 0-5%
      cutarray[ 1] = "6120001032092970023220000000"; mesonCutArray[ 1] = "01523015009000"; // 5-10%
      cutarray[ 2] = "5010001032092970023220000000"; mesonCutArray[ 2] = "01523015009000"; // 0-10%
      cutarray[ 3] = "5120001032092970023220000000"; mesonCutArray[ 3] = "01523015009000"; // 10-20%
      cutarray[ 4] = "5020001032092970023220000000"; mesonCutArray[ 4] = "01523015009000"; // 0-20%
   } else if (trainConfig == 30) {  // eta 0.65 (new standard), y = 0.6 (new Standard) pt dependent alpha
      cutarray[ 0] = "5240001032092970023220000000"; mesonCutArray[ 0] = "01523015009000"; // 20-40%
      cutarray[ 1] = "5460001032092970023220000000"; mesonCutArray[ 1] = "01523025009000"; // 40-60%
      cutarray[ 2] = "5680001032092970023220000000"; mesonCutArray[ 2] = "01523025009000"; // 60-80%
      cutarray[ 3] = "5480001032092970023220000000"; mesonCutArray[ 3] = "01523025009000"; // 40-80%
      cutarray[ 4] = "5490001032092970023220000000"; mesonCutArray[ 4] = "01523025009000"; // 40-90%
   } else if (trainConfig == 31) { // eta 0.65 (new standard), y = 0.6 (new Standard) cuts only added signals, pt dependent alpha
      cutarray[ 0] = "6010002032092970023220000000"; mesonCutArray[ 0] = "01523015009000"; // 0-5%
      cutarray[ 1] = "6120002032092970023220000000"; mesonCutArray[ 1] = "01523015009000"; // 5-10%
      cutarray[ 2] = "5010002032092970023220000000"; mesonCutArray[ 2] = "01523015009000"; // 0-10%
      cutarray[ 3] = "5120002032092970023220000000"; mesonCutArray[ 3] = "01523015009000"; // 10-20%
      cutarray[ 4] = "5020002032092970023220000000"; mesonCutArray[ 4] = "01523015009000"; // 0-20%
   } else if (trainConfig == 32) { // eta 0.65 (new standard), y = 0.6 (new Standard) cuts only added signals, pt dependent alpha
      cutarray[ 0] = "5240002032092970023220000000"; mesonCutArray[ 0] = "01523015009000"; // 20-40%
      cutarray[ 1] = "5460002032092970023220000000"; mesonCutArray[ 1] = "01523025009000"; // 40-60%
      cutarray[ 2] = "5680002032092970023220000000"; mesonCutArray[ 2] = "01523025009000"; // 60-80%
      cutarray[ 3] = "5480002032092970023220000000"; mesonCutArray[ 3] = "01523025009000"; // 20-40% 
      cutarray[ 4] = "5490002032092970023220000000"; mesonCutArray[ 4] = "01523025009000"; // 40-90%
   } else if (trainConfig == 33){ // Standard cuts, eta 0.9, only to be run on data
      cutarray[ 0] = "6010001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 0-5%
      cutarray[ 1] = "6120001002092970023220000000"; mesonCutArray[ 1] = "01525045009000"; // 5-10%
      cutarray[ 2] = "5010001002092970023220000000"; mesonCutArray[ 2] = "01525045009000"; // 0-10%
      cutarray[ 3] = "5120001002092970023220000000"; mesonCutArray[ 3] = "01525045009000"; // 10-20%
      cutarray[ 4] = "5020001002092970023220000000"; mesonCutArray[ 4] = "01525045009000"; // 0-20%
   } else if (trainConfig == 34) { // Standard cuts, eta 0.9, only to be run on data
      cutarray[ 0] = "5240001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 20-40%
      cutarray[ 1] = "5460001002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680001002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480001002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350001002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; // 40-90%   
   } else if (trainConfig == 35){ // Standard cuts, eta 1.4, only to be run on data
      cutarray[ 0] = "6010001022092970023220000000"; mesonCutArray[ 0] = "01520045009000"; // 0-5%
      cutarray[ 1] = "6120001022092970023220000000"; mesonCutArray[ 1] = "01520045009000"; // 5-10%
      cutarray[ 2] = "5010001022092970023220000000"; mesonCutArray[ 2] = "01520045009000"; // 0-10%
      cutarray[ 3] = "5120001022092970023220000000"; mesonCutArray[ 3] = "01520045009000"; // 10-20%
      cutarray[ 4] = "5020001022092970023220000000"; mesonCutArray[ 4] = "01520045009000"; // 0-20%
   } else if (trainConfig == 36) { // Standard cuts, eta 1.4, only to be run on data
      cutarray[ 0] = "5240001022092970023220000000"; mesonCutArray[ 0] = "01520045009000"; // 20-40%
      cutarray[ 1] = "5460001022092970023220000000"; mesonCutArray[ 1] = "01520065009000"; // 40-60%
      cutarray[ 2] = "5680001022092970023220000000"; mesonCutArray[ 2] = "01520065009000"; // 60-80%
      cutarray[ 3] = "5480001022092970023220000000"; mesonCutArray[ 3] = "01520065009000"; // 40-80%
      cutarray[ 4] = "5350001022092970023220000000"; mesonCutArray[ 4] = "01520065009000"; // 40-90%   
   } else if (trainConfig == 37){ // Standard cuts, eta 0.9, only to be run on data : kMB
      cutarray[ 0] = "6014001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 0-5%
      cutarray[ 1] = "6124001002092970023220000000"; mesonCutArray[ 1] = "01525045009000"; // 5-10%
      cutarray[ 2] = "5014001002092970023220000000"; mesonCutArray[ 2] = "01525045009000"; // 0-10%
      cutarray[ 3] = "5124001002092970023220000000"; mesonCutArray[ 3] = "01525045009000"; // 10-20%
      cutarray[ 4] = "5024001002092970023220000000"; mesonCutArray[ 4] = "01525045009000"; // 0-20%
   } else if (trainConfig == 38) { // Standard cuts, eta 0.9, only to be run on data : kMB
      cutarray[ 0] = "5244001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 20-40%
      cutarray[ 1] = "5464001002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5684001002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5484001002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5354001002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; // 40-90%   
   } else if (trainConfig == 39){ // Standard cuts, eta 0.9, only to be run on data :kCentral
      cutarray[ 0] = "6015001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 0-5%
      cutarray[ 1] = "6125001002092970023220000000"; mesonCutArray[ 1] = "01525045009000"; // 5-10%
      cutarray[ 2] = "5015001002092970023220000000"; mesonCutArray[ 2] = "01525045009000"; // 0-10%
      cutarray[ 3] = "5125001002092970023220000000"; mesonCutArray[ 3] = "01525045009000"; // 10-20%
      cutarray[ 4] = "5025001002092970023220000000"; mesonCutArray[ 4] = "01525045009000"; // 0-20%
   } else if (trainConfig == 40) { // Standard cuts, eta 0.9, only to be run on data : kCentral
      cutarray[ 0] = "5245001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 20-40%
      cutarray[ 1] = "5465001002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5685001002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5485001002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5355001002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; // 40-90%   
   } else if (trainConfig == 41){ // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
      cutarray[ 0] = "6016001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; // 0-5%
      cutarray[ 1] = "6126001002092970023220000000"; mesonCutArray[ 1] = "01525045009000"; // 5-10%
      cutarray[ 2] = "5016001002092970023220000000"; mesonCutArray[ 2] = "01525045009000"; // 0-10%
      cutarray[ 3] = "5126001002092970023220000000"; mesonCutArray[ 3] = "01525045009000"; // 10-20%
      cutarray[ 4] = "5026001002092970023220000000"; mesonCutArray[ 4] = "01525045009000"; // 0-20%
   } else if (trainConfig == 42) { // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
      cutarray[ 0] = "5236001002092970023220000000"; mesonCutArray[ 0] = "01525065009000"; 
      cutarray[ 1] = "5346001002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; 
      cutarray[ 2] = "5456001002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; 
      cutarray[ 3] = "5566001002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; 
      cutarray[ 4] = "5676001002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; 
   } else if (trainConfig == 43){ // Standard cuts, eta 0.9, only to be run on data
      cutarray[ 0] = "5230001002092970023220000000"; mesonCutArray[ 0] = "01525045009000"; 
      cutarray[ 1] = "5340001002092970023220000000"; mesonCutArray[ 1] = "01525065009000"; 
      cutarray[ 2] = "5450001002092970023220000000"; mesonCutArray[ 2] = "01525065009000"; 
      cutarray[ 3] = "5560001002092970023220000000"; mesonCutArray[ 3] = "01525065009000"; 
      cutarray[ 4] = "5250001002092970023220000000"; mesonCutArray[ 4] = "01525065009000"; 
   } else if ( trainConfig == 44){ // qt elipse cut 0.05
      cutarray[ 0] = "6010001002092970028220000000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
      cutarray[ 1] = "6120001002092970028220000000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
      cutarray[ 2] = "5010001002092970028220000000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
      cutarray[ 3] = "5120001002092970028220000000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
      cutarray[ 4] = "5020001002092970028220000000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
   } else if ( trainConfig == 45) { // qt elipse cut 0.05
      cutarray[ 0] = "5240001002092970028220000000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
      cutarray[ 1] = "5460001002092970028220000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680001002092970028220000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480001002092970028220000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350001002092970028220000000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
   } else if ( trainConfig == 46){ // qt elipse cut 0.05
      cutarray[ 0] = "5230001002092970028220000000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
      cutarray[ 1] = "5340001002092970028220000000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
      cutarray[ 2] = "5450001002092970028220000000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
      cutarray[ 3] = "5560001002092970028220000000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
      cutarray[ 4] = "5250001002092970028220000000"; mesonCutArray[ 4] = "01525065009000"; // 60-70% 
   } else if ( trainConfig == 47){ // cos(theta_point) cut 0.85
      cutarray[ 0] = "6010001002092970023220400000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
      cutarray[ 1] = "6120001002092970023220400000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
      cutarray[ 2] = "5010001002092970023220400000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
      cutarray[ 3] = "5120001002092970023220400000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
      cutarray[ 4] = "5020001002092970023220400000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
   } else if ( trainConfig == 48) { // cos(theta_point) cut 0.85
      cutarray[ 0] = "5240001002092970023220400000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
      cutarray[ 1] = "5460001002092970023220400000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680001002092970023220400000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480001002092970023220400000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350001002092970023220400000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
   } else if ( trainConfig == 49){ // cos(theta_point) cut 0.85
      cutarray[ 0] = "5230001002092970023220400000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
      cutarray[ 1] = "5340001002092970023220400000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
      cutarray[ 2] = "5450001002092970023220400000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
      cutarray[ 3] = "5560001002092970023220400000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
      cutarray[ 4] = "5250001002092970023220400000"; mesonCutArray[ 4] = "01525065009000"; // 60-70% 
   } else if ( trainConfig == 50){ // psi pair 2D 0.05
      cutarray[ 0] = "6010001002092970023260000000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
      cutarray[ 1] = "6120001002092970023260000000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
      cutarray[ 2] = "5010001002092970023260000000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
      cutarray[ 3] = "5120001002092970023260000000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
      cutarray[ 4] = "5020001002092970023260000000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
   } else if ( trainConfig == 51) { // psi pair 2D 0.05
      cutarray[ 0] = "5240001002092970023260000000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
      cutarray[ 1] = "5460001002092970023260000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680001002092970023260000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480001002092970023260000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350001002092970023260000000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
   } else if ( trainConfig == 52){ // psi pair 2D 0.05
      cutarray[ 0] = "5230001002092970023260000000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
      cutarray[ 1] = "5340001002092970023260000000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
      cutarray[ 2] = "5450001002092970023260000000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
      cutarray[ 3] = "5560001002092970023260000000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
      cutarray[ 4] = "5250001002092970023260000000"; mesonCutArray[ 4] = "01525065009000"; // 60-70%       
   } else if ( trainConfig == 53){ // psi pair 2D 0.1
      cutarray[ 0] = "6010001002092970023250000000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
      cutarray[ 1] = "6120001002092970023250000000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
      cutarray[ 2] = "5010001002092970023250000000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
      cutarray[ 3] = "5120001002092970023250000000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
      cutarray[ 4] = "5020001002092970023250000000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
   } else if ( trainConfig == 54) { // psi pair 2D 0.1
      cutarray[ 0] = "5240001002092970023250000000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
      cutarray[ 1] = "5460001002092970023250000000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680001002092970023250000000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480001002092970023250000000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350001002092970023250000000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
   } else if ( trainConfig == 55){ // psi pair 2D 0.1
      cutarray[ 0] = "5230001002092970023250000000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
      cutarray[ 1] = "5340001002092970023250000000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
      cutarray[ 2] = "5450001002092970023250000000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
      cutarray[ 3] = "5560001002092970023250000000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
      cutarray[ 4] = "5250001002092970023250000000"; mesonCutArray[ 4] = "01525065009000"; // 60-70%          
   } else if ( trainConfig == 56){ // cleaner cuts
      cutarray[ 0] = "6010001002092970028250400000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
      cutarray[ 1] = "6120001002092970028250400000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
      cutarray[ 2] = "5010001002092970028250400000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
      cutarray[ 3] = "5120001002092970028250400000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
      cutarray[ 4] = "5020001002092970028250400000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
   } else if ( trainConfig == 57) { // cleaner cuts
      cutarray[ 0] = "5240001002092970028250400000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
      cutarray[ 1] = "5460001002092970028250400000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680001002092970028250400000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480001002092970028250400000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350001002092970028250400000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
   } else if ( trainConfig == 58){ // cleaner cuts
      cutarray[ 0] = "5230001002092970028250400000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
      cutarray[ 1] = "5340001002092970028250400000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
      cutarray[ 2] = "5450001002092970028250400000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
      cutarray[ 3] = "5560001002092970028250400000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
      cutarray[ 4] = "5250001002092970028250400000"; mesonCutArray[ 4] = "01525065009000"; // 60-70%                
   } else if ( trainConfig == 59){ // cleaner cuts
      cutarray[ 0] = "6010002002092970028250400000"; mesonCutArray[ 0] = "01525065009000"; // 0-5%
      cutarray[ 1] = "6120002002092970028250400000"; mesonCutArray[ 1] = "01525065009000"; // 5-10%
      cutarray[ 2] = "5010002002092970028250400000"; mesonCutArray[ 2] = "01525065009000"; // 0-10%
      cutarray[ 3] = "5120002002092970028250400000"; mesonCutArray[ 3] = "01525065009000"; // 10-20%
      cutarray[ 4] = "5020002002092970028250400000"; mesonCutArray[ 4] = "01525065009000"; // 0-20%
   } else if ( trainConfig == 60) { // cleaner cuts
      cutarray[ 0] = "5240002002092970028250400000"; mesonCutArray[ 0] = "01525065009000"; // 20-40%
      cutarray[ 1] = "5460002002092970028250400000"; mesonCutArray[ 1] = "01525065009000"; // 40-60%
      cutarray[ 2] = "5680002002092970028250400000"; mesonCutArray[ 2] = "01525065009000"; // 60-80%
      cutarray[ 3] = "5480002002092970028250400000"; mesonCutArray[ 3] = "01525065009000"; // 40-80%
      cutarray[ 4] = "5350002002092970028250400000"; mesonCutArray[ 4] = "01525065009000"; // 30-50%  
   } else if ( trainConfig == 61){ // cleaner cuts
      cutarray[ 0] = "5230002002092970028250400000"; mesonCutArray[ 0] = "01525065009000"; // 20-30% 
      cutarray[ 1] = "5340002002092970028250400000"; mesonCutArray[ 1] = "01525065009000"; // 30-40% 
      cutarray[ 2] = "5450002002092970028250400000"; mesonCutArray[ 2] = "01525065009000"; // 40-50% 
      cutarray[ 3] = "5560002002092970028250400000"; mesonCutArray[ 3] = "01525065009000"; // 50-60% 
      cutarray[ 4] = "5670002002092970028250400000"; mesonCutArray[ 4] = "01525065009000"; // 60-70%                   
   } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
   }

   TList *ConvCutList = new TList();
   TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   if (periodName.CompareTo("LHC13d2")==0){
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
   //    TObjString *Header3 = new TObjString("eta_2");
   //    HeaderList->Add(Header3);

   } else if (periodName.CompareTo("LHC12a17x_fix")==0){
      TObjString *Header1 = new TObjString("PARAM");
      HeaderList->Add(Header1);
   } else if (periodName.CompareTo("LHC14a1x")==0){
      TObjString *Header1 = new TObjString("BOX");
      HeaderList->Add(Header1);
   }
   
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

   for(Int_t i = 0; i<numberOfCuts; i++){
      analysisCuts[i] = new AliConversionCuts();
      if (trainConfig == 1 ||trainConfig == 5 || trainConfig == 9 || trainConfig == 13 || trainConfig == 17 || trainConfig == 21 || trainConfig == 25 || trainConfig == 29 ){ //|| trainConfig == 33 || trainConfig == 37 || trainConfig == 41 
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
         if (i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
      } else if (trainConfig == 2 ||trainConfig == 6 || trainConfig == 10 || trainConfig == 14 || trainConfig == 18 || trainConfig == 22 || trainConfig == 26 || trainConfig == 30  ){ //|| trainConfig == 34 || trainConfig == 38 || trainConfig == 42
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      } else if (trainConfig == 3 ||trainConfig == 7 || trainConfig == 11 || trainConfig == 15 || trainConfig == 19 || trainConfig == 23 || trainConfig == 27 || trainConfig == 31 ){ //|| trainConfig == 35 || trainConfig == 39 || trainConfig == 43 
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
          if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
      } else if (trainConfig == 4 ||trainConfig == 8 || trainConfig == 12 || trainConfig == 16 || trainConfig == 20 || trainConfig == 24 || trainConfig == 28 || trainConfig == 32 ){ //|| trainConfig == 36 || trainConfig == 40 || trainConfig == 44 
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

      if (trainConfig == 37 || trainConfig == 38){
         analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kMB, "AliVEvent::kMB" );
      }
      if (trainConfig == 39 || trainConfig == 40){   
         analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kCentral,"AliVEvent::kCentral" );
      }   
      if (trainConfig == 41 || trainConfig == 42){   
         analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kSemiCentral,"AliVEvent::kSemiCentral" );
      }   
      
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
