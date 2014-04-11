void AddTask_GammaConvV1_PbPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Int_t doWeightingInt = 0,  //enable Weighting
                              TString cutnumberAODBranch = "1000000060084000001500000",
                              TString periodName = "LHC13d2",
							  Bool_t doWeighting = kFALSE
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
      cutarray[ 0] = "6010001002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "6120001002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "5010001002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "5120001002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "5020001002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
   } else if ( trainConfig == 57){ // cleaner cuts added signal
      cutarray[ 0] = "6010002002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "6120002002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "5010002002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "5120002002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "5020002002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%   
	} else if ( trainConfig == 58) { // cleaner cuts
      cutarray[ 0] = "5240001002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-40%
      cutarray[ 1] = "5460001002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 40-60%
      cutarray[ 2] = "5680001002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 60-80%
      cutarray[ 3] = "5480001002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 40-80%
      cutarray[ 4] = "5350001002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 30-50%  
	} else if ( trainConfig == 59) { // cleaner cuts added signal
      cutarray[ 0] = "5240002002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-40%
      cutarray[ 1] = "5460002002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 40-60%
      cutarray[ 2] = "5680002002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 60-80%
      cutarray[ 3] = "5480002002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 40-80%
      cutarray[ 4] = "5350002002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 30-50%  		
	} else if ( trainConfig == 60){ // cleaner cuts
      cutarray[ 0] = "5230001002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-30% 
      cutarray[ 1] = "5340001002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 30-40% 
      cutarray[ 2] = "5450001002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 40-50% 
      cutarray[ 3] = "5560001002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 50-60% 
      cutarray[ 4] = "5250001002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 60-70%                
   } else if ( trainConfig == 61){ // cleaner cuts added signal
      cutarray[ 0] = "5230002002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 20-30% 
      cutarray[ 1] = "5340002002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 30-40% 
      cutarray[ 2] = "5450002002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 40-50% 
      cutarray[ 3] = "5560002002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 50-60% 
      cutarray[ 4] = "5670002002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 60-70%                   
   } else if ( trainConfig == 62){ // cleaner cuts
      cutarray[ 0] = "6230001002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "6340001002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "6450001002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "6560001002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "6670001002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
   } else if ( trainConfig == 63){ // cleaner cuts added signal
      cutarray[ 0] = "6230002002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "6340002002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "6450002002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "6560002002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "6670002002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
   } else if ( trainConfig == 64){ // cleaner cuts
      cutarray[ 0] = "6780001002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "6890001002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "5670001002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "5780001002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "5890001002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
   } else if ( trainConfig == 65){ // cleaner cuts added signal
      cutarray[ 0] = "6780002002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "6890002002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "5670002002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "5780002002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "5890002002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 66){ // cleaner cuts
      cutarray[ 0] = "7010001002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "7120001002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "7230001002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "7340001002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "7450001002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 67){ // cleaner cuts added signal
      cutarray[ 0] = "7010002002092970028250400000"; mesonCutArray[ 0] = "01525065000000"; // 0-5%
      cutarray[ 1] = "7120002002092970028250400000"; mesonCutArray[ 1] = "01525065000000"; // 5-10%
      cutarray[ 2] = "7230002002092970028250400000"; mesonCutArray[ 2] = "01525065000000"; // 0-10%
      cutarray[ 3] = "7340002002092970028250400000"; mesonCutArray[ 3] = "01525065000000"; // 10-20%
      cutarray[ 4] = "7450002002092970028250400000"; mesonCutArray[ 4] = "01525065000000"; // 0-20%
	} else if ( trainConfig == 68){ // cleaner cuts
		cutarray[ 0] = "7560001002092970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "7670001002092970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "7780001002092970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "7890001002092970028250400000"; mesonCutArray[3]= "01525065000000"; // 10-20%
		cutarray[ 4] = "7090001002092970028250400000"; mesonCutArray[4]= "01525065000000"; // 0-20%	
   } else if ( trainConfig == 69){ // cleaner cuts added signal
		cutarray[ 0] = "7560002002092970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "7670002002092970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "7780002002092970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "7890002002092970028250400000"; mesonCutArray[3]= "01525065000000"; // 10-20%
		cutarray[ 4] = "7090002002092970028250400000"; mesonCutArray[4]= "01525065000000"; // 0-20%			
	} else if ( trainConfig == 70){ // variation eta  0.65
		cutarray[ 0] = "6010001032092970028250400000"; mesonCutArray[0]= "01523065000000"; // 0-5%
		cutarray[ 1] = "6120001032092970028250400000"; mesonCutArray[1]= "01523065000000"; // 5-10%
		cutarray[ 2] = "5010001032092970028250400000"; mesonCutArray[2]= "01523065000000"; // 0-10%
		cutarray[ 3] = "5240001032092970028250400000"; mesonCutArray[3]= "01523065000000"; // 20-40%
		cutarray[ 4] = "5250001032092970028250400000"; mesonCutArray[4]= "01523065000000"; // 20-50% 
	} else if ( trainConfig == 71){ // variation eta  0.65 added signal
		cutarray[ 0] = "6010002032092970028250400000"; mesonCutArray[0]= "01523065000000"; // 0-5%
		cutarray[ 1] = "6120002032092970028250400000"; mesonCutArray[1]= "01523065000000"; // 5-10%
		cutarray[ 2] = "5010002032092970028250400000"; mesonCutArray[2]= "01523065000000"; // 0-10%
		cutarray[ 3] = "5240002032092970028250400000"; mesonCutArray[3]= "01523065000000"; // 20-40%
		cutarray[ 4] = "5250002032092970028250400000"; mesonCutArray[4]= "01523065000000"; // 20-50% 		
	} else if ( trainConfig == 72){ // variation eta  0.75
		cutarray[ 0] = "6010001072092970028250400000"; mesonCutArray[0]= "01522065000000"; // 0-5%
		cutarray[ 1] = "6120001072092970028250400000"; mesonCutArray[1]= "01522065000000"; // 5-10%
		cutarray[ 2] = "5010001072092970028250400000"; mesonCutArray[2]= "01522065000000"; // 0-10%
		cutarray[ 3] = "5240001072092970028250400000"; mesonCutArray[3]= "01522065000000"; // 20-40%
		cutarray[ 4] = "5250001072092970028250400000"; mesonCutArray[4]= "01522065000000"; // 20-50% 
	} else if ( trainConfig == 73){ // variation eta  0.75 added signal
		cutarray[ 0] = "6010002072092970028250400000"; mesonCutArray[0]= "01522065000000"; // 0-5%
		cutarray[ 1] = "6120002072092970028250400000"; mesonCutArray[1]= "01522065000000"; // 5-10%
		cutarray[ 2] = "5010002072092970028250400000"; mesonCutArray[2]= "01522065000000"; // 0-10%
		cutarray[ 3] = "5240002072092970028250400000"; mesonCutArray[3]= "01522065000000"; // 20-40%
		cutarray[ 4] = "5250002072092970028250400000"; mesonCutArray[4]= "01522065000000"; // 20-50% 
	} else if ( trainConfig == 74){ // single pt 0.075
		cutarray[ 0] = "6010001002492970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002492970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002492970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002492970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002492970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 75){ // single pt 0.075 added signal
		cutarray[ 0] = "6010002002492970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002492970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002492970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002492970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002492970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 76){ // single pt 0.1
		cutarray[ 0] = "6010001002192970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002192970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002192970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002192970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002192970028250400000"; mesonCutArray[4]= "01525065000000"; //20-50%
	} else if ( trainConfig == 77){ // single pt 0.1 added signal
		cutarray[ 0] = "6010002002192970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002192970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002192970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002192970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002192970028250400000"; mesonCutArray[4]= "01525065000000"; //20-50%
	} else if ( trainConfig == 78){ // variation TPC cls 0.7
		cutarray[ 0] = "6010001002062970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002062970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002062970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002062970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002062970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 79){ // variation TPC cls 0.7 added signal
		cutarray[ 0] = "6010001002062970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002062970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002062970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002062970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002062970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 80){ // variation TPC cls 0.35
		cutarray[ 0] = "6010001002082970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002082970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002082970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 0] = "5240001002082970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002082970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 81){ // variation TPC cls 0.35 added signal
		cutarray[ 0] = "6010002002082970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002082970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002082970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 0] = "5240002002082970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002082970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 82){ // variation edEdx  -4,5
		cutarray[ 0] = "6010001002093970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002093970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002093970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002093970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002093970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 83){ // variation edEdx  -4,5  added signal
		cutarray[ 0] = "6010002002093970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002093970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002093970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002093970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002093970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 84){ // variation edEdx  -2.5,4
		cutarray[ 0] = "6010001002096970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002096970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002096970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002096970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002096970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 85){ // variation edEdx  -2.5,4  added signal
		cutarray[ 0] = "6010002002096970028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002096970028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002096970028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002096970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002096970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 86){ //variation pion p dEdx 0.3-5.
		cutarray[ 0] = "6010001002092951028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092951028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092951028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092951028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092951028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 87){ //variation pion p dEdx 0.3-5.  added signal
		cutarray[ 0] = "6010002002092951028250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092951028250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092951028250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092951028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092951028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 88){ // TOF el. PID -3,5
		cutarray[ 0] = "6010001002092970038250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970038250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970038250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970038250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970038250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 89){ // TOF el. PID -3,5  added signal
		cutarray[ 0] = "6010002002092970038250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970038250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970038250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970038250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970038250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%		
	} else if ( trainConfig == 90){ // TOF el. PID -2,3
		cutarray[ 0] = "6010001002092970048250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970048250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970048250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970048250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970048250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 91){ // TOF el. PID -2,3  added signal
		cutarray[ 0] = "6010002002092970048250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970048250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970048250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970048250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970048250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 92){ // qt 0.03
		cutarray[ 0] = "6010001002092970029250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970029250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970029250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970029250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970029250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 93){ // qt 0.03  added signal
		cutarray[ 0] = "6010002002092970029250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970029250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970029250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970029250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970029250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 94){ // qt 0.07 no2D
		cutarray[ 0] = "6010001002092970022250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970022250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970022250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970022250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970022250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 95){ // qt 0.07 no2D  added signal
		cutarray[ 0] = "6010002002092970022250400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970022250400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970022250400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970022250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970022250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 96){ // chi2  50.
		cutarray[ 0] = "6010001002092970028150400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028150400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028150400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028150400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028150400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 97){ // chi2  50.  added signal
		cutarray[ 0] = "6010002002092970028150400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970028150400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970028150400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970028150400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970028150400000"; mesonCutArray[4]= "01525065000000"; // 20-50%		
	} else if ( trainConfig == 98){ // chi2  20.
		cutarray[ 0] = "6010001002092970028850400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028850400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028850400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028850400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028850400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 99){ // chi2  20.  added signal
		cutarray[ 0] = "6010002002092970028850400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970028850400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970028850400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970028850400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970028850400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 100){ // psi pair 0.05
		cutarray[ 0] = "6010001002092970028260400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028260400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028260400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028260400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028260400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 101){ // psi pair 0.05  added signal
		cutarray[ 0] = "6010002002092970028260400000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970028260400000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970028260400000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970028260400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970028260400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 102){ // cosPA -1
		cutarray[ 0] = "6010001002092970028250000000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028250000000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028250000000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028250000000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028250000000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 103){ // cosPA -1  added signal
		cutarray[ 0] = "6010002002092970028250000000"; mesonCutArray[0]= "01525065000000"; // 0-5%
		cutarray[ 1] = "6120002002092970028250000000"; mesonCutArray[1]= "01525065000000"; // 5-10%
		cutarray[ 2] = "5010002002092970028250000000"; mesonCutArray[2]= "01525065000000"; // 0-10%
		cutarray[ 3] = "5240002002092970028250000000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970028250000000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 104){ // variation alpha 0.6&0.8
		cutarray[ 0] = "6010001002092970028250400000"; mesonCutArray[0]= "01525085000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028250400000"; mesonCutArray[1]= "01525085000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028250400000"; mesonCutArray[2]= "01525085000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%
	} else if ( trainConfig == 105){ // variation alpha 0.6&0.8  added signal
		cutarray[ 0] = "6010002002092970028250400000"; mesonCutArray[0]= "01525085000000"; // 0-5%
		cutarray[ 1] = "6120002002092970028250400000"; mesonCutArray[1]= "01525085000000"; // 5-10%
		cutarray[ 2] = "5010002002092970028250400000"; mesonCutArray[2]= "01525085000000"; // 0-10%
		cutarray[ 3] = "5240002002092970028250400000"; mesonCutArray[3]= "01525065000000"; // 20-40%
		cutarray[ 4] = "5250002002092970028250400000"; mesonCutArray[4]= "01525065000000"; // 20-50%	
	} else if ( trainConfig == 106){ // variation alpha 0.65&0.75
		cutarray[ 0] = "6010001002092970028250400000"; mesonCutArray[0]= "01525045000000"; // 0-5%
		cutarray[ 1] = "6120001002092970028250400000"; mesonCutArray[1]= "01525045000000"; // 5-10%
		cutarray[ 2] = "5010001002092970028250400000"; mesonCutArray[2]= "01525045000000"; // 0-10%
		cutarray[ 3] = "5240001002092970028250400000"; mesonCutArray[3]= "01525055000000"; // 20-40%
		cutarray[ 4] = "5250001002092970028250400000"; mesonCutArray[4]= "01525055000000"; // 20-50%
	} else if ( trainConfig == 107){ // variation alpha 0.65&0.75  added signal
		cutarray[ 0] = "6010002002092970028250400000"; mesonCutArray[0]= "01525045000000"; // 0-5%
		cutarray[ 1] = "6120002002092970028250400000"; mesonCutArray[1]= "01525045000000"; // 5-10%
		cutarray[ 2] = "5010002002092970028250400000"; mesonCutArray[2]= "01525045000000"; // 0-10%
		cutarray[ 3] = "5240002002092970028250400000"; mesonCutArray[3]= "01525055000000"; // 20-40%
		cutarray[ 4] = "5250002002092970028250400000"; mesonCutArray[4]= "01525055000000"; // 20-50%
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
	} else if (periodName.CompareTo("LHC14a1a")==0){
// 		if (doWeightingInt == 1){ 
// 			TObjString *Header1 = new TObjString("pi0_1");
// 			HeaderList->Add(Header1);
// 		} else if (doWeightingInt == 2){
// 			TObjString *Header1 = new TObjString("eta_2");
// 			HeaderList->Add(Header1);
// 		} else {
			TObjString *Header1 = new TObjString("pi0_1");
			HeaderList->Add(Header1);
			TObjString *Header2 = new TObjString("eta_2");
			HeaderList->Add(Header2);
// 		}  
	} else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
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

<<<<<<< HEAD
		if (trainConfig == 56 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
			}	
		}	  
		if (trainConfig == 57 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
			}	
		}	  
		if (trainConfig == 58 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
			}	
		}	  
		if (trainConfig == 59 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
			}	
		}	  
	
		if (trainConfig == 60 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		}	  
		if (trainConfig == 61 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		}	  
		
		if ( trainConfig == 70 || trainConfig == 72  || trainConfig == 74  || trainConfig == 76  || trainConfig == 78  || trainConfig == 80  || trainConfig == 82  || trainConfig == 84 || trainConfig == 86  || trainConfig == 88  || trainConfig == 90 || trainConfig == 92 || trainConfig == 94 || trainConfig == 96  || trainConfig == 98  || trainConfig == 100 || trainConfig == 102  || trainConfig == 104 || trainConfig == 106 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		} 
		if ( trainConfig == 71 || trainConfig == 73  || trainConfig == 75  || trainConfig == 77  || trainConfig == 79  || trainConfig == 81  || trainConfig == 83  || trainConfig == 85 || trainConfig == 87  || trainConfig == 89  || trainConfig == 91 || trainConfig == 93 || trainConfig == 95 || trainConfig == 97  || trainConfig == 99  || trainConfig == 101 || trainConfig == 103  || trainConfig == 105 || trainConfig == 107 ){
			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
			}	
		} 
=======
// 		if (trainConfig == 56 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
// 			}	
// 		}	  
// 		if (trainConfig == 57 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
// 			}	
// 		}	  
// 		if (trainConfig == 58 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
// 			}	
// 		}	  
// 		if (trainConfig == 59 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
// 			}	
// 		}	  
// 	
// 		if (trainConfig == 60 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
// 			}	
// 		}	  
// 		if (trainConfig == 61 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
// 			}	
// 		}	  
// 		
// 		if ( trainConfig == 70 || trainConfig == 72  || trainConfig == 74  || trainConfig == 76  || trainConfig == 78  || trainConfig == 80  || trainConfig == 82  || trainConfig == 84 || trainConfig == 86  || trainConfig == 88  || trainConfig == 90 || trainConfig == 92 || trainConfig == 94 || trainConfig == 96  || trainConfig == 98  || trainConfig == 100 || trainConfig == 102  || trainConfig == 104 || trainConfig == 106 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
// 			}	
// 		} 
// 		if ( trainConfig == 71 || trainConfig == 73  || trainConfig == 75  || trainConfig == 77  || trainConfig == 79  || trainConfig == 81  || trainConfig == 83  || trainConfig == 85 || trainConfig == 87  || trainConfig == 89  || trainConfig == 91 || trainConfig == 93 || trainConfig == 95 || trainConfig == 97  || trainConfig == 99  || trainConfig == 101 || trainConfig == 103  || trainConfig == 105 || trainConfig == 107 ){
// 			if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
// 				if ( i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
// 				if ( i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
// 				if ( i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
// 				if ( i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
// 				if ( i == 4 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
// 			}	
// 		} 
>>>>>>> 014acd7... - changes for the reweighting part
		
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
