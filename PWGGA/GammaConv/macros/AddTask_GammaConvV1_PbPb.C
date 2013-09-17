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

   if (trainConfig == 1){ 
      cutarray[ 0] = "601000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 2) { 
      cutarray[ 0] = "524000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 3) { 
      cutarray[ 0] = "502000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 4) { 
      cutarray[ 0] = "601000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 5) { 
      cutarray[ 0] = "524000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 6) { 
      cutarray[ 0] = "502000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 7){ // Psi pair 10000
      cutarray[ 0] = "601000104209297002320000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002320000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002320000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002320000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 8) { // Psi pair 10000
      cutarray[ 0] = "524000104209297002320000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002320000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002320000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209297002320000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 9) { // Psi pair 10000
      cutarray[ 0] = "502000104209297002320000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209297002320000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209297002320000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209297002320000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 10) { // Psi pair 10000
      cutarray[ 0] = "601000204209297002320000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002320000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002320000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002320000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 11) { // Psi pair 10000
      cutarray[ 0] = "524000204209297002320000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002320000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002320000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209297002320000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 12) { // Psi pair 10000
      cutarray[ 0] = "502000204209297002320000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209297002320000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209297002320000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209297002320000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%   
   } else if (trainConfig == 13){ // Psi pair 0.1
      cutarray[ 0] = "601000104209297002321000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002321000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002321000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002321000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 14) { // Psi pair 0.1
      cutarray[ 0] = "524000104209297002321000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002321000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002321000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209297002321000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 15) { // Psi pair 0.1
      cutarray[ 0] = "502000104209297002321000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209297002321000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209297002321000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209297002321000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 16) { // Psi pair 0.1
      cutarray[ 0] = "601000204209297002321000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002321000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002321000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002321000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 17) { // Psi pair 0.1
      cutarray[ 0] = "524000204209297002321000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002321000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002321000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209297002321000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 18) { // Psi pair 0.1
      cutarray[ 0] = "502000204209297002321000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209297002321000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209297002321000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209297002321000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%      
   } else if (trainConfig == 19){ // Psi pair 0.035
      cutarray[ 0] = "601000104209297002323000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002323000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002323000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002323000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 20) { // Psi pair 0.035
      cutarray[ 0] = "524000104209297002323000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002323000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002323000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209297002323000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 21) { // Psi pair 0.035
      cutarray[ 0] = "502000104209297002323000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209297002323000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209297002323000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209297002323000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 22) { // Psi pair 0.035
      cutarray[ 0] = "601000204209297002323000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002323000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002323000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002323000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 23) { // Psi pair 0.035
      cutarray[ 0] = "524000204209297002323000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002323000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002323000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209297002323000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 24) { // Psi pair 0.035
      cutarray[ 0] = "502000204209297002323000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209297002323000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209297002323000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209297002323000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%         
   } else if (trainConfig == 25){ // Eta 0.75
      cutarray[ 0] = "601000103209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000103209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000103209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000103209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 26) { // Eta 0.75
      cutarray[ 0] = "524000103209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000103209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000103209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000103209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 27) { // Eta 0.75
      cutarray[ 0] = "502000103209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000103209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000103209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000103209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 28) { // Eta 0.75
      cutarray[ 0] = "601000203209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000203209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000203209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000203209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 29) { // Eta 0.75
      cutarray[ 0] = "524000203209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000203209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000203209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000203209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 30) { // Eta 0.75
      cutarray[ 0] = "502000203209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000203209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000203209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000203209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%         
   } else if (trainConfig == 31){ // Eta 0.8
      cutarray[ 0] = "601000100209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000100209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000100209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000100209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 32) { // Eta 0.8
      cutarray[ 0] = "524000100209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000100209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000100209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000100209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 33) { // Eta 0.8
      cutarray[ 0] = "502000100209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000100209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000100209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000100209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 34) { // Eta 0.8
      cutarray[ 0] = "601000200209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000200209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000200209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000200209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 35) { // Eta 0.8
      cutarray[ 0] = "524000200209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000200209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000200209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000200209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 36) { // Eta 0.8
      cutarray[ 0] = "502000200209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000200209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000200209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000200209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%         
   } else if (trainConfig == 37){ // Eta 0.6
      cutarray[ 0] = "601000101209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000101209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000101209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000101209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 38) { // Eta 0.6
      cutarray[ 0] = "524000101209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000101209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000101209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000101209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 39) { // Eta 0.6
      cutarray[ 0] = "502000101209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000101209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000101209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000101209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 40) { // Eta 0.6
      cutarray[ 0] = "601000201209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000201209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000201209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000201209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 41) { // Eta 0.6
      cutarray[ 0] = "524000201209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000201209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000201209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000201209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 42) { // Eta 0.6
      cutarray[ 0] = "502000201209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000201209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000201209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000201209297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 43){ // R 10-180
      cutarray[ 0] = "601000104509297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104509297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104509297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104509297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 44) {  // R 10-180
      cutarray[ 0] = "524000104509297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104509297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104509297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104509297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 45) {  // R 10-180
      cutarray[ 0] = "502000104509297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104509297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104509297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104509297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 46) {  // R 10-180
      cutarray[ 0] = "601000204509297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204509297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204509297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204509297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 47) {  // R 10-180
      cutarray[ 0] = "524000204509297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204509297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204509297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204509297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 48) {  // R 10-180
      cutarray[ 0] = "502000204509297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204509297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204509297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204509297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%  
    } else if (trainConfig == 49){ // R 2.8-180
      cutarray[ 0] = "601000104109297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104109297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104109297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104109297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 50) {  // R 2.8-180
      cutarray[ 0] = "524000104109297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104109297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104109297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104109297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 51) {  // R 2.8-180
      cutarray[ 0] = "502000104109297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104109297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104109297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104109297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 52) {  // R 2.8-180
      cutarray[ 0] = "601000204109297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204109297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204109297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204109297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 53) {  // R 2.8-180
      cutarray[ 0] = "524000204109297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204109297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204109297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204109297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 54) {  // R 2.8-180
      cutarray[ 0] = "502000204109297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204109297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204109297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204109297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%  
    } else if (trainConfig == 55){ // R 2.8-180
      cutarray[ 0] = "601000104609297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104609297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104609297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104609297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 56) {  // R 2.8-180
      cutarray[ 0] = "524000104609297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104609297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104609297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104609297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 57) {  // R 2.8-180
      cutarray[ 0] = "502000104609297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104609297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104609297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104609297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 58) {  // R 2.8-180
      cutarray[ 0] = "601000204609297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204609297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204609297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204609297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 59) {  // R 2.8-180
      cutarray[ 0] = "524000204609297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204609297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204609297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204609297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 60) {  // R 2.8-180
      cutarray[ 0] = "502000204609297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204609297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204609297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204609297002322000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%  
     
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
      if (trainConfig == 1 ||trainConfig == 7 || trainConfig == 13 || trainConfig == 19 || trainConfig == 25 || trainConfig == 31 || trainConfig == 37 || trainConfig == 43 || trainConfig == 49 || trainConfig == 55 || trainConfig == 61 || trainConfig == 67 || trainConfig == 73 || trainConfig == 79 || trainConfig == 85){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 2 ||trainConfig == 8 || trainConfig == 14 || trainConfig == 20 || trainConfig == 26 || trainConfig == 32 || trainConfig == 38 || trainConfig == 44 || trainConfig == 50 || trainConfig == 56 || trainConfig == 62 || trainConfig == 68 || trainConfig == 74 || trainConfig == 80 || trainConfig == 86){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 3 ||trainConfig == 9 || trainConfig == 15 || trainConfig == 21 || trainConfig == 27 || trainConfig == 33 || trainConfig == 39 || trainConfig == 45 || trainConfig == 51 || trainConfig == 57 || trainConfig == 63 || trainConfig == 69 || trainConfig == 75 || trainConfig == 81 || trainConfig == 87){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      } else if (trainConfig == 4 ||trainConfig == 10 || trainConfig == 16 || trainConfig == 22 || trainConfig == 28 || trainConfig == 34 || trainConfig == 40 || trainConfig == 46 || trainConfig == 52 || trainConfig == 58 || trainConfig == 64 || trainConfig == 70 || trainConfig == 76 || trainConfig == 82 || trainConfig == 88){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 5 ||trainConfig == 11 || trainConfig == 17 || trainConfig == 23 || trainConfig == 29 || trainConfig == 35 || trainConfig == 41 || trainConfig == 47 || trainConfig == 53 || trainConfig == 59 || trainConfig == 65 || trainConfig == 71 || trainConfig == 77 || trainConfig == 83 || trainConfig == 89){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 6 ||trainConfig == 12 || trainConfig == 18 || trainConfig == 24 || trainConfig == 30 || trainConfig == 36 || trainConfig == 42 || trainConfig == 48 || trainConfig == 54 || trainConfig == 60 || trainConfig == 66 || trainConfig == 72 || trainConfig == 78 || trainConfig == 84 || trainConfig == 90){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
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
