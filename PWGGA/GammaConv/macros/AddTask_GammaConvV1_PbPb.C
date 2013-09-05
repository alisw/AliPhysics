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
      cutarray[ 3] = "124000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% direct photon cut
   } else if (trainConfig == 3) { // Standard direct photon cuts
      cutarray[ 0] = "102000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "101000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10% 
      cutarray[ 2] = "112000104209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "148000104209237002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 4) { // Standard neutral pion cuts only added signals
      cutarray[ 0] = "301000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "312000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "101000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "112000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 5) { // Standard neutral pion cuts only added signals
      cutarray[ 0] = "124000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "146000204209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "168000204209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "124000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon cut
   } else if (trainConfig == 6) { // Standard direct photon cuts
      cutarray[ 0] = "102000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "101000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "112000204209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20
      cutarray[ 3] = "148000204209237002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 7){ // Standard neutral pion cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "601000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 8) { // Standard neutral pion cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "524000104209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 9) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209237002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 10) { // Standard neutral pion cuts only added signals, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "601000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209297002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209297002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209297002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 11) { // Standard neutral pion cuts only added signals, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "524000204209297002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209297002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209297002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 12) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209237002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 13) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209295102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209295102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209295102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209295102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 14) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209295102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209295102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209295102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209295102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 15) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209295102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209295102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209295102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209295102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 16) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209295102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209295102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209295102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209295102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 17) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209295102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209295102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209295102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209295102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 18) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209295102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209295102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209295102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209295102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 19) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209255102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209255102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209255102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209255102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 20) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209255102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209255102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209255102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209255102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 21) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209255102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209255102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209255102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209255102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 22) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209255102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209255102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209255102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209255102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 23) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209255102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209255102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209255102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209255102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 24) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209255102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209255102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209255102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209255102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%   
   } else if (trainConfig == 25) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000104209237002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209237002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209237002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209237002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 26) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000104209237002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209237002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209237002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 27) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000104209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209237002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 28) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000204209237002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209237002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209237002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209237002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 29) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000204209237002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209237002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209237002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209237002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 30) {// dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000204209237002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209237002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209237002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209237002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%         
   } else if (trainConfig == 31) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000104209277002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209277002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209277002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209277002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 32) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000104209277002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209277002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209277002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209277002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 33) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000104209277002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209277002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209277002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209277002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 34) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000204209277002322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209277002322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209277002322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209277002322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 35) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000204209277002322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209277002322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209277002322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209277002370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 36) {// dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000204209277002370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209277002370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209277002370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209277002370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%            
   } else if (trainConfig == 37) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209267102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209267102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209267102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209267102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 38) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209267102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209267102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209267102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209267102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 39) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209267102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209267102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209267102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209267102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 40) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209267102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209267102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209267102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209267102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 41) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209267102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209267102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209267102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209267102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 42) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209267102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209267102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209267102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209267102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%            
   } else if (trainConfig == 43) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209287102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000104209287102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000104209287102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000104209287102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 44) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209287102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000104209287102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000104209287102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000104209287102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 45) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209287102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000104209287102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000104209287102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000104209287102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%
   } else if (trainConfig == 46) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209287102322000000"; mesonCutArray[ 0] = "01522045009000"; // 0-5%
      cutarray[ 1] = "612000204209287102322000000"; mesonCutArray[ 1] = "01522045009000"; // 5-10%
      cutarray[ 2] = "501000204209287102322000000"; mesonCutArray[ 2] = "01522045009000"; // 0-10%
      cutarray[ 3] = "512000204209287102322000000"; mesonCutArray[ 3] = "01522045009000"; // 10-20%
   } else if (trainConfig == 47) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209287102322000000"; mesonCutArray[ 0] = "01522045009000"; // 20-40%
      cutarray[ 1] = "546000204209287102322000000"; mesonCutArray[ 1] = "01522065009000"; // 40-60%
      cutarray[ 2] = "568000204209287102322000000"; mesonCutArray[ 2] = "01522065009000"; // 60-80%
      cutarray[ 3] = "524000204209287102370000000"; mesonCutArray[ 3] = "01522045009000"; // 20-40% //direct photon
   } else if (trainConfig == 48) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209287102370000000"; mesonCutArray[ 0] = "01522045009000"; // 0-20%
      cutarray[ 1] = "501000204209287102370000000"; mesonCutArray[ 1] = "01522045009000"; // 0-10%
      cutarray[ 2] = "512000204209287102370000000"; mesonCutArray[ 2] = "01522045009000"; // 10-20%
      cutarray[ 3] = "548000204209287102370000000"; mesonCutArray[ 3] = "01522065009000"; // 40-80%            
    } else if (trainConfig == 49){ // Standard neutral pion cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "601000104209297002322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209297002322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209297002322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209297002322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 50) { // Standard neutral pion cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "524000104209297002322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209297002322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209297002322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209237002370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 51) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000104209237002370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209237002370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209237002370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209237002370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 52) { // Standard neutral pion cuts only added signals, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "601000204209297002322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209297002322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209297002322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209297002322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 53) { // Standard neutral pion cuts only added signals, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "524000204209297002322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209297002322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209297002322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209237002370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 54) { // Standard direct photon cuts, Centrality V0M for data and same TPC mult for MC
      cutarray[ 0] = "502000204209237002370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209237002370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209237002370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209237002370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 55) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209295102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209295102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209295102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209295102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 56) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209295102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209295102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209295102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209295102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 57) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209295102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209295102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209295102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209295102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 58) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209295102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209295102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209295102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209295102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 59) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209295102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209295102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209295102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209295102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 60) { // dE/dx E varied to \sigma > 3, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209295102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209295102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209295102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209295102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 61) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209255102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209255102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209255102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209255102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 62) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209255102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209255102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209255102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209255102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 63) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209255102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209255102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209255102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209255102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 64) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209255102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209255102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209255102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209255102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 65) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209255102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209255102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209255102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209255102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 66) { // dE/dx E varied to \sigma > 2, \sigma > -10, min p\pi = 0.3 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209255102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209255102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209255102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209255102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%   
   } else if (trainConfig == 67) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000104209237002322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209237002322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209237002322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209237002322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 68) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000104209237002322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209237002322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209237002322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209237002370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 69) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000104209237002370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209237002370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209237002370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209237002370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 70) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000204209237002322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209237002322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209237002322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209237002322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 71) { // dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000204209237002322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209237002322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209237002322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209237002370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 72) {// dE/dx E varied to \sigma > 2.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000204209237002370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209237002370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209237002370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209237002370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%         
   } else if (trainConfig == 73) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000104209277002322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209277002322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209277002322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209277002322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 74) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000104209277002322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209277002322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209277002322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209277002370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 75) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000104209277002370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209277002370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209277002370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209277002370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 76) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "601000204209277002322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209277002322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209277002322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209277002322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 77) { // dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "524000204209277002322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209277002322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209277002322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209277002370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 78) {// dE/dx E varied to \sigma > 3.5, \sigma > -10, min p\pi = 0.4 GeV/c, max p\pi = 100 GeV/c
      cutarray[ 0] = "502000204209277002370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209277002370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209277002370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209277002370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%            
   } else if (trainConfig == 79) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209267102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209267102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209267102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209267102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 80) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209267102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209267102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209267102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209267102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 81) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209267102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209267102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209267102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209267102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 82) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209267102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209267102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209267102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209267102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 83) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209267102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209267102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209267102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209267102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 84) { // dE/dx E varied to \sigma > 2, \sigma > 0.5, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209267102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209267102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209267102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209267102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%            
   } else if (trainConfig == 85) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000104209287102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000104209287102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000104209287102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000104209287102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 86) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000104209287102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000104209287102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000104209287102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000104209287102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 87) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000104209287102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000104209287102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000104209287102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000104209287102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%
   } else if (trainConfig == 88) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "601000204209287102322000000"; mesonCutArray[ 0] = "01522045009300"; // 0-5%
      cutarray[ 1] = "612000204209287102322000000"; mesonCutArray[ 1] = "01522045009300"; // 5-10%
      cutarray[ 2] = "501000204209287102322000000"; mesonCutArray[ 2] = "01522045009300"; // 0-10%
      cutarray[ 3] = "512000204209287102322000000"; mesonCutArray[ 3] = "01522045009300"; // 10-20%
   } else if (trainConfig == 89) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "524000204209287102322000000"; mesonCutArray[ 0] = "01522045009300"; // 20-40%
      cutarray[ 1] = "546000204209287102322000000"; mesonCutArray[ 1] = "01522065009300"; // 40-60%
      cutarray[ 2] = "568000204209287102322000000"; mesonCutArray[ 2] = "01522065009300"; // 60-80%
      cutarray[ 3] = "524000204209287102370000000"; mesonCutArray[ 3] = "01522045009300"; // 20-40% //direct photon
   } else if (trainConfig == 90) { // dE/dx E varied to \sigma > 2, \sigma > 1, min p\pi = 0.4 GeV/c, max p\pi = 5 GeV/c
      cutarray[ 0] = "502000204209287102370000000"; mesonCutArray[ 0] = "01522045009300"; // 0-20%
      cutarray[ 1] = "501000204209287102370000000"; mesonCutArray[ 1] = "01522045009300"; // 0-10%
      cutarray[ 2] = "512000204209287102370000000"; mesonCutArray[ 2] = "01522045009300"; // 10-20%
      cutarray[ 3] = "548000204209287102370000000"; mesonCutArray[ 3] = "01522065009300"; // 40-80%              
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
      if (trainConfig == 1 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M"); 
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 2){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 3){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      } else if (trainConfig == 4 ){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M"); 
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 5){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 6){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      } else if (trainConfig == 7 || trainConfig == 13 || trainConfig == 19 || trainConfig == 25 || trainConfig == 31 || trainConfig == 37 || trainConfig == 43 || trainConfig == 49 || trainConfig == 55 || trainConfig == 61 || trainConfig == 67 || trainConfig == 73 || trainConfig == 79 || trainConfig == 85){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 8 || trainConfig == 14 || trainConfig == 20 || trainConfig == 26 || trainConfig == 32 || trainConfig == 38 || trainConfig == 44 || trainConfig == 50 || trainConfig == 56 || trainConfig == 62 || trainConfig == 68 || trainConfig == 74 || trainConfig == 80 || trainConfig == 86){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2 && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 9 || trainConfig == 15 || trainConfig == 21 || trainConfig == 27 || trainConfig == 33 || trainConfig == 39 || trainConfig == 45 || trainConfig == 51 || trainConfig == 57 || trainConfig == 63 || trainConfig == 69 || trainConfig == 75 || trainConfig == 81 || trainConfig == 87){
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
        if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
      } else if (trainConfig == 10 || trainConfig == 16 || trainConfig == 22 || trainConfig == 28 || trainConfig == 34 || trainConfig == 40 || trainConfig == 46 || trainConfig == 52 || trainConfig == 58 || trainConfig == 64 || trainConfig == 70 || trainConfig == 76 || trainConfig == 82 || trainConfig == 88){
         if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
         if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
         if (i == 2 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
         if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      } else if (trainConfig == 11 || trainConfig == 17 || trainConfig == 23 || trainConfig == 29 || trainConfig == 35 || trainConfig == 41 || trainConfig == 47 || trainConfig == 53 || trainConfig == 59 || trainConfig == 65 || trainConfig == 71 || trainConfig == 77 || trainConfig == 83 || trainConfig == 89){ 
        if (i == 0 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
        if (i == 1 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
        if (i == 2  && doWeighting) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
        if (i == 3 && doWeighting)  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      } else if (trainConfig == 12 || trainConfig == 18 || trainConfig == 24 || trainConfig == 30 || trainConfig == 36 || trainConfig == 42 || trainConfig == 48 || trainConfig == 54 || trainConfig == 60 || trainConfig == 66 || trainConfig == 72 || trainConfig == 78 || trainConfig == 84 || trainConfig == 90){
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
