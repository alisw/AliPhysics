void AddTask_GammaConvV1_pPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Bool_t enableQAMesonTask = kFALSE, //enable QA in AliAnalysisTaskGammaConvV1
                              Bool_t enableQAPhotonTask = kFALSE, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root" // path to file for weigting input
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
   TString cutnumber = "8000000060084001001500000"; 
   Bool_t doEtaShift = kFALSE;
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
            fCuts->DoEtaShift(doEtaShift);
            fV0ReaderV1->SetConversionCuts(fCuts);
            fCuts->SetFillCutHistograms("",kTRUE);
         }
      }
      fV0ReaderV1->Init();

      AliLog::SetGlobalLogLevel(AliLog::kInfo);

      //================================================
      //              data containers for V0Reader
      //================================================
      AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
      //connect input V0Reader
      mgr->AddTask(fV0ReaderV1);
      mgr->ConnectInput(fV0ReaderV1,0,cinput);

   }

   //================================================
   //========= Add task to the ANALYSIS manager =====
   //================================================
   //              data containers
   //================================================
   //            find input container
   AliAnalysisTaskGammaConvV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
   task->SetIsHeavyIon(2);
   task->SetIsMC(isMC);
   // Cut Numbers to use in Analysis
   Int_t numberOfCuts = 6;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];
   Bool_t doEtaShiftIndCuts = kFALSE;
   TString stringShift = "";
   if(trainConfig === 1){ 
      // Shifting in pPb direction
      doEtaShiftIndCuts = kTRUE;
      stringShift = "pPb";
      cutarray[ 0] = "8000000082093172003290000"; mesonCutArray[ 0] = "01629045000";  //standard cut Pi0 pPb 00-100 Midrapidity
      cutarray[ 1] = "8020000082093172003290000"; mesonCutArray[ 1] = "01629045000";  //standard cut Pi0 pPb 00-20 Midrapidity
      cutarray[ 2] = "8240000082093172003290000"; mesonCutArray[ 2] = "01629045000";  //standard cut Pi0 pPb 20-40 Midrapidity
      cutarray[ 3] = "8460000082093172003290000"; mesonCutArray[ 3] = "01629045000";  //standard cut Pi0 pPb 40-60 Midrapidity
      cutarray[ 4] = "8600000082093172003290000"; mesonCutArray[ 4] = "01629045000";  //standard cut Pi0 pPb 60-100 Midrapidity
      cutarray[ 5] = "8000000082093172003290000"; mesonCutArray[ 5] = "01620045000";  //standard cut Pi0 pPb0-100 Midrapidity open y cut
   } else if (trainConfig == 2) { 
      // no shifting
      cutarray[ 0] = "8000000072093172003290000"; mesonCutArray[ 0] = "01627045000";  //standard cut Pi0 pPb 00-100 Fwd (PHOS) rapidity
      cutarray[ 1] = "8020000072093172003290000"; mesonCutArray[ 1] = "01627045000";  //standard cut Pi0 pPb 00-20 Fwd (PHOS) rapidity
      cutarray[ 2] = "8240000072093172003290000"; mesonCutArray[ 2] = "01627045000";  //standard cut Pi0 pPb 20-40 Fwd (PHOS) rapidity
      cutarray[ 3] = "8460000072093172003290000"; mesonCutArray[ 3] = "01627045000";  //standard cut Pi0 pPb 40-60 Fwd (PHOS) rapidity
      cutarray[ 4] = "8600000072093172003290000"; mesonCutArray[ 4] = "01627045000";  //standard cut Pi0 pPb 60-100 Fwd (PHOS) rapidity
      cutarray[ 5] = "8000000072093172003290000"; mesonCutArray[ 5] = "01620045000";  //standard cut Pi0 pPb 00-100 Fwd (PHOS) rapidity open y cut
   } else if (trainConfig == 3) { 
      // no shifting
      cutarray[ 0] = "8000000002093172003290000"; mesonCutArray[ 0] = "01621045000";  //standard cut Pi0 pPb 00-100 Full Eta Range 
      cutarray[ 1] = "8020000002093172003290000"; mesonCutArray[ 1] = "01621045000";  //standard cut Pi0 pPb 00-20 Full Eta Range 
      cutarray[ 2] = "8240000002093172003290000"; mesonCutArray[ 2] = "01621045000";  //standard cut Pi0 pPb 20-40 Full Eta Range 
      cutarray[ 3] = "8460000002093172003290000"; mesonCutArray[ 3] = "01621045000";  //standard cut Pi0 pPb 40-60 Full Eta Range 
      cutarray[ 4] = "8600000002093172003290000"; mesonCutArray[ 4] = "01621045000";  //standard cut Pi0 pPb 60-100 Full Eta Range 
      cutarray[ 5] = "8000000002093172003290000"; mesonCutArray[ 5] = "01620045000";  //standard cut Pi0 pPb 00-100 Full Eta Range open y cut

   } else if (trainConfig == 4) { 
      // shifting in Pbp direction
      doEtaShiftIndCuts = kTRUE;
      stringShift = "Pbp";
      cutarray[ 0] = "8000000082093172003290000"; mesonCutArray[ 0] = "01629045000";  //standard cut Pi0 pPb 00-100 Pbp shift
      cutarray[ 1] = "8020000082093172003290000"; mesonCutArray[ 1] = "01629045000";  //standard cut Pi0 pPb 00-20 Pbp shift
      cutarray[ 2] = "8240000082093172003290000"; mesonCutArray[ 2] = "01629045000";  //standard cut Pi0 pPb 20-40 Pbp shift
      cutarray[ 3] = "8460000082093172003290000"; mesonCutArray[ 3] = "01629045000";  //standard cut Pi0 pPb 40-60 Pbp shift
      cutarray[ 4] = "8600000082093172003290000"; mesonCutArray[ 4] = "01629045000";  //standard cut Pi0 pPb 60-100 Pbp shift
      cutarray[ 5] = "8000000082093172003290000"; mesonCutArray[ 5] = "01620045000";  //standard cut Pi0 pPb0-100 Pbp shift open y cut
   } else if (trainConfig == 5) { 
      // Special trigger setup
      cutarray[ 0] = "8000000002093172003290000"; mesonCutArray[ 0] = "01620045000";  //standard cut Pi0 pPb open eta and y cut
      cutarray[ 1] = "8004000002093172003290000"; mesonCutArray[ 1] = "01620045000";  //Special trigger kTRD
      cutarray[ 2] = "8005000002093172003290000"; mesonCutArray[ 2] = "01620045000";  //Special trigger kPHI7
      cutarray[ 3] = "8006000002093172003290000"; mesonCutArray[ 3] = "01620045000"; // Special trigger kEMCEJE
      cutarray[ 4] = "8007000002093172003290000"; mesonCutArray[ 4] = "01620045000"; // Special trigger kEMCEGA
      cutarray[ 5] = "8008000002093172003290000"; mesonCutArray[ 5] = "01620045000"; // Special trigger kHighMult
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
      analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
      if (doEtaShiftIndCuts) {
         analysisCuts[i]->DoEtaShift(doEtaShiftIndCuts);
         analysisCuts[i]->SetEtaShift(stringShift);
      }
      
      if (trainConfig == 5) { 
         if (i == 1) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kTRD, "AliVEvent::kTRD" );
         if (i == 2) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kPHI7,"AliVEvent::kPHI7" );
         if (i == 3) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kEMCEJE,"AliVEvent::kEMCEJE" );
         if (i == 4) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kEMCEGA,"AliVEvent::kEMCEGA" );
         if (i == 5) analysisCuts[i]->SelectSpecialTrigger(AliVEvent::kHighMult,"AliVEvent::kHighMult" );
      }
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
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return;

}
