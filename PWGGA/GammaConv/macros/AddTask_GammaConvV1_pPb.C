void AddTask_GammaConvV1_pPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Bool_t enableQAMesonTask = kFALSE, //enable QA in AliAnalysisTaskGammaConvV1
                              Bool_t enableQAPhotonTask = kFALSE, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              TString cutnumberAODBranch = "8000000060084000001500000" // cutnumber for AOD branch
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
   TString cutnumber = "800000006008400100150000000"; 
   Bool_t doEtaShift = kFALSE;
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
            fCuts->DoEtaShift(doEtaShift);
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
   if(trainConfig == 1){ 
      // Shifting in pPb direction
      doEtaShiftIndCuts = kTRUE;
      stringShift = "pPb";
      cutarray[ 0] = "800000008009317200329000000"; mesonCutArray[ 0] = "01629045000000";  //standard cut Pi0 pPb RCut Changed
      cutarray[ 1] = "800000008109317200329000000"; mesonCutArray[ 1] = "01629045000000";  //standard cut Pi0 pPb RCut Changed
      cutarray[ 2] = "800000008309317200329000000"; mesonCutArray[ 2] = "01629045000000";  //standard cut Pi0 pPb RCut Changed
      cutarray[ 3] = "800000008409317200329000000"; mesonCutArray[ 3] = "01629045000000";  //standard cut Pi0 pPb RCut Changed
      cutarray[ 4] = "800000008509317200329000000"; mesonCutArray[ 4] = "01629045000000";  //standard cut Pi0 pPb RCut Changed
      cutarray[ 5] = "800000008609317200329000000"; mesonCutArray[ 5] = "01629045000000";  //standard cut Pi0 pPb RCut Changed
   } else if (trainConfig == 2) { 
      doEtaShiftIndCuts = kTRUE;
      stringShift = "pPb";
      cutarray[ 0] = "800000008219317200329000000"; mesonCutArray[ 0] = "01629045000000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 1] = "800000008229317200329000000"; mesonCutArray[ 1] = "01629045000000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 2] = "800000008249317200329000000"; mesonCutArray[ 2] = "01629045000000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 3] = "800000008269317200329000000"; mesonCutArray[ 3] = "01629045000000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 4] = "800000008279317200329000000"; mesonCutArray[ 4] = "01629045000000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 5] = "800000008259317200329000000"; mesonCutArray[ 5] = "01629045000000";  //standard cut Pi0 pPb Single pT Cut changed
   } else if (trainConfig == 3) { 
      doEtaShiftIndCuts = kTRUE;
      stringShift = "pPb";
      cutarray[ 0] = "800000008209017200329000000"; mesonCutArray[ 0] = "01629045000000";  //standard cut Pi0 pPb dedx electron Cut changed
      cutarray[ 1] = "800000008209117200329000000"; mesonCutArray[ 1] = "01629045000000";  //standard cut Pi0 pPb  dedx electron Cut changed
      cutarray[ 2] = "800000008209217200329000000"; mesonCutArray[ 2] = "01629045000000";  //standard cut Pi0 pPb  dedx electron Cut changed
      cutarray[ 3] = "800000008209417200329000000"; mesonCutArray[ 3] = "01629045000000";  //standard cut Pi0 pPb  dedx electron Cut changed
      cutarray[ 4] = "800000008209517200329000000"; mesonCutArray[ 4] = "01629045000000";  //standard cut Pi0 pPb  dedx electron Cut changed
      cutarray[ 5] = "800000008209617200329000000"; mesonCutArray[ 5] = "01629045000000";  //standard cut Pi0 pPb  dedx electron Cut changed
   } else if (trainConfig == 4) { 
      doEtaShiftIndCuts = kTRUE;
      stringShift = "pPb";
      cutarray[ 0] = "800000008209316200329000000"; mesonCutArray[ 0] = "01629045000000";  //standard cut Pi0 pPb dedx Pion Line changed
      cutarray[ 1] = "800000008209310200329000000"; mesonCutArray[ 1] = "01629045000000";  //standard cut Pi0 pPb dedx Pion Line changed
      cutarray[ 2] = "800000008209317400329000000"; mesonCutArray[ 2] = "01629045000000";  //standard cut Pi0 pPb dedx Pion Line changed
      cutarray[ 3] = "800000008209317300329000000"; mesonCutArray[ 3] = "01629045000000";  //standard cut Pi0 pPb dedx Pion Line changed
      cutarray[ 4] = "800000008209357200329000000"; mesonCutArray[ 4] = "01629045000000";  //standard cut Pi0 pPb dedx Pion Line changed
      cutarray[ 5] = "800000008209307200329000000"; mesonCutArray[ 5] = "01629045000000";  //standard cut Pi0 pPb dedx Pion Line changed
   } else if (trainConfig == 5) { 
      doEtaShiftIndCuts = kFALSE;
      stringShift = "pPb"; //GammaCuts
      cutarray[ 0] = "800001100209326000380000000"; mesonCutArray[ 0] = "01631031009000"; 
      cutarray[ 1] = "802001100209326000380000000"; mesonCutArray[ 1] = "01631031009000"; 
      cutarray[ 2] = "804001100209326000380000000"; mesonCutArray[ 2] = "01631031009000"; 
      cutarray[ 3] = "828001100209326000380000000"; mesonCutArray[ 3] = "01631031009000"; 
      cutarray[ 4] = "848001100209326000380000000"; mesonCutArray[ 4] = "01631031009000"; 
      cutarray[ 5] = "800001100209226000380000000"; mesonCutArray[ 5] = "01631031009000"; 
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
