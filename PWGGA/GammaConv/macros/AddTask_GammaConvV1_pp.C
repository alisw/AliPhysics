void AddTask_GammaConvV1_pp(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC 
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              TString cutnumberAODBranch = "0000000060084001001500000" // cutnumber for AOD branch
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
   TString cutnumber = "0000000002084000002200000000"; 
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
   //            find input container
   AliAnalysisTaskGammaConvV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
   task->SetIsHeavyIon(0);
   task->SetIsMC(isMC);
   // Cut Numbers to use in Analysis
   Int_t numberOfCuts = 4;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if(trainConfig == 1){
      cutarray[ 0] = "0000012002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , only boxes
      cutarray[ 1] = "0001012002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND , only boxes
      cutarray[ 2] = "0000012002093260003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Gamma pp 2-76TeV , only boxes
      cutarray[ 3] = "0000012002093260003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Gamma pp 2-76TeV , only boxes
   } else if (trainConfig == 2) {
      cutarray[ 0] = "0000011002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , only Minbias MC
      cutarray[ 1] = "0001011002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND
      cutarray[ 2] = "0000011002093260003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Gamma pp 2-76TeV
      cutarray[ 3] = "0000011002093260003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Gamma pp 2-76TeV , only boxes
   } else if (trainConfig == 3) {
      cutarray[ 0] = "0002011002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
      cutarray[ 1] = "0003011002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
      cutarray[ 2] = "0002012002093663003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
      cutarray[ 3] = "0003012002093663003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
   } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
   }

   TList *ConvCutList = new TList();
   TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   TObjString *Header2 = new TObjString("BOX");
   HeaderList->Add(Header2);

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
