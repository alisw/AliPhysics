
void AddTask_GammaConvV1(TString trainConfig = "pp",   Bool_t isMC	= kFALSE, UInt_t triggerMaskpPb = AliVEvent::kINT7 ){

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
      

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_GammaConvV1", "No analysis manager found.");
      return 0;
   }
      
   Int_t IsHeavyIon=0;
   if (trainConfig.Contains("PbPb")) IsHeavyIon=1;
   else if (trainConfig.Contains("pPb")) IsHeavyIon=2;

   Bool_t doEtaShift = kFALSE;
   
   TString cutnumber = "";
   if(IsHeavyIon == 1){
      cutnumber = "1000000002084001001500000";
    } else if (IsHeavyIon==2){
      cutnumber = "8000000062084001001500000";
      doEtaShift = kFALSE; // Only needed if Eta Range is small, but default is now +- 5.0
   } else{
      cutnumber = "0000000002084000002200000";
   }

   //========= Add PID Reponse to ANALYSIS manager ====
   if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AddTaskPIDResponse();
   }
      
   //========= Add V0 Reader to  ANALYSIS manager =====
   AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
   
   fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
   fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
   fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

   if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
   }
   AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

   if(inputHandler->IsA()==AliESDInputHandler::Class()){
      // ESD mode
   }

   if(inputHandler->IsA()==AliAODInputHandler::Class()){
      // AOD mode
      //   task->SetUseSatelliteAODs(kTRUE);
   }

   // Set AnalysisCut Number
   AliConversionCuts *fCuts=NULL;
   if(cutnumber!=""){
      fCuts= new AliConversionCuts(cutnumber.Data(),cutnumber.Data());
      if(fCuts->InitializeCutsFromCutString(cutnumber.Data())){
         if (IsHeavyIon==2){
            fCuts->SelectCollisionCandidates(triggerMaskpPb);
            fCuts->DoEtaShift(doEtaShift);
         }
         fV0ReaderV1->SetConversionCuts(fCuts);
         fCuts->SetFillCutHistograms("",kTRUE);
      }
   }
   fV0ReaderV1->Init();



   AliLog::SetGlobalLogLevel(AliLog::kInfo);

   //================================================
   //              data containers
   //================================================
   //            find input container
   //below the trunk version
   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   //connect input V0Reader
   mgr->AddTask(fV0ReaderV1);
   mgr->ConnectInput(fV0ReaderV1,0,cinput);

   //================================================
   //========= Add task to the ANALYSIS manager =====
   //================================================
   //              data containers
   //================================================
   //            find input container
   AliAnalysisTaskGammaConvV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvV1("GammaConvV1");
   task->SetIsHeavyIon(IsHeavyIon);
   // Cut Numbers to use in Analysis
   Int_t numberOfCuts = 2;
   if(trainConfig.Contains("PbPb")) numberOfCuts = 7;
   else if(trainConfig.Contains("pPb")) numberOfCuts = 5;
   else numberOfCuts = 3;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if(trainConfig.Contains("PbPb")){
      cutarray[ 0] = "3010001042092970023220000"; mesonCutArray[ 0] = "01522045000";
      cutarray[ 1] = "1010001042092970023220000"; mesonCutArray[ 1] = "01522045000";
      cutarray[ 2] = "3120001042092970023220000"; mesonCutArray[ 2] = "01522045000";
      cutarray[ 3] = "1120001042092970023220000"; mesonCutArray[ 3] = "01522045000";
      cutarray[ 4] = "1240001042092970023220000"; mesonCutArray[ 4] = "01522045000";
      cutarray[ 5] = "1460001042092970023220000"; mesonCutArray[ 5] = "01522065000";
      cutarray[ 6] = "1680001042092970023220000"; mesonCutArray[ 6] = "01522065000";
   } else if(trainConfig.Contains("pPb")){ //pA needs thighter rapidity cut y < 0.5
      cutarray[ 0] = "8000000082093172023290000"; mesonCutArray[ 0] = "01629045000"; //shifted Pbp
      cutarray[ 1] = "8020000082093172023290000"; mesonCutArray[ 1] = "01629045000";
      cutarray[ 2] = "8240000082093172023290000"; mesonCutArray[ 2] = "01629045000";
      cutarray[ 3] = "8460000082093172023290000"; mesonCutArray[ 3] = "01629045000";
      cutarray[ 4] = "8600000082093172023290000"; mesonCutArray[ 4] = "01629045000";

   } else {
      cutarray[ 0] = "0000012002093663003800000"; mesonCutArray[0] = "01631031009"; //standard cut Pi0 pp 2.76TeV without SDD , only boxes
      cutarray[ 1] = "0001012002093663003800000"; mesonCutArray[1] = "01631031009"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND , only boxes
      cutarray[ 2] = "0000012002093260003800000"; mesonCutArray[2] = "01631031009"; //standard cut Gamma pp 2-76TeV , only boxes

   }

   TList *ConvCutList = new TList();
   TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   // TObjString *Header1 = new TObjString("PARAM");
   // HeaderList->Add(Header1);
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
//   TObjString *Header3 = new TObjString("Pythia");
//   HeaderList->Add(Header3);
//   TObjString *Header4 = new TObjString("Hijing");
//   HeaderList->Add(Header4);


   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];


   for(Int_t i = 0; i<numberOfCuts; i++){
      analysisCuts[i] = new AliConversionCuts();
      analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
      if (trainConfig.Contains("pPb")){
         analysisCuts[i]->SelectCollisionCandidates(triggerMaskpPb);
         if (i<5){
            analysisCuts[i]->DoEtaShift(kTRUE);
            analysisCuts[i]->SetEtaShift("Pbp");
         }
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
   if (trainConfig.Contains("pPb") || trainConfig.Contains("pp") )task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
   if (trainConfig.Contains("pPb") || trainConfig.Contains("pp") )task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"GammaConvV1.root");

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return task;

}
