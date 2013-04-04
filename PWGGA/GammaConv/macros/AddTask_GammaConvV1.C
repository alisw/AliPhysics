
void AddTask_GammaConvV1(TString trainConfig = "pp",   Bool_t isMC	= kFALSE){

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

   

   TString cutnumber = "";
   if(IsHeavyIon == 1){
      cutnumber = "1000000002084001001500000";
    } else if (IsHeavyIon==2){
     cutnumber = "8000000002084001001500000";   
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
   Int_t numberOfCuts = 1;
   if(trainConfig.Contains("PbPb")) numberOfCuts = 3;
   else numberOfCuts = 1;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if(trainConfig.Contains("PbPb")){
      cutarray[ 0] = "1000003042092970723220000"; mesonCutArray[ 0] = "01022045000";  // all centralities
      cutarray[ 1] = "3010003042092970723220000"; mesonCutArray[ 1] = "01022045000";  // most central
      cutarray[ 2] = "1680003042092970723220000"; mesonCutArray[ 2] = "01022065000";  // peripheral

   }
   else if(trainConfig.Contains("pPb")){ //pA needs thighter rapidity cut y < 0.5
      cutarray[ 0] = "8000000042092172023290000"; mesonCutArray[0] = "01024045000";  //standard cut Pi0 PbPb 00-100
   } else {
      cutarray[ 0] = "0000011002093663003800000"; mesonCutArray[0] = "01631031009";
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
      if (trainConfig.Contains("pPb")) analysisCuts[i]->SelectCollisionCandidates(AliVEvent::kINT7|AliVEvent::kTRD);
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

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"GammaConvV1.root");

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return task;

}
