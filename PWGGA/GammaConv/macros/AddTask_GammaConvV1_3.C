
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

   Bool_t doEtaShift = kFALSE;
   Int_t forceEtaShift = 0; // Carefull !!! Should be zero otherwise will shift eta cut for all periods
                            // Use doEtaShift flag for pPb or Pbp instead (1: shift +0.465, 2: shift -0.465)

   TString cutnumber = "";
   if(IsHeavyIon == 1){
      cutnumber = "1000000002084001001500000";
    } else if (IsHeavyIon==2){
     cutnumber = "8000000062084001001500000";
     doEtaShift = kTRUE;
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
            fCuts->SelectCollisionCandidates(AliVEvent::kINT7);
            fCuts->DoEtaShift(doEtaShift);
         }
         fCuts->ForceEtaShift(forceEtaShift);
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
   if(trainConfig.Contains("PbPb")) numberOfCuts = 3;
   else if(trainConfig.Contains("pPb")) numberOfCuts = 4;
   else numberOfCuts = 4;

   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   if(trainConfig.Contains("PbPb")){
      cutarray[ 0] = "1460001042092970023220000"; mesonCutArray[ 0] = "01522065000"; // Standard cut 40-60
      cutarray[ 1] = "1020003042092370023750000"; mesonCutArray[ 1] = "01522045009";// Standard cut gamma 0-20%
      cutarray[ 2] = "1480003042092370023750000"; mesonCutArray[ 2] = "01522065009"; // Standard cut gamma 40-80
   } else if(trainConfig.Contains("pPb")){ //pA needs thighter rapidity cut y < 0.5
     cutarray[ 0] = "8020000072093172023290000"; mesonCutArray[0] = "01627045000";  //standard cut Pi0 Pb 00-20  wo shifted Eta 0.3
     cutarray[ 1] = "8240000072093172023290000"; mesonCutArray[1] = "01627045000";  //standard cut Pi0 Pb 20-40 wo shifted Eta 0.3
     cutarray[ 2] = "8460000072093172023290000"; mesonCutArray[2] = "01627045000";  //standard cut Pi0 Pb 40-60 wo shifted Eta 0.3
     cutarray[ 3] = "8680000072093172023290000"; mesonCutArray[3] = "01627045000";  //standard cut Pi0 Pb 60-80 wo shifted Eta 0.3
   } else {
      cutarray[ 0] = "0002011002093663003800000"; mesonCutArray[0] = "01631031009"; //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
      cutarray[ 1] = "0003011002093663003800000"; mesonCutArray[1] = "01631031009"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
      cutarray[ 2] = "0002011002093663003800000"; mesonCutArray[2] = "01631031009"; //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
      cutarray[ 3] = "0003012002093663003800000"; mesonCutArray[3] = "01631031009"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
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
         analysisCuts[i]->SelectCollisionCandidates(AliVEvent::kINT7);
//          if (i<4) analysisCuts[i]->DoEtaShift(1);
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
//    if (trainConfig.Contains("pPb") || trainConfig.Contains("pp") )task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
//    if (trainConfig.Contains("pPb") || trainConfig.Contains("pp") )task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"GammaConvV1.root");

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return task;

}
