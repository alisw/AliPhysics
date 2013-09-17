AliAnalysisTask *AddTask_gamconv_GammaConvV1_PbPb(){
   
   gSystem->Load("libANALYSISalice.so");

   //get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_GammaConvV1", "No analysis manager found.");
      return 0;
   }

   TString trainConfig=gSystem->Getenv("CONFIG_FILE");
   Int_t IsHeavyIon=1;
   
   Bool_t IsMC = kFALSE;
   if (trainConfig.Contains("MC")) IsMC=kTRUE;


   TString cutnumber = "108000000008400100150000000";
   Bool_t doEtaShift = kFALSE;
   
   //========= Add PID Reponse to ANALYSIS manager ====
   if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AddTaskPIDResponse(IsMC);
   }
   
   //========= Add V0 Reader to  ANALYSIS manager =====
   AliV0ReaderV1 *fV0ReaderV1=new AliV0ReaderV1("V0ReaderV1");
   ConfigV0ReaderV1(fV0ReaderV1,cutnumber,IsHeavyIon,doEtaShift);
   mgr->AddTask(fV0ReaderV1);

   AliLog::SetGlobalLogLevel(AliLog::kInfo);

   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
   //connect input V0Reader
   mgr->ConnectInput (fV0ReaderV1,0,cinput);

   AliAnalysisTaskGammaConvV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvV1("GammaConvV1");
   task->SetIsHeavyIon(IsHeavyIon);
   task->SetIsMC(IsMC);
   // Cut Numbers to use in Analysis

   Int_t numberOfCuts = 4;
   
   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];

   cutarray[ 0] = "301000204209297002322000000"; mesonCutArray[ 0] = "01522045009000";
   cutarray[ 1] = "124000204209297002322000000"; mesonCutArray[ 1] = "01522045009000";
   cutarray[ 2] = "146000204209297002322000000"; mesonCutArray[ 2] = "01522065009000";
   cutarray[ 3] = "168000204209297002322000000"; mesonCutArray[ 3] = "01522065009000";

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
//       if (i==0) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
//       if (i==1) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
//       if (i==2) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
//       if (i==3) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");     
      if (i==0) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i==1) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i==2) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i==3) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080V0M", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      
//       if (i==1) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0010");
//        if (i==2) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0510"); 
//       if (i==3) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE,"MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_1020");      
//       if (i==7) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0020");
//       if (i==8) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0040", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0040");
//       if (i==9) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0080", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0080");
//       if (i==10) analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kTRUE, "MCSpectraInput.root", "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080", "", "K0s_RatioDataToMC_Hijing_PbPb_2760GeV_4080");

      analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
      ConvCutList->Add(analysisCuts[i]);
      
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
      analysisCuts[i]->SetAcceptedHeader(HeaderList);
      analysisMesonCuts[i] = new AliConversionMesonCuts();
      analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");

//      if ( i < 5){
//         AliConversionCuts *QACuts = new AliConversionCuts();
//         QACuts->InitializeCutsFromCutString(cutarray[i].Data());
//         QACuts->SetFillCutHistograms("",kTRUE);
//         QACuts->SetAcceptedHeader(HeaderList);
//         AddQATaskV1(QACuts,IsHeavyIon,kTRUE,kFALSE,IsMC);
//       }

   }
   
   task->SetConversionCutList(numberOfCuts,ConvCutList);
   task->SetMesonCutList(numberOfCuts,MesonCutList);
   task->SetMoveParticleAccordingToVertex(kTRUE);
   task->SetDoMesonAnalysis(kTRUE);
   task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
   task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA
   mgr->AddTask(task);

   //connect containers
   AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"GammaConvV1.root");
   mgr->ConnectInput  (task,  0, cinput );
   mgr->ConnectOutput (task,  1, coutput1);

   return task;
}

void ConfigV0ReaderV1(AliV0ReaderV1 *fV0ReaderV1,TString analysiscut="",Int_t IsHeavyIon=0,Bool_t doEtaShift = kFALSE){

   fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
   fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
   fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
   }
   AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

   if(inputHandler->IsA()==AliESDInputHandler::Class()){
      // ESD mode
   }

   // Set AnalysisCut Number
   AliConversionCuts *fCuts=NULL;
   if(analysiscut!=""){
      fCuts= new AliConversionCuts(analysiscut.Data(),analysiscut.Data());
      if(fCuts->InitializeCutsFromCutString(analysiscut.Data())){
         if (IsHeavyIon==2){
            fCuts->DoEtaShift(doEtaShift);
         }
         fCuts->SetPreSelectionCutFlag(kTRUE);
         fV0ReaderV1->SetConversionCuts(fCuts);
         fCuts->SetFillCutHistograms("",kTRUE);
      }
   }

   fV0ReaderV1->Init();
}

void AddQATaskV1(AliConversionCuts *ConversionCuts, Int_t IsHeavyIon, Bool_t tree, Bool_t histograms, Bool_t IsMCQA){

   //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
   Error("AddTask_V0ReaderV1QA", "No analysis manager found.");
   return 0;
    }
    AliAnalysisTaskConversionQA *fQA =
       new AliAnalysisTaskConversionQA(Form("%s_QA",(ConversionCuts->GetCutNumber()).Data()));
    fQA->SetConversionCuts(ConversionCuts,IsHeavyIon);
    fQA->FillType(tree,histograms);
    fQA->SetIsMC(IsMCQA);
    mgr->AddTask(fQA);

    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput =
       mgr->CreateContainer(Form("GammaConv_V1QA_%s",(ConversionCuts->GetCutNumber()).Data()), TList::Class(),
                            AliAnalysisManager::kOutputContainer,Form("GammaConvV1_QAHist_%s.root",(ConversionCuts->GetCutNumber()).Data()));


    mgr->ConnectInput(fQA,0,cinput);
    if(histograms) mgr->ConnectOutput(fQA,  1, coutput);

    if(tree){
       TString addoutput=gSystem->Getenv("ADD_OUTPUT_FILES");
       if (addoutput.Length()) addoutput+=",";
       addoutput+=Form("GammaConvV1_QATree_%s.root",(ConversionCuts->GetCutNumber()).Data());
       gSystem->Setenv("ADD_OUTPUT_FILES",addoutput.Data());
       cout<<"Adding addoutput.Data(): "<<Form("GammaConvV1_QATree_%s.root",(ConversionCuts->GetCutNumber()).Data())<<endl;
    }

}
