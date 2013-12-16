void AddTask_GammaConvV1_pPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Bool_t doWeighting = kFALSE,  //enable Weighting
                              TString generatorName = "DPMJET",
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
   TString cutnumber = "8000000060084001001500000000";
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
         cout << "AOD handler: adding " << cutnumberAODBranch.Data() << " as conversion branch" << endl;
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
      cutarray[ 0] = "8000011082093172003290000000"; mesonCutArray[ 0] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 1] = "8020011082093172003290000000"; mesonCutArray[ 1] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 2] = "8240011082093172003290000000"; mesonCutArray[ 2] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 3] = "8460011082093172003290000000"; mesonCutArray[ 3] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 4] = "8680011082093172003290000000"; mesonCutArray[ 4] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 5] = "8600011082093172003290000000"; mesonCutArray[ 5] = "01629045009000";  //standard cut Pi0 pPb
   } else if (trainConfig == 2) {
     // doEtaShiftIndCuts = kTRUE;
     // stringShift = "pPb";
      cutarray[ 0] = "8000011072093172003290000000"; mesonCutArray[ 0] = "01627045009000"; 
      cutarray[ 1] = "8020011072093172003290000000"; mesonCutArray[ 1] = "01627045009000"; 
      cutarray[ 2] = "8240011072093172003290000000"; mesonCutArray[ 2] = "01627045009000"; 
      cutarray[ 3] = "8460011072093172003290000000"; mesonCutArray[ 3] = "01627045009000"; 
      cutarray[ 4] = "8680011072093172003290000000"; mesonCutArray[ 4] = "01627045009000"; 
      cutarray[ 5] = "8600011072093172003290000000"; mesonCutArray[ 5] = "01627045009000"; 
   } else if (trainConfig == 3) {
     // doEtaShiftIndCuts = kTRUE;
     //  stringShift = "pPb";
      cutarray[ 0] = "8000011002093172003290000000"; mesonCutArray[ 0] = "01621035009000"; 
      cutarray[ 1] = "8020011002093172003290000000"; mesonCutArray[ 1] = "01621035009000"; 
      cutarray[ 2] = "8240011002093172003290000000"; mesonCutArray[ 2] = "01621035009000"; 
      cutarray[ 3] = "8460011002093172003290000000"; mesonCutArray[ 3] = "01621035009000"; 
      cutarray[ 4] = "8680011002093172003290000000"; mesonCutArray[ 4] = "01621035009000"; 
      cutarray[ 5] = "8600011002093172003290000000"; mesonCutArray[ 5] = "01621035009000"; 
   } else if(trainConfig == 4){
      // Shifting in pPb direction
      doEtaShiftIndCuts = kTRUE;
      stringShift = "pPb";
      cutarray[ 0] = "8000012082093172003290000000"; mesonCutArray[ 0] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 1] = "8020012082093172003290000000"; mesonCutArray[ 1] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 2] = "8240012082093172003290000000"; mesonCutArray[ 2] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 3] = "8460012082093172003290000000"; mesonCutArray[ 3] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 4] = "8680012082093172003290000000"; mesonCutArray[ 4] = "01629045009000";  //standard cut Pi0 pPb
      cutarray[ 5] = "8600012082093172003290000000"; mesonCutArray[ 5] = "01629045009000";  //standard cut Pi0 pPb
   } else if (trainConfig == 5) {
     // doEtaShiftIndCuts = kTRUE;
     // stringShift = "pPb";
      cutarray[ 0] = "8000012072093172003290000000"; mesonCutArray[ 0] = "01627045009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 1] = "8020012072093172003290000000"; mesonCutArray[ 1] = "01627045009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 2] = "8240012072093172003290000000"; mesonCutArray[ 2] = "01627045009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 3] = "8460012072093172003290000000"; mesonCutArray[ 3] = "01627045009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 4] = "8680012072093172003290000000"; mesonCutArray[ 4] = "01627045009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 5] = "8600012072093172003290000000"; mesonCutArray[ 5] = "01627045009000";  //standard cut Pi0 pPb Single pT Cut changed
   } else if (trainConfig == 6) {
     // doEtaShiftIndCuts = kTRUE;
     //  stringShift = "pPb";
      cutarray[ 0] = "8000012002093172003290000000"; mesonCutArray[ 0] = "01621035009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 1] = "8020012002093172003290000000"; mesonCutArray[ 1] = "01621035009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 2] = "8240012002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 3] = "8460012002093172003290000000"; mesonCutArray[ 3] = "01621035009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 4] = "8680012002093172003290000000"; mesonCutArray[ 4] = "01621035009000";  //standard cut Pi0 pPb Single pT Cut changed
      cutarray[ 5] = "8600012002093172003290000000"; mesonCutArray[ 5] = "01621035009000";  //standard cut Pi0 pPb Single pT Cut changed
   } else if (trainConfig == 7) {
     // doEtaShiftIndCuts = kTRUE;
     //  stringShift = "pPb";
      cutarray[ 0] = "8000011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
      cutarray[ 1] = "8000011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
      cutarray[ 2] = "8000011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
      cutarray[ 3] = "8000011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
      cutarray[ 4] = "8000011002092770023220000000"; mesonCutArray[ 4] = "01621035009000";                       
      cutarray[ 5] = "8000011002092551023220000000"; mesonCutArray[ 5] = "01621035009000";                       
   } else if (trainConfig == 8) {
     // doEtaShiftIndCuts = kTRUE;
     //  stringShift = "pPb";
      cutarray[ 0] = "8000012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
      cutarray[ 1] = "8000012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
      cutarray[ 2] = "8000012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
      cutarray[ 3] = "8000012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
      cutarray[ 4] = "8000012002092770023220000000"; mesonCutArray[ 4] = "01621035009000";                       
      cutarray[ 5] = "8000012002092551023220000000"; mesonCutArray[ 5] = "01621035009000";                          
   } else if (trainConfig == 9) {
     // doEtaShiftIndCuts = kTRUE;
     //  stringShift = "pPb";
      cutarray[ 0] = "8000011002092170003220000000"; mesonCutArray[ 0] = "01621035009000"; //just tighten Psi pair
      cutarray[ 1] = "8000011002092170003260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten Psi pair and chi2 in 2D
      cutarray[ 2] = "8000011002092170008220000000"; mesonCutArray[ 2] = "01621035009000"; //tighten psi pair and qt in 2D
      cutarray[ 3] = "8000011002092170008260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
      cutarray[ 4] = "8000011002092170008220400000"; mesonCutArray[ 4] = "01621035009000"; //clean cuts
      cutarray[ 5] = "8000011002092170008260400000"; mesonCutArray[ 5] = "01621035009000"; //clean cuts
      
   } else if (trainConfig == 10) {
     // doEtaShiftIndCuts = kTRUE;
     //  stringShift = "pPb";
      cutarray[ 0] = "8000012002092170003220000000"; mesonCutArray[ 0] = "01621035009000"; //just tighten Psi pair
      cutarray[ 1] = "8000012002092170003260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten Psi pair and chi2 in 2D
      cutarray[ 2] = "8000012002092170008220000000"; mesonCutArray[ 2] = "01621035009000"; //tighten psi pair and qt in 2D
      cutarray[ 3] = "8000012002092170008260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
      cutarray[ 4] = "8000012002092170008220400000"; mesonCutArray[ 4] = "01621035009000"; //clean cuts
      cutarray[ 5] = "8000012002092170008260400000"; mesonCutArray[ 5] = "01621035009000"; //clean cuts
            
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
      if (trainConfig == 1 ||trainConfig == 2 || trainConfig == 3 ){
         if (i == 0 && doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
            }
         }
         if (i == 1 && doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_0020V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_0020V0A", "","Pi0_Fit_Data_pPb_5023GeV_0020V0A","Eta_Fit_Data_pPb_5023GeV_0020V0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_0020V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_0020V0A", "","Pi0_Fit_Data_pPb_5023GeV_0020V0A","Eta_Fit_Data_pPb_5023GeV_0020V0A");
            }
         }
         if (i == 2 && doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_2040V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_2040V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
            }
         }
         if (i == 3 && doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_4060V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_4060V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
            }
         }
         if (i == 4 && doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_6080V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_6080V0A", "","Pi0_Fit_Data_pPb_5023GeV_6080V0A","Eta_Fit_Data_pPb_5023GeV_6080V0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_6080V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_6080V0A", "","Pi0_Fit_Data_pPb_5023GeV_6080V0A","Eta_Fit_Data_pPb_5023GeV_6080V0A");
            }
         }
         if (i == 5 && doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_60100V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_60100V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
            }
         }
      }
      if (trainConfig == 4 ||trainConfig == 5 || trainConfig == 6 ){
         if (i == 0 && doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
           
         }
         if (i == 1 && doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_0020V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_0020V0A", "","Pi0_Fit_Data_pPb_5023GeV_0020V0A","Eta_Fit_Data_pPb_5023GeV_0020V0A");
            
         }
         if (i == 2 && doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
            
         }
         if (i == 3 && doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
         }
         if (i == 4 && doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A", "","Pi0_Fit_Data_pPb_5023GeV_6080V0A","Eta_Fit_Data_pPb_5023GeV_6080V0A");
         }
         if (i == 5 && doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
         }
      }
      if (trainConfig == 7 || trainConfig == 9 ){
         if (doWeighting){
            if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
            } else if (generatorName.CompareTo("HIJING")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
            }
         }
      }   
      if (trainConfig == 8 || trainConfig == 10){
         if (doWeighting){
            analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
         }
         
      }   
         
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
   task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
   task->SetDoPhotonQA(enableQAMesonTask);  //Attention new switch small for Photon QA
 
   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));
 
   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);
 
   return;
 
}
