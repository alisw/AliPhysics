void AddTask_GammaConvV1_pPb(  Int_t trainConfig = 1,  //change different set of cuts
                              Bool_t isMC   = kFALSE, //run MC
                              Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvV1
                              Int_t enableQAPhotonTask = 0, // enable additional QA task
                              TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                              Int_t doWeightingPart = 0,  //enable Weighting
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
   Int_t numberOfCuts = 4;
 
   TString *cutarray = new TString[numberOfCuts];
   TString *mesonCutArray = new TString[numberOfCuts];
   Bool_t doEtaShiftIndCuts = kFALSE;
   TString stringShift = "";
	if(trainConfig == 1){
		cutarray[ 0] = "8000011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8000011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8000011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8000011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 2) {
		cutarray[ 0] = "8000012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8000012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8000012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8000012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 3) {
		cutarray[ 0] = "8000011002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8000011002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                       
		cutarray[ 2] = "8000011002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8000011002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 4) {
		cutarray[ 0] = "8000012002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8000012002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                          
		cutarray[ 2] = "8000012002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8000012002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 5) {	  
		cutarray[ 0] = "8000011002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8000011002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8000011002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8000011002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 6) {
		cutarray[ 0] = "8000012002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8000012002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8000012002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8000012002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 7) {
		cutarray[ 0] = "8000011001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8000011009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8000011002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8000011002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 8) {
		cutarray[ 0] = "8000012001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8000012009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8000012002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8000012002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 9) {
		cutarray[ 0] = "8000011002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8000011002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8000011002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8000011002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 10) {
		cutarray[ 0] = "8000012002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8000012002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8000012002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8000012002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 11) {
		cutarray[ 0] = "8000011002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8000011002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8000011002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8000011002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 12) {
		cutarray[ 0] = "8000012002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8000012002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8000012002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8000012002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 13) {
		cutarray[ 0] = "8000011002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8000011002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8000011002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8000011002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 14) {
		cutarray[ 0] = "8000012002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8000012002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8000012002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8000012002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 15) {
		cutarray[ 0] = "8000011002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8000011002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8000011002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8000011002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 16) {
		cutarray[ 0] = "8000012002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8000012002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8000012002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8000012002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 17) {	
		cutarray[ 0] = "8000011002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8000011002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8000011002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8000011002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 18) {
		cutarray[ 0] = "8000012002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8000012002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8000012002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8000012002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 19) {
		cutarray[ 0] = "8000011002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8000011002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8000011002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8000011002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 20) {
		cutarray[ 0] = "8000012002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8000012002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8000012002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8000012002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 21) {
		cutarray[ 0] = "8000011002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8000011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8000011002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8000011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 22) {
		cutarray[ 0] = "8000012002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8000012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8000012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8000012002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 23){
		cutarray[ 0] = "8020011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8020011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8020011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8020011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 24) {
		cutarray[ 0] = "8020012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8020012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8020012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8020012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 25) {
		cutarray[ 0] = "8020011002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8020011002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                       
		cutarray[ 2] = "8020011002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8020011002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 26) {
		cutarray[ 0] = "8020012002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8020012002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                          
		cutarray[ 2] = "8020012002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8020012002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 27) {	  
		cutarray[ 0] = "8020011002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8020011002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8020011002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8020011002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 28) {
		cutarray[ 0] = "8020012002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8020012002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8020012002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8020012002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 29) {
		cutarray[ 0] = "8020011001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8020011009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8020011002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8020011002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 30) {
		cutarray[ 0] = "8020012001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8020012009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8020012002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8020012002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 31) {
		cutarray[ 0] = "8020011002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8020011002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8020011002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8020011002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 32) {
		cutarray[ 0] = "8020012002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8020012002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8020012002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8020012002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 33) {
		cutarray[ 0] = "8020011002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8020011002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8020011002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8020011002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 34) {
		cutarray[ 0] = "8020012002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8020012002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8020012002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8020012002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 35) {
		cutarray[ 0] = "8020011002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8020011002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8020011002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8020011002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 36) {
		cutarray[ 0] = "8020012002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8020012002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8020012002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8020012002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 37) {
		cutarray[ 0] = "8020011002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8020011002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8020011002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8020011002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 38) {
		cutarray[ 0] = "8020012002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8020012002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8020012002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8020012002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 39) {	
		cutarray[ 0] = "8020011002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8020011002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8020011002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8020011002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 40) {
		cutarray[ 0] = "8020012002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8020012002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8020012002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8020012002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 41) {
		cutarray[ 0] = "8020011002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8020011002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8020011002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8020011002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 42) {
		cutarray[ 0] = "8020012002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8020012002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8020012002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8020012002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 43) {
		cutarray[ 0] = "8020011002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8020011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8020011002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8020011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 44) {
		cutarray[ 0] = "8020012002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8020012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8020012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8020012002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 45){
		cutarray[ 0] = "8240011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8240011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8240011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8240011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 46) {
		cutarray[ 0] = "8240012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8240012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8240012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8240012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 47) {
		cutarray[ 0] = "8240011002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8240011002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                       
		cutarray[ 2] = "8240011002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8240011002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 48) {
		cutarray[ 0] = "8240012002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8240012002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                          
		cutarray[ 2] = "8240012002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8240012002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 49) {	  
		cutarray[ 0] = "8240011002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8240011002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8240011002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8240011002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 50) {
		cutarray[ 0] = "8240012002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8240012002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8240012002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8240012002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 51) {
		cutarray[ 0] = "8240011001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8240011009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8240011002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8240011002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 52) {
		cutarray[ 0] = "8240012001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8240012009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8240012002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8240012002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 53) {
		cutarray[ 0] = "8240011002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8240011002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8240011002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8240011002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 54) {
		cutarray[ 0] = "8240012002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8240012002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8240012002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8240012002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 55) {
		cutarray[ 0] = "8240011002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8240011002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8240011002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8240011002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 56) {
		cutarray[ 0] = "8240012002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8240012002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8240012002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8240012002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 57) {
		cutarray[ 0] = "8240011002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8240011002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8240011002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8240011002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 58) {
		cutarray[ 0] = "8240012002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8240012002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8240012002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8240012002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 59) {
		cutarray[ 0] = "8240011002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8240011002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8240011002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8240011002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 60) {
		cutarray[ 0] = "8240012002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8240012002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8240012002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8240012002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 61) {	
		cutarray[ 0] = "8240011002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8240011002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8240011002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8240011002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 62) {
		cutarray[ 0] = "8240012002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8240012002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8240012002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8240012002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 63) {
		cutarray[ 0] = "8240011002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8240011002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8240011002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8240011002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 64) {
		cutarray[ 0] = "8240012002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8240012002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8240012002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8240012002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 65) {
		cutarray[ 0] = "8240011002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8240011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8240011002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8240011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 66) {
		cutarray[ 0] = "8240012002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8240012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8240012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8240012002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 67){
		cutarray[ 0] = "8460011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8460011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8460011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8460011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 68) {
		cutarray[ 0] = "8460012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8460012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8460012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8460012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 69) {
		cutarray[ 0] = "8460011002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8460011002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                       
		cutarray[ 2] = "8460011002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8460011002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 70) {
		cutarray[ 0] = "8460012002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8460012002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                          
		cutarray[ 2] = "8460012002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8460012002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 71) {	  
		cutarray[ 0] = "8460011002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8460011002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8460011002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8460011002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 72) {
		cutarray[ 0] = "8460012002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8460012002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8460012002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8460012002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 73) {
		cutarray[ 0] = "8460011001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8460011009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8460011002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8460011002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 74) {
		cutarray[ 0] = "8460012001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8460012009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8460012002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8460012002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 75) {
		cutarray[ 0] = "8460011002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8460011002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8460011002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8460011002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 76) {
		cutarray[ 0] = "8460012002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8460012002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8460012002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8460012002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 77) {
		cutarray[ 0] = "8460011002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8460011002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8460011002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8460011002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 78) {
		cutarray[ 0] = "8460012002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8460012002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8460012002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8460012002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 79) {
		cutarray[ 0] = "8460011002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8460011002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8460011002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8460011002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 80) {
		cutarray[ 0] = "8460012002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8460012002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8460012002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8460012002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 81) {
		cutarray[ 0] = "8460011002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8460011002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8460011002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8460011002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 82) {
		cutarray[ 0] = "8460012002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8460012002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8460012002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8460012002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 83) {	
		cutarray[ 0] = "8460011002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8460011002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8460011002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8460011002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 84) {
		cutarray[ 0] = "8460012002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8460012002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8460012002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8460012002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 85) {
		cutarray[ 0] = "8460011002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8460011002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8460011002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8460011002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 86) {
		cutarray[ 0] = "8460012002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8460012002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8460012002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8460012002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 87) {
		cutarray[ 0] = "8460011002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8460011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8460011002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8460011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 88) {
		cutarray[ 0] = "8460012002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8460012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8460012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8460012002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;	
	} else if (trainConfig == 89){
		cutarray[ 0] = "8680011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8680011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8680011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8680011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 90) {
		cutarray[ 0] = "8680012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8680012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8680012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8680012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 91) {
		cutarray[ 0] = "8680011002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8680011002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                       
		cutarray[ 2] = "8680011002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8680011002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 92) {
		cutarray[ 0] = "8680012002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8680012002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                          
		cutarray[ 2] = "8680012002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8680012002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 93) {	  
		cutarray[ 0] = "8680011002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8680011002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8680011002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8680011002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 94) {
		cutarray[ 0] = "8680012002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8680012002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8680012002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8680012002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 95) {
		cutarray[ 0] = "8680011001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8680011009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8680011002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8680011002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 96) {
		cutarray[ 0] = "8680012001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8680012009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8680012002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8680012002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 97) {
		cutarray[ 0] = "8680011002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8680011002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8680011002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8680011002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 98) {
		cutarray[ 0] = "8680012002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8680012002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8680012002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8680012002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 99) {
		cutarray[ 0] = "8680011002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8680011002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8680011002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8680011002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 100) {
		cutarray[ 0] = "8680012002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8680012002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8680012002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8680012002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 101) {
		cutarray[ 0] = "8680011002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8680011002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8680011002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8680011002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 102) {
		cutarray[ 0] = "8680012002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8680012002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8680012002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8680012002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 103) {
		cutarray[ 0] = "8680011002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8680011002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8680011002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8680011002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 104) {
		cutarray[ 0] = "8680012002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8680012002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8680012002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8680012002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 105) {	
		cutarray[ 0] = "8680011002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8680011002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8680011002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8680011002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 106) {
		cutarray[ 0] = "8680012002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8680012002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8680012002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8680012002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 107) {
		cutarray[ 0] = "8680011002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8680011002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8680011002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8680011002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 108) {
		cutarray[ 0] = "8680012002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8680012002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8680012002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8680012002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 109) {
		cutarray[ 0] = "8680011002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8680011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8680011002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8680011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 110) {
		cutarray[ 0] = "8680012002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8680012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8680012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8680012002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;	
	} else if (trainConfig == 111){
		cutarray[ 0] = "8600011002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8600011002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8600011002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8600011002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 112) {
		cutarray[ 0] = "8600012002092970023220000000"; mesonCutArray[ 0] = "01621035009000"; 
		cutarray[ 1] = "8600012002093663003800000000"; mesonCutArray[ 1] = "01621035009000"; 
		cutarray[ 2] = "8600012002093260003800000000"; mesonCutArray[ 2] = "01621035009000"; 
		cutarray[ 3] = "8600012002092370023220000000"; mesonCutArray[ 3] = "01621035009000";                       
	} else if (trainConfig == 113) {
		cutarray[ 0] = "8600011002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8600011002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                       
		cutarray[ 2] = "8600011002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8600011002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 114) {
		cutarray[ 0] = "8600012002092770023220000000"; mesonCutArray[ 0] = "01621035009000";                       
		cutarray[ 1] = "8600012002092551023220000000"; mesonCutArray[ 1] = "01621035009000";                          
		cutarray[ 2] = "8600012002092170003220000000"; mesonCutArray[ 2] = "01621035009000"; //just tighten Psi pair
		cutarray[ 3] = "8600012002092170003260000000"; mesonCutArray[ 3] = "01621035009000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 115) {	  
		cutarray[ 0] = "8600011002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8600011002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8600011002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8600011002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 116) {
		cutarray[ 0] = "8600012002092170008220000000"; mesonCutArray[ 0] = "01621035009000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8600012002092170008260000000"; mesonCutArray[ 1] = "01621035009000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8600012002092170008220400000"; mesonCutArray[ 2] = "01621035009000"; //clean cuts
		cutarray[ 3] = "8600012002092170008260400000"; mesonCutArray[ 3] = "01621035009000"; //clean cuts
	} else if (trainConfig == 117) {
		cutarray[ 0] = "8600011001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8600011009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8600011002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8600011002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 118) {
		cutarray[ 0] = "8600012001092170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8600012009092170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8600012002792170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8600012002012170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 119) {
		cutarray[ 0] = "8600011002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8600011002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8600011002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8600011002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 120) {
		cutarray[ 0] = "8600012002082170008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8600012002062170008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8600012002093170008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8600012002096170008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 121) {
		cutarray[ 0] = "8600011002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8600011002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8600011002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8600011002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 122) {
		cutarray[ 0] = "8600012002092270008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8600012002092570008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8600012002092160008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8600012002092150008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 123) {
		cutarray[ 0] = "8600011002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8600011002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8600011002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8600011002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 124) {
		cutarray[ 0] = "8600012002092172008260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8600012002092162008260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8600012002092260008260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8600012002092262008260400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 125) {
		cutarray[ 0] = "8600011002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8600011002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8600011002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8600011002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 126) {
		cutarray[ 0] = "8600012002092170003260400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8600012002092170009260400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8600012002092170002260400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8600012002092170008220400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 127) {	
		cutarray[ 0] = "8600011002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8600011002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8600011002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8600011002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 128) {
		cutarray[ 0] = "8600012002092170008160400000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8600012002092170008860400000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8600012002092170008250400000"; mesonCutArray[ 2] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8600012002092170008270400000"; mesonCutArray[ 3] = "01621035009000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 129) {
		cutarray[ 0] = "8600011002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8600011002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8600011002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8600011002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 130) {
		cutarray[ 0] = "8600012002092170008260300000"; mesonCutArray[ 0] = "01621035009000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8600012002092170008260600000"; mesonCutArray[ 1] = "01621035009000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8600012002092170008260400000"; mesonCutArray[ 2] = "01621065009000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8600012002092170008260400000"; mesonCutArray[ 3] = "01621034009000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 131) {
		cutarray[ 0] = "8600011002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8600011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8600011002093172003290000000"; mesonCutArray[ 2] = "01621035009000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8600011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 132) {
		cutarray[ 0] = "8600012002092170008260400000"; mesonCutArray[ 0] = "02621035009000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8600012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8600012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8600012002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;	
	} else if (trainConfig == 133) {
		cutarray[ 0] = "8000011002092770023220000000"; mesonCutArray[ 0] = "01621035000000";                       
		cutarray[ 1] = "8000011002092551023220000000"; mesonCutArray[ 1] = "01621035000000";                       
		cutarray[ 2] = "8000011002092170003220000000"; mesonCutArray[ 2] = "01621035000000"; //just tighten Psi pair
		cutarray[ 3] = "8000011002092170003260000000"; mesonCutArray[ 3] = "01621035000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 134) {
		cutarray[ 0] = "8000012002092770023220000000"; mesonCutArray[ 0] = "01621035000000";                       
		cutarray[ 1] = "8000012002092551023220000000"; mesonCutArray[ 1] = "01621035000000";                          
		cutarray[ 2] = "8000012002092170003220000000"; mesonCutArray[ 2] = "01621035000000"; //just tighten Psi pair
		cutarray[ 3] = "8000012002092170003260000000"; mesonCutArray[ 3] = "01621035000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 135) {	  
		cutarray[ 0] = "8000011002092170008220000000"; mesonCutArray[ 0] = "01621035000000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8000011002092170008260000000"; mesonCutArray[ 1] = "01621035000000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8000011002092170008220400000"; mesonCutArray[ 2] = "01621035000000"; //clean cuts
		cutarray[ 3] = "8000011002092170008260400000"; mesonCutArray[ 3] = "01621035000000"; //clean cuts
	} else if (trainConfig == 136) {
		cutarray[ 0] = "8000012002092170008220000000"; mesonCutArray[ 0] = "01621035000000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8000012002092170008260000000"; mesonCutArray[ 1] = "01621035000000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8000012002092170008220400000"; mesonCutArray[ 2] = "01621035000000"; //clean cuts
		cutarray[ 3] = "8000012002092170008260400000"; mesonCutArray[ 3] = "01621035000000"; //clean cuts
	} else if (trainConfig == 137) {
		cutarray[ 0] = "8000011001092170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8000011009092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8000011002792170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8000011002012170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 138) {
		cutarray[ 0] = "8000012001092170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8000012009092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8000012002792170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8000012002012170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 139) {
		cutarray[ 0] = "8000011002082170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8000011002062170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8000011002093170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8000011002096170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 140) {
		cutarray[ 0] = "8000012002082170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8000012002062170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8000012002093170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8000012002096170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 141) {
		cutarray[ 0] = "8000011002092270008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8000011002092570008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8000011002092160008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8000011002092150008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 142) {
		cutarray[ 0] = "8000012002092270008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8000012002092570008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8000012002092160008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8000012002092150008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 143) {
		cutarray[ 0] = "8000011002092172008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8000011002092162008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8000011002092260008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8000011002092262008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 144) {
		cutarray[ 0] = "8000012002092172008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8000012002092162008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8000012002092260008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8000012002092262008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 145) {
		cutarray[ 0] = "8000011002092170003260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8000011002092170009260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8000011002092170002260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8000011002092170008220400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 146) {
		cutarray[ 0] = "8000012002092170003260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8000012002092170009260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8000012002092170002260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8000012002092170008220400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 147) {	
		cutarray[ 0] = "8000011002092170008160400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8000011002092170008860400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8000011002092170008250400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8000011002092170008270400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 148) {
		cutarray[ 0] = "8000012002092170008160400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8000012002092170008860400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8000012002092170008250400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8000012002092170008270400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 149) {
		cutarray[ 0] = "8000011002092170008260300000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8000011002092170008260600000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8000011002092170008260400000"; mesonCutArray[ 2] = "01621065000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8000011002092170008260400000"; mesonCutArray[ 3] = "01621034000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 150) {
		cutarray[ 0] = "8000012002092170008260300000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8000012002092170008260600000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8000012002092170008260400000"; mesonCutArray[ 2] = "01621065000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8000012002092170008260400000"; mesonCutArray[ 3] = "01621034000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 151) {
		cutarray[ 0] = "8000011002092170008260400000"; mesonCutArray[ 0] = "02621035000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8000011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8000011002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8000011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 152) {
		cutarray[ 0] = "8000012002092170008260400000"; mesonCutArray[ 0] = "02621035000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8000012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8000012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8000012002092170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 153) {
		cutarray[ 0] = "8020011002092770023220000000"; mesonCutArray[ 0] = "01621035000000";                       
		cutarray[ 1] = "8020011002092551023220000000"; mesonCutArray[ 1] = "01621035000000";                       
		cutarray[ 2] = "8020011002092170003220000000"; mesonCutArray[ 2] = "01621035000000"; //just tighten Psi pair
		cutarray[ 3] = "8020011002092170003260000000"; mesonCutArray[ 3] = "01621035000000"; //tighten Psi pair and chi2 in 2D	
	} else if (trainConfig == 154) {
		cutarray[ 0] = "8020012002092770023220000000"; mesonCutArray[ 0] = "01621035000000";                       
		cutarray[ 1] = "8020012002092551023220000000"; mesonCutArray[ 1] = "01621035000000";                          
		cutarray[ 2] = "8020012002092170003220000000"; mesonCutArray[ 2] = "01621035000000"; //just tighten Psi pair
		cutarray[ 3] = "8020012002092170003260000000"; mesonCutArray[ 3] = "01621035000000"; //tighten Psi pair and chi2 in 2D
	} else if (trainConfig == 155) {	  
		cutarray[ 0] = "8020011002092170008220000000"; mesonCutArray[ 0] = "01621035000000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8020011002092170008260000000"; mesonCutArray[ 1] = "01621035000000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8020011002092170008220400000"; mesonCutArray[ 2] = "01621035000000"; //clean cuts
		cutarray[ 3] = "8020011002092170008260400000"; mesonCutArray[ 3] = "01621035000000"; //clean cuts
	} else if (trainConfig == 156) {
		cutarray[ 0] = "8020012002092170008220000000"; mesonCutArray[ 0] = "01621035000000"; //tighten psi pair and qt in 2D
		cutarray[ 1] = "8020012002092170008260000000"; mesonCutArray[ 1] = "01621035000000"; //tighten psi pair and chi2 in 2D and qt in 2D                      
		cutarray[ 2] = "8020012002092170008220400000"; mesonCutArray[ 2] = "01621035000000"; //clean cuts
		cutarray[ 3] = "8020012002092170008260400000"; mesonCutArray[ 3] = "01621035000000"; //clean cuts
	} else if (trainConfig == 157) {
		cutarray[ 0] = "8020011001092170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8020011009092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8020011002792170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8020011002012170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 158) {
		cutarray[ 0] = "8020012001092170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8  // minR 2.8
		cutarray[ 1] = "8020012009092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8  //minR 7.5
		cutarray[ 2] = "8020012002792170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
		cutarray[ 3] = "8020012002012170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
	} else if (trainConfig == 159) {
		cutarray[ 0] = "8020011002082170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8020011002062170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7      
		cutarray[ 2] = "8020011002093170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8020011002096170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 160) {
		cutarray[ 0] = "8020012002082170008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
		cutarray[ 1] = "8020012002062170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
		cutarray[ 2] = "8020012002093170008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -4,5
		cutarray[ 3] = "8020012002096170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
	} else if (trainConfig == 161) {
		cutarray[ 0] = "8020011002092270008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8020011002092570008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8020011002092160008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8020011002092150008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
	} else if (trainConfig == 162) {
		cutarray[ 0] = "8020012002092270008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
		cutarray[ 1] = "8020012002092570008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
		cutarray[ 2] = "8020012002092160008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi minMom 0.25
		cutarray[ 3] = "8020012002092150008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi min Mom 0.3		
	} else if (trainConfig == 163) {
		cutarray[ 0] = "8020011002092172008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8020011002092162008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8020011002092260008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8020011002092262008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 164) {
		cutarray[ 0] = "8020012002092172008260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
		cutarray[ 1] = "8020012002092162008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
		cutarray[ 2] = "8020012002092260008260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
		cutarray[ 3] = "8020012002092262008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
	} else if (trainConfig == 165) {
		cutarray[ 0] = "8020011002092170003260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8020011002092170009260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8020011002092170002260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8020011002092170008220400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 166) {
		cutarray[ 0] = "8020012002092170003260400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.05 1D
		cutarray[ 1] = "8020012002092170009260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.03 2D
		cutarray[ 2] = "8020012002092170002260400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //qT 0.07 1D
		cutarray[ 3] = "8020012002092170008220400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D 
	} else if (trainConfig == 167) {	
		cutarray[ 0] = "8020011002092170008160400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8020011002092170008860400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8020011002092170008250400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8020011002092170008270400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035	
	} else if (trainConfig == 168) {
		cutarray[ 0] = "8020012002092170008160400000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 50  2D
		cutarray[ 1] = "8020012002092170008860400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //chi2 20  2D
		cutarray[ 2] = "8020012002092170008250400000"; mesonCutArray[ 2] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.1
		cutarray[ 3] = "8020012002092170008270400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 //psi pair 0.035
	} else if (trainConfig == 169) {
		cutarray[ 0] = "8020011002092170008260300000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8020011002092170008260600000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8020011002092170008260400000"; mesonCutArray[ 2] = "01621065000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8020011002092170008260400000"; mesonCutArray[ 3] = "01621034000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500
	} else if (trainConfig == 170) {
		cutarray[ 0] = "8020012002092170008260300000"; mesonCutArray[ 0] = "01621035000000";  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
		cutarray[ 1] = "8020012002092170008260600000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
		cutarray[ 2] = "8020012002092170008260400000"; mesonCutArray[ 2] = "01621065000000";  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8 
		cutarray[ 3] = "8020012002092170008260400000"; mesonCutArray[ 3] = "01621034000000";  //new standard eta=0.9 y=0.8 //chi2 meson 500		
	} else if (trainConfig == 171) {
		cutarray[ 0] = "8020011002092170008260400000"; mesonCutArray[ 0] = "02621035000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8020011002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing
		cutarray[ 2] = "8020011002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8020011002092170008260400000"; mesonCutArray[ 3] = "01621035008000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
	} else if (trainConfig == 172) {
		cutarray[ 0] = "8020012002092170008260400000"; mesonCutArray[ 0] = "02621035000000";  //new standard eta=0.9 y=0.8 // BG track multiplicity
		cutarray[ 1] = "8020012002092170008260400000"; mesonCutArray[ 1] = "01621035000000";  //new standard eta=0.9 y=0.8 //no MCP smearing      
		cutarray[ 2] = "8020012002093172003290000000"; mesonCutArray[ 2] = "01621035000000";  //old standard eta=0.9 y=0.8 
		cutarray[ 3] = "8020012002092170008260400000"; mesonCutArray[ 3] = "01621035000000";  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;		
	} else {
		Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
		return;
	}
 
	TList *ConvCutList = new TList();
	TList *MesonCutList = new TList();
	
	TList *HeaderList = new TList();
	if (doWeightingPart==1) {
		TObjString *Header1 = new TObjString("pi0_1");
		HeaderList->Add(Header1);
	}
	if (doWeightingPart==2){
		TObjString *Header3 = new TObjString("eta_2");
		HeaderList->Add(Header3);
	}
	if (doWeightingPart==3) {
		TObjString *Header1 = new TObjString("pi0_1");
		HeaderList->Add(Header1);
		TObjString *Header3 = new TObjString("eta_2");
		HeaderList->Add(Header3);
	}
	
	
	Bool_t doWeighting = kFALSE;
	if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;
	
	ConvCutList->SetOwner(kTRUE);
	AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
	
	
	for(Int_t i = 0; i<numberOfCuts; i++){
		
		analysisCuts[i] = new AliConversionCuts();
		if ( trainConfig == 1 || trainConfig == 3 || trainConfig == 5 || trainConfig == 7 || trainConfig == 9 || trainConfig == 11 || trainConfig == 13 || trainConfig == 15|| trainConfig == 17|| trainConfig == 19 || trainConfig == 21 || trainConfig == 133 || trainConfig == 135 || trainConfig == 137 || trainConfig == 139 || trainConfig == 141 || trainConfig == 143 || trainConfig == 145 || trainConfig == 147 || trainConfig == 149 || trainConfig == 151){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				}
			}
		}   
		if ( trainConfig == 2 || trainConfig == 4 || trainConfig == 6 || trainConfig == 8 || trainConfig == 10 || trainConfig == 12 || trainConfig == 14 || trainConfig == 16|| trainConfig == 18|| trainConfig == 20|| trainConfig == 22 || trainConfig == 134 || trainConfig == 136 || trainConfig == 138 || trainConfig == 140 || trainConfig == 142 || trainConfig == 144 || trainConfig == 146 || trainConfig == 148 || trainConfig == 150 || trainConfig == 152){
			if (doWeighting){
				analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
			}
			
		}   
		if ( trainConfig == 23 || trainConfig == 25 || trainConfig == 27 || trainConfig == 29 || trainConfig == 31 || trainConfig == 33 || trainConfig == 35 || trainConfig == 37|| trainConfig == 39|| trainConfig == 41 || trainConfig == 43 || trainConfig == 153 || trainConfig == 155 || trainConfig == 157 || trainConfig == 159 || trainConfig == 161 || trainConfig == 163 || trainConfig == 165 || trainConfig == 167 || trainConfig == 169 || trainConfig == 171){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_0020V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_0020V0A", "","Pi0_Fit_Data_pPb_5023GeV_0020V0A","Eta_Fit_Data_pPb_5023GeV_0020V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_0020V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_0020V0A", "","Pi0_Fit_Data_pPb_5023GeV_0020V0A","Eta_Fit_Data_pPb_5023GeV_0020V0A");
				}
			}
		}   
		if ( trainConfig == 24 || trainConfig == 26 || trainConfig == 28 || trainConfig == 30 || trainConfig == 32 || trainConfig == 34 || trainConfig == 36 || trainConfig == 38|| trainConfig == 40|| trainConfig == 42|| trainConfig == 44 || trainConfig == 154 || trainConfig == 156 || trainConfig == 158 || trainConfig == 160 || trainConfig == 162 || trainConfig == 164 || trainConfig == 166 || trainConfig == 168 || trainConfig == 170 || trainConfig == 172){
			if (doWeighting){
				analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_0020V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_0020V0A", "","Pi0_Fit_Data_pPb_5023GeV_0020V0A","Eta_Fit_Data_pPb_5023GeV_0020V0A");
			}
		}   
		if ( trainConfig == 45 || trainConfig == 47 || trainConfig == 49 || trainConfig == 51 || trainConfig == 53 || trainConfig == 55 || trainConfig == 57 || trainConfig == 59|| trainConfig == 61|| trainConfig == 63 || trainConfig == 65 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_2040V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_2040V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
				}
			}
		}   
		if ( trainConfig == 46 || trainConfig == 48 || trainConfig == 50 || trainConfig == 52 || trainConfig == 54 || trainConfig == 56 || trainConfig == 58 || trainConfig == 60|| trainConfig == 62|| trainConfig == 64|| trainConfig == 66){
			if (doWeighting){
				analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A", "","Pi0_Fit_Data_pPb_5023GeV_2040V0A","Eta_Fit_Data_pPb_5023GeV_2040V0A");
			}
		}   
		if ( trainConfig == 67 || trainConfig == 69 || trainConfig == 71 || trainConfig == 73 || trainConfig == 75 || trainConfig == 77 || trainConfig == 79 || trainConfig == 81 || trainConfig == 83 || trainConfig == 85 || trainConfig == 87 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_4060V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_4060V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
				}
			}
		}   
		if ( trainConfig == 68 || trainConfig == 70 || trainConfig == 72 || trainConfig == 74 || trainConfig == 76 || trainConfig == 78 || trainConfig == 80 || trainConfig == 82 || trainConfig == 84 || trainConfig == 86 || trainConfig == 88){
			if (doWeighting){
				analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A", "","Pi0_Fit_Data_pPb_5023GeV_4060V0A","Eta_Fit_Data_pPb_5023GeV_4060V0A");
			}
		}
		if ( trainConfig == 89 || trainConfig == 91 || trainConfig == 93 || trainConfig == 95 || trainConfig == 97 || trainConfig == 99 || trainConfig == 101 || trainConfig == 103 || trainConfig == 105 || trainConfig == 107 || trainConfig == 109 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_6080V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_6080V0A", "","Pi0_Fit_Data_pPb_5023GeV_6080V0A","Eta_Fit_Data_pPb_5023GeV_6080V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_6080V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_6080V0A", "","Pi0_Fit_Data_pPb_5023GeV_6080V0A","Eta_Fit_Data_pPb_5023GeV_6080V0A");
				}
			}
		}   
		if ( trainConfig == 90 || trainConfig == 92 || trainConfig == 94 || trainConfig == 96 || trainConfig == 98 || trainConfig == 100 || trainConfig == 102 || trainConfig == 104 || trainConfig == 106 || trainConfig == 108 || trainConfig == 108 ){
			if (doWeighting){
				analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A", "","Pi0_Fit_Data_pPb_5023GeV_6080V0A","Eta_Fit_Data_pPb_5023GeV_6080V0A");
			}
		}
		if ( trainConfig == 111 || trainConfig == 113 || trainConfig == 115 || trainConfig == 117 || trainConfig == 119 || trainConfig == 121 || trainConfig == 123 || trainConfig == 125 || trainConfig == 127 || trainConfig == 129 || trainConfig == 131 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_60100V0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_60100V0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
				}
			}
		}   
		if ( trainConfig == 112 || trainConfig == 114 || trainConfig == 116 || trainConfig == 118 || trainConfig == 120 || trainConfig == 122 || trainConfig == 124 || trainConfig == 126 || trainConfig == 128 || trainConfig == 130 || trainConfig == 130){
			if (doWeighting){
				analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A", "","Pi0_Fit_Data_pPb_5023GeV_60100V0A","Eta_Fit_Data_pPb_5023GeV_60100V0A");
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
