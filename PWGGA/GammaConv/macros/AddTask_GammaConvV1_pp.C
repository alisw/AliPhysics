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
      cutarray[ 3] = "0000011002093260003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Gamma pp 2-76TeV
   } else if (trainConfig == 3) {
      cutarray[ 0] = "0002011002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
      cutarray[ 1] = "0003011002093663003800000000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
      cutarray[ 2] = "0002012002093663003800000000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
      cutarray[ 3] = "0003012002093663003800000000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
   } else if (trainConfig == 4) {
      cutarray[ 0] = "0000011002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities
      cutarray[ 1] = "0000011002093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1
      cutarray[ 2] = "0000011002093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2
      cutarray[ 3] = "0000011002093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3
   } else if (trainConfig == 5) {
      cutarray[ 0] = "0000011007093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities, min R = 35 cm
      cutarray[ 1] = "0000011007093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1, min R = 35 cm
      cutarray[ 2] = "0000011007093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2, min R = 35 cm
	  cutarray[ 3] = "0000011007093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3, min R = 35 cm
   } else if (trainConfig == 6) {
      cutarray[ 0] = "0000011002083663003200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities
      cutarray[ 1] = "0000011002083663003200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1
      cutarray[ 2] = "0000011002083663003200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2
	  cutarray[ 3] = "0000011002083663003200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3   
   } else if (trainConfig == 7) {
      cutarray[ 0] = "0000011007083663003200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities, min R = 35 cm
      cutarray[ 1] = "0000011007083663003200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1, min R = 35 cm
      cutarray[ 2] = "0000011007083663003200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2, min R = 35 cm
	  cutarray[ 3] = "0000011007083663003200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3, min R = 35 cm
   } else if (trainConfig == 8) {
      cutarray[ 0] = "0000011002083663000200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, all photon qualities
      cutarray[ 1] = "0000011002083663000200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 1
      cutarray[ 2] = "0000011002083663000200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 2
	  cutarray[ 3] = "0000011002083663000200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 3   
   } else if (trainConfig == 9) {
      cutarray[ 0] = "0000011007083663000200000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 7TeV, all photon qualities, min R = 35 cm
      cutarray[ 1] = "0000011007083663000200020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 1, min R = 35 cm
      cutarray[ 2] = "0000011007083663000200030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 2, min R = 35 cm
	  cutarray[ 3] = "0000011007083663000200040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 7TeV, photon quality 3, min R = 35 cm	   
   } else if (trainConfig == 10) {
      cutarray[ 0] = "0002011002093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities
      cutarray[ 1] = "0002011002093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1
      cutarray[ 2] = "0002011002093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2
      cutarray[ 3] = "0002011002093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3
   } else if (trainConfig == 11) {
      cutarray[ 0] = "0002011007093663003800000000"; mesonCutArray[0] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities, min R = 35 cm
      cutarray[ 1] = "0002011007093663003800020000"; mesonCutArray[1] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1, min R = 35 cm
      cutarray[ 2] = "0002011007093663003800030000"; mesonCutArray[2] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2, min R = 35 cm
      cutarray[ 3] = "0002011007093663003800040000"; mesonCutArray[3] = "01631031009000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3, min R = 35 cm
   } else if (trainConfig == 12) {
      cutarray[ 0] = "0000011002092970028250400000"; mesonCutArray[0] = "01525065000000"; //standard cut LHC11h pp 2.76TeV 
      cutarray[ 1] = "0000011032092970028250400000"; mesonCutArray[1] = "01525065000000"; //variation eta 0.65
      cutarray[ 2] = "0000011042092970028250400000"; mesonCutArray[2] = "01525065000000"; //variation eta 0.75
      cutarray[ 3] = "0000011002092950028250400000"; mesonCutArray[3] = "01525065000000"; //variation pion p dEdx 0.3-5.
   } else if (trainConfig == 13) { //added signals
      cutarray[ 0] = "0000012002092970028250400000"; mesonCutArray[0] = "01525065000000"; //standard cut LHC11h pp 2.76TeV 
      cutarray[ 1] = "0000012032092970028250400000"; mesonCutArray[1] = "01525065000000"; //variation eta 0.65
      cutarray[ 2] = "0000012042092970028250400000"; mesonCutArray[2] = "01525065000000"; //variation eta 0.75
      cutarray[ 3] = "0000012002092950028250400000"; mesonCutArray[3] = "01525065000000"; //variation pion p dEdx 0.3-5.
	} else if (trainConfig == 14) {
      cutarray[ 0] = "0000011002492970028250400000"; mesonCutArray[0] = "01525065000000"; //variation pt 0.075 
      cutarray[ 1] = "0000011002192970028250400000"; mesonCutArray[1] = "01525065000000"; //variation pt 0.1
      cutarray[ 2] = "0000011002062970028250400000"; mesonCutArray[2] = "01525065000000"; //variation TPC cls 0.7
      cutarray[ 3] = "0000011002082970028250400000"; mesonCutArray[3] = "01525065000000"; //variation TPC cls 0.35 
	} else if (trainConfig == 15) { //added signals
      cutarray[ 0] = "0000012002492970028250400000"; mesonCutArray[0] = "01525065000000"; //variation pt 0.075 
      cutarray[ 1] = "0000012002192970028250400000"; mesonCutArray[1] = "01525065000000"; //variation pt 0.1
      cutarray[ 2] = "0000012002062970028250400000"; mesonCutArray[2] = "01525065000000"; //variation TPC cls 0.7
      cutarray[ 3] = "0000012002082970028250400000"; mesonCutArray[3] = "01525065000000"; //variation TPC cls 0.35 
	} else if (trainConfig == 16) {
      cutarray[ 0] = "0000011002093970028250400000"; mesonCutArray[0] = "01525065000000"; //variation edEdx -4,5
      cutarray[ 1] = "0000011002096970028250400000"; mesonCutArray[1] = "01525065000000"; //variation edEdx -2.5,4
      cutarray[ 2] = "0000011002092970038250400000"; mesonCutArray[2] = "01525065000000"; //variation TOF el. PID -3,5
      cutarray[ 3] = "0000011002092970048250400000"; mesonCutArray[3] = "01525065000000"; //variation TOF el. PID -2,3
	} else if (trainConfig == 17) { //added signals
      cutarray[ 0] = "0000012002093970028250400000"; mesonCutArray[0] = "01525065000000"; //variation edEdx -4,5
      cutarray[ 1] = "0000012002096970028250400000"; mesonCutArray[1] = "01525065000000"; //variation edEdx -2.5,4
      cutarray[ 2] = "0000012002092970038250400000"; mesonCutArray[2] = "01525065000000"; //variation TOF el. PID -3,5
      cutarray[ 3] = "0000012002092970048250400000"; mesonCutArray[3] = "01525065000000"; //variation TOF el. PID -2,3
	} else if (trainConfig == 18) {
      cutarray[ 0] = "0000011002092970029250400000"; mesonCutArray[0] = "01525065000000"; //variation qt 0.03
      cutarray[ 1] = "0000011002092970022250400000"; mesonCutArray[1] = "01525065000000"; //variation qt 0.07 no2D
      cutarray[ 2] = "0000011002092970028150400000"; mesonCutArray[2] = "01525065000000"; //variation chi2 50.
      cutarray[ 3] = "0000011002092970028850400000"; mesonCutArray[3] = "01525065000000"; //variation chi2 20.
	} else if (trainConfig == 19) { //added signals
      cutarray[ 0] = "0000012002092970029250400000"; mesonCutArray[0] = "01525065000000"; //variation qt 0.03
      cutarray[ 1] = "0000012002092970022250400000"; mesonCutArray[1] = "01525065000000"; //variation qt 0.07 no2D
      cutarray[ 2] = "0000012002092970028150400000"; mesonCutArray[2] = "01525065000000"; //variation chi2 50.
      cutarray[ 3] = "0000012002092970028850400000"; mesonCutArray[3] = "01525065000000"; //variation chi2 20.
	} else if (trainConfig == 20) {
      cutarray[ 0] = "0000011002092970028260400000"; mesonCutArray[0] = "01525065000000"; //variation psi pair 0.05
      cutarray[ 1] = "0000011002092970028280400000"; mesonCutArray[1] = "01525065000000"; //variation psi pair 0.2
      cutarray[ 2] = "0000011002092970028250000000"; mesonCutArray[2] = "01525065000000"; //variation cosPA -1
      cutarray[ 3] = "0000011002092970028250400000"; mesonCutArray[3] = "01525055000000"; //variation alpha 0.75
	} else if (trainConfig == 21) { //added signals
      cutarray[ 0] = "0000012002092970028260400000"; mesonCutArray[0] = "01525065000000"; //variation psi pair 0.05
      cutarray[ 1] = "0000012002092970028280400000"; mesonCutArray[1] = "01525065000000"; //variation psi pair 0.2
      cutarray[ 2] = "0000012002092970028250000000"; mesonCutArray[2] = "01525065000000"; //variation cosPA -1
      cutarray[ 3] = "0000012002092970028250400000"; mesonCutArray[3] = "01525055000000"; //variation alpha 0.75
	} else if (trainConfig == 22) {
      cutarray[ 0] = "0004011002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kTRD
      cutarray[ 1] = "0005011002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMC
      cutarray[ 2] = "0006011002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kPHI
      cutarray[ 3] = "0007011002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kHighMult
	} else if (trainConfig == 23) {
      cutarray[ 0] = "0008011002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEGA
      cutarray[ 1] = "0009011002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEJE
      cutarray[ 2] = "0000011002092970028250400000"; mesonCutArray[2] = "01525065000000"; // minimum bias
      cutarray[ 3] = "0001111002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kINT8
	} else if (trainConfig == 24) {
      cutarray[ 0] = "0004211002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kTRD CINT8 HEE
      cutarray[ 1] = "0004411002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kTRD CINT8 HSE
      cutarray[ 2] = "0004611002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kTRD CINT8 HJE
      cutarray[ 3] = "0004811002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kTRD CINT8 HQU
	} else if (trainConfig == 25) {
      cutarray[ 0] = "0004111002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kTRD CINT7 HEE
      cutarray[ 1] = "0004311002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kTRD CINT7 HSE
      cutarray[ 2] = "0004511002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kTRD CINT7 HJE
      cutarray[ 3] = "0004711002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kTRD CINT7 HQU
	} else if (trainConfig == 26) {
      cutarray[ 0] = "0005211002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMC7
      cutarray[ 1] = "0005311002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMC8
      cutarray[ 2] = "0006211002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kPHI7
      cutarray[ 3] = "0006311002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kPHI8
	} else if (trainConfig == 27) {
      cutarray[ 0] = "0005111002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMC1
      cutarray[ 1] = "0007111002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kSHM1
      cutarray[ 2] = "0007211002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kSHM7
      cutarray[ 3] = "0007311002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kSHM8
	} else if (trainConfig == 28) {
      cutarray[ 0] = "0008111002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEGA + CINT7
      cutarray[ 1] = "0008211002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEGA + CINT8
      cutarray[ 2] = "0008311002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kEMCEG1 + CINT7
      cutarray[ 3] = "0008411002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kEMCEG1 + CINT8
	} else if (trainConfig == 29) {
      cutarray[ 0] = "0008511002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEG2 + CINT7
      cutarray[ 1] = "0008611002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEG2 + CINT8
      cutarray[ 2] = "0009111002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kEMCEJE + CINT7
      cutarray[ 3] = "0009211002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kEMCEJE + CINT8
	} else if (trainConfig == 30) {
      cutarray[ 0] = "0009311002092970028250400000"; mesonCutArray[0] = "01525065000000"; // trigger kEMCEJ1 + CINT7
      cutarray[ 1] = "0009411002092970028250400000"; mesonCutArray[1] = "01525065000000"; // trigger kEMCEJ1 + CINT8
      cutarray[ 2] = "0009511002092970028250400000"; mesonCutArray[2] = "01525065000000"; // trigger kEMCEJ2 + CINT7
      cutarray[ 3] = "0009611002092970028250400000"; mesonCutArray[3] = "01525065000000"; // trigger kEMCEJ2 + CINT8		
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
