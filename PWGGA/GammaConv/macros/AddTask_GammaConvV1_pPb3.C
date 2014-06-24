void AddTask_GammaConvV1_pPb3(  Int_t trainConfig = 1,  //change different set of cuts
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
   TString cutnumber = "8000000160084001001500000000"; //offline V0Finder
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
	
     if (trainConfig == 1) {
           cutarray[ 0] = "8000011102092170008260400000"; mesonCutArray[ 0] = "01621035009000";  // all Photon Qualities  //offline V0Finder
           cutarray[ 1] = "8000011102092170008260420000"; mesonCutArray[ 1] = "01621035009000";  // only photon quality 1  //offline V0Finder
           cutarray[ 2] = "8000011102092170008260430000"; mesonCutArray[ 2] = "01621035009000";  // only photon quality 2  //offline V0Finder
           cutarray[ 3] = "8000011102092170008260440000"; mesonCutArray[ 3] = "01621035009000";  // only photon quality 3  //offline V0Finder
     } else if (trainConfig == 2) {
           cutarray[ 0] = "8000012102092170008260400000"; mesonCutArray[ 0] = "01621035009000";  // all Photon Qualities  //offline V0Finder
           cutarray[ 1] = "8000012102092170008260420000"; mesonCutArray[ 1] = "01621035009000";  // only photon quality 1  //offline V0Finder
           cutarray[ 2] = "8000012102092170008260430000"; mesonCutArray[ 2] = "01621035009000";  // only photon quality 2  //offline V0Finder
           cutarray[ 3] = "8000012102092170008260440000"; mesonCutArray[ 3] = "01621035009000";  // only photon quality 3  //offline V0Finder
      
     
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
		if ( trainConfig == 1  ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				} else if (generatorName.CompareTo("HIJING")==0){
					analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				}
			}
		}   
		if ( trainConfig == 2 ){
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
