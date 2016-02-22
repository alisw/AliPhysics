void AddTask_GammaConvDalitzV1_pPb(    	Int_t trainConfig = 1,
										Bool_t isMC       = kFALSE, //run MC 
										Int_t enableQAMesonTask = 0, //enable QA in AliAnalysisTaskGammaConvDalitzV1
										Bool_t enableDoMesonChic = kFALSE, // enable additional Chic analysis
										TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
										Bool_t doWeighting = kFALSE,  //enable Weighting
										TString generatorName = "DPMJET",				
										TString cutnumberAODBranch = "0000000060084001001500000",
										Bool_t enableV0findingEffi = kFALSE
                                  ) {
	cout<<"*********Parameters*******"<<endl;
	cout<<"trainConfig: "<<trainConfig<<endl;
	cout<<"isMC: "<<isMC<<endl;
	cout<<"enableQAMesonTask: "<<enableQAMesonTask<<endl;
	cout<<"enableDoMesonChic: "<<enableDoMesonChic<<endl;
	cout<<"fileNameInputForWeighting: "<<fileNameInputForWeighting.Data()<<endl;
	cout<<"doWeighting: "<<doWeighting<<endl;
	cout<<"generatorName: "<<generatorName.Data()<<endl;
	cout<<"cutnumberAODBranch: "<<cutnumberAODBranch.Data()<<endl;

	// ================= Load Librariers =================================
	gSystem->Load("libCore");  
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libPhysics");
	gSystem->Load("libMinuit");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");  
	gSystem->Load("libPWGGAGammaConv");
	gSystem->Load("libCDB");
	gSystem->Load("libSTEER");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libTender");
	gSystem->Load("libTenderSupplies");
	
	Int_t isHeavyIon = 2;

	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error(Form("AddTask_GammaConvDalitzV1_pPb_%i",trainConfig), "No analysis manager found.");
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
	TString cutnumberEvent  = "80000103";
	TString cutnumberPhoton = "06000008000100007500000000";   //Online  V0 finder                         
	TString ElecCuts        = "30105400000003300000";            //Electron Cuts

       
				 
	Bool_t doEtaShift = kFALSE;
	
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	
	//========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
	if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
		
	        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
		fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
		fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
		fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
		fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);
		
		if (!mgr) {
			Error("AddTask_V0ReaderV1", "No analysis manager found.");
			return;
		}
	
		AliConvEventCuts *fEventCuts=NULL;
		if(cutnumberEvent!=""){
			fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
			fEventCuts->SetPreSelectionCutFlag(kTRUE);
			if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
				fEventCuts->DoEtaShift(doEtaShift);
				fV0ReaderV1->SetEventCuts(fEventCuts);
				fEventCuts->SetFillCutHistograms("",kTRUE);
			}
		}

		// Set AnalysisCut Number
		AliConversionPhotonCuts *fCuts=NULL;
		if(cutnumberPhoton!=""){
			fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
			fCuts->SetPreSelectionCutFlag(kTRUE);
			fCuts->SetIsHeavyIon(isHeavyIon);
			if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
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
	//========= Add Electron Selector ================


	if( !(AliDalitzElectronSelector*)mgr->GetTask("ElectronSelector") ){

		AliDalitzElectronSelector *fElectronSelector = new AliDalitzElectronSelector("ElectronSelector");
		// Set AnalysisCut Number
		AliDalitzElectronCuts *fElecCuts=0;	
		if( ElecCuts!=""){
			fElecCuts= new AliDalitzElectronCuts(ElecCuts.Data(),ElecCuts.Data());
			fElecCuts->SetUseCrossedRows(kTRUE);
			if(fElecCuts->InitializeCutsFromCutString(ElecCuts.Data())){
				fElectronSelector->SetDalitzElectronCuts(fElecCuts);
				fElecCuts->SetFillCutHistograms("",kTRUE);
			}
		}

		fElectronSelector->Init();
		mgr->AddTask(fElectronSelector);
		
		AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

		//connect input V0Reader
		mgr->ConnectInput (fElectronSelector,0,cinput1);
	}

	//================================================
	//========= Add task to the ANALYSIS manager =====
	//================================================
	//            find input container
	AliAnalysisTaskGammaConvDalitzV1 *task=NULL;

	task= new AliAnalysisTaskGammaConvDalitzV1(Form("GammaConvDalitzV1_%i",trainConfig));
	task->SetIsHeavyIon(isHeavyIon);
	task->SetIsMC(isMC);

	// Cut Numbers to use in Analysis
	Int_t numberOfCuts = 4;

	TString *eventCutArray   = new TString[numberOfCuts];
	TString *photonCutArray  = new TString[numberOfCuts];
	TString *ElecCutarray    = new TString[numberOfCuts];
	TString *MesonCutarray   = new TString[numberOfCuts];

	Bool_t doEtaShiftIndCuts = kFALSE;
	TString stringShift = "";

	// Shifting in pPb direction
	doEtaShiftIndCuts = kFALSE;
	stringShift = "pPb";
	
	if( trainConfig == 1 ) {  // No eta shift |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007900000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 15
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007800000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 20
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007100000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 50
	}  else if( trainConfig == 2 ) {  // No eta shift |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300002200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.7
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300003200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.5
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009365300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360100007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
	}  else if( trainConfig == 3 ) {  // No eta shift |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009380300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec  sec  2.0sigmas Low and  1 High momentum
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20435400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 2.0sigmas Low and 0 High momentum
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20477400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475200233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum
	}  else if( trainConfig == 4 ) {  // No eta shift  |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20425400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx pion rejec  primary  2.0sigmas Low and -1 High momentum 
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400133102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + SPD first layer
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233302223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + PsiPair cut 0.52
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263100500900000"; //standard cut Pi0 pPb 00-100 Standard cut + Alpha cut < 0.7	
	}  else if( trainConfig == 5 ) { // No eta shift  |Y| < 0.8	
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "90375400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -5,5
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "90575400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -3,5
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009160300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -5,5
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009260300007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -3,5
	} else if ( trainConfig == 6 ) { //No eta shift   |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "04200009360300007200000000"; ElecCutarray[0] = "20475400235102223710"; MesonCutarray[0] = "0263203500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Y < 0.70  and prim and sec e |eta| < 0.75 //NOTE revisar
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233102233710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single prim Pt cut > 0.150 
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233102231710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + Single prim Pt cut > 0.100 
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102222710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + DCAxy < 1 cm
	} else if ( trainConfig == 7 ) {  // No eta shift |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200049360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.075
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200019360300007200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.100
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200008360300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls sec  > 0.35
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400273102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls prim > 0.60
	} else if ( trainConfig == 8 ) {  //No eta shift |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223810"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.015 < InvMass(e+,e-) < 0.050
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233102223910"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.025 < InvMass(e+,e-) < 0.035
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300001200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + qT < 0.1
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	} else if ( trainConfig == 9 ) {  //No eta shift |Y| < 0.8
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0273103500900000"; //standard cut Pi0 pPb 00-100  New Stardad cut +100 events background 
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0163103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Background method V0 multiplicity 
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500000000"; //standard cut Pi0 PbPb 00-100 + No extra smearing
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400253102221710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100 + Old Standard
	} else if( trainConfig == 10 ) {  // No eta shift |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007900000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 15
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300007800000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 20
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007100000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 50
	}  else if( trainConfig == 11 ) {  // No eta shift |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300002200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.7
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300003200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.5
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009365300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360100007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
	}  else if( trainConfig == 12 ) {  // No eta shift |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009380300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec  sec  2.0sigmas Low and  1 High momentum
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20435400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 2.0sigmas Low and 0 High momentum
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20477400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475200233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum
	}  else if( trainConfig == 13 ) {  // No eta shift  |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20425400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx pion rejec  primary  2.0sigmas Low and -1 High momentum 
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400133102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + SPD first layer
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233302223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + PsiPair cut 0.52
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263100500900000"; //standard cut Pi0 pPb 00-100 Standard cut + Alpha cut < 0.7	
	}  else if( trainConfig == 14 ) { // No eta shift  |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20375400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -5,5
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20575400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -3,5
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009160300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -5,5
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009260300007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -3,5
	} else if ( trainConfig == 15 ) { //No eta shift   |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "04200009360300007200000000"; ElecCutarray[0] = "20475400235102223710"; MesonCutarray[0] = "0103203500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Y < 0.70  and prim and sec e |eta| < 0.75 //NOTE revisar
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233102233710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single prim Pt cut > 0.150 
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233102231710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + Single prim Pt cut > 0.100 
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102222710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + DCAxy < 1 cm
	} else if ( trainConfig == 16 ) {  // No eta shift |Y| < 0.8 
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200049360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.075
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200019360300007200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.100
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200008360300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls sec  > 0.35
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400273102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls prim > 0.60
	} else if ( trainConfig == 17 ) {  //No eta shift |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223810"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.015 < InvMass(e+,e-) < 0.050
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233102223910"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.025 < InvMass(e+,e-) < 0.035
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300001200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100 Standard cut + qT < 0.1
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	} else if ( trainConfig == 18 ) {  //No eta shift |Y| < 0.8 + AddedSignals
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0273103500900000"; //standard cut Pi0 pPb 00-100  New Stardad cut +100 events background 
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0163103500900000"; //standard cut Pi0 pPb 00-100  New Standard cut + Background method V0 multiplicity 
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233102223710"; MesonCutarray[2] = "0263103500000000"; //standard cut Pi0 PbPb 00-100 + No extra smearing
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400253102221710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100 + Old Standard
	} else if ( trainConfig == 19 ) {
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
		eventCutArray[1]="80000113"; photonCutArray[1] = "03200009360300007200000000"; ElecCutarray[1] = "20475400239102223710"; MesonCutarray[1] = "0103303500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
		eventCutArray[2]="80000113"; photonCutArray[2] = "04200009360300007200000000"; ElecCutarray[2] = "20475400235102223710"; MesonCutarray[2] = "0103203500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
		eventCutArray[3]="80000113"; photonCutArray[3] = "01200009360300007200000000"; ElecCutarray[3] = "20475400236102223710"; MesonCutarray[3] = "0103403500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
	} else if ( trainConfig == 20 ) {
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
		eventCutArray[1]="80000123"; photonCutArray[1] = "03200009360300007200000000"; ElecCutarray[1] = "20475400239102223710"; MesonCutarray[1] = "0103303500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
		eventCutArray[2]="80000123"; photonCutArray[2] = "04200009360300007200000000"; ElecCutarray[2] = "20475400235102223710"; MesonCutarray[2] = "0103203500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
		eventCutArray[3]="80000123"; photonCutArray[3] = "01200009360300007200000000"; ElecCutarray[3] = "20475400236102223710"; MesonCutarray[3] = "0103403500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
	} else if ( trainConfig == 21 ) {
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400433102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  3Cls ITS
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400533102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  4Cls ITS
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400633102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  5Cls ITS
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400733102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  4Cls ITS no Any
	} else if ( trainConfig == 22 ) {
		eventCutArray[0]="80000123"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400433102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 3 ITScls
		eventCutArray[1]="80000123"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400533102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls
		eventCutArray[2]="80000123"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400633102223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 5 ITScls
		eventCutArray[3]="80000123"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400733102223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls no Any
	} else if ( trainConfig == 23 ) {
		eventCutArray[0]="80000113"; photonCutArray[0] = "00200049360300007200000000"; ElecCutarray[0] = "20475400233102223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt > 0.075
		eventCutArray[1]="80000113"; photonCutArray[1] = "00200019360300007200000000"; ElecCutarray[1] = "20475400233102223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt > 0.100
		eventCutArray[2]="80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233102233710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt{e} > 0.150 
		eventCutArray[3]="80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233102253710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt{e} > 0.175 
	} else if ( trainConfig == 24 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400833202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth  + new psiPair cut   0.60, 0.0 0.12
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400133202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kFirst + new psiPair cut + 0.60  0.0 0.12
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut   0.60, 0.0 0.12
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00500009360300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  10cm  0.60, 0.0 0.12
	} else if ( trainConfig == 25 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00800009360300007200000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  12.5cm  0.60, 0.0 0.12
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00600009360300007200000000"; ElecCutarray[1] = "20475400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  20 cm  0.60, 0.0 0.12
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00700009360300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  35 cm  0.60, 0.0 0.12
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00900009360300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  7.5 cm  0.60, 0.0 0.12
	} else if ( trainConfig == 26 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00800009360300007200000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  12.5cm  0.60, 0.0 0.12
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00600009360300007200000000"; ElecCutarray[1] = "20475400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  20 cm  0.60, 0.0 0.12
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00700009360300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  35 cm  0.60, 0.0 0.12
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00900009360300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + gammaR >  7.5 cm  0.60, 0.0 0.12
	} else if ( trainConfig == 27 ) {
		eventCutArray[0]= "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233202213710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Ptprim > 100 MeV 
		eventCutArray[1]= "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233202233710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Ptprim > 150 MeV 
		eventCutArray[2]= "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233202203710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Ptprim > 75 MeV  
		eventCutArray[3]= "80000113"; photonCutArray[3] = "00200049360300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Ptsec  > 75 MeV
	} else if ( trainConfig == 28 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200019360300007200000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Pt sec  > 100 MeV
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233202653710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Pt prim > 175 GeV
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300001200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Qt > 0.1
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300002200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12 + Qt > 0.07
	} else if ( trainConfig == 29 ) {
		eventCutArray[0]= "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20375400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx primary electron -5,5                            
		eventCutArray[1]= "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20575400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx primary electron -3,5
		eventCutArray[2]= "80000113"; photonCutArray[2] = "00200009160300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx secondary electron -5,5
		eventCutArray[3]= "80000113"; photonCutArray[3] = "00200009260300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx secondary electron -3,5
	} else if ( trainConfig == 30 ) {
		eventCutArray[0]= "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20435400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx pion rejec primary 2.0 sigmas Low and 0 High momentum
		eventCutArray[1]= "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20425400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx pion rejec primary 2.0 sigmas Low and -1 High momentum 
		eventCutArray[2]= "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20477400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum 
		eventCutArray[3]= "80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475200233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum*/
	} else if ( trainConfig == 31 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009365300007200000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360100007200000000"; ElecCutarray[1] = "20475400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009380300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + dEdx pion rejec sec  2.0 sigmas Low and  1 High momentum
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007900000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + Chi2  < 15
	} else if ( trainConfig == 32 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007800000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + Chi2  < 20
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007100000000"; ElecCutarray[1] = "20475400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + Chi2  < 50
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263100500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair cut + Alpha < 0.7
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233002223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + no psi pair +  weights
	} else if ( trainConfig == 33 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233202222710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + DCAxy < 1 cm
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400233202223810"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + 0.015 < InvMass(e+,e-) < 0.050
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400233202223910"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + 0.025 < InvMass(e+,e-) < 0.035
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0273103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  
	} else if ( trainConfig == 34 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0163103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + Background method V0 multiplicity 
		eventCutArray[1] = "80000113"; photonCutArray[1] = "03200009360300007200000000"; ElecCutarray[1] = "20475400239202223710"; MesonCutarray[1] = "0263303500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
		eventCutArray[2] = "80000113"; photonCutArray[2] = "04200009360300007200000000"; ElecCutarray[2] = "20475400235202223710"; MesonCutarray[2] = "0263203500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
		eventCutArray[3] = "80000113"; photonCutArray[3] = "01200009360300007200000000"; ElecCutarray[3] = "20475400236202223710"; MesonCutarray[3] = "0263403500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60
	} else if ( trainConfig == 35 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + new psiPair Cut    0.60, 0.0 0.12
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300003200000000"; ElecCutarray[1] = "20475400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  +  Qt < 0.05
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200008360300007200000000"; ElecCutarray[2] = "20475400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + Findable Cls > 0.35 Secondary
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400273202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + Findable Cls > 0.60 primary
	} else if ( trainConfig == 36 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400433202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + ITScls >= 3
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400533202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + ITScls >= 4
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400633202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + ITScls >= 5
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233202223710"; MesonCutarray[3] = "0263103500000000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny  + No extra smearing
	} else if ( trainConfig == 37 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200000000"; ElecCutarray[0] = "20475400233202223700"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny    + new psiPair Cut    0.60, 0.0 0.12
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200000000"; ElecCutarray[1] = "20475400833202223700"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth   + new psiPair cut    0.60, 0.0 0.12
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200000000"; ElecCutarray[2] = "20475400133202223700"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kFirst  + new psiPair Cut    0.60, 0.0 0.12
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200000000"; ElecCutarray[3] = "20475400233002223700"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny    + no psi pair + no weights
	} else if ( trainConfig == 38 ){
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20405400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	New standard
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20405400233202223710"; MesonCutarray[1] = "0163103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   Background gammas
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263605500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Y < 0.75
		eventCutArray[3] = "80000113"; photonCutArray[3] = "04200009360300007200004000"; ElecCutarray[3] = "20405400235202223710"; MesonCutarray[3] = "0263107500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	n < 0.75
	} else if ( trainConfig == 39 ){
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007800004000"; ElecCutarray[0] = "20405400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Gamma Chi2 < 20
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007100004000"; ElecCutarray[1] = "20405400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Gamma Chi2 < 50
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300001200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	qT < 1.0 (1D)
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300003200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	qT < 0.5 (1D)
	} else if ( trainConfig == 40 ){
                eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20415400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   dEdx prim pion reje  3 low, -10 high
                eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009310300007200004000"; ElecCutarray[1] = "20405400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   dEdx sec pion  reje  0 low, -10 high
                eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009320300007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   dEdx sec pion  reje  1 low, -10 high
                eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300001200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	qT < 1.0 (1D)	
       } else if ( trainConfig == 41 ){
                eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20415400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   
                eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009310300007200004000"; ElecCutarray[1] = "20405400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   
                eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009320300007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   
                eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300001200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection
       } else  if ( trainConfig == 42 ){
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20705400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx prim electron   -3.0, 5.0
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20505400233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx prim electron   -2.5, 4.0
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009260300007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx sec  electron   -3.0, 5.0
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009660300007200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx sec  electron   -2.5, 4.0
       } else  if ( trainConfig == 43 ){
	 	eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20455400233202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx prim pion reje   0.5 < p < 3.5  2 low, -10 high                      
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20487200233202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx prim pion reje   0.3 < p < 5.0  1.5 low,  -1 high
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009315500007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx sec  pion reje   0.3 < p < 7.0  0.0 low, -10 high
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009340300007200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	dEdx sec  pion reje   0.5 < p < 3.5  3 low, 1 high
       } else  if ( trainConfig == 44 ){
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20405400223202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Findable/NCls  crossroad 0.9                                                                      
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20405400213202223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Findable/NCls  crossroad 0.7
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200008360300007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Findable/NCls  sec 0.35
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200006360300007200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Findable/NCls  sec 0.70
       } else if  ( trainConfig == 45 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20405400233202213710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Pt prim > 100 MeV
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20405400233202233710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Pt prim > 150 MeV
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200049360300007200004000"; ElecCutarray[2] = "20405400233202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Pt sec  > 75 MeV
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200019360300007200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Pt sec  > 100 MeV
       } else if ( trainConfig == 46 ) {	 
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20405400233302223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection   PsiPairVsDelthaPhi    	0.52, 0 - 0.12
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20405400233602223710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	PsiPairVsDelthaPhi      0.65  0 - 0.14
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200004000"; ElecCutarray[2] = "20405400233202223810"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	InvMass cut  0.015 < 1GeV < 0.050
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200004000"; ElecCutarray[3] = "20405400233202223710"; MesonCutarray[3] = "0273103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Background mult 100		
       }  else if ( trainConfig == 47 ) {
		eventCutArray[0] = "80000113"; photonCutArray[0] = "00200009360300007200004000"; ElecCutarray[0] = "20405400243202223710"; MesonCutarray[0] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Prim Crossroad > 90
		eventCutArray[1] = "80000113"; photonCutArray[1] = "00200009360300007200004000"; ElecCutarray[1] = "20405400233202222710"; MesonCutarray[1] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Prim DCAxy > 1 cm
		eventCutArray[2] = "80000113"; photonCutArray[2] = "00200009360300007200004000"; ElecCutarray[2] = "20405400133202223710"; MesonCutarray[2] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Prim ITS cluster kFirst     
		eventCutArray[3] = "80000113"; photonCutArray[3] = "00200009360300007200004000"; ElecCutarray[3] = "20405400833202223710"; MesonCutarray[3] = "0263103500900000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny   + new psiPair Cut + double counting rejection	Prim ITS cluster kBoth
       }
       

		
	TList *EventCutList = new TList();
	TList *ConvCutList  = new TList();
	TList *MesonCutList = new TList();
	TList *ElecCutList  = new TList();

	TList *HeaderList = new TList();
	TObjString *Header1 = new TObjString("pi0_1");
	HeaderList->Add(Header1);
	TObjString *Header3 = new TObjString("eta_2");
	HeaderList->Add(Header3);
	
	EventCutList->SetOwner(kTRUE);
	AliConvEventCuts **analysisEventCuts 		= new AliConvEventCuts*[numberOfCuts];
	ConvCutList->SetOwner(kTRUE);
	AliConversionPhotonCuts **analysisCuts       	= new AliConversionPhotonCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
	AliConversionMesonCuts **analysisMesonCuts   	= new AliConversionMesonCuts*[numberOfCuts];
	ElecCutList->SetOwner(kTRUE);
	AliDalitzElectronCuts **analysisElecCuts     	= new AliDalitzElectronCuts*[numberOfCuts];

	for(Int_t i = 0; i<numberOfCuts; i++){

		analysisEventCuts[i] = new AliConvEventCuts();
		if(  ( trainConfig >= 1 && trainConfig <= 9 ) || trainConfig == 19  || trainConfig == 21 || trainConfig == 23 || ( trainConfig >= 24 && trainConfig <=36 )  || trainConfig == 38 
		     || trainConfig == 39 || trainConfig == 40 || trainConfig == 41  || trainConfig == 42 || trainConfig == 43 || trainConfig == 44 || trainConfig == 45    || trainConfig == 46
		     || trainConfig == 47 ){
			if (doWeighting){
				if (generatorName.CompareTo("DPMJET")==0){
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				} else if (generatorName.CompareTo("HIJING")==0){   
					analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
				}
			}
		} else if (  ( trainConfig >= 10 && trainConfig <= 18 ) || trainConfig == 20 || trainConfig == 22 ){
			if (doWeighting){
				analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
			}	
		}
	
	
		analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
		if (doEtaShiftIndCuts) {
			analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
			analysisEventCuts[i]->SetEtaShift(stringShift);
		}
		EventCutList->Add(analysisEventCuts[i]);
		analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
		analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
		
		analysisCuts[i] = new AliConversionPhotonCuts();
		
		if( ! analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data()) ) {
			cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
			return 0;
		}
		
		
					
		
		
		analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
		ConvCutList->Add(analysisCuts[i]);
		analysisCuts[i]->SetFillCutHistograms("",kFALSE);
			
		analysisMesonCuts[i] = new AliConversionMesonCuts();
		if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
				cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
				return 0;
		}
		MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("");
		
		analysisElecCuts[i] = new AliDalitzElectronCuts();
		analysisElecCuts[i]->SetUseCrossedRows(kTRUE);
		if( !analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data())) {

				cout<< "ERROR:  analysisElecCuts [ " <<i<<" ] "<<endl;
				return 0;
		}
		ElecCutList->Add(analysisElecCuts[i]);
		analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,ElecCutarray[i].Data()); 
		
	}

	task->SetEventCutList(numberOfCuts,EventCutList);
	task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetMesonCutList(MesonCutList);
	task->SetElectronCutList(ElecCutList);

	task->SetMoveParticleAccordingToVertex(kTRUE);
	task->SetProductionVertextoVGamma(kTRUE);


	task->SetDoMesonQA(enableQAMesonTask);
	if(enableDoMesonChic) task->SetDoChicAnalysis(kTRUE);

	//connect containers
	AliAnalysisDataContainer *coutput =
	mgr->CreateContainer(Form("GammaConvDalitzV1_%i",trainConfig), TList::Class(),
							AliAnalysisManager::kOutputContainer,Form("GammaConvV1Dalitz_%i.root",trainConfig));

	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return;

}
