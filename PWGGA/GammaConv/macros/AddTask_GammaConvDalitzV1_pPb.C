void AddTask_GammaConvDalitzV1_pPb(    Int_t trainConfig = 1,
                                       Bool_t isMC       = kFALSE, //run MC 
                                       Bool_t enableQAMesonTask = kTRUE, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                       Bool_t enableDoMesonChic = kFALSE, // enable additional Chic analysis
				       TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                       Bool_t doWeighting = kFALSE,  //enable Weighting
                                       TString generatorName = "DPMJET",				
                                       TString cutnumberAODBranch = "0000000060084001001500000"
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


   cout<<"Entro 0"<<endl;
   
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
   TString cutnumberEvent = "8000000";
   
   TString cutnumberPhoton="060084001001500000000";   //Online  V0 finder
   
   TString ElecCuts      = "90005400000002000000";            //Electron Cuts
   
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



    cout<<"Entro"<<endl;
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
	
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007900000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 15
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007800000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 20
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007100000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 50
	
}  else if( trainConfig == 2 ) {  // No eta shift |Y| < 0.8
							 
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603002200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.7
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603003200000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.5
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093653007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093601007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 3 ) {  // No eta shift |Y| < 0.8

	eventCutArray[0]="8000011"; photonCutArray[0] = "002093803007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec  sec  2.0sigmas Low and  1 High momentum
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90435400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 2.0sigmas Low and 0 High momentum
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90477400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475200233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 4 ) {  // No eta shift  |Y| < 0.8
  
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90425400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx pion rejec  primary  2.0sigmas Low and -1 High momentum 
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400133102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + SPD first layer
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233302623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + PsiPair cut 0.52
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031005009000"; //standard cut Pi0 pPb 00-100 Standard cut + Alpha cut < 0.7	
	
}  else if( trainConfig == 5 ) { // No eta shift  |Y| < 0.8	
  
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90375400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -5,5
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90575400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -3,5
	eventCutArray[2]="8000011"; photonCutArray[2] = "002091603007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -5,5
	eventCutArray[3]="8000011"; photonCutArray[3] = "002092603007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -3,5
	
} else if ( trainConfig == 6 ) { //No eta shift   |Y| < 0.8
	
	eventCutArray[0]="8000011"; photonCutArray[0] = "042093603007200000000"; ElecCutarray[0] = "90475400235102623710"; MesonCutarray[0] = "01032035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Y < 0.70  and prim and sec e |eta| < 0.75 //NOTE revisar
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400233102633710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single prim Pt cut > 0.150
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233102631710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + Single prim Pt cut > 0.100
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102622710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + DCAxy < 1 cm
	
	
} else if ( trainConfig == 7 ) {  // No eta shift |Y| < 0.8

	eventCutArray[0]="8000011"; photonCutArray[0] = "002493603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.075
	eventCutArray[1]="8000011"; photonCutArray[1] = "002193603007200000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.100
	eventCutArray[2]="8000011"; photonCutArray[2] = "002083603007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls sec  > 0.35
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400273102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls prim > 0.60
	
	
} else if ( trainConfig == 8 ) {  //No eta shift |Y| < 0.8
	
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102623810"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.015 < InvMass(e+,e-) < 0.050
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400233102623910"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.025 < InvMass(e+,e-) < 0.035
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603001200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + qT < 0.1
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 9 ) {  //No eta shift |Y| < 0.8

  
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102723710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Stardad cut +100 events background
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400233101623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Background method V0 multiplicity
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035000000"; //standard cut Pi0 PbPb 00-100 + No extra smearing
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400253102621710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard

} else if( trainConfig == 10 ) {  // No eta shift |Y| < 0.8 + AddedSignals
	
	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007900000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 15
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603007800000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 20
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007100000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 50
	
}  else if( trainConfig == 11 ) {  // No eta shift |Y| < 0.8 + AddedSignals

	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603002200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.7
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603003200000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.5
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093653007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093601007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 12 ) {  // No eta shift |Y| < 0.8 + AddedSignals

	eventCutArray[0]="8000012"; photonCutArray[0] = "002093803007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec  sec  2.0sigmas Low and  1 High momentum
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90435400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 2.0sigmas Low and 0 High momentum
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90477400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475200233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 13 ) {  // No eta shift  |Y| < 0.8 + AddedSignals
  
	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90425400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx pion rejec  primary  2.0sigmas Low and -1 High momentum 
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400133102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + SPD first layer
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233302623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + PsiPair cut 0.52
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031005009000"; //standard cut Pi0 pPb 00-100 Standard cut + Alpha cut < 0.7	
	
}  else if( trainConfig == 14 ) { // No eta shift  |Y| < 0.8 + AddedSignals
  
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90375400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -5,5
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90575400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -3,5
	eventCutArray[2]="8000011"; photonCutArray[2] = "002091603007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -5,5
	eventCutArray[3]="8000011"; photonCutArray[3] = "002092603007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -3,5
	
} else if ( trainConfig == 15 ) { //No eta shift   |Y| < 0.8 + AddedSignals
	
	eventCutArray[0]="8000012"; photonCutArray[0] = "042093603007200000000"; ElecCutarray[0] = "90475400235102623710"; MesonCutarray[0] = "01032035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Y < 0.70  and prim and sec e |eta| < 0.75 //NOTE revisar
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400233102633710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single prim Pt cut > 0.150
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233102631710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + Single prim Pt cut > 0.100
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102622710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + DCAxy < 1 cm
	
	
} else if ( trainConfig == 16 ) {  // No eta shift |Y| < 0.8 

	eventCutArray[0]="8000011"; photonCutArray[0] = "002493603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.075
	eventCutArray[1]="8000011"; photonCutArray[1] = "002193603007200000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.100
	eventCutArray[2]="8000011"; photonCutArray[2] = "002083603007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls sec  > 0.35
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400273102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls prim > 0.60
	
	
} else if ( trainConfig == 17 ) {  //No eta shift |Y| < 0.8 + AddedSignals
	
	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102623810"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.015 < InvMass(e+,e-) < 0.050
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400233102623910"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.025 < InvMass(e+,e-) < 0.035
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603001200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + qT < 0.1
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 18 ) {  //No eta shift |Y| < 0.8 + AddedSignals
  
	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102723710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Stardad cut +100 events background
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400233101623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Background method V0 multiplicity
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233102623710"; MesonCutarray[2] = "01031035000000"; //standard cut Pi0 PbPb 00-100 + No extra smearing
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400253102621710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard

} else if ( trainConfig == 19 ) {
	
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	eventCutArray[1]="8000011"; photonCutArray[1] = "032093603007200000000"; ElecCutarray[1] = "90475400239102623710"; MesonCutarray[1] = "01033035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
	eventCutArray[2]="8000011"; photonCutArray[2] = "042093603007200000000"; ElecCutarray[2] = "90475400235102623710"; MesonCutarray[2] = "01032035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
	eventCutArray[3]="8000011"; photonCutArray[3] = "012093603007200000000"; ElecCutarray[3] = "90475400236102623710"; MesonCutarray[3] = "01034035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
	
} else if ( trainConfig == 20 ) {
	
	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	eventCutArray[1]="8000012"; photonCutArray[1] = "032093603007200000000"; ElecCutarray[1] = "90475400239102623710"; MesonCutarray[1] = "01033035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
	eventCutArray[2]="8000012"; photonCutArray[2] = "042093603007200000000"; ElecCutarray[2] = "90475400235102623710"; MesonCutarray[2] = "01032035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
	eventCutArray[3]="8000012"; photonCutArray[3] = "012093603007200000000"; ElecCutarray[3] = "90475400236102623710"; MesonCutarray[3] = "01034035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
	
} else if ( trainConfig == 21 ) {
	
	eventCutArray[0]="8000011"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400433102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  3Cls ITS
	eventCutArray[1]="8000011"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400533102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  4Cls ITS
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400633102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  5Cls ITS
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400733102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  +  4Cls ITS no Any
	
	
} else if ( trainConfig == 22 ) {
	
	eventCutArray[0]="8000012"; photonCutArray[0] = "002093603007200000000"; ElecCutarray[0] = "90475400433102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 3 ITScls
	eventCutArray[1]="8000012"; photonCutArray[1] = "002093603007200000000"; ElecCutarray[1] = "90475400533102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls
	eventCutArray[2]="8000012"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400633102623710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 5 ITScls
	eventCutArray[3]="8000012"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400733102623710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls no Any
	
} else if ( trainConfig == 23 ) {
  
	eventCutArray[0]="8000011"; photonCutArray[0] = "002493603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt > 0.075
	eventCutArray[1]="8000011"; photonCutArray[1] = "002193603007200000000"; ElecCutarray[1] = "90475400233102623710"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt > 0.100
	eventCutArray[2]="8000011"; photonCutArray[2] = "002093603007200000000"; ElecCutarray[2] = "90475400233102633710"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt{e} > 0.150
	eventCutArray[3]="8000011"; photonCutArray[3] = "002093603007200000000"; ElecCutarray[3] = "90475400233102653710"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt{e} > 0.175
	
  
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
   AliConvEventCuts **analysisEventCuts 	= new AliConvEventCuts*[numberOfCuts];
   ConvCutList->SetOwner(kTRUE);
   AliConversionPhotonCuts **analysisCuts       = new AliConversionPhotonCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
   ElecCutList->SetOwner(kTRUE);
   AliDalitzElectronCuts **analysisElecCuts     = new AliDalitzElectronCuts*[numberOfCuts];



   for(Int_t i = 0; i<numberOfCuts; i++){

      analysisEventCuts[i] = new AliConvEventCuts();

	  if (  ( trainConfig >= 1 && trainConfig <= 9 ) || trainConfig == 19  || trainConfig == 21 || trainConfig == 23 ){
	    
	    if (doWeighting){
	      if (generatorName.CompareTo("DPMJET")==0){
               analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      } else if (generatorName.CompareTo("HIJING")==0){   
               analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      }
	    }
	  }
	  else if (  ( trainConfig >= 10 && trainConfig <= 18 ) || trainConfig == 20 || trainConfig == 22 ){
	    
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
      


      TString cutName( Form("%s_%s_%s",photonCutArray[i].Data(),ElecCutarray[i].Data(),MesonCutarray[i].Data() ) );


       analysisElecCuts[i] = new AliDalitzElectronCuts();
       if( !analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data())) {

            cout<< "ERROR:  analysisElecCuts [ " <<i<<" ] "<<endl;
            return 0;
       }
       ElecCutList->Add(analysisElecCuts[i]);
       analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
       

   }

   task->SetEventCutList(numberOfCuts,EventCutList);
   task->SetConversionCutList(numberOfCuts,ConvCutList);
   task->SetMesonCutList(MesonCutList);
   task->SetElectronCutList(ElecCutList);

   task->SetMoveParticleAccordingToVertex(kTRUE);
   task->SetProductionVertextoVGamma(kTRUE);


   if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);
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
