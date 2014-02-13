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
   TString ConvCutnumber = "8000000060084001001500000000";   //Online  V0 finder
   TString ElecCuts      = "9000540000000200000";            //Electron Cuts
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
      if( ConvCutnumber !=""){
         fCuts= new AliConversionCuts(ConvCutnumber.Data(),ConvCutnumber.Data());
         fCuts->SetPreSelectionCutFlag(kTRUE);
         if(fCuts->InitializeCutsFromCutString(ConvCutnumber.Data())){
            fCuts->DoEtaShift(doEtaShift);
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
   //========= Add Electron Selector ================


   if( !(AliDalitzElectronSelector*)mgr->GetTask("ElectronSelector") ){

   AliDalitzElectronSelector *fElectronSelector = new AliDalitzElectronSelector("ElectronSelector");

   // Set AnalysisCut Number

   AliDalitzElectronCuts *fElecCuts=0;

   //ElecCuts = "900054000000020000";

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

   task->SetIsHeavyIon(2);
   task->SetIsMC(isMC);



   // Cut Numbers to use in Analysis
   Int_t numberOfCuts = 4;

   TString *ConvCutarray    = new TString[numberOfCuts];

   TString *ElecCutarray    = new TString[numberOfCuts];

   TString *MesonCutarray   = new TString[numberOfCuts];

   Bool_t doEtaShiftIndCuts = kFALSE;
   TString stringShift = "";

   // Shifting in pPb direction

   doEtaShiftIndCuts = kFALSE;
   stringShift = "pPb";


   
   
if( trainConfig == 1 ) {  // No eta shift |Y| < 0.8
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	ConvCutarray[1] = "8000011002093603007900000000"; ElecCutarray[1] = "9047540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 15
	ConvCutarray[2] = "8000011002093603007800000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 20
	ConvCutarray[3] = "8000011002093603007100000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 50
	
}  else if( trainConfig == 2 ) {  // No eta shift |Y| < 0.8

	ConvCutarray[0] = "8000011002093603002200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.7
	ConvCutarray[1] = "8000011002093603003200000000"; ElecCutarray[1] = "9047540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.5
	ConvCutarray[2] = "8000011002093653007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
	ConvCutarray[3] = "8000011002093601007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 3 ) {  // No eta shift |Y| < 0.8

	ConvCutarray[0] = "8000011002093803007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec  sec  2.0sigmas Low and  1 High momentum
	ConvCutarray[1] = "8000011002093603007200000000"; ElecCutarray[1] = "9043540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 2.0sigmas Low and 0 High momentum
	ConvCutarray[2] = "8000011002093603007200000000"; ElecCutarray[2] = "9047740023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047520023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 4 ) {  // No eta shift  |Y| < 0.8
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "9042540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx pion rejec  primary  2.0sigmas Low and -1 High momentum 
	ConvCutarray[1] = "8000011002093603007200000000"; ElecCutarray[1] = "9047540013310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + SPD first layer
	ConvCutarray[2] = "8000011002093603007200000000"; ElecCutarray[2] = "9047540023330262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + PsiPair cut 0.52
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031005009000"; //standard cut Pi0 pPb 00-100 Standard cut + Alpha cut < 0.7	
	
}  else if( trainConfig == 5 ) { // No eta shift  |Y| < 0.8	
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "9037540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -5,5
	ConvCutarray[1] = "8000011002093603007200000000"; ElecCutarray[1] = "9057540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -3,5
	ConvCutarray[2] = "8000011002091603007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -5,5
	ConvCutarray[3] = "8000011002092603007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -3,5
	
} else if ( trainConfig == 6 ) { //No eta shift   |Y| < 0.8
	
	ConvCutarray[0] = "8000011042093603007200000000"; ElecCutarray[0] = "9047540023510262371"; MesonCutarray[0] = "01032035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Y < 0.70  and prim and sec e |eta| < 0.75 //NOTE revisar
	ConvCutarray[1] = "8000011002093603007200000000"; ElecCutarray[1] = "9047540023310263371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single prim Pt cut > 0.150
	ConvCutarray[2] = "8000011002093603007200000000"; ElecCutarray[2] = "9047540023310263171"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + Single prim Pt cut > 0.100
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047540023310262271"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + DCAxy < 1 cm
	
	
} else if ( trainConfig == 7 ) {  // No eta shift |Y| < 0.8

	ConvCutarray[0] = "8000011002493603007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.075
	ConvCutarray[1] = "8000011002193603007200000000"; ElecCutarray[1] = "9047540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.100
	ConvCutarray[2] = "8000011002083603007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls sec  > 0.35
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047540027310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls prim > 0.60
	
	
} else if ( trainConfig == 8 ) {  //No eta shift |Y| < 0.8
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "9047540023310262381"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.015 < InvMass(e+,e-) < 0.050
	ConvCutarray[1] = "8000011002093603007200000000"; ElecCutarray[1] = "9047540023310262391"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.025 < InvMass(e+,e-) < 0.035
	ConvCutarray[2] = "8000011002093603001200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + qT < 0.1
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 9 ) {  //No eta shift |Y| < 0.8

  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "9047540023310272371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Stardad cut +100 events background
	ConvCutarray[1] = "8000011002093603007200000000"; ElecCutarray[1] = "9047540023310162371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Background method V0 multiplicity
	ConvCutarray[2] = "8000011002093603007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035000000"; //standard cut Pi0 PbPb 00-100 + No extra smearing
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047540025310262171"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard

} else if( trainConfig == 10 ) {  // No eta shift |Y| < 0.8 + AddedSignals
	
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	ConvCutarray[1] = "8000012002093603007900000000"; ElecCutarray[1] = "9047540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 15
	ConvCutarray[2] = "8000012002093603007800000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 20
	ConvCutarray[3] = "8000012002093603007100000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Chi2 < 50
	
}  else if( trainConfig == 11 ) {  // No eta shift |Y| < 0.8 + AddedSignals

	ConvCutarray[0] = "8000012002093603002200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.7
	ConvCutarray[1] = "8000012002093603003200000000"; ElecCutarray[1] = "9047540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Qt < 0.5
	ConvCutarray[2] = "8000012002093653007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.3 GeV Low and 3.5 High momentum
	ConvCutarray[3] = "8000012002093601007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec sec 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 12 ) {  // No eta shift |Y| < 0.8 + AddedSignals

	ConvCutarray[0] = "8000012002093803007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec  sec  2.0sigmas Low and  1 High momentum
	ConvCutarray[1] = "8000012002093603007200000000"; ElecCutarray[1] = "9043540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 2.0sigmas Low and 0 High momentum
	ConvCutarray[2] = "8000012002093603007200000000"; ElecCutarray[2] = "9047740023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.3 GeV Low and 3.5 High momentum
	ConvCutarray[3] = "8000012002093603007200000000"; ElecCutarray[3] = "9047520023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + dEdx pion rejec primary 0.5 GeV Low and 5.0 High momentum
	
}  else if( trainConfig == 13 ) {  // No eta shift  |Y| < 0.8 + AddedSignals
  
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "9042540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx pion rejec  primary  2.0sigmas Low and -1 High momentum 
	ConvCutarray[1] = "8000012002093603007200000000"; ElecCutarray[1] = "9047540013310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + SPD first layer
	ConvCutarray[2] = "8000012002093603007200000000"; ElecCutarray[2] = "9047540023330262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + PsiPair cut 0.52
	ConvCutarray[3] = "8000012002093603007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031005009000"; //standard cut Pi0 pPb 00-100 Standard cut + Alpha cut < 0.7	
	
}  else if( trainConfig == 14 ) { // No eta shift  |Y| < 0.8 + AddedSignals
  
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "9037540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -5,5
	ConvCutarray[1] = "8000012002093603007200000000"; ElecCutarray[1] = "9057540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx primary   electron -3,5
	ConvCutarray[2] = "8000012002091603007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -5,5
	ConvCutarray[3] = "8000012002092603007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 New Standard cut + dEdx secondary electron -3,5
	
} else if ( trainConfig == 15 ) { //No eta shift   |Y| < 0.8 + AddedSignals
	
	ConvCutarray[0] = "8000012042093603007200000000"; ElecCutarray[0] = "9047540023510262371"; MesonCutarray[0] = "01032035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Y < 0.70  and prim and sec e |eta| < 0.75 //NOTE revisar
	ConvCutarray[1] = "8000012002093603007200000000"; ElecCutarray[1] = "9047540023310263371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single prim Pt cut > 0.150
	ConvCutarray[2] = "8000012002093603007200000000"; ElecCutarray[2] = "9047540023310263171"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + Single prim Pt cut > 0.100
	ConvCutarray[3] = "8000012002093603007200000000"; ElecCutarray[3] = "9047540023310262271"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + DCAxy < 1 cm
	
	
} else if ( trainConfig == 16 ) {  // No eta shift |Y| < 0.8 + AddedSignals

	ConvCutarray[0] = "8000011002493603007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.075
	ConvCutarray[1] = "8000011002193603007200000000"; ElecCutarray[1] = "9047540023310262371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Single sec  Pt cut > 0.100
	ConvCutarray[2] = "8000011002083603007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls sec  > 0.35
	ConvCutarray[3] = "8000011002093603007200000000"; ElecCutarray[3] = "9047540027310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Findable Cls prim > 0.60
	
	
} else if ( trainConfig == 17 ) {  //No eta shift |Y| < 0.8 + AddedSignals
	
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "9047540023310262381"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.015 < InvMass(e+,e-) < 0.050
	ConvCutarray[1] = "8000012002093603007200000000"; ElecCutarray[1] = "9047540023310262391"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + 0.025 < InvMass(e+,e-) < 0.035
	ConvCutarray[2] = "8000012002093603001200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035009000"; //standard cut Pi0 pPb 00-100 Standard cut + qT < 0.1
	ConvCutarray[3] = "8000012002093603007200000000"; ElecCutarray[3] = "9047540023310262371"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 18 ) {  //No eta shift |Y| < 0.8 + AddedSignals
  
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "9047540023310272371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Stardad cut +100 events background
	ConvCutarray[1] = "8000012002093603007200000000"; ElecCutarray[1] = "9047540023310162371"; MesonCutarray[1] = "01031035009000"; //standard cut Pi0 pPb 00-100  New Standard cut + Background method V0 multiplicity
	ConvCutarray[2] = "8000012002093603007200000000"; ElecCutarray[2] = "9047540023310262371"; MesonCutarray[2] = "01031035000000"; //standard cut Pi0 PbPb 00-100 + No extra smearing
	ConvCutarray[3] = "8000012002093603007200000000"; ElecCutarray[3] = "9047540025310262171"; MesonCutarray[3] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard

} else if ( trainConfig == 19 ) {
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	ConvCutarray[1] = "8000011032093603007200000000"; ElecCutarray[1] = "9047540023910262371"; MesonCutarray[1] = "01033035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
	ConvCutarray[2] = "8000011042093603007200000000"; ElecCutarray[2] = "9047540023510262371"; MesonCutarray[2] = "01032035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
	ConvCutarray[3] = "8000011012093603007200000000"; ElecCutarray[3] = "9047540023610262371"; MesonCutarray[3] = "01034035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
	
} else if ( trainConfig == 20 ) {
	
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "9047540023310262371"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	ConvCutarray[1] = "8000012032093603007200000000"; ElecCutarray[1] = "9047540023910262371"; MesonCutarray[1] = "01033035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
	ConvCutarray[2] = "8000012042093603007200000000"; ElecCutarray[2] = "9047540023510262371"; MesonCutarray[2] = "01032035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75
	ConvCutarray[3] = "8000012012093603007200000000"; ElecCutarray[3] = "9047540023610262371"; MesonCutarray[3] = "01034035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
}




   TList *ConvCutList  = new TList();
   TList *MesonCutList = new TList();
   TList *ElecCutList  = new TList();

   TList *HeaderList = new TList();
   TObjString *Header1 = new TObjString("pi0_1");
   HeaderList->Add(Header1);
   TObjString *Header3 = new TObjString("eta_2");
   HeaderList->Add(Header3);
   
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts             = new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
   ElecCutList->SetOwner(kTRUE);
   AliDalitzElectronCuts **analysisElecCuts     = new AliDalitzElectronCuts*[numberOfCuts];



   for(Int_t i = 0; i<numberOfCuts; i++){


      analysisCuts[i] = new AliConversionCuts();
      if( ! analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data()) ) {
            cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
            return 0;
      } else {

	  if (  ( trainConfig >= 1 && trainConfig <= 9 ) || trainConfig == 19   ){
	    
	    if (doWeighting){
	      if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      } else if (generatorName.CompareTo("HIJING")==0){   
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      }
	    }
	  }
	  else if (  ( trainConfig >= 10 && trainConfig <= 18 ) || trainConfig == 20 ){
	    
    	    if (doWeighting){
	      analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	    }
	    
	  }
   

      
	  if (doEtaShiftIndCuts) {
	  
	      analysisCuts[i]->DoEtaShift(doEtaShiftIndCuts);
	      analysisCuts[i]->SetEtaShift(stringShift);
	      
	  }
	  ConvCutList->Add(analysisCuts[i]);
	  analysisCuts[i]->SetFillCutHistograms("",kFALSE);
	  analysisCuts[i]->SetAcceptedHeader(HeaderList);
	}



      analysisMesonCuts[i] = new AliConversionMesonCuts();
    
      if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
            cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
            return 0;
      }
      else {
            MesonCutList->Add(analysisMesonCuts[i]);
            analysisMesonCuts[i]->SetFillCutHistograms("");
      }


       TString cutName( Form("%s_%s_%s",ConvCutarray[i].Data(),ElecCutarray[i].Data(),MesonCutarray[i].Data() ) );


       analysisElecCuts[i] = new AliDalitzElectronCuts();
       if( !analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data())) {

            cout<< "ERROR:  analysisElecCuts [ " <<i<<" ] "<<endl;
            return 0;
       }
       else { 
        ElecCutList->Add(analysisElecCuts[i]);
        analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
       }
     

   }


   task->SetConversionCutList(numberOfCuts,ConvCutList);
   task->SetMesonCutList(MesonCutList);
   task->SetElectronCutList(ElecCutList);

   task->SetMoveParticleAccordingToVertex(kTRUE);


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
