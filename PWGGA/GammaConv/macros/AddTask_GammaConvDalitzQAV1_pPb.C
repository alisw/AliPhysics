void AddTask_GammaConvDalitzQAV1_pPb(  Int_t trainConfig = 1,
                                       Bool_t isMC       = kFALSE, //run MC 
                                       Bool_t enableQAMesonTask = kTRUE, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                       Bool_t enableDoMesonChic = kFALSE, // enable additional Chic analysis
				       Bool_t enableSetProdVtxVGamma = kTRUE,
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
   cout<<"enableSetProdVtxVGamma: "<<enableSetProdVtxVGamma<<endl;
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
   TString ConvCutnumber="";
   
   if(trainConfig == 9 || trainConfig == 10 ){
   ConvCutnumber = "8000000160084001001500000000";   //Offline  V0 finder 
   }
   else {
   ConvCutnumber = "8000000060084001001500000000";   //Online  V0 finder
   }
   
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
   Int_t numberOfCuts = 1;

   TString *ConvCutarray    = new TString[numberOfCuts];

   TString *ElecCutarray    = new TString[numberOfCuts];

   TString *MesonCutarray   = new TString[numberOfCuts];

   Bool_t doEtaShiftIndCuts = kFALSE;
   TString stringShift = "";

   // Shifting in pPb direction

   doEtaShiftIndCuts = kFALSE;
   stringShift = "pPb";


   
   
if( trainConfig == 1 ) {  // No eta shift |Y| < 0.8
										     
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
}  else if( trainConfig == 2 ) {  // No eta shift |Y| < 0.8
  
	ConvCutarray[0] = "8000011032093603007200000000"; ElecCutarray[0] = "90475400239102623710"; MesonCutarray[0] = "01033035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.6 and |Gamma_eta| < 0.65 and |e+_eta| < 0.65 and |e-_eta| < 0.65 
	
}  else if( trainConfig == 3 ) {  // No eta shift |Y| < 0.8
  
	ConvCutarray[0] = "8000011042093603007200000000"; ElecCutarray[0] = "90475400235102623710"; MesonCutarray[0] = "01032035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.7 and |Gamma_eta| < 0.75 and |e+_eta| < 0.75 and |e-_eta| < 0.75

}  else if( trainConfig == 4 ) {  // No eta shift  |Y| < 0.8
  
	ConvCutarray[0] = "8000011012093603007200000000"; ElecCutarray[0] = "90475400236102623710"; MesonCutarray[0] = "01034035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  |Y| < 0.5 and |Gamma_eta| < 0.60 and |e+_eta| < 0.60 and |e-_eta| < 0.60  
    
} else if ( trainConfig == 5 ) {

	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102623310"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011

} else if ( trainConfig == 6 ) {  // No eta shift |Y| < 0.8
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102643710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 7 ) {
  
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 8 ) {  // No eta shift |Y| < 0.8
	
	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "90475400233102643710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 9  ) {
	
	ConvCutarray[0] = "8000011102093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 10 ) {

        ConvCutarray[0] = "8000012102093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011

} else if ( trainConfig == 11  ) {
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102623010"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011
	
} else if ( trainConfig == 12 ) {

        ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "90475400233102623010"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011

}  else if ( trainConfig == 13 ) {
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400533102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls
	
}  else if ( trainConfig == 14 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400733102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls no Any
	
}  else if ( trainConfig == 15 ) {

	ConvCutarray[0] = "8000012002093603007200000000"; ElecCutarray[0] = "90475400533102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + 4 ITScls
	
}  else if ( trainConfig == 16 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833002623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth  +  No psipair
		
}  else if ( trainConfig == 17 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth
	
} else if ( trainConfig  == 18 ) {

        ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400933102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + 4ITS cls

} else if ( trainConfig  == 19 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400133102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kFirts
	
} else if ( trainConfig  == 20 ) {

	ConvCutarray[0] = "8000011002092170008260400000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01621035009000"; // standard cut Annika analysis:
	
} else if ( trainConfig  == 21 ) {

	ConvCutarray[0] = "8000011002092170008260400000"; ElecCutarray[0] = "90475400133102623710"; MesonCutarray[0] = "01621035009000"; // standard cut Annika analysis: + kFirst
	
} else if ( trainConfig  == 22 ){
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400153102621710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard 2010 + kFirtst

} else if ( trainConfig  == 23 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400853102621710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard 2010 + kBoth
	
} else if ( trainConfig  == 24 ){
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400153102621700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard 2010 + kFirtst No weights

} else if ( trainConfig  == 25 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400853102621700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100 + Old Standard 2010 + kBoth No weights
	
} else if ( trainConfig  == 26 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400133102623700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kFirts + No weights
	
} else if ( trainConfig  == 27 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833102623700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth  +  No weights
	
} else if ( trainConfig  == 28 ) {
  
	ConvCutarray[0] = "8000011002493603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt > 0.075
	
} else if ( trainConfig  == 29 ) {
  
	ConvCutarray[0] = "8000011002193603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt > 0.100
	
} else if ( trainConfig  == 30 ) {
      
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102633710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt{e} > 0.150
	
} else if ( trainConfig  == 31 ) {
	
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102653710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + Pt{e} > 0.175
	
} else if ( trainConfig  == 32  ) {
  
	ConvCutarray[0] = "8000011007093603007200000000"; ElecCutarray[0] = "90475400233102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + Photon R > 35 cm
  
} else if ( trainConfig  == 33  ) {
  
	ConvCutarray[0] = "8000011007093603007200000000"; ElecCutarray[0] = "90475400833102623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + Photon R > 35 cm 
	
} else if ( trainConfig  == 34  ) {
  
	ConvCutarray[0] = "8000011007093603007200000000"; ElecCutarray[0] = "90475400833102623700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + Photon R > 35 cm + No weights 
	
} else if ( trainConfig  == 35 ) {						
									     
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833002623700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth  + NoPsiPair + No weights
	
} else if ( trainConfig  == 36 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102623700"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kAny no Weights

} else if ( trainConfig  == 37 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833102623711"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + smearing photon virtual

} else if ( trainConfig  == 38 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400133102623711"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kFirts + smearing photon virtual
	
} else if( trainConfig   == 39 ) {  
					     
        ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102623711"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + kAny  + smearing photon virtual 

} else if ( trainConfig  == 40 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833102623712"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + smearing photon virtual  electrons

} else if ( trainConfig  == 41 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400133102623712"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kFirts + smearing photon virtual electrons
	
} else if( trainConfig   == 42 ) {  
					     
        ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233102623712"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + kAny  + smearing photon virtual electrons 

} else if( trainConfig  == 43 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833502623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 +  kBoth + New psi pair cut  fPsiPairCut = 0.60;    fDeltaPhiCutMin = 0.0; fDeltaPhiCutMax = 0.06;

} else if( trainConfig  == 44 ) {
  
	ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400833502623712"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011 + kBoth + New psi pair cut  fPsiPairCut = 0.60;     fDeltaPhiCutMin = 0.0; fDeltaPhiCutMax = 0.06; +  Electron Smearing
	
} else if( trainConfig   == 45 ) {  
					     
        ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233502623710"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + kAny  + New psi pair cut + New psi pair cut  fPsiPairCut = 0.60    fDeltaPhiCutMin = 0.0 fDeltaPhiCutMax = 0.06

} else if( trainConfig   == 46 ) {  

        ConvCutarray[0] = "8000011002093603007200000000"; ElecCutarray[0] = "90475400233502623712"; MesonCutarray[0] = "01031035009000"; //standard cut Pi0 pPb 00-100  //Tracks 2011  + kAny  + New psi pair cut + New psi pair cut  fPsiPairCut = 0.60;    fDeltaPhiCutMin = 0.0; fDeltaPhiCutMax = 0.06; + photon virtual electrons 
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
      

	  if (  ( trainConfig >= 1 && trainConfig <= 6 ) || trainConfig == 9  ||  trainConfig == 11  || trainConfig == 13 || trainConfig == 14 || trainConfig == 16 || trainConfig == 17 || trainConfig == 18 || trainConfig == 19 || trainConfig == 20 || trainConfig == 21 || trainConfig == 22 || trainConfig == 23 ||
		  trainConfig == 28 || trainConfig == 29 || trainConfig == 30 ||  trainConfig == 31  || trainConfig == 32 || trainConfig == 33 || trainConfig == 37 || trainConfig == 38 || trainConfig == 39 || trainConfig == 40 || trainConfig == 41 || trainConfig == 41 || trainConfig == 43 || trainConfig == 44 ||
		  trainConfig == 45 || trainConfig == 46 )
	  ){
	    
	    if (doWeighting){
	      if (generatorName.CompareTo("DPMJET")==0){
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      } else if (generatorName.CompareTo("HIJING")==0){   
               analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      }
	    }
	  } else if ( trainConfig == 7 || trainConfig == 8 || trainConfig == 10 || trainConfig == 12  || trainConfig == 15 ){
	    
	      if (doWeighting){
		  analysisCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
	      }
	  }
	  if( ! analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data()) ) {
            cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
            return 0;
	  } else { 
	  	  	  
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
   
   if(enableSetProdVtxVGamma) task->SetProductionVertextoVGamma(kTRUE);
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
