void AddTask_GammaConvDalitzV1_pPb(    Int_t trainConfig = 1,
                                       Bool_t isMC       = kFALSE, //run MC 
                                       Bool_t enableQAMesonTask = kTRUE, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                       Bool_t enableDoMesonChic = kFALSE, // enable additional Chic analysis
				       TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                       Bool_t doWeighting = kFALSE,  //enable Weighting
                                       TString generatorName = "DPMJET",				
                                       TString cutnumberAODBranch = "0000000060084001001500000"
                                  ) {


   
   cout<<"Entro -1"<<endl;

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
   Int_t numberOfCuts = 6;

   TString *ConvCutarray    = new TString[numberOfCuts];

   TString *ElecCutarray    = new TString[numberOfCuts];

   TString *MesonCutarray   = new TString[numberOfCuts];

   Bool_t doEtaShiftIndCuts = kFALSE;
   TString stringShift = "";

   // Shifting in pPb direction

   doEtaShiftIndCuts = kTRUE;
   stringShift = "pPb";


   if( trainConfig == 1 ) {

        ConvCutarray[0] = "8000011082093603007200000000"; ElecCutarray[0] = "9047540025810262170"; MesonCutarray[0] = "01039035009000"; //standard cut Pi0 PbPb 00-100
        ConvCutarray[1] = "8000011082093603007200000000"; ElecCutarray[1] = "9047540025810261170"; MesonCutarray[1] = "01039035009000"; //standard cut Pi0 PbPb 00-100 + Single Pt primary > 0.100 GeV
        ConvCutarray[2] = "8000011082094603007200000000"; ElecCutarray[2] = "9047540025810262170"; MesonCutarray[2] = "01039035009000"; //standard cut Pi0 PbPb 00-100 + dEdx electron gamma   -6 ,7 sigmas
        ConvCutarray[3] = "8000011082093603007203000000"; ElecCutarray[3] = "9047540025810262170"; MesonCutarray[3] = "01039035009000"; //standard cut Pi0 PbPb 00-100  do Aysemtri cut
        ConvCutarray[4] = "8000011082093603007200000000"; ElecCutarray[4] = "9051540025810262170"; MesonCutarray[4] = "01039035009000"; //standard cut Pi0 PbPb 00-100
        ConvCutarray[5] = "8000011082093603007200000000"; ElecCutarray[5] = "9051540025810262170"; MesonCutarray[5] = "01039035009000"; //standard cut Pi0 PbPb 00-100 Standard cut + dEdx primary -3, 5 and  3.0  , -10 pion rejection

   } 
  else if( trainConfig == 2 ) {

        ConvCutarray[0] = "8000011082093603007200000000"; ElecCutarray[0] = "9047540025810262170"; MesonCutarray[0] = "01039035009000"; //standard cut Pi0 PbPb 00-100
        ConvCutarray[1] = "8020011082093603007200000000"; ElecCutarray[1] = "9047540025810262170"; MesonCutarray[1] = "01039035009000"; //standard cut Pi0 PbPb 00-20
        ConvCutarray[2] = "8240011082093603007200000000"; ElecCutarray[2] = "9047540025810262170"; MesonCutarray[2] = "01039035009000"; //standard cut Pi0 PbPb 20-40
        ConvCutarray[3] = "8460011082093603007200000000"; ElecCutarray[3] = "9047540025810262170"; MesonCutarray[3] = "01039035009000"; //standard cut Pi0 PbPb 40-60
        ConvCutarray[4] = "8680011082093603007200000000"; ElecCutarray[4] = "9047540025810262170"; MesonCutarray[4] = "01039035009000"; //standard cut Pi0 PbPb 60-80        
        ConvCutarray[5] = "8600011082093603007200000000"; ElecCutarray[5] = "9047540025810262170"; MesonCutarray[5] = "01039035009000"; //standard cut Pi0 PbPb 60-100

 }

else if( trainConfig == 3 ) {

        ConvCutarray[0] = "8000011082093603007200000000"; ElecCutarray[0] = "9047540025810262171"; MesonCutarray[0] = "01039035009000"; //standard cut Pi0 PbPb 00-100
        ConvCutarray[1] = "8020011082093603007200000000"; ElecCutarray[1] = "9047540025810262171"; MesonCutarray[1] = "01039035009000"; //standard cut Pi0 PbPb 00-20
        ConvCutarray[2] = "8240011082093603007200000000"; ElecCutarray[2] = "9047540025810262171"; MesonCutarray[2] = "01039035009000"; //standard cut Pi0 PbPb 20-40
        ConvCutarray[3] = "8460011082093603007200000000"; ElecCutarray[3] = "9047540025810262171"; MesonCutarray[3] = "01039035009000"; //standard cut Pi0 PbPb 40-60
        ConvCutarray[4] = "8680011082093603007200000000"; ElecCutarray[4] = "9047540025810262171"; MesonCutarray[4] = "01039035009000"; //standard cut Pi0 PbPb 60-80        
        ConvCutarray[5] = "8600011082093603007200000000"; ElecCutarray[5] = "9047540025810262171"; MesonCutarray[5] = "01039035009000"; //standard cut Pi0 PbPb 60-100

 } else if( trainConfig == 4 ) {

        ConvCutarray[0] = "8000012082093603007200000000"; ElecCutarray[0] = "9047540025810262171"; MesonCutarray[0] = "01039035009000"; //standard cut Pi0 PbPb 00-100
        ConvCutarray[1] = "8020012082093603007200000000"; ElecCutarray[1] = "9047540025810262171"; MesonCutarray[1] = "01039035009000"; //standard cut Pi0 PbPb 00-20
        ConvCutarray[2] = "8240012082093603007200000000"; ElecCutarray[2] = "9047540025810262171"; MesonCutarray[2] = "01039035009000"; //standard cut Pi0 PbPb 20-40
        ConvCutarray[3] = "8460012082093603007200000000"; ElecCutarray[3] = "9047540025810262171"; MesonCutarray[3] = "01039035009000"; //standard cut Pi0 PbPb 40-60
        ConvCutarray[4] = "8680012082093603007200000000"; ElecCutarray[4] = "9047540025810262171"; MesonCutarray[4] = "01039035009000"; //standard cut Pi0 PbPb 60-80        
        ConvCutarray[5] = "8600012082093603007200000000"; ElecCutarray[5] = "9047540025810262171"; MesonCutarray[5] = "01039035009000"; //standard cut Pi0 PbPb 60-100
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
      }
   else {

   if ( trainConfig == 3 ){

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
      else if (trainConfig == 4 ){

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
