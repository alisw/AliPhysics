void AddTask_GammaConvDalitzV1_pp(  Int_t trainConfig = 1,  //change different set of cuts
				    Bool_t isMC   = kFALSE, //run MC 
				    TString fileNameInputForWeighting = "MCSpectraInput.root" // path to file for weigting input
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
      Error(Form("AddTask_GammaConvDalitzV1_%i",trainConfig), "No analysis manager found.");
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
   
		        
   TString cutnumber = "00000000000840010015000000"; 
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
      fV0ReaderV1->Init();

      AliLog::SetGlobalLogLevel(AliLog::kInfo);

      //connect input V0Reader
      mgr->AddTask(fV0ReaderV1);
      mgr->ConnectInput(fV0ReaderV1,0,cinput);

   }

   
   if( !(AliDalitzElectronSelector*)mgr->GetTask("ElectronSelector") ){



   AliDalitzElectronSelector *fElectronSelector = new AliDalitzElectronSelector("ElectronSelector");

   //ConfigV0ReaderV1(fV0ReaderV1,ConvCutnumber,IsHeavyIon);



   // Set AnalysisCut Number

   AliDalitzElectronCuts *fElecCuts=0;

   TString ElecCuts = "900054000000020000";



   if( ElecCuts!=""){

       fElecCuts= new AliDalitzElectronCuts(ElecCuts.Data(),ElecCuts.Data());

      if(fElecCuts->InitializeCutsFromCutString(ElecCuts.Data())){

         fElectronSelector->SetDalitzElectronCuts(fElecCuts);

         fElecCuts->SetFillCutHistograms("",kTRUE);

      }

   }

   fElectronSelector->Init();
 }



   mgr->AddTask(fElectronSelector);





   //================================================
   //========= Add task to the ANALYSIS manager =====
   //            find input container
   AliAnalysisTaskGammaConvDalitzV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvDalitzV1(Form("GammaConvDalitzV1_%i",trainConfig));
   task->SetIsHeavyIon(0);
   task->SetIsMC(isMC);
   
   
   // Cut Numbers to use in Analysis
   Int_t numberOfCuts = 2;

   TString *ConvCutarray 		= new TString[numberOfCuts];
   TString *MesonCutarray 		= new TString[numberOfCuts];
   TString *ElecCutarray    		= new TString[numberOfCuts];
   
            
   
   

   if(trainConfig == 1){
     //TOF PID
     ConvCutarray[0] = "00000110020936630278000000"; MesonCutarray[0] = "01631031009";ElecCutarray[0] = "904784032531026210";  //TOF[-3,5] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0
     ConvCutarray[1] = "00000110020936630278000000"; MesonCutarray[1] = "01631031009";ElecCutarray[1] = "904784042531026210";  //TOF[-2,3] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0     
    } else if (trainConfig == 2) {
     //TOF PID
     ConvCutarray[0] = "00000110020936630278000000"; MesonCutarray[0] = "01631031009";ElecCutarray[0] = "904784032531026210";  //TOF[-3,5] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0
     ConvCutarray[1] = "00000110020936630278000000"; MesonCutarray[1] = "01631031009";ElecCutarray[1] = "904784042531026210";  //TOF[-2,3] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0     
    } else if (trainConfig == 3) {
     //TOF PID
     ConvCutarray[0] = "00000110020936630278000000"; MesonCutarray[0] = "01631031009";ElecCutarray[0] = "904784032531026210";  //TOF[-3,5] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0
     ConvCutarray[1] = "00000110020936630278000000"; MesonCutarray[1] = "01631031009";ElecCutarray[1] = "904784042531026210";  //TOF[-2,3] 0.0 sigmas at low Pt for pion rejection, Pt 0.125 cut,  DCAxy Pt Dep, No Mass(e+,e-)  FindCluster > 0.0     
    } else {
      Error(Form("GammaConvDalitzV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
   }

   TList  *ConvCutList  = new TList();
   TList  *MesonCutList = new TList();
   TList  *ElecCutList  = new TList();

   TList *HeaderList = new TList();
   TObjString *Header2 = new TObjString("BOX");
   HeaderList->Add(Header2);

   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts 		= new AliConversionCuts*[numberOfCuts];
   MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts 	= new AliConversionMesonCuts*[numberOfCuts];
   ElecCutList->SetOwner(kTRUE);
   AliDalitzElectronCuts **analysisElecCuts 	= new AliDalitzElectronCuts*[numberOfCuts];
   


   for(Int_t i = 0; i<numberOfCuts; i++){
     
     
      TString cutName( Form("%s_%s_%s",ConvCutarray[i].Data(),ElecCutarray[i].Data(),MesonCutarray[i].Data() ) );
    
      analysisCuts[i] = new AliConversionCuts();
      analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data());
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
   
      analysisMesonCuts[i] = new AliConversionMesonCuts();
      analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data());
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
      
      
      analysisElecCuts[i] = new AliDalitzElectronCuts();
      analysisElecCuts[i]->InitializeCutsFromCutString(ElecCutarray[i].Data());
      ElecCutList->Add(analysisElecCuts[i]);
      analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
      //analysisElecCuts[i]->PrintCuts();





 
      
      analysisCuts[i]->SetAcceptedHeader(HeaderList);
     
   }
   
 
   task->SetConversionCutList(numberOfCuts,ConvCutList);
   task->SetMesonCutList(MesonCutList);
   task->SetElectronCutList(ElecCutList);
   task->SetMoveParticleAccordingToVertex(kTRUE);
   //task->SetDoMesonAnalysis(kTRUE);
   //if (enableQAMesonTask) task->SetDoMesonQA(kTRUE); //Attention new switch for Pi0 QA
   //if (enableQAMesonTask) task->SetDoPhotonQA(kTRUE);  //Attention new switch small for Photon QA

   //connect containers
   AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(Form("GammaConvDalitzV1_%i",trainConfig), TList::Class(),
                           AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

   mgr->AddTask(task);
   mgr->ConnectInput(task,0,cinput);
   mgr->ConnectOutput(task,1,coutput);

   return;

}
