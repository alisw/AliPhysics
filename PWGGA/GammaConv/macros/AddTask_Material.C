void AddTask_Material(TString V0ReaderCutNumber = "0000000060084001001500000000",
                      TString TaskCutnumber = "0000000090092663743800000000",
                      Bool_t IsMC = kFALSE, 
                      Int_t IsHeavyIon = 0, 
                      TString cutnumberAODBranch = "0000000060084001001500000",
                      Bool_t doEtaShiftV0Reader = kFALSE 
                     ){
  
   // Suitable Cutnumbers for the V0 Reader for
   // PbPb: V0ReaderCutNumber =  "1000000060084001001500000000"; (V0Mult MC)
   //  or   V0ReaderCutNumber =  "5000000060084001001500000000" (TPC mult MC)
   // pPb: V0ReaderCutNumber =   "8000000060084001001500000000";
   // pp: V0ReaderCutNumber =    "0000000060084001001500000000";


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
   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
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
      if(V0ReaderCutNumber!=""){
         fCuts= new AliConversionCuts(V0ReaderCutNumber.Data(),V0ReaderCutNumber.Data());
         fCuts->SetPreSelectionCutFlag(kTRUE);
         if(fCuts->InitializeCutsFromCutString(V0ReaderCutNumber.Data())){
            if (IsHeavyIon==2){
               fCuts->SelectCollisionCandidates(AliVEvent::kINT7);
               fCuts->DoEtaShift(doEtaShiftV0Reader);
            }
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

   } else {
      Error("AddTask_V0ReaderV1", "Cannot execute AddTask, V0ReaderV1 already exists.");
   }   


   // suitable cuts for the material Task:
   // PbPb:  TaskCutnumber = "5680001060092663044803000000"; TPC mult in MC - 60-80% central
   //   or:  TaskCutnumber = "1680001060092663044803000000"; V0 mult in MC  - 60-80% central
   //  pPb:  TaskCutnumber = "8000000090092663743800000000";
   //   pp:  TaskCutnumber = "0000000090092663743800000000";

   AliConversionCuts *analysisCuts = new AliConversionCuts();
   analysisCuts->InitializeCutsFromCutString(TaskCutnumber.Data());
   analysisCuts->SetFillCutHistograms("",kFALSE);
   
   AliAnalysisTaskMaterial *fMaterial= new AliAnalysisTaskMaterial(Form("%s_Material",(analysisCuts->GetCutNumber()).Data()));
   fMaterial->SetConversionCuts(analysisCuts,IsHeavyIon);
   fMaterial->SetIsMC(IsMC);
   mgr->AddTask(fMaterial);
   
   AliAnalysisDataContainer *coutput1 =
   mgr->CreateContainer(Form("GammaConvMaterial_%s",TaskCutnumber.Data()), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("GammaConv_Material_%s.root",TaskCutnumber.Data()));

   AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
   mgr->ConnectInput(fMaterial,  0, cinput1 );
   mgr->ConnectOutput (fMaterial,  1, coutput1);
   //connect containers
   return ;
}
