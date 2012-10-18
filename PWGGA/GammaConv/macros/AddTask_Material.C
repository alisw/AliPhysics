AliAnalysisTask *AddTask_Material(){
                       
   //get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_Material", "No analysis manager found.");
      return 0;
   }
   
   TString trainConfig=gSystem->Getenv("CONFIG_FILE");
	cout << trainConfig.Data() << endl;
   Bool_t IsHeavyIon=trainConfig.Contains("PbPb");

   TString cutnumber = "";
   if(IsHeavyIon) cutnumber = "1080000020084001001500000";
   else cutnumber = "0000000020084001001500000"; 
	
   //========= Add PID Reponse to ANALYSIS manager ====
   if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AddTaskPIDResponse();
   }
   
   //========= Add V0 Reader to  ANALYSIS manager =====
   AliV0ReaderV1 *fV0ReaderV1=new AliV0ReaderV1("V0ReaderV1");
   ConfigV0ReaderV1(fV0ReaderV1,cutnumber,IsHeavyIon);
   mgr->AddTask(fV0ReaderV1);
  
   AliLog::SetGlobalLogLevel(AliLog::kInfo);

   //================================================
   //              data containers
   //================================================
   //            find input container
   //below the trunk version
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
   //connect input V0Reader
   mgr->ConnectInput (fV0ReaderV1,0,cinput);
   
//   TString cutarray = "0000010020093663003800000"; //Standard Cut Pi0 with PileUp Rejection
//   TString cutarray = "0000010020092663003800000"; //Standard Cut Pi0 with PileUp Rejection
    TString cutarray = "0000010020092663043800000"; //Standard Cut Pi0 with PileUp Rejection

   TList *ConvCutList = new TList();
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts *analysisCuts = new AliConversionCuts();
	analysisCuts->InitializeCutsFromCutString(cutarray.Data());
	ConvCutList->Add(analysisCuts);
	analysisCuts->SetFillCutHistograms("",kFALSE);
	AliAnalysisTaskMaterial *fMaterial= new AliAnalysisTaskMaterial(Form("%s_Material",(analysisCuts->GetCutNumber()).Data()));
	fMaterial->SetConversionCuts(analysisCuts,IsHeavyIon);
	TString addoutput=gSystem->Getenv("ADD_OUTPUT_FILES");
	if (addoutput.Length()) addoutput+=",";
	addoutput+=Form("GammaConvV1_Material_%s.root",(analysisCuts->GetCutNumber()).Data());
	if (trainConfig.Contains("MC"){
		if (addoutput.Length()) addoutput+=",";
		addoutput+=Form("GammaConvV1_Resolution_%s.root",(analysisCuts->GetCutNumber()).Data());
	}
	AddQATaskV1(analysisCuts,IsHeavyIon);
	if (addoutput.Length()) addoutput+=",";
	addoutput+=Form("GammaConvV1_QA_%s.root",(analysisCuts->GetCutNumber()).Data());
	gSystem->Setenv("ADD_OUTPUT_FILES",addoutput.Data());
	
	cout<<"Adding addoutput.Data()"<<endl; 

	mgr->AddTask(fMaterial);
	mgr->ConnectInput(fMaterial,  0, cinput );
	
   //connect containers
   return fMaterial;
}

void ConfigV0ReaderV1(AliV0ReaderV1 *fV0ReaderV1,TString analysiscut="",Bool_t IsHeavyIon=kFALSE){

   fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
   fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
   fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
   }
   AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
   
   if(inputHandler->IsA()==AliESDInputHandler::Class()){
      // ESD mode
   }
   
   if(inputHandler->IsA()==AliAODInputHandler::Class()){
      // AOD mode
      //   task->SetUseSatelliteAODs(kTRUE);
   }

   // Set AnalysisCut Number
   AliConversionCuts *fCuts=NULL;
   if(analysiscut!=""){
      fCuts= new AliConversionCuts(analysiscut.Data(),analysiscut.Data());
      if(fCuts->InitializeCutsFromCutString(analysiscut.Data())){
         fV0ReaderV1->SetConversionCuts(fCuts);
         fCuts->SetFillCutHistograms("",kTRUE);
      }
   }
   fV0ReaderV1->Init();
   
}

void AddQATaskV1(AliConversionCuts *ConversionCuts, Bool_t IsHeavyIon){
   
   //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
	Error("AddTask_V0ReaderV1QA", "No analysis manager found.");
	return 0;
    }
    AliAnalysisTaskConversionQA *fQA = new AliAnalysisTaskConversionQA(Form("%s_QA",(ConversionCuts->GetCutNumber()).Data()));
    fQA->SetConversionCuts(ConversionCuts,IsHeavyIon);
    mgr->AddTask(fQA);
    
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    mgr->ConnectInput (fQA,0,cinput);
}
