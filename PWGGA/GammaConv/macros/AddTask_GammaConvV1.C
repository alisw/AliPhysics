AliAnalysisTask *AddTask_GammaConvV1(TString collisionSystem = "pp"){
                       
   //get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_GammaConvV1", "No analysis manager found.");
      return 0;
   }
   
   Bool_t IsHeavyIon=collisionSystem.Contains("PbPb");

   TString cutnumber = "";
   if(IsHeavyIon) cutnumber = "900400508050113211200001080000000";
   else cutnumber = "900400508050113211360000000000000"; 
	
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
   
   //================================================
   //========= Add task to the ANALYSIS manager =====
   //================================================
   //              data containers
   //================================================
   //            find input container
   AliAnalysisTaskGammaConvV1 *task=NULL;
   task= new AliAnalysisTaskGammaConvV1("GammaConvV1");
   task->SetIsHeavyIon(IsHeavyIon); 
   // Cut Numbers to use in Analysis
   
   Int_t numberOfCuts = 1;
   if(IsHeavyIon) numberOfCuts = 8;
   else numberOfCuts = 1;

   TString *cutarray = new TString[numberOfCuts];
   if(IsHeavyIon){
		// Standard Cuts 
			cutarray[0] = "900297209450304221200003012002000"; //standard cut Pi0 PbPb 00-05
			cutarray[1] = "900297209450304221200003122002000"; //standard cut Pi0 PbPb 05-10
			cutarray[2] = "900297209450304221200001012002000"; //standard cut Pi0 PbPb 00-10
			cutarray[3] = "900297209450304221200001122002000"; //standard cut Pi0 PbPb 10-20
			cutarray[4] = "900297209450304221200001022002000"; //standard cut Pi0 PbPb 00-20
			cutarray[5] = "900297209450304221200001242002000"; //standard cut Pi0 PbPb 20-40     
			cutarray[6] = "900297209450306221200001462002000"; //standard cut Pi0 PbPb 40-60
			cutarray[7] = "900297209450306221200001682002000"; //standard cut Pi0 PbPb 60-80      

	// 	CutStudies more BG events 20 events
//			cutarray[0] = "900297209450304221230003012002000"; 
//			cutarray[1] = "900297209450304221230003122002000"; 
//			cutarray[2] = "900297209450304221230001012002000"; 
//			cutarray[3] = "900297209450304221230001122002000"; 
//			cutarray[4] = "900297209450304221230001022002000"; 
//			cutarray[5] = "900297209450304221230001242002000"; 
//			cutarray[6] = "900297209450306221230001462002000"; 
//			cutarray[7] = "900297209450306221230001682002000"; 

	// 	CutStudies more BG events 50 events
//			cutarray[0] = "900297209450304221250003012002000"; 
//			cutarray[1] = "900297209450304221250003122002000"; 
//			cutarray[2] = "900297209450304221250001012002000"; 
//			cutarray[3] = "900297209450304221250001122002000"; 
//			cutarray[4] = "900297209450304221250001022002000"; 
//			cutarray[5] = "900297209450304221250001242002000"; 
//			cutarray[6] = "900297209450306221250001462002000"; 
//			cutarray[7] = "900297209450306221250001682002000"; 

		// Cutstudies 0-5% dEdx
//			cutarray[0] = "900397209450304221200003012002000"; 
//			cutarray[1] = "900697209450304221200003012002000"; 
//			cutarray[2] = "900247209450304221200003012002000"; 
//			cutarray[3] = "900277209450304221200003012002000"; 
//			cutarray[4] = "900295209450304221200003012002000"; 
		// Cutstudies 5-10% dEdx
//			cutarray[0] = "900397209450304221200003122002000"; 
//			cutarray[1] = "900697209450304221200003122002000"; 
//			cutarray[2] = "900247209450304221200003122002000"; 
//			cutarray[3] = "900277209450304221200003122002000"; 
//			cutarray[4] = "900295209450304221200003122002000"; 
		// Cutstudies 0-10% dEdx
//			cutarray[0] = "900397209450304221200001012002000"; 
//			cutarray[1] = "900697209450304221200001012002000"; 
//			cutarray[2] = "900247209450304221200001012002000"; 
//			cutarray[3] = "900277209450304221200001012002000"; 
//			cutarray[4] = "900295209450304221200001012002000"; 
		// Cutstudies 10-20% dEdx
//			cutarray[0] = "900397209450304221200001122002000"; 
//			cutarray[1] = "900697209450304221200001122002000"; 
//			cutarray[2] = "900247209450304221200001122002000"; 
//			cutarray[3] = "900277209450304221200001122002000"; 
//			cutarray[4] = "900295209450304221200001122002000"; 
		// Cutstudies 40-60% dEdx
//			cutarray[0] = "900397209450306221200001462002000"; 
//			cutarray[1] = "900697209450306221200001462002000"; 
//			cutarray[2] = "900247209450306221200001462002000"; 
//			cutarray[3] = "900277209450306221200001462002000"; 
//			cutarray[4] = "900295209450306221200001462002000"; 
		// Cutstudies 60-80% dEdx
//			cutarray[0] = "900397209450306221200001682002000"; 
//			cutarray[1] = "900697209450306221200001682002000"; 
//			cutarray[2] = "900247209450306221200001682002000"; 
//			cutarray[3] = "900277209450306221200001682002000"; 
//			cutarray[4] = "900295209450306221200001682002000"; 

		// Cutstudies 0-5% TOF
//			cutarray[0] = "900297209450304221200003013002000"; 
//			cutarray[1] = "900297209450304221200003014002000"; 
		// Cutstudies 5-10% TOF
//			cutarray[0] = "900297209450304221200003123002000"; 
//			cutarray[1] = "900297209450304221200003124002000"; 
		// Cutstudies 0-10% TOF
//			cutarray[0] = "900297209450304221200001013002000"; 
//			cutarray[1] = "900297209450304221200001014002000"; 
		// Cutstudies 10-20% TOF
//			cutarray[0] = "900297209450304221200001123002000"; 
//			cutarray[1] = "900297209450304221200001124002000"; 
		// Cutstudies 20-40% TOF
//			cutarray[0] = "900297209450304221200001243002000"; 
//			cutarray[1] = "900297209450304221200001244002000"; 
		// Cutstudies 40-60% TOF
//			cutarray[0] = "900297209450306221200001463002000"; 
//			cutarray[1] = "900297209450306221200001464002000"; 
		// Cutstudies 60-80% TOF
//			cutarray[0] = "900297209450306221200001683002000"; 
//			cutarray[1] = "900297209450306221200001684002000"; 

		// Cutstudies 0-5% Alpha
//			cutarray[0] = "900297209450308221200003012002000"; 
//			cutarray[1] = "900297209450300221200003012002000"; 
		// Cutstudies 5-10% Alpha
//			cutarray[0] = "900297209450308221200003122002000"; 
//			cutarray[1] = "900297209450300221200003122002000"; 
		// Cutstudies 0-10% Alpha
//			cutarray[0] = "900297209450308221200001012002000"; 
//			cutarray[1] = "900297209450300221200001012002000"; 
		// Cutstudies 10-20% Alpha
//			cutarray[0] = "900297209450308221200001122002000"; 
//			cutarray[1] = "900297209450300221200001122002000"; 
		// Cutstudies 20-40% Alpha
//       cutarray[0] = "900297209450308221200001242002000";
//       cutarray[1] = "900297209450300221200001242002000";
		// Cutstudies 40-60% Alpha
//       cutarray[0] = "900297209450307221200001682002000";
//       cutarray[1] = "900297209450305221200001682002000";
		// Cutstudies 60-80% Alpha
//       cutarray[8] = "900297209450307221200001462002000";
//       cutarray[9] = "900297209450305221200001462002000";

		// Cutstudies 0-5% Qt
//			cutarray[0] = "900297209450404221200003012002000"; 
//			cutarray[1] = "900297209450204221200003012002000"; 
		// Cutstudies 5-10% Qt
//			cutarray[0] = "900297209450404221200003122002000"; 
//			cutarray[1] = "900297209450204221200003122002000"; 
		// Cutstudies 0-10% Qt
//			cutarray[0] = "900297209450404221200001012002000"; 
//			cutarray[1] = "900297209450204221200001012002000"; 
		// Cutstudies 10-20% Qt
//			cutarray[0] = "900297209450404221200001122002000"; 
//			cutarray[1] = "900297209450204221200001122002000"; 
		// Cutstudies 20-40% Qt
//			cutarray[0] = "900297209450404221200001242002000"; 
//			cutarray[1] = "900297209450204221200001242002000"; 
		// Cutstudies 40-60% Qt
//			cutarray[0] = "900297209450406221200001462002000"; 
//			cutarray[1] = "900297209450206221200001462002000"; 
		// Cutstudies 60-80% Qt
//			cutarray[0] = "900297209450406221200001682002000"; 
//			cutarray[1] = "900297209450206221200001682002000"; 

		// Cutstudies 0-5% Single Pt
//			cutarray[0] = "900297249450304221200003012002000"; 
//			cutarray[1] = "900297219450304221200003012002000"; 
		// Cutstudies 5-10% Single Pt
//			cutarray[0] = "900297249450304221200003122002000"; 
//			cutarray[1] = "900297219450304221200003122002000"; 
		// Cutstudies 0-10% Single Pt
//			cutarray[0] = "900297249450304221200001012002000"; 
//			cutarray[1] = "900297219450304221200001012002000"; 
		// Cutstudies 10-20% Single Pt
//			cutarray[0] = "900297249450304221200001122002000"; 
//			cutarray[1] = "900297219450304221200001122002000"; 
		// Cutstudies 20-40% Single Pt
//			cutarray[0] = "900297249450304221200001242002000"; 
//			cutarray[1] = "900297219450304221200001242002000"; 
		// Cutstudies 40-60% Single Pt
//			cutarray[0] = "900297249450306221200001462002000"; 
//			cutarray[1] = "900297219450306221200001462002000"; 
		// Cutstudies 60-80% Single Pt
//			cutarray[0] = "900297249450306221200001682002000"; 
//			cutarray[1] = "900297219450306221200001682002000"; 

		// Cutstudies 0-5% Chi2 Gamma
//			cutarray[0] = "900297109450304221200003012002000"; 
//			cutarray[1] = "900297809450304221200004012002000"; 
		// Cutstudies 5-10% Chi2 Gamma
//			cutarray[0] = "900297109450304221200003122002000"; 
//			cutarray[1] = "900297809450304221200003122002000"; 
		// Cutstudies 0-10% Chi2 Gamma
//			cutarray[0] = "900297109450304221200001012002000"; 
//			cutarray[1] = "900297809450304221200001012002000"; 
		// Cutstudies 10-20% Chi2 Gamma
//			cutarray[0] = "900297109450304221200001122002000"; 
//			cutarray[1] = "900297809450304221200001122002000"; 
		// Cutstudies 20-40% Chi2 Gamma
//			cutarray[0] = "900297109450304221200001242002000"; 
//			cutarray[1] = "900297809450304221200001242002000"; 
		// Cutstudies 40-60% Chi2 Gamma
//			cutarray[0] = "900297109450306221200001462002000"; 
//			cutarray[1] = "900297809450306221200001462002000"; 
		// Cutstudies 60-80% Chi2 Gamma
//			cutarray[0] = "900297109450306221200001682002000"; 
//			cutarray[1] = "900297809450306221200001682002000"; 

		// Cutstudies 0-5% TPC Cluster
//			cutarray[0] = "900297206450304221200003012002000"; 
//			cutarray[1] = "900297208450304221200003012002000"; 
		// Cutstudies 5-10% TPC Cluster
//			cutarray[0] = "900297206450304221200003122002000"; 
//			cutarray[1] = "900297208450304221200003122002000"; 
		// Cutstudies 0-10% TPC Cluster
//			cutarray[0] = "900297206450304221200001012002000"; 
//			cutarray[1] = "900297208450304221200001012002000"; 
		// Cutstudies 10-20% TPC Cluster
//			cutarray[0] = "900297206450304221200001122002000"; 
//			cutarray[1] = "900297208450304221200001122002000"; 
		// Cutstudies 20-40% TPC Cluster
//			cutarray[0] = "900297206450304221200001242002000"; 
//			cutarray[1] = "900297208450304221200001242002000"; 
		// Cutstudies 40-60% TPC Cluster
//			cutarray[0] = "900297206450306221200001462002000"; 
//			cutarray[1] = "900297208450306221200001462002000"; 
		// Cutstudies 60-80% TPC Cluster
//			cutarray[0] = "900297206450306221200001682002000"; 
//			cutarray[1] = "900297208450306221200001682002000"; 

		// Cutstudies 0-5% Psi Pair
//			cutarray[0] = "900297209450304221200003012001000"; 
//			cutarray[1] = "900297209450304221200003012003000"; 
		// Cutstudies 5-10% Psi Pair
//			cutarray[0] = "900297209450304221200003122001000"; 
//			cutarray[1] = "900297209450304221200003122003000"; 
		// Cutstudies 0-10% Psi Pair
//			cutarray[0] = "900297209450304221200001012001000"; 
//			cutarray[1] = "900297209450304221200001012003000"; 
		// Cutstudies 10-20% Psi Pair
//			cutarray[0] = "900297209450304221200001122001000"; 
//			cutarray[1] = "900297209450304221200001122003000"; 
		// Cutstudies 20-40% Psi Pair
//			cutarray[0] = "900297209450304221200001242001000"; 
//			cutarray[1] = "900297209450304221200001242003000"; 
		// Cutstudies 40-60% Psi Pair
//			cutarray[0] = "900297209450306221200001462001000"; 
//			cutarray[1] = "900297209450306221200001462003000"; 
		// Cutstudies 60-80% Psi Pair
//			cutarray[0] = "900297209450306221200001682001000"; 
//			cutarray[1] = "900297209450306221200001682003000"; 

		// Cutstudies 0-5% R Cut
//			cutarray[0] = "900297209450304121200003012002000"; 
//			cutarray[1] = "900297209450304421200003012002000"; 
		// Cutstudies 5-10% R Cut
//			cutarray[0] = "900297209450304121200003122002000"; 
//			cutarray[1] = "900297209450304421200003122002000"; 
		// Cutstudies 0-10% R Cut
//			cutarray[0] = "900297209450304121200001012002000"; 
//			cutarray[1] = "900297209450304421200001012002000"; 
		// Cutstudies 10-20% R Cut
//			cutarray[0] = "900297209450304121200001122002000"; 
//			cutarray[1] = "900297209450304421200001122002000"; 
		// Cutstudies 20-40% R Cut
//			cutarray[0] = "900297209450304121200001242002000"; 
//			cutarray[1] = "900297209450304421200001242002000"; 
		// Cutstudies 40-60% R Cut
//			cutarray[0] = "900297209450306121200001462002000"; 
//			cutarray[1] = "900297209450306421200001462002000"; 
		// Cutstudies 60-80% R Cut
//			cutarray[0] = "900297209450306121200001682002000"; 
//			cutarray[1] = "900297209450306421200001682002000"; 

   } else {
       cutarray[0] = "900366809010333211361000000900000"; //Standard Cut Pi0 with PileUp Rejection
//       cutarray[1] = "900226609010303211361000004900000"; //Standard Cut Gamma with PileUp Rejection
// 		// V0 finder      
//		cutarray[2] = "910366809010333211361000000900000"; //Standard Cut Pi0 with PileUp Rejection V0 offline
// 		// dEdx e-line
//       cutarray[2] = "900166809010333211361000000900000"; 
//       cutarray[3] = "900266809010333211361000000900000"; 
// 		// dEdx pi-line
// 		cutarray[4] = "900386809010333211361000000900000";
// 		cutarray[5] = "900355809010333211361000000900000";
// 		cutarray[6] = "900310809010303211361000000900000";
// 		cutarray[7] = "900367809010343211361000000900000";
// 		// Alpha 
// 		cutarray[8] = "900366809010330211361000000900000";
// 		// chi2 gamma
//       cutarray[9] = "900366909010333211361000000900000";
//       cutarray[10] = "900366209010333211361000000900000";
// 		// pt single e+/e-
//       cutarray[11] = "900366849010333211361000000900000";
//       cutarray[12] = "900366819010333211361000000900000";
// 		// findable cluster TPC
// 		cutarray[0] = "900366800010333211361000000900000";
// 		cutarray[1] = "900366808010333211361000000900000";
// 		// qt variation
// 		cutarray[2] = "900366809010433211361000000900000";
// 		cutarray[3] = "900366809010233211361000000900000";
// 		// TOF
// 		cutarray[4] = "900366809010333211361000002900000";
// 		cutarray[5] = "900366809010333211361000003900000";
// 		cutarray[6] = "900366809010333211361000004900000";
// 		// track sharing
// 		cutarray[7] = "900366809010333211361000000900013";
// 		cutarray[8] = "900366809010333211361000000900012";

      // cutarray[1] = "9006266090103032113600000049000";
      // cutarray[2] = "9007266090103032113600000049000";
      // cutarray[3] = "9006267090103032113600000049000";
      // cutarray[4] = "9006268090103032113600000049000";
      // cutarray[5] = "9006265090103032113600000049000";
      // cutarray[6] = "9006265090103032113600000049010";
      // cutarray[7] = "9006265090103032113600000049020";
      // cutarray[8] = "9006265090103032113600000049030";
      // cutarray[9] = "9003662080100332113600000009000";
      // cutarray[0] = "9006166090103032113600000049000";
      // cutarray[1] = "9006266090103032113600000049000";
      // cutarray[2] = "9006366090103032113600000049000";
      // cutarray[3] = "9006206090103032113600000049000";
      // cutarray[4] = "9006666090103032113600000049000";
      // cutarray[5] = "9006566090103032113600000049000";
      // cutarray[6] = "9006266090104032113600000049000";
      // cutarray[7] = "9006266090102032113600000049000";
      // cutarray[8] = "9006266090103032113600000049020";
      // cutarray[9] = "9003662080100332113600000009000";
      // cutarray[0] = "9006266090103032113600000049000";
      // cutarray[1] = "9006267090103032113600000049000";
      // cutarray[2] = "9006268090103032113600000049000";
      // cutarray[3] = "9006266190103032113600000049000";
      // cutarray[4] = "9006266290103032113600000049000";
      // cutarray[5] = "9006266080103032113600000049000";
      // cutarray[6] = "9006266091103032113600000049000";
      // cutarray[7] = "9006266091403032113600000049000";
      // cutarray[8] = "9002265090103032113600000049000";
      // cutarray[9] = "9003662080100332113600000009000";
      
      // cutarray[1] = "9003662080100332113600000009000";
      // cutarray[2] = "9006266090103031113600000049000";
      // cutarray[3] = "9006266090103035113600000049000";
      // cutarray[4] = "9006266090103032113600000029000";
      // cutarray[5] = "9006266090103032113600000039000";
      // cutarray[6] = "9002266090103032113600000049000";
      // cutarray[7] = "9003266090103032113600000049000";
      // cutarray[8] = "9002267090103032113600000049000";
      // cutarray[9] = "9002268090103032113600000049000";
   }
  
   TList *ConvCutList = new TList();
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
   for(Int_t i = 0; i<numberOfCuts; i++){
      analysisCuts[i] = new AliConversionCuts();
      analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
//      if(i == 0 ){ //|| i == 1
//       if(!IsHeavyIon){
//          AddQATaskV1(analysisCuts[i],IsHeavyIon);
//          TString addoutput=gSystem->Getenv("ADD_OUTPUT_FILES");
//          if (addoutput.Length()) addoutput+=",";
//          addoutput+=Form("GammaConvV1_QA_%s.root",((analysisCuts[i])->GetCutNumber()).Data());
//          gSystem->Setenv("ADD_OUTPUT_FILES",addoutput.Data());
//          cout<<"Adding addoutput.Data()"<<endl;
 //      }
//		}
   }
   task->SetConversionCutList(numberOfCuts,ConvCutList);
   task->SetMoveParticleAccordingToVertex(kTRUE);
   task->SetDoMesonAnalysis(kTRUE);
   mgr->AddTask(task);
  
   //connect containers
   AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"GammaConvV1.root");
   AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("gammaConv_tree",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "gammaConv_default");
   mgr->ConnectInput  (task,  0, cinput );
   mgr->ConnectOutput (task,  0, coutput0);
   mgr->ConnectOutput (task,  1, coutput1);
   return task;
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
