AliAnalysisTask *AddTask_GammaConvV1(TString collisionSystem = "pp"){
                       
   //get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_GammaConvV1", "No analysis manager found.");
      return 0;
   }
   
   Bool_t IsHeavyIon=collisionSystem.Contains("PbPb");

   TString cutnumber = "";
	TString cutnumberMeson = "";
   if(IsHeavyIon){ 
		cutnumber = "1080000002084001001500000";cutnumberMeson = "01021035000"; 
	} else{
		cutnumber = "0000000002084001001500000";cutnumberMeson = "01631035000"; 
	}
	
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
   if(trainConfig.Contains("PbPb")) numberOfCuts = 21;
   else if(trainConfig.Contains("pPb")) numberOfCuts = 1;
	else numberOfCuts = 9;
	
   TString *cutarray = new TString[numberOfCuts];
	TString *mesonCutArray = new TString[numberOfCuts];
   if(trainConfig.Contains("PbPb")){
		// Standard Cuts 
//       cutarray[0] = "1000000042092970023220000"; mesonCutArray[0] = "01022045000";  //standard cut Pi0 PbPb 00-100
// 			cutarray[0] = "1000001042092970023220000"; mesonCutArray[0] = "01022045000";  //standard cut Pi0 PbPb 00-100
//       cutarray[2] = "1000002042092970023220000"; mesonCutArray[2] = "01022045000";  //standard cut Pi0 PbPb 00-100
// 			cutarray[0] = "3010001042092970023220000"; mesonCutArray[0] = "01022045000";  //standard cut Pi0 PbPb 00-05
// 			cutarray[1] = "3120001042092970023220000"; mesonCutArray[1] = "01022045000";  //standard cut Pi0 PbPb 05-10
// 			cutarray[2] = "1010001042092970023220000"; mesonCutArray[2] = "01022045000";  //standard cut Pi0 PbPb 00-10
// 			cutarray[3] = "1120001042092970023220000"; mesonCutArray[3] = "01022045000";  //standard cut Pi0 PbPb 10-20
// // 			cutarray[4] = "1020001042092970023220000"; mesonCutArray[4] = "01022045000";  //standard cut Pi0 PbPb 00-20
// 			cutarray[4] = "1240001042092970023220000"; mesonCutArray[4] = "01022045000";  //standard cut Pi0 PbPb 20-40     
// 			cutarray[5] = "1460001042092970023220000"; mesonCutArray[5] = "01022065000";  //standard cut Pi0 PbPb 40-60
// 			cutarray[6] = "1680001042092970023220000"; mesonCutArray[6] = "01022065000";  //standard cut Pi0 PbPb 60-80      
// 
// 	// 	CutStudies more BG events 20 events
// 			cutarray[7] = "3010001042092970023220000"; mesonCutArray[7] = "01322045000";  
// 			cutarray[8] = "3120001042092970023220000"; mesonCutArray[8] = "01322045000"; 
// 			cutarray[9] = "1010001042092970023220000"; mesonCutArray[9] = "01322045000"; 
// 			cutarray[10] = "1120001042092970023220000"; mesonCutArray[10] = "01322045000"; 
// // 			cutarray[12] = "1020001042092970023220000"; mesonCutArray[12] = "01322045000"; 
// 			cutarray[11] = "1240001042092970023220000"; mesonCutArray[11] = "01322045000"; 
// 			cutarray[12] = "1460001042092970023220000"; mesonCutArray[12] = "01322065000";
// 			cutarray[13] = "1680001042092970023220000"; mesonCutArray[13] = "01322065000"; 

	// 	CutStudies more BG events 50 events
// 			cutarray[0] = "3010001042092970023220000"; mesonCutArray[0] = "01522045000"; 
// 			cutarray[1] = "3120001042092970023220000"; mesonCutArray[1] = "01522045000"; 
// 			cutarray[2] = "1010001042092970023220000"; mesonCutArray[2] = "01522045000"; 
// 			cutarray[3] = "1120001042092970023220000"; mesonCutArray[3] = "01522045000"; 
// // 			cutarray[4] = "1020001042092970023220000"; mesonCutArray[4] = "01522045000"; 
// 			cutarray[4] = "1240001042092970023220000"; mesonCutArray[4] = "01522045000"; 
// 			cutarray[5] = "1460001042092970023220000"; mesonCutArray[5] = "01522065000"; 
// 			cutarray[6] = "1680001042092970023220000"; mesonCutArray[6] = "01522065000"; 

		// Cutstudies 0-5% dEdx
		//	cutarray[7] = "3010001042093970023220000"; mesonCutArray[7] = "01022045000"; 
		//	cutarray[8] = "3010001042096970023220000"; mesonCutArray[8] = "01022045000"; 
		//	cutarray[9] = "3010001042092470023220000"; mesonCutArray[9] = "01022045000"; 
		//	cutarray[10] = "3010001042092770023220000"; mesonCutArray[10] = "01022045000"; 
		//	cutarray[11] = "3010001042092950023220000"; mesonCutArray[11] = "01022045000"; 
		// // Cutstudies 5-10% dEdx
		//	cutarray[12] = "3120001042093970023220000"; mesonCutArray[12] = "01022045000"; 
		//	cutarray[13] = "3120001042096970023220000"; mesonCutArray[13] = "01022045000"; 
// 			cutarray[7] = "3120001042092470023220000"; mesonCutArray[7] = "01022045000"; 
// 			cutarray[8] = "3120001042092770023220000"; mesonCutArray[8] = "01022045000"; 
// 			cutarray[9] = "3120001042092950023220000"; mesonCutArray[9] = "01022045000"; 
// 		// Cutstudies 0-10% dEdx
// 			cutarray[10] = "1010001042093970023220000"; mesonCutArray[10] = "01022045000"; 
// 			cutarray[11] = "1010001042096970023220000"; mesonCutArray[11] = "01022045000"; 
// 			cutarray[12] = "1010001042092470023220000"; mesonCutArray[12] = "01022045000"; 
// 			cutarray[13] = "1010001042092770023220000"; mesonCutArray[13] = "01022045000"; 
			cutarray[0] = "1010001042092950023220000"; mesonCutArray[0] = "01022045000"; 
		// Cutstudies 10-20% dEdx
			cutarray[1] = "1120001042093970023220000"; mesonCutArray[1] = "01022045000"; 
			cutarray[2] = "1120001042096970023220000"; mesonCutArray[2] = "01022045000";
			cutarray[3] = "1120001042092470023220000"; mesonCutArray[3] = "01022045000"; 
			cutarray[4] = "1120001042092770023220000"; mesonCutArray[4] = "01022045000"; 
			cutarray[5] = "1120001042092950023220000"; mesonCutArray[5] = "01022045000"; 
// 		Cutstudies 40-60% dEdx
			cutarray[6] = "1460001042093970023220000"; mesonCutArray[6] = "01022065000"; 
			cutarray[7] = "1460001042096970023220000"; mesonCutArray[7] = "01022065000"; 
			cutarray[8] = "1460001042092470023220000"; mesonCutArray[8] = "01022065000"; 
			cutarray[9] = "1460001042092770023220000"; mesonCutArray[9] = "01022065000"; 
			cutarray[10] = "1460001042092950023220000"; mesonCutArray[10] = "01022065000"; 
		// Cutstudies 60-80% dEdx
			cutarray[11] = "1680001042093970023220000"; mesonCutArray[11] = "01022065000"; 
			cutarray[12] = "1680001042096970023220000"; mesonCutArray[12] = "01022065000"; 
			cutarray[13] = "1680001042092470023220000"; mesonCutArray[13] = "01022065000";
			cutarray[14] = "1680001042092770023220000"; mesonCutArray[14] = "01022065000"; 
			cutarray[15] = "1680001042092950023220000"; mesonCutArray[15] = "01022065000"; 
		// Cutstudies 20-40% dEdx
			cutarray[16] = "1240001042093970023220000"; mesonCutArray[16] = "01022045000"; 
			cutarray[17] = "1240001042096970023220000"; mesonCutArray[17] = "01022045000";
			cutarray[18] = "1240001042092470023220000"; mesonCutArray[18] = "01022045000"; 
			cutarray[19] = "1240001042092770023220000"; mesonCutArray[19] = "01022045000"; 
			cutarray[20] = "1240001042092950023220000"; mesonCutArray[20] = "01022045000"; 

		// Cutstudies 0-5% TOF
//			cutarray[5] = "3010001042092970033220000"; mesonCutArray[5] = "01022045000"; 
//			cutarray[6] = "3010001042092970043220000"; mesonCutArray[6] = "01022045000"; 
		// Cutstudies 5-10% TOF
//			cutarray[7] = "3120001042092970033220000"; mesonCutArray[7] = "01022045000"; 
//			cutarray[8] = "3120001042092970043220000"; mesonCutArray[8] = "01022045000"; 
		// Cutstudies 0-10% TOF
//			cutarray[9] = "1010001042092970033220000"; mesonCutArray[9] = "01022045000"; 
//			cutarray[10] = "1010001042092970043220000"; mesonCutArray[10] = "01022045000"; 
		// Cutstudies 10-20% TOF
//			cutarray[11] = "1120001042092970033220000"; mesonCutArray[11] = "01022045000"; 
//			cutarray[12] = "1120001042092970043220000"; mesonCutArray[12] = "01022045000"; 
		// Cutstudies 20-40% TOF
//			cutarray[13] = "1240001042092970033220000"; mesonCutArray[13] = "01022045000"; 
//			cutarray[14] = "1240001042092970043220000"; mesonCutArray[14] = "01022045000"; 
		// Cutstudies 40-60% TOF
//			cutarray[15] = "1460001042092970033220000"; mesonCutArray[15] = "01022065000"; 
//			cutarray[16] = "1460001042092970043220000"; mesonCutArray[16] = "01022065000"; 
		// Cutstudies 60-80% TOF
//			cutarray[17] = "1680001042092970033220000"; mesonCutArray[17] = "01022065000"; 
//			cutarray[18] = "1680001042092970043220000"; mesonCutArray[18] = "01022065000"; 

		// Cutstudies 0-5% Alpha
//			cutarray[0] = "3010001042092970023220000"; mesonCutArray[0] = "01022085000"; 
//			cutarray[1] = "3010001042092970023220000"; mesonCutArray[1] = "01022005000"; 
		// Cutstudies 5-10% Alpha
//			cutarray[2] = "3120001042092970023220000"; mesonCutArray[2] = "01022085000"; 
//			cutarray[3] = "3120001042092970023220000"; mesonCutArray[3] = "01022005000"; 
		// Cutstudies 0-10% Alpha
//			cutarray[4] = "1010001042092970023220000"; mesonCutArray[4] = "01022085000"; 
//			cutarray[5] = "1010001042092970023220000"; mesonCutArray[5] = "01022005000"; 
		// Cutstudies 10-20% Alpha
//			cutarray[6] = "1120001042092970023220000"; mesonCutArray[6] = "01022085000"; 
//			cutarray[7] = "1120001042092970023220000"; mesonCutArray[7] = "01022005000"; 
		// Cutstudies 20-40% Alpha
//       cutarray[8] = "1240001042092970023220000"; mesonCutArray[8] = "01022085000";
//       cutarray[9] = "1240001042092970023220000"; mesonCutArray[9] = "01022005000"; 
		// Cutstudies 60-80% Alpha
//       cutarray[10] = "1680001042092970023220000"; mesonCutArray[10] = "01022075000"; 
//       cutarray[11] = "1680001042092970023220000"; mesonCutArray[11] = "01022055000"; 
		// Cutstudies 40-60% Alpha
//       cutarray[12] = "1460001042092970023220000"; mesonCutArray[12] = "01022075000"; 
//       cutarray[13] = "1460001042092970023220000"; mesonCutArray[13] = "01022055000"; 

		// Cutstudies 0-5% Qt
//			cutarray[14] = "3010001042092970024220000"; mesonCutArray[14] = "01022045000"; 
//			cutarray[15] = "3010001042092970022220000"; mesonCutArray[15] = "01022045000"; 
		// Cutstudies 5-10% Qt
//			cutarray[16] = "3120001042092970024220000"; mesonCutArray[16] = "01022045000"; 
//			cutarray[17] = "3120001042092970022220000"; mesonCutArray[17] = "01022045000"; 
		// Cutstudies 0-10% Qt
//			cutarray[0] = "1010001042092970024220000"; mesonCutArray[0] = "01022045000"; 
//			cutarray[1] = "1010001042092970022220000"; mesonCutArray[1] = "01022045000"; 
		// Cutstudies 10-20% Qt
//			cutarray[2] = "1120001042092970024220000"; mesonCutArray[2] = "01022045000"; 
//			cutarray[3] = "1120001042092970022220000"; mesonCutArray[3] = "01022045000"; 
		// Cutstudies 20-40% Qt
//			cutarray[4] = "1240001042092970024220000"; mesonCutArray[4] = "01022045000"; 
//			cutarray[5] = "1240001042092970022220000"; mesonCutArray[5] = "01022045000"; 
		// Cutstudies 40-60% Qt
//			cutarray[6] = "1460001042092970024220000"; mesonCutArray[6] = "01022065000"; 
//			cutarray[7] = "1460001042092970022220000"; mesonCutArray[7] = "01022065000"; 
		// Cutstudies 60-80% Qt
//			cutarray[8] = "1680001042092970024220000"; mesonCutArray[8] = "01022065000"; 
//			cutarray[9] = "1680001042092970022220000"; mesonCutArray[9] = "01022065000"; 

		// Cutstudies 0-5% Single Pt
//			cutarray[10] = "3010001042492970023220000"; mesonCutArray[10] = "01022045000";
//			cutarray[11] = "3010001042192970023220000"; mesonCutArray[11] = "01022045000";
		// Cutstudies 5-10% Single Pt
//			cutarray[12] = "3120001042492970023220000"; mesonCutArray[12] = "01022045000"; 
//			cutarray[13] = "3120001042192970023220000"; mesonCutArray[13] = "01022045000"; 
		// Cutstudies 0-10% Single Pt
//			cutarray[14] = "1010001042492970023220000"; mesonCutArray[14] = "01022045000"; 
//			cutarray[15] = "1010001042192970023220000"; mesonCutArray[15] = "01022045000"; 
		// Cutstudies 10-20% Single Pt
//			cutarray[16] = "1120001042492970023220000"; mesonCutArray[16] = "01022045000"; 
//			cutarray[17] = "1120001042192970023220000"; mesonCutArray[17] = "01022045000"; 
		// Cutstudies 20-40% Single Pt
//			cutarray[0] = "1240001042492970023220000"; mesonCutArray[0] = "01022045000"; 
//			cutarray[1] = "1240001042192970023220000"; mesonCutArray[1] = "01022045000"; 
		// Cutstudies 40-60% Single Pt
//			cutarray[2] = "1460001042492970023220000"; mesonCutArray[2] = "01022065000"; 
//			cutarray[3] = "1460001042192970023220000"; mesonCutArray[3] = "01022065000"; 
		// Cutstudies 60-80% Single Pt
//			cutarray[4] = "1680001042492970023220000"; mesonCutArray[4] = "01022065000"; 
//			cutarray[5] = "1680001042192970023220000"; mesonCutArray[5] = "01022065000"; 

		// Cutstudies 0-5% Chi2 Gamma
//			cutarray[6] = "3010001042092970023120000"; mesonCutArray[6] = "01022045000"; 
//			cutarray[7] = "3010001042092970023820000"; mesonCutArray[7] = "01022045000"; 
		// Cutstudies 5-10% Chi2 Gamma
//			cutarray[8] = "3120001042092970023120000"; mesonCutArray[8] = "01022045000"; 
//			cutarray[9] = "3120001042092970023820000"; mesonCutArray[9] = "01022045000"; 
		// Cutstudies 0-10% Chi2 Gamma
//			cutarray[10] = "1010001042092970023120000"; mesonCutArray[10] = "01022045000"; 
//			cutarray[11] = "1010001042092970023820000"; mesonCutArray[11] = "01022045000"; 
		// Cutstudies 10-20% Chi2 Gamma
//			cutarray[12] = "1120001042092970023120000"; mesonCutArray[12] = "01022045000"; 
//			cutarray[13] = "1120001042092970023820000"; mesonCutArray[13] = "01022045000"; 
		// Cutstudies 20-40% Chi2 Gamma
//			cutarray[14] = "1240001042092970023120000"; mesonCutArray[14] = "01022045000"; 
//			cutarray[15] = "1240001042092970023820000"; mesonCutArray[15] = "01022045000"; 
		// Cutstudies 40-60% Chi2 Gamma
//			cutarray[16] = "1460001042092970023120000"; mesonCutArray[16] = "01022065000"; 
//			cutarray[17] = "1460001042092970023820000"; mesonCutArray[17] = "01022065000"; 
		// Cutstudies 60-80% Chi2 Gamma
//			cutarray[18] = "1680001042092970023120000"; mesonCutArray[18] = "01022065000"; 
//			cutarray[19] = "1680001042092970023820000"; mesonCutArray[19] = "01022065000"; 

		// Cutstudies 0-5% TPC Cluster
//			cutarray[0] = "3010001042062970023220000"; mesonCutArray[0] = "01022045000"; 
//			cutarray[1] = "3010001042082970023220000"; mesonCutArray[1] = "01022045000";
		// Cutstudies 5-10% TPC Cluster
//			cutarray[2] = "3120001042062970023220000"; mesonCutArray[2] = "01022045000"; 
//			cutarray[3] = "3120001042082970023220000"; mesonCutArray[3] = "01022045000"; 
		// Cutstudies 0-10% TPC Cluster
//			cutarray[4] = "1010001042062970023220000"; mesonCutArray[4] = "01022045000"; 
//			cutarray[5] = "1010001042082970023220000"; mesonCutArray[5] = "01022045000"; 
		// Cutstudies 10-20% TPC Cluster
//			cutarray[6] = "1120001042062970023220000"; mesonCutArray[6] = "01022045000"; 
//			cutarray[7] = "1120001042082970023220000"; mesonCutArray[7] = "01022045000"; 
		// Cutstudies 20-40% TPC Cluster
//			cutarray[8] = "1240001042062970023220000"; mesonCutArray[8] = "01022045000"; 
//			cutarray[9] = "1240001042082970023220000"; mesonCutArray[9] = "01022045000"; 
		// Cutstudies 40-60% TPC Cluster
//			cutarray[10] = "1460001042062970023220000"; mesonCutArray[10] = "01022065000"; 
//			cutarray[11] = "1460001042082970023220000"; mesonCutArray[11] = "01022065000"; 
		// Cutstudies 60-80% TPC Cluster
//			cutarray[12] = "1680001042062970023220000"; mesonCutArray[12] = "01022065000"; 
//			cutarray[13] = "1680001042082970023220000"; mesonCutArray[13] = "01022065000"; 

		// Cutstudies 0-5% Psi Pair
//			cutarray[14] = "3010001042092970023210000"; mesonCutArray[14] = "01022045000"; 
//			cutarray[15] = "3010001042092970023230000"; mesonCutArray[15] = "01022045000"; 
		// Cutstudies 5-10% Psi Pair
//			cutarray[16] = "3120001042092970023210000"; mesonCutArray[16] = "01022045000"; 
//			cutarray[17] = "3120001042092970023230000"; mesonCutArray[17] = "01022045000"; 
		// Cutstudies 0-10% Psi Pair
//			cutarray[0] = "1010001042092970023210000"; mesonCutArray[0] = "01022045000"; 
//			cutarray[1] = "1010001042092970023230000"; mesonCutArray[1] = "01022045000"; 
		// Cutstudies 10-20% Psi Pair
//			cutarray[2] = "1120001042092970023210000"; mesonCutArray[2] = "01022045000"; 
//			cutarray[3] = "1120001042092970023230000"; mesonCutArray[3] = "01022045000";
		// Cutstudies 20-40% Psi Pair
//			cutarray[4] = "1240001042092970023210000"; mesonCutArray[4] = "01022045000"; 
//			cutarray[5] = "1240001042092970023230000"; mesonCutArray[5] = "01022045000"; 
		// Cutstudies 40-60% Psi Pair
//			cutarray[6] = "1460001042092970023210000"; mesonCutArray[6] = "01022065000"; 
//			cutarray[7] = "1460001042092970023230000"; mesonCutArray[7] = "01022065000"; 
		// Cutstudies 60-80% Psi Pair
//			cutarray[8] = "1680001042092970023210000"; mesonCutArray[8] = "01022065000"; 
//			cutarray[9] = "1680001042092970023230000"; mesonCutArray[9] = "01022065000"; 

		// Cutstudies 0-5% R Cut
//			cutarray[10] = "3010001041092970023220000"; mesonCutArray[10] = "01022045000"; 
//			cutarray[11] = "3010001044092970023220000"; mesonCutArray[11] = "01022045000"; 
		// Cutstudies 5-10% R Cut
//			cutarray[12] = "3120001041092970023220000"; mesonCutArray[12] = "01022045000"; 
//			cutarray[13] = "3120001044092970023220000"; mesonCutArray[13] = "01022045000";
		// Cutstudies 0-10% R Cut
//			cutarray[14] = "1010001041092970023220000"; mesonCutArray[14] = "01022045000"; 
//			cutarray[15] = "1010001044092970023220000"; mesonCutArray[15] = "01022045000"; 
		// Cutstudies 10-20% R Cut
//			cutarray[16] = "1120001041092970023220000"; mesonCutArray[16] = "01022045000"; 
//			cutarray[17] = "1120001044092970023220000"; mesonCutArray[17] = "01022045000"; 
		// Cutstudies 20-40% R Cut
//			cutarray[0] = "1240001041092970023220000"; mesonCutArray[0] = "01022045000"; 
//			cutarray[1] = "1240001044092970023220000"; mesonCutArray[1] = "01022045000"; 
		// Cutstudies 40-60% R Cut
//			cutarray[2] = "1460001041092970023220000"; mesonCutArray[2] = "01022065000"; 
//			cutarray[3] = "1460001044092970023220000"; mesonCutArray[3] = "01022065000"; 
		// Cutstudies 60-80% R Cut
//			cutarray[4] = "1680001041092970023220000"; mesonCutArray[4] = "01022065000"; 
//			cutarray[5] = "1680001044092970023220000"; mesonCutArray[5] = "01022065000"; 
	} else if(trainConfig.Contains("pPb")){ //pA needs thighter rapidity cut y < 0.5
		 cutarray[0] = "1000000042092172023220000"; mesonCutArray[0] = "01024045000";  //standard cut Pi0 PbPb 00-100
   } else {
//        cutarray[0] = "0000010002093663003800000"; mesonCutArray[0] = "01631031009"; //Standard Cut Pi0 with PileUp Rejection
//  		// V0 finder      
// // 		 cutarray[1] = "0000010102093663003800000"; mesonCutArray[1] = "01631031009"; //Standard Cut Pi0 with PileUp Rejection
//        cutarray[1] = "0000011002093663003800000"; mesonCutArray[1] = "01631031009"; //Standard Cut only MinBias in MC
// 		 cutarray[2] = "0000012002093663003800000"; mesonCutArray[2] = "01631031009"; //Standard Cut only added Signals in MC		 
// //       cutarray[1] = "0000010002093260043800000"; mesonCutArray[0] = "01631031009";  //Standard Cut Gamma with PileUp Rejection
// 
// 		// dEdx e-line
//       cutarray[3] = "0000011002091663003800000"; mesonCutArray[3] = "01631031009"; 
//       cutarray[4] = "0000011002092663003800000"; mesonCutArray[4] = "01631031009";
// 		// dEdx pi-line
// 		cutarray[5] = "0000011002093863003800000"; mesonCutArray[5] = "01631031009"; 
// 		cutarray[6] = "0000011002093553003800000"; mesonCutArray[6] = "01631031009"; 
// 		cutarray[7] = "0000011002093100003800000"; mesonCutArray[7] = "01631031009"; 
// 		cutarray[8] = "0000011002093674003800000"; mesonCutArray[8] = "01631031009"; 
// 		// Alpha 
// 		cutarray[9] = "0000011002093663003800000"; mesonCutArray[9] = "01631001009"; 
// 		// chi2 gamma
//       cutarray[10] = "0000011002093663003900000"; mesonCutArray[10] = "01631031009"; 
//       cutarray[11] = "0000011002093663003200000"; mesonCutArray[11] = "01631031009"; 
// 		// pt single e+/e-
//       cutarray[12] = "0000011002493663003800000"; mesonCutArray[12] = "01631031009"; 
//       cutarray[13] = "0000011002193663003800000"; mesonCutArray[13] = "01631031009"; 
		// findable cluster TPC
		cutarray[0] = "0000010002003663003800000"; mesonCutArray[0] = "01631031009"; 
		cutarray[1] = "0000010002083663003800000"; mesonCutArray[1] = "01631031009"; 
		// qt variation
		cutarray[2] = "0000010002093663004800000"; mesonCutArray[2] = "01631031009"; 
		cutarray[3] = "0000010002093663002800000"; mesonCutArray[3] = "01631031009"; 
		// TOF
		cutarray[4] = "0000010002093663023800000"; mesonCutArray[4] = "01631031009"; 
		cutarray[5] = "0000010002093663033800000"; mesonCutArray[5] = "01631031009"; 
		cutarray[6] = "0000010002093663043800000"; mesonCutArray[6] = "01631031009"; 
		// track sharing
		cutarray[7] = "0000010002093663003800013"; mesonCutArray[7] = "01631031009"; 
		cutarray[8] = "0000010002093663003800012"; mesonCutArray[8] = "01631031009"; 
  }
  
   TList *ConvCutList = new TList();
	TList *MesonCutList = new TList();

   TList *HeaderList = new TList();
   TObjString *Header1 = new TObjString("PARAM");
   HeaderList->Add(Header1);
   TObjString *Header2 = new TObjString("BOX");
   HeaderList->Add(Header2);
   TObjString *Header3 = new TObjString("Pythia");
   HeaderList->Add(Header3);
   
   ConvCutList->SetOwner(kTRUE);
   AliConversionCuts **analysisCuts = new AliConversionCuts*[numberOfCuts];
	MesonCutList->SetOwner(kTRUE);
   AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

   for(Int_t i = 0; i<numberOfCuts; i++){
      analysisCuts[i] = new AliConversionCuts();
      analysisCuts[i]->InitializeCutsFromCutString(cutarray[i].Data());
 		if (trainConfig.Contains("pPb")) analysisCuts[i]->SelectCollisionCandidates(AliVEvent::kCINT5);
      ConvCutList->Add(analysisCuts[i]);
		
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
//      if(i == 2 ){
 //        analysisCuts[i]->SetAcceptedHeader(HeaderList);
//      }
      analysisMesonCuts[i] = new AliConversionMesonCuts();
      analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
      MesonCutList->Add(analysisMesonCuts[i]);
		analysisMesonCuts[i]->SetFillCutHistograms("",kFALSE);
 		if (i < 7 && trainConfig.Contains("PbPb")){
			AddQATaskV1(analysisCuts[i],IsHeavyIon);
 			TString addoutput=gSystem->Getenv("ADD_OUTPUT_FILES");
 			if (addoutput.Length()) addoutput+=",";
 			addoutput+=Form("GammaConvV1_QA_%s.root",((analysisCuts[i])->GetCutNumber()).Data());
 			gSystem->Setenv("ADD_OUTPUT_FILES",addoutput.Data());
     	}
   }
   task->SetConversionCutList(numberOfCuts,ConvCutList);
	task->SetMesonCutList(numberOfCuts,MesonCutList);
   task->SetMoveParticleAccordingToVertex(kTRUE);
   task->SetDoMesonAnalysis(kTRUE);
   mgr->AddTask(task);
  
   //connect containers
   AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("GammaConvV1", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"GammaConvV1.root");
   AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("kkoch_tree",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "kkoch_default");
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
