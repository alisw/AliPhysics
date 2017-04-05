
AliAnalysisTaskEMCalHFEpA2* ConfigEMCalHFEpA(
											

Bool_t isMC=kFALSE, 
Int_t triggerIndex=0, 
Int_t configIndex=0, 
Int_t centralityIndex=0, 
Bool_t isAOD = kFALSE,
Bool_t isEMCal = kFALSE,
Bool_t isTrigger = kFALSE,
Int_t EMCalThreshould = 0, //0 == EG1, 1 == EG2
Bool_t isTender = kFALSE,
char* period = "b",
Int_t   centralityEstimator = 0,
Bool_t isCentralitySys 		= kFALSE,
Bool_t isTOFdet 		= kFALSE
											 
)

{
///_______________________________________________________________________________________________________________
///Track selection: Cuts used to ensure a minimum quality level of the tracks selected to perform the analysis
	AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsMinBias","HFE Cuts");
	hfecuts->CreateStandardCuts();
	
	//TPC Cuts
	hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);	
	if(configIndex==1) 	hfecuts->SetMinNClustersTPC(90);			                //Minimum number of clusters on TPC
	else if(configIndex==2) 	hfecuts->SetMinNClustersTPC(110);
	else if(configIndex==3) 	hfecuts->SetMinNClustersTPC(80);
	else if(configIndex==4) 	hfecuts->SetMinNClustersTPC(85);
	else if(configIndex==5) 	hfecuts->SetMinNClustersTPC(115);
	else if(configIndex==6) 	hfecuts->SetMinNClustersTPC(120);
	else if(configIndex==80) 	hfecuts->SetMinNClustersTPC(130);
	else if(configIndex==81) 	hfecuts->SetMinNClustersTPC(140);//Minimum number of clusters on TPC
	else hfecuts->SetMinNClustersTPC(100);							                //Minimum number of clusters on TPC
	
	if(configIndex==7) hfecuts->SetMinNClustersTPCPID(70); 
	else if (configIndex==8) hfecuts->SetMinNClustersTPCPID(90);
	else if (configIndex==9) hfecuts->SetMinNClustersTPCPID(60);
	else if (configIndex==10) hfecuts->SetMinNClustersTPCPID(65);
	else if (configIndex==11) hfecuts->SetMinNClustersTPCPID(100);
	else if (configIndex==12) hfecuts->SetMinNClustersTPCPID(95);
	else hfecuts->SetMinNClustersTPCPID(80);										//Minimum number of clusters for dE/dx
	
	hfecuts->SetMinRatioTPCclusters(0.6);						                    //Number of clusters (Found/Findable)
	
	//ITS
	if(configIndex==13) hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);			//Require at least one cluster on SPD
	else if(configIndex==82) hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
	else hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);							//Require at least one cluster on SPD
	//hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 			                    //Require at least one cluster on SDD
	hfecuts->SetCheckITSLayerStatus(kFALSE); 
	
	
	
	
	if(configIndex==14) hfecuts->SetMinNClustersITS(2);								//Minimum number of clusters on ITS
	else if(configIndex==15) hfecuts->SetMinNClustersITS(4);	
	else if(configIndex==16) hfecuts->SetMinNClustersITS(1);
	else if(configIndex==17) hfecuts->SetMinNClustersITS(5);
	else hfecuts->SetMinNClustersITS(3);	
		
	if(!isEMCal){
		if(configIndex==13) hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);			//Require at least one cluster on SPD
		else if(configIndex==82) hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
		else hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);	
		
		if(configIndex==14) hfecuts->SetMinNClustersITS(3);								//Minimum number of clusters on ITS
		else if(configIndex==15) hfecuts->SetMinNClustersITS(4);	
		else if(configIndex==16) hfecuts->SetMinNClustersITS(1);
		else if(configIndex==17) hfecuts->SetMinNClustersITS(5);
		else hfecuts->SetMinNClustersITS(4);
	}//Minimum number of clusters on ITS
	
	//Additional Cuts
	if(!isEMCal){
		hfecuts->SetPtRange(0.5, 1e6);
	}
	else{	
		hfecuts->SetPtRange(2, 1e6);//Transversal momentum range in GeV/c
	}
	
	//testing this line for the DCA cut
	//hfecuts->SetRequireDCAToVertex();
	//DCA cut included in the analysis 12 March 2014
	if(configIndex==18) hfecuts->SetMaxImpactParam(2,3); //changed z to 3
	else if(configIndex==19) hfecuts->SetMaxImpactParam(0.5,1);//r,z   default on AOD 2.4, 3.2
	else if(configIndex==83) hfecuts->SetMaxImpactParam(0.1,0.2);
	else hfecuts->SetMaxImpactParam(1,2); 							                //DCA to vertex
	
	//Event Selection
	hfecuts->SetVertexRange(10.);													//
	//hfecuts->SetProductionVertex(0,0.3,0,0.3);									//
///_______________________________________________________________________________________________________________
	// new cuts for event selection
	
		// done inside the task, by hand
		// hfecuts->SetUseCorrelationVertex();
		// hfecuts->SetSPDVtxResolutionCut();
		// hfecuts->SetpApileupCut();

///_________________________________________________________________________________________________________________________
///Task config
	AliAnalysisTaskEMCalHFEpA2 *task = new AliAnalysisTaskEMCalHFEpA2(Form("HFECuts%d_%d_%d",triggerIndex,configIndex,centralityIndex));
	printf("task ------------------------ %p\n ", task);
	task->SetHFECuts(hfecuts);
	task->SetCorrelationAnalysis(kFALSE);
	task->SetAODanalysis(isAOD);
	task->SetEventMixing(kFALSE);
	
		//centrality
	task->SetCentralityEstimator(centralityEstimator);
	task->SetCentralitySys(isCentralitySys);
	//in order to use same task in pp analysis
	task->SetPPanalysis(kFALSE);
	
	//to separate trigger threshold
	if(EMCalThreshould==0 && triggerIndex==2) task->SetEMCalTriggerEG1();
	if(EMCalThreshould==1 && triggerIndex==2) task->SetEMCalTriggerEG2();
	
	if(isEMCal) task->SetUseEMCal();
		
	
	if(isTrigger){
		task->SetUseTrigger();
		task->SetUseShowerShapeCut(kTRUE);
				
		if(configIndex==120)task->SetM20Cut(0.0,0.2);
		else if(configIndex==121)task->SetM20Cut(0.0,0.4);
		else if(configIndex==122)task->SetM20Cut(0.01,0.3);
		else if(configIndex==123)task->SetM20Cut(0.0,2);
		else task->SetM20Cut(0.0,0.3);
		
		if(configIndex==124)task->SetM02Cut(0.0,0.2);
		else if(configIndex==125)task->SetM02Cut(0.0,0.4);
		else if(configIndex==126)task->SetM02Cut(0.0,0.5);
		else if(configIndex==127)task->SetM02Cut(0.0,0.6);
		else task->SetM02Cut(0.0,2);
			
		
	}
	
	if(period == "d"){
		
			printf("======================================================================================\n ");
			printf("\n\n Running on 13d period!!! WITH CALIBRATION \n\n  \n\n ");
			printf("======================================================================================\n ");
		
			task->SetTPCCalibration();
	        task->SetTPC_mean_sigma(1.02, 1.68);
			// task->SetTPC_mean_sigma(0.63, 1.17);

		    task->SetTPCcal_cut_min(-1);
		    task->SetTPCcal_cut_max(3);
			
	}
	
	if((period == "b" || period == "c") && !isEMCal && !isTOFdet){
		
		printf("======================================================================================\n ");
		printf("\n\n Running on 13b or 13c period!!! WITH CALIBRATION for TPCnsigma mean and width      \n\n  \n\n ");
		printf("======================================================================================\n ");
		
		task->SetTPCCalibration();
		task->SetTPC_mean_sigma(0.09, 1.03);
			// task->SetTPC_mean_sigma(0.63, 1.17);
		
		task->SetTPCcal_cut_min(0);
		task->SetTPCcal_cut_max(3);
	}
	
	
	if(period == "e" || period == "f"){
		task->SetTPCCalibration_eta(kTRUE);
			//task->SetTPCCalibration_eta(kFALSE);
	}
	 
	
	
	if(configIndex==300){
		task->SetTPCCalibration_eta(kFALSE);
	}
	
	
			
		//if(period == "d3"){
		//task->SetUseShowerShapeCut(kTRUE);
		//printf("Shower Shape on for d3 because is for the trigger data correction!!! Cannot run as trigger because the corrections with weight\n");
		//}
	
	if(isTender) task->SetUseTender();
		
		//not used anymore  - using hfe package for this cut
	if(configIndex==104)  task->SetdcaCut(2,3);//r,z
	else if(configIndex==105)  task->SetdcaCut(1.5,2.5);//r,z
	else if(configIndex==106)  task->SetdcaCut(0.5,1);//r,z
	else if(configIndex==107)  task->SetdcaCut(0.5,1.5);//r,z
	else if(configIndex==108)  task->SetdcaCut(0.1,0.2);//r,z
	//else task->SetdcaCut(1,2);//r,z
	

	
	task->SetBackground(kTRUE);
	
	//nonHFE cuts
	if(configIndex==20) task->SetNonHFEmassCut(0.05);
	else if(configIndex==21) task->SetNonHFEmassCut(0.10);
	else if(configIndex==22) task->SetNonHFEmassCut(0.03);
	else if(configIndex==23) task->SetNonHFEmassCut(0.18);
	else if(configIndex==24) task->SetNonHFEmassCut(0.01);
	else if(configIndex==25) task->SetNonHFEmassCut(0.2);
	else task->SetNonHFEmassCut(0.15);
	
	if(configIndex==26) task->SetNonHFEangleCut(0.1);
	if(configIndex==27) task->SetNonHFEangleCut(0.15);
	if(configIndex==28) task->SetNonHFEangleCut(0.05);
	
	//partner cuts
	
	if(configIndex==29) task->SetAdditionalCuts(0.1,80);
    else if(configIndex==30) task->SetAdditionalCuts(0.15,80);
    else if(configIndex==31) task->SetAdditionalCuts(0.2,80);
    else if(configIndex==32) task->SetAdditionalCuts(0.25,80);
	else if(configIndex==84) task->SetAdditionalCuts(0.3,80);
	else if(configIndex==85) task->SetAdditionalCuts(0.35,80);
	else if(configIndex==86) task->SetAdditionalCuts(0.4,80);
	else if(configIndex==87) task->SetAdditionalCuts(0.45,80);
    
	
	
    
	else if(configIndex==33) task->SetAdditionalCuts(0,60);
	else if(configIndex==34) task->SetAdditionalCuts(0,70);
	else if(configIndex==35) task->SetAdditionalCuts(0,90);
	else if(configIndex==36) task->SetAdditionalCuts(0,100);
	 
	else task->SetAdditionalCuts(0,80);
	 
	
	//eta cuts
	if(configIndex==40) task->SetEtaCut(-0.6,0);
	else if (configIndex==41) task->SetEtaCut(-0.5,0.1);
	else if (configIndex==42) task->SetEtaCut(0,0.6);
	else if (configIndex==43) task->SetEtaCut(-0.1,0.5);
	else if (configIndex==44) task->SetEtaCut(-0.5,0.5);
	else if (configIndex==45) task->SetEtaCut(-0.4,0.4);
	else if (configIndex==46) task->SetEtaCut(-0.3,0.3);
	else task->SetEtaCut(-0.6,0.6);
	
	//track matching cuts
	if(configIndex==50) task->SetdPhidEtaCut(0.02,0.02);
	else if (configIndex==51) task->SetdPhidEtaCut(0.03,0.03);
	else if (configIndex==52) task->SetdPhidEtaCut(0.04,0.04);
	else task->SetdPhidEtaCut(0.05,0.05);

	//E/p Cuts
	if (configIndex==60) task->SetEoverPCut(0.85,1.2);
	else if (configIndex==61) task->SetEoverPCut(0.70,1.2);
	else if (configIndex==62) task->SetEoverPCut(0.75,1.2);
	else if (configIndex==63) task->SetEoverPCut(0.76,1.2);
	else if (configIndex==64) task->SetEoverPCut(0.77,1.2);
	else if (configIndex==65) task->SetEoverPCut(0.78,1.2);
	else if (configIndex==66) task->SetEoverPCut(0.79,1.2);
	
	else if (configIndex==67) task->SetEoverPCut(0.81,1.2);
	else if (configIndex==68) task->SetEoverPCut(0.82,1.2);
	else if (configIndex==69) task->SetEoverPCut(0.83,1.2);
	
	else if (configIndex==88) task->SetEoverPCut(0.84,1.2);
	else if (configIndex==89) task->SetEoverPCut(0.86,1.2);
	
	else if (configIndex==90) task->SetEoverPCut(0.85,1.3);
	else if (configIndex==91) task->SetEoverPCut(0.80,1.3);
	else if (configIndex==92) task->SetEoverPCut(0.81,1.3);
	else if (configIndex==93) task->SetEoverPCut(0.78,1.3);
	
	
	else task->SetEoverPCut(0.80,1.2);
	
		
	//this line is to set the change on the E/p cut used in the efficiency calculations.
    //task->SetEoverPnsigma(kTRUE);
	
	
	if(centralityIndex==0) task->SetCentrality(0,20);
	if(centralityIndex==1) task->SetCentrality(20,40);
	if(centralityIndex==2) task->SetCentrality(40,60);
	if(centralityIndex==3) task->SetCentrality(60,80);
	if(centralityIndex==4) task->SetCentrality(80,100);
	if(centralityIndex==5) task->SetCentrality(0,100);
	
	if(centralityIndex==6) task->SetCentrality(0,10);
	if(centralityIndex==7) task->SetCentrality(10,20);
	if(centralityIndex==8) task->SetCentrality(20,30);
	if(centralityIndex==9) task->SetCentrality(30,40);
	if(centralityIndex==10) task->SetCentrality(40,50);
	if(centralityIndex==11) task->SetCentrality(50,60);
	if(centralityIndex==12) task->SetCentrality(60,70);
	if(centralityIndex==13) task->SetCentrality(70,80);
	if(centralityIndex==14) task->SetCentrality(80,90);
	if(centralityIndex==15) task->SetCentrality(90,100);
	
///_______________________________________________________________________________________________________________

///_______________________________________________________________________________________________________________
///Particle identification
	AliHFEpid *pid = task->GetPID();

//______________________________________
//In the case of a simulation
	if(isMC)
	{
	  pid->SetHasMCData(kTRUE);
	  task->SetMCanalysis();
	}
//______________________________________

//______________________________________________________
//Configure PID
	
	
	if(isTOFdet){	
		pid->AddDetector("TOF", 0); //Add TOF PID
		pid->ConfigureTOF(3.0); //Configure TOF cut: Defaut = 3 sigmas
	}
	
	//_________________________
	//TPC PID
	pid->AddDetector("TPC", 1);				//Add TPC PID
	
	
	
	//_________________________
	//Configure TPC cut
	//Defaul = -1 to 3 sigmas
	//Note that it is also possible to define a model instead of a constant
	//--------->For this change the "cut model"
	
	Double_t params[4];
	char *cutmodel;
	cutmodel = "pol0";
	
	if(configIndex==70){ 
		params[0] = -1.5;
		if(period == "d"){
			params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
			task->SetTPCcal_cut_min(-1.5);
		}
	}
	else if (configIndex==71){ 
		params[0] = -0.5;
		if(period == "d" || period == "b" || period == "c"){
			params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
			task->SetTPCcal_cut_min(-0.5);
		}
		
	}
	else if (configIndex==72){
		params[0] = 0;
		if(period == "d"){
			params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
			task->SetTPCcal_cut_min(0);
			
		}
	}	
	else if (configIndex==73){ 
		params[0] = 0.25;
		if(period == "d" || period == "b" || period == "c"){
			params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
			task->SetTPCcal_cut_min(0.25);
			
		}
	}
	else if (configIndex==74){ 
		params[0] = 0.5;
		if(period == "d"){
			params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
			task->SetTPCcal_cut_min(0.5);
			
		}
	}
	else{ 
		params[0] = -1;
	
		if(period =="e"||period =="f"){
				if(configIndex==200){
					params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
					task->SetTPCcal_cut_min(-1);
				}
		}
		if(period == "d"){
				params[0] = -3; //looser cut. Real cut applied inside the task define by line below:
				task->SetTPCcal_cut_min(-1);
			
		}
		
		if(!isEMCal)params[0] = 0;
		
		
		if((period == "b" || period =="c") && !isEMCal){
			params[0] = -1; //looser cut. Real cut applied inside the task define by line below:
			task->SetTPCcal_cut_min(0);
		}
		
		
	
	}
	if(isTOFdet){	
		params[0] = -1;

	}

	
		//=========================================
	
	
	if(configIndex==75)Double_t max=1.5;
	else if(configIndex==76)Double_t max=2.0;
	else if(configIndex==77){
		Double_t max=2.5;
		if(period =="d"){
			Double_t max=5;//looser cut for hfe package. Real cut inside the task.
			task->SetTPCcal_cut_max(2.5);
		}
	}
	else if(configIndex==78){
		Double_t max=3.5;
		if(period =="d"){
			Double_t max=6;//looser cut for hfe package. Real cut inside the task.
			task->SetTPCcal_cut_max(3.5);
		}
	}
	else if(configIndex==79){
		Double_t max=4.0;
		if(period =="d"){
			Double_t max=6;//looser cut for hfe package. Real cut inside the task.
			task->SetTPCcal_cut_max(4.0);
		}
	}
	else{
		Double_t max=3.0;
		if(period =="e"||period =="f"){
			
			if(configIndex==200){
				Double_t max=5;//looser cut for hfe package. Real cut inside the task.
				task->SetTPCcal_cut_max(3);
			}
			
		}
		if(period =="d"){
				Double_t max=5;//looser cut for hfe package. Real cut inside the task.
				task->SetTPCcal_cut_max(3);
		}
		
		if((period == "b" || period =="c") && !isEMCal && !isTOFdet){
			Double_t max=5;//looser cut for hfe package. Real cut inside the task.
			task->SetTPCcal_cut_max(3);
		}
	}
	
	
		
   //testing hadron contamination in E/p cut
	if(configIndex==94){
		params[0] = 0;
		task->SetEoverPCut(0.76,1.2);
	}
	if(configIndex==95){
		params[0] = 0;
		task->SetEoverPCut(0.78,1.2);
	}
	if(configIndex==96){
		params[0] = 0;
		task->SetEoverPCut(0.80,1.2);
	}
	if(configIndex==97){
		params[0] = 0;
		task->SetEoverPCut(0.82,1.2);
	}
	if(configIndex==98){
		params[0] = 0;
		task->SetEoverPCut(0.83,1.2);
	}
	if(configIndex==99){
		params[0] = 0;
		task->SetEoverPCut(0.84,1.2);
	}
	if(configIndex==100){
		params[0] = 0;
		task->SetEoverPCut(0.85,1.2);
	}
	if(configIndex==101){
		params[0] = 0;
		task->SetEoverPCut(0.86,1.2);
	}
	if(configIndex==102){
		params[0] = 0;
		task->SetEoverPCut(0.85,1.3);
	}
	if(configIndex==103){
		params[0] = 0;
		task->SetEoverPCut(0.86,1.3);
	}
	
	//testing hadron contamination with other TPCsignal : -0.5
	
	if(configIndex==109){
		params[0] = -0.5;
		task->SetEoverPCut(0.76,1.2);
	}
	if(configIndex==110){
		params[0] =  -0.5;
		task->SetEoverPCut(0.78,1.2);
	}
	if(configIndex==111){
		params[0] =  -0.5;
		task->SetEoverPCut(0.80,1.2);
	}
	if(configIndex==112){
		params[0] =  -0.5;
		task->SetEoverPCut(0.82,1.2);
	}
	if(configIndex==113){
		params[0] =  -0.5;
		task->SetEoverPCut(0.83,1.2);
	}
	if(configIndex==114){
		params[0] =  -0.5;
		task->SetEoverPCut(0.84,1.2);
	}
	if(configIndex==115){
		params[0] =  -0.5;
		task->SetEoverPCut(0.85,1.2);
	}
	if(configIndex==116){
		params[0] =  -0.5;
		task->SetEoverPCut(0.86,1.2);
	}
	if(configIndex==117){
		params[0] =  -0.5;
		task->SetEoverPCut(0.85,1.3);
	}
	if(configIndex==118){
		params[0] =  -0.5;
		task->SetEoverPCut(0.86,1.3);
	}
	

	
	
	
	
	
	pid->ConfigureTPCdefaultCut(cutmodel,params,max); 
//_______________________________________________________
	
///_______________________________________________________________________________________________________________
/// New configurations for random cuts -- March, 05, 2014 -- Values in the macro "Random_configurations.C"
	/*
	if (configIndex==80){
		hfecuts->SetMinNClustersTPC(86);
		hfecuts->SetMinNClustersTPCPID(76);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.087);
		task->SetNonHFEangleCut(0.069);
		task->SetAdditionalCuts(0.152, 91);
		task->SetdPhidEtaCut(0.019, 0.044);
		task->SetEoverPCut(0.798, 1.225);
		Double_t params[0]=-0.61;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==81){
		hfecuts->SetMinNClustersTPC(95);
		hfecuts->SetMinNClustersTPCPID(99);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.057);
		task->SetNonHFEangleCut(0.135);
		task->SetAdditionalCuts(0.351, 61);
		task->SetdPhidEtaCut(0.012, 0.044);
		task->SetEoverPCut(0.792, 1.235);
		Double_t params[0]=-1.17;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==82){
		hfecuts->SetMinNClustersTPC(117);
		hfecuts->SetMinNClustersTPCPID(69);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.054);
		task->SetNonHFEangleCut(0.062);
		task->SetAdditionalCuts(0.842, 91);
		task->SetdPhidEtaCut(0.018, 0.033);
		task->SetEoverPCut(0.818, 1.212);
		Double_t params[0]=-1.15;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==83){
		hfecuts->SetMinNClustersTPC(98);
		hfecuts->SetMinNClustersTPCPID(93);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.083);
		task->SetNonHFEangleCut(0.051);
		task->SetAdditionalCuts(0.415, 83);
		task->SetdPhidEtaCut(0.047, 0.016);
		task->SetEoverPCut(0.826, 1.225);
		Double_t params[0]=-1.06;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==84){
		hfecuts->SetMinNClustersTPC(99);
		hfecuts->SetMinNClustersTPCPID(99);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.058);
		task->SetNonHFEangleCut(0.145);
		task->SetAdditionalCuts(0.654, 99);
		task->SetdPhidEtaCut(0.025, 0.014);
		task->SetEoverPCut(0.757, 1.228);
		Double_t params[0]=-1.29;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==85){
		hfecuts->SetMinNClustersTPC(85);
		hfecuts->SetMinNClustersTPCPID(91);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.167);
		task->SetNonHFEangleCut(0.144);
		task->SetAdditionalCuts(0.897, 78);
		task->SetdPhidEtaCut(0.046, 0.043);
		task->SetEoverPCut(0.771, 1.238);
		Double_t params[0]=-1.16;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==86){
		hfecuts->SetMinNClustersTPC(104);
		hfecuts->SetMinNClustersTPCPID(75);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.078);
		task->SetNonHFEangleCut(0.112);
		task->SetAdditionalCuts(0.036, 93);
		task->SetdPhidEtaCut(0.019, 0.013);
		task->SetEoverPCut(0.824, 1.211);
		Double_t params[0]=-1.24;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==87){
		hfecuts->SetMinNClustersTPC(108);
		hfecuts->SetMinNClustersTPCPID(93);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.128);
		task->SetNonHFEangleCut(0.140);
		task->SetAdditionalCuts(0.814, 89);
		task->SetdPhidEtaCut(0.041, 0.022);
		task->SetEoverPCut(0.762, 1.205);
		Double_t params[0]=-0.93;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==88){
		hfecuts->SetMinNClustersTPC(80);
		hfecuts->SetMinNClustersTPCPID(82);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.064);
		task->SetNonHFEangleCut(0.102);
		task->SetAdditionalCuts(0.092, 97);
		task->SetdPhidEtaCut(0.014, 0.031);
		task->SetEoverPCut(0.784, 1.216);
		Double_t params[0]=-1.10;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==89){
		hfecuts->SetMinNClustersTPC(100);
		hfecuts->SetMinNClustersTPCPID(66);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.100);
		task->SetNonHFEangleCut(0.082);
		task->SetAdditionalCuts(0.339, 76);
		task->SetdPhidEtaCut(0.040, 0.011);
		task->SetEoverPCut(0.811, 1.223);
		Double_t params[0]=-0.61;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==90){
		hfecuts->SetMinNClustersTPC(106);
		hfecuts->SetMinNClustersTPCPID(90);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.175);
		task->SetNonHFEangleCut(0.098);
		task->SetAdditionalCuts(0.630, 91);
		task->SetdPhidEtaCut(0.034, 0.026);
		task->SetEoverPCut(0.771, 1.249);
		Double_t params[0]=-0.64;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==91){
		hfecuts->SetMinNClustersTPC(84);
		hfecuts->SetMinNClustersTPCPID(80);
		hfecuts->SetMinNClustersITS(2);
		task->SetNonHFEmassCut(0.083);
		task->SetNonHFEangleCut(0.098);
		task->SetAdditionalCuts(0.594, 66);
		task->SetdPhidEtaCut(0.013, 0.043);
		task->SetEoverPCut(0.809, 1.246);
		Double_t params[0]=-0.57;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==92){
		hfecuts->SetMinNClustersTPC(112);
		hfecuts->SetMinNClustersTPCPID(80);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.090);
		task->SetNonHFEangleCut(0.069);
		task->SetAdditionalCuts(0.803, 91);
		task->SetdPhidEtaCut(0.043, 0.031);
		task->SetEoverPCut(0.805, 1.227);
		Double_t params[0]=-1.29;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==93){
		hfecuts->SetMinNClustersTPC(115);
		hfecuts->SetMinNClustersTPCPID(82);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.075);
		task->SetNonHFEangleCut(0.113);
		task->SetAdditionalCuts(0.637, 86);
		task->SetdPhidEtaCut(0.043, 0.032);
		task->SetEoverPCut(0.760, 1.204);
		Double_t params[0]=-0.84;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	if (configIndex==94){
		hfecuts->SetMinNClustersTPC(81);
		hfecuts->SetMinNClustersTPCPID(63);
		hfecuts->SetMinNClustersITS(3);
		task->SetNonHFEmassCut(0.105);
		task->SetNonHFEangleCut(0.050);
		task->SetAdditionalCuts(0.529, 87);
		task->SetdPhidEtaCut(0.026, 0.047);
		task->SetEoverPCut(0.789, 1.200);
		Double_t params[0]=-1.15;
		pid->ConfigureTPCdefaultCut(cutmodel,params,3.0);
	}
	 */
	
	
///_______________________________________________________________________________________________________________

	printf("*************************************\n");
	printf("Configuring standard Task:\n");
	pid->PrintStatus();
	printf("*************************************\n");

	return task;
}
