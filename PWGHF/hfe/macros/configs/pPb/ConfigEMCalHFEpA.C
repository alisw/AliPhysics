///*******************************************************
///Config Description

///*******************************************************

AliAnalysisTaskEMCalHFEpA* ConfigEMCalHFEpA(
											
										

Bool_t isMC=kFALSE, 
Int_t triggerIndex=0, 
Int_t configIndex=0, 
Int_t centralityIndex=0, 
Bool_t isAOD = kFALSE,
Bool_t isEMCal = kFALSE,
Int_t EMCalThreshould = 0 //0 == EG1, 1 == EG2
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
	else if(configIndex==6) 	hfecuts->SetMinNClustersTPC(120);					//Minimum number of clusters on TPC
	else hfecuts->SetMinNClustersTPC(100);							                //Minimum number of clusters on TPC
	
	if(configIndex==7) hfecuts->SetMinNClustersTPCPID(70); 
	else if (configIndex==8) hfecuts->SetMinNClustersTPCPID(90);
	else if (configIndex==9) hfecuts->SetMinNClustersTPCPID(60)
	else if (configIndex==10) hfecuts->SetMinNClustersTPCPID(65);
	else if (configIndex==11) hfecuts->SetMinNClustersTPCPID(100);
	else if (configIndex==12) hfecuts->SetMinNClustersTPCPID(95);
	else hfecuts->SetMinNClustersTPCPID(80);										//Minimum number of clusters for dE/dx
	
	hfecuts->SetMinRatioTPCclusters(0.6);						                    //Number of clusters (Found/Findable)
	
	//ITS
	if(configIndex==13) hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);				//Require at least one cluster on SPD
	else hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);							//Require at least one cluster on SPD
	//hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 			                    //Require at least one cluster on SDD
	hfecuts->SetCheckITSLayerStatus(kFALSE); 
	
	if(configIndex==14) hfecuts->SetMinNClustersITS(2);								//Minimum number of clusters on ITS
	else if(configIndex==15) hfecuts->SetMinNClustersITS(4);	
	else if(configIndex==16) hfecuts->SetMinNClustersITS(1);
	else if(configIndex==17) hfecuts->SetMinNClustersITS(5);
	else hfecuts->SetMinNClustersITS(3);								            //Minimum number of clusters on ITS
	
	//Additional Cuts
	hfecuts->SetPtRange(2, 1e6);								                    //Transversal momentum range in GeV/c
	//hfecuts->SetMaxImpactParam(1,2); 							                    //DCA to vertex
	
	//Event Selection
	hfecuts->SetVertexRange(10.);													//
	//hfecuts->SetProductionVertex(0,0.3,0,0.3);									//
///_______________________________________________________________________________________________________________
	// new cuts for event selection
	
	//hfecuts->SetUseCorrelationVertex();
	//hfecuts->SetSPDVtxResolutionCut();
	//hfecuts->SetpApileupCut();

///_________________________________________________________________________________________________________________________
///Task config
	AliAnalysisTaskEMCalHFEpA *task = new AliAnalysisTaskEMCalHFEpA(Form("HFECuts%d_%d_%d",triggerIndex,configIndex,centralityIndex));
	printf("task ------------------------ %p\n ", task);
	task->SetHFECuts(hfecuts);
	task->SetCorrelationAnalysis(kFALSE);
	task->SetAODanalysis(isAOD);
	task->SetEventMixing(kTRUE);
	
	//to separate trigger threshould
	if(EMCalThreshould==0 && triggerIndex==2) task->SetEMCalTriggerEG1();
	if(EMCalThreshould==1 && triggerIndex==2) task->SetEMCalTriggerEG2();
	
	if(isEMCal) task->SetUseEMCal();
	
	
	if(configIndex==100){
		task->SetUseShowerShapeCut(kTRUE);
		//task->SetM02Cut(0.0,0.3);
		task->SetM20Cut(0.0,0.3);
	}
	task->SetBackground(kTRUE);
	
	//nonHFE cuts
	if(configIndex==20) task->SetNonHFEmassCut(0.05);
	else if(configIndex==21) task->SetNonHFEmassCut(0.15);
	else if(configIndex==22) task->SetNonHFEmassCut(0.03);
	else if(configIndex==23) task->SetNonHFEmassCut(0.18);
	else if(configIndex==24) task->SetNonHFEmassCut(0.01);
	else if(configIndex==25) task->SetNonHFEmassCut(0.2);
	else task->SetNonHFEmassCut(0.1);
	
	if(configIndex==26) task->SetNonHFEangleCut(0.1);
	if(configIndex==27) task->SetNonHFEangleCut(0.15);
	if(configIndex==28) task->SetNonHFEangleCut(0.05);
	
	//partner cuts
	
	if(configIndex==29) task->SetAdditionalCuts(0,80);
	if(configIndex==30) task->SetAdditionalCuts(0.5,80);
	if(configIndex==31) task->SetAdditionalCuts(0.7,80);
	if(configIndex==32) task->SetAdditionalCuts(0.9,80);
	
	if(configIndex==33) task->SetAdditionalCuts(0.3,60);
	if(configIndex==34) task->SetAdditionalCuts(0.3,70);
	if(configIndex==35) task->SetAdditionalCuts(0.3,90);
	if(configIndex==36) task->SetAdditionalCuts(0.3,100);
	 
	 task->SetAdditionalCuts(0.3,80);
	 
	
	//eta cuts
	if(configIndex==40) task->SetEtaCut(-0.6,-0.2);
	else if (configIndex==41) task->SetEtaCut(-0.5,-0.1);
	else if (configIndex==42) task->SetEtaCut(-0.4,0);
	else if (configIndex==43) task->SetEtaCut(-0.3,0.1);
	else if (configIndex==44) task->SetEtaCut(-0.2,0.2);
	else if (configIndex==45) task->SetEtaCut(-0.1,0.3);
	else if (configIndex==46) task->SetEtaCut(0,0.4);
	else if (configIndex==47) task->SetEtaCut(-0.4,0.4);
	else if (configIndex==48) task->SetEtaCut(-0.3,0.3);
	else task->SetEtaCut(-0.6,0.6);
	
	//track matching cuts
	if(configIndex==50) task->SetdPhidEtaCut(0.02,0.02);
	else if (configIndex==51) task->SetdPhidEtaCut(0.03,0.03);
	else if (configIndex==52) task->SetdPhidEtaCut(0.04,0.04);
	else task->SetdPhidEtaCut(0.05,0.05);

	//E/p Cuts
	if (configIndex==60) task->SetEoverPCut(0.85,1.2);
	else if (configIndex==61) task->SetEoverPCut(0.75,1.25);
	else if (configIndex==62) task->SetEoverPCut(0.70,1.2);
	else if (configIndex==63) task->SetEoverPCut(0.80,1.25);
	else if (configIndex==64) task->SetEoverPCut(0.9,1.3);
	else if (configIndex==65) task->SetEoverPCut(0.95,1.3);
	else if (configIndex==66) task->SetEoverPCut(0.75,1.2);
	else task->SetEoverPCut(0.8,1.2);

	
	
	if(centralityIndex==0) task->SetCentrality(0,20);
	if(centralityIndex==1) task->SetCentrality(20,40);
	if(centralityIndex==2) task->SetCentrality(40,60);
	if(centralityIndex==3) task->SetCentrality(60,80);
	if(centralityIndex==4) task->SetCentrality(80,100);
	if(centralityIndex==5) task->SetCentrality(0,100);
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
	
	if(configIndex==70) params[0] = -1.5;
	else if (configIndex==71) params[0] = -0.5;
	else if (configIndex==72) params[0] = 0;
	else if (configIndex==73) params[0] = 0.25;
	else if (configIndex==74) params[0] = -1.75;
	else params[0] = -1;
	
	pid->ConfigureTPCdefaultCut(cutmodel,params,3.0); 
//_______________________________________________________
///_______________________________________________________________________________________________________________

	printf("*************************************\n");
	printf("Configuring standard Task:\n");
	pid->PrintStatus();
	printf("*************************************\n");

	return task;
}
