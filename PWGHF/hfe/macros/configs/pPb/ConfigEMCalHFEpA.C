///*******************************************************
///Config Description
//configIndex = 0 ---> Default cuts and PID
//configIndex = 1 ---> Non-HFE, Op Angle < 0.1
//configIndex = 2 ---> 
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
	if(configIndex==2) 	hfecuts->SetMinNClustersTPC(90);			                //Minimum number of clusters on TPC
	else if(configIndex==3) 	hfecuts->SetMinNClustersTPC(110);			        //Minimum number of clusters on TPC
	else hfecuts->SetMinNClustersTPC(100);							                //Minimum number of clusters on TPC
	
	if(configIndex==22) hfecuts->SetMinNClustersTPCPID(70); 
	else if (configIndex==23) hfecuts->SetMinNClustersTPCPID(80);					//Minimum number of clusters for dE/dx
	else hfecuts->SetMinNClustersTPCPID(90);					    //Minimum number of clusters for dE/dx
	
	hfecuts->SetMinRatioTPCclusters(0.6);						                    //Number of clusters (Found/Findable)
	
	//ITS
	if(configIndex==25) hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);			//Require at least one cluster on SPD
	else hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);							//Require at least one cluster on SPD
	//hfecuts->SetCutITSdrift(AliHFEextraCuts::kAny); 			                    //Require at least one cluster on SDD
	hfecuts->SetCheckITSLayerStatus(kFALSE); 
	
	if(configIndex==4) hfecuts->SetMinNClustersITS(2);								//Minimum number of clusters on ITS
	else if(configIndex==5) hfecuts->SetMinNClustersITS(4);								//Minimum number of clusters on ITS
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
	
	
	if(configIndex==26){
		task->SetUseShowerShapeCut(kTRUE);
		//task->SetM02Cut(0.0,0.3);
		task->SetM20Cut(0.0,0.3);
	}
	task->SetBackground(kTRUE);
	
	if(configIndex==6) task->SetNonHFEmassCut(0.05);
	else task->SetNonHFEmassCut(0.1);
	
	//if(isEMCal) task->SetEtaCut(-0.6,0.6);
	//else task->SetEtaCut(-0.9,0.9);
	
	
	if(configIndex==12) task->SetEtaCut(-0.6,-0.2);
	else if (configIndex==13) task->SetEtaCut(-0.5,-0.1);
	else if (configIndex==14) task->SetEtaCut(-0.4,0);
	else if (configIndex==15) task->SetEtaCut(-0.3,0.1);
	else if (configIndex==16) task->SetEtaCut(-0.2,0.2);
	else if (configIndex==17) task->SetEtaCut(-0.1,0.3);
	else if (configIndex==18) task->SetEtaCut(0,0.4);
	else task->SetEtaCut(-0.6,0.6);
	
	if(configIndex==19) task->SetdPhidEtaCut(0.02,0.02);
	else if (configIndex==20) task->SetdPhidEtaCut(0.03,0.03);
	else if (configIndex==21) task->SetdPhidEtaCut(0.04,0.04);
	else task->SetdPhidEtaCut(0.05,0.05);

	
	if (configIndex==7) task->SetEoverPCut(0.85,1.2);
	else if (configIndex==8) task->SetEoverPCut(0.75,1.25);
	else task->SetEoverPCut(0.8,1.2);

	if(configIndex==1) task->SetNonHFEangleCut(0.1);
	
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
	
	if(configIndex==9) params[0] = -1.5;
	else if (configIndex==10) params[0] = -0.5;
	else if (configIndex==11) params[0] = 0;
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
