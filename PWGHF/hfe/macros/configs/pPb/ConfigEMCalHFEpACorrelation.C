///*******************************************************
///Config Description
//configIndex = 0 ---> Default cuts and PID
//configIndex = 1 ---> TPC Ncls = 100
//configIndex = 2 ---> TPC Ncls = 60
//configIndex = 3 ---> SPD kBoth + 3 ITS cls
//configIndex = 4 ---> SPD kBoth + 4 ITS cls
//configIndex = 5 ---> SPD kAny + 3 ITS cls
//configIndex = 6 ---> Mass < 0.05
//configIndex = 7 ---> Mass < 0.15
//configIndex = 8 ---> Op Angle < 0.1
//configIndex = 9 ---> TPC PID: -0.5 to 3.0
//configIndex = 10 ---> V0A -> other
//configIndex = 11 ---> Associated hadron with SPD::kAny cut
///*******************************************************

AliAnalysisTaskEMCalHFEpA* ConfigEMCalHFEpACorrelation(
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
	if(configIndex==1) 	hfecuts->SetMinNClustersTPC(100);			                //Minimum number of clusters on TPC = 100
	else if(configIndex==2) hfecuts->SetMinNClustersTPC(60);			        	//Minimum number of clusters on TPC = 60
	else hfecuts->SetMinNClustersTPC(80);							                //Minimum number of clusters on TPC = 80
	
	hfecuts->SetMinNClustersTPCPID(80); 						                    //Minimum number of clusters for dE/dx
	hfecuts->SetMinRatioTPCclusters(0.6);						                    //Number of clusters (Found/Findable)
	
	//ITS
	if(configIndex==3) 
	{
		hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);							//Require 2 cluster on SPD
		hfecuts->SetMinNClustersITS(3);												//Minimum number of clusters on ITS
	}
	else if(configIndex==4) 
	{
		hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);							//Require 2 cluster on SPD
		hfecuts->SetMinNClustersITS(4);												//Minimum number of clusters on ITS
	}
	else if(configIndex==5) 
	{
		hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);				            	//Require at least one cluster on SPD
		hfecuts->SetMinNClustersITS(3);												//Minimum number of clusters on ITS
	}
	else
	{
		hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);				            	//Require at least one cluster on SPD
		hfecuts->SetMinNClustersITS(2);												//Minimum number of clusters on ITS
	}
	
	hfecuts->SetCheckITSLayerStatus(kFALSE); 
	
	//Additional Cuts
	hfecuts->SetPtRange(0.5, 1e6);								                    //Transversal momentum range in GeV/c
	//hfecuts->SetMaxImpactParam(1,2); 							                    //DCA to vertex
	
	//Event Selection
	hfecuts->SetVertexRange(10.);													//
	//hfecuts->SetProductionVertex(0,0.3,0,0.3);									//
///_______________________________________________________________________________________________________________

///_________________________________________________________________________________________________________________________
///Task config
	AliAnalysisTaskEMCalHFEpA *task = new AliAnalysisTaskEMCalHFEpA(Form("HFECuts%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould));
	printf("task ------------------------ %p\n ", task);
	
	task->SetHFECuts(hfecuts);
	task->SetCorrelationAnalysis();
	task->SetAODanalysis(isAOD);
	task->SetEventMixing(kTRUE);
	
	task->SetAssHadronPtRange(0.5,2.0);
	
	task->SetAdditionalCuts(0.0,80);
	if(configIndex==20) task->SetAdditionalCuts(0.0,80);
	if(configIndex==21) task->SetAdditionalCuts(0.3,80);
	if(configIndex==22) task->SetAdditionalCuts(0.5,80);
	if(configIndex==23) task->SetAdditionalCuts(0.7,80);
	
	if(configIndex==11) task->SetSPDCutForHadrons();
	
	if(configIndex==10) task->SetCentralityEstimator(1);
	else task->SetCentralityEstimator(0);
	
	if(EMCalThreshould==0 && triggerIndex==2) task->SetEMCalTriggerEG1();
	if(EMCalThreshould==1 && triggerIndex==2) task->SetEMCalTriggerEG2();
	
	if(isEMCal) task->SetUseEMCal();
	
	if(configIndex==6) task->SetNonHFEmassCut(0.05);
	else if(configIndex==7) task->SetNonHFEmassCut(0.15);
	else task->SetNonHFEmassCut(0.1);
	
	if(isEMCal) task->SetEtaCut(-0.6,0.6);
	else task->SetEtaCut(-0.9,0.9);
	
	task->SetEoverPCut(0.8,1.2);	//Will work only in case isEMCal = kTRUE

	if(configIndex==8) task->SetNonHFEangleCut(0.1);
	
	if(centralityIndex==0) task->SetCentrality(0,20);
	if(centralityIndex==1) task->SetCentrality(20,60);
	if(centralityIndex==2) task->SetCentrality(60,100);
	if(centralityIndex==3) task->SetCentrality(0,10);
	if(centralityIndex==4) task->SetCentrality(10,20);
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
	//TPC+TOF PID
	pid->AddDetector("TOF", 0);				//Add TOF PID
	pid->AddDetector("TPC", 1);				//Add TPC PID
	
	//_________________________
	//Configure TPC cut
	//Defaul = -1 to 3 sigmas
	//Note that it is also possible to define a model instead of a constant
	//--------->For this change the "cut model"
	
	Double_t params[4];
	char *cutmodel;
	cutmodel = "pol0";
	
	if(configIndex==9) params[0] = 0.0;
	else params[0] = -0.5;
	
	pid->ConfigureTPCdefaultCut(cutmodel,params,3.0); 
//_______________________________________________________
///_______________________________________________________________________________________________________________

	printf("*************************************\n");
	printf("Configuring standard Task:\n");
	pid->PrintStatus();
	printf("*************************************\n");

	return task;
}
