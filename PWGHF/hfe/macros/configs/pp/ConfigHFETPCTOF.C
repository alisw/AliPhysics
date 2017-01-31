///*******************************************************
///Config Description
/// 05th March 2015 - Cristiane Jahnke
///*******************************************************

AliAnalysisHFETPCTOF* ConfigHFETPCTOF(
											
										

Bool_t isMC=kFALSE, 
Bool_t isAOD = kTRUE,
Bool_t isPP = kTRUE

)

{
///_______________________________________________________________________________________________________________
///Track selection: Cuts used to ensure a minimum quality level of the tracks selected to perform the analysis
	AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsMinBias","HFE Cuts");
	hfecuts->CreateStandardCuts();
	
	//TPC Cuts
	hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);	
    hfecuts->SetMinNClustersTPC(100);							                //Minimum number of clusters on TPC
	
    hfecuts->SetMinNClustersTPCPID(80);										//Minimum number of clusters for dE/dx
	
	hfecuts->SetMinRatioTPCclusters(0.6);						                    //Number of clusters (Found/Findable)
	
	//ITS
    //hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);							//Require at least one cluster on SPD
	hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);  //Require cluster in at least two layers of SPD: to reduce backg of conversion electrons.
	hfecuts->SetCheckITSLayerStatus(kFALSE); 
	
	hfecuts->SetMinNClustersITS(3);								            //Minimum number of clusters on ITS
	
	//Additional Cuts
	//hfecuts->SetPtRange(2, 1e6);								                    //Transversal momentum range in GeV/c
	hfecuts->SetPtRange(0.1, 1e6);
	
	//testing this line for the DCA cut
	//hfecuts->SetRequireDCAToVertex();
	hfecuts->SetMaxImpactParam(1,2); 							                //DCA to vertex
	
	//Event Selection
	hfecuts->SetVertexRange(10.);													//
	//hfecuts->SetProductionVertex(0,0.3,0,0.3);									//
///_______________________________________________________________________________________________________________
	// new cuts for event selection
	
	//hfecuts->SetUseCorrelationVertex();
	//hfecuts->SetSPDVtxResolutionCut();
	//hfecuts->SetpApileupCut();

//___________________________________________________________________________________________________________
///Task config
	AliAnalysisHFETPCTOF *task = new AliAnalysisHFETPCTOF();
	printf("task ------------------------ %p\n ", task);
	task->SetHFECuts(hfecuts);
	task->SetAODanalysis(isAOD);
	task->SetPPanalysis(isPP);

		
//______________________________________
///Created by the user
    
	
//______________________________________
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
	
    params[0] = -1;
	
    Double_t max=3.0;
	
	//params[0] = 3;
	
    //Double_t max=7.0;
	
	pid->ConfigureTPCdefaultCut(cutmodel,params,max);

	
///_______________________________________________________________________________________________________________

	printf("*************************************\n");
	printf("Configuring standard Task:\n");
	pid->PrintStatus();
	printf("*************************************\n");

	return task;
}
