AliAnalysisHFEppTPCTOFBeauty* AddTaskHFEppTPCTOFBeauty( ///-> to run locally
			TString uniqueID        = "",
			Bool_t 	isMC 			= kFALSE, 
			Bool_t 	isAOD 			= kFALSE,
			Bool_t 	isPP 			= kFALSE,
			Double_t tpcPIDmincut,
			Double_t tpcPIDmaxcut,
			Double_t tofPIDcut,
			Int_t MinNClustersTPC =	100,
		    Int_t MinNClustersTPCPID =	80,
		    Float_t MinRatioTPCclusters =	0.6,
		    Int_t  MinNClustersITS = 3,
            AliHFEextraCuts::ITSPixel_t pixel,
            Float_t EtaMin,
            Float_t EtaMax,
            AliVEvent::EOfflineTriggerTypes trigger,
            Float_t DCAxy,
            Float_t DCAz,
            Int_t IsBcorr = 0,
            Int_t IsDcorr = 0
            
)           
{
    
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
	::Error("AddTaskHFEppTPCTOFBeauty", "No analysis manager to connect to.");
	return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
	::Error("AddTaskHFEppTPCTOFBeauty", "This task requires an input event handler");
	return NULL;
	}
	
	//_______________________
    AliAnalysisHFEppTPCTOFBeauty *task = ConfigHFEppTPCTOFBeauty(isMC,isAOD,isPP,tpcPIDmincut,tpcPIDmaxcut,tofPIDcut,MinNClustersTPC,MinNClustersTPCPID,MinRatioTPCclusters,MinNClustersITS,pixel,EtaMin,EtaMax,DCAxy,DCAz,IsBcorr,IsDcorr);
    //_____________________________________________________
	//Trigger
		//if(!isMC){
			task->SelectCollisionCandidates(trigger); //Selecting Minumum Bias events (selected randomlly)
	//	}
	//_____________________________________________________
	
	
	
	mgr->AddTask(task);
	//added to run on train
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":ElectroID_";
    fileName += uniqueID;
	
	//Create containers for input/output -> to run locally
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("ccontainer0_%s",uniqueID.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
	//Connect input/output
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
    
}


AliAnalysisHFEppTPCTOFBeauty* ConfigHFEppTPCTOFBeauty(Bool_t isMCc, Bool_t isAODc, Bool_t isPPc,Double_t tpcPIDmincut,Double_t tpcPIDmaxcut, Double_t tofPID, Int_t minNClustersTPC, Int_t minNClustersTPCPID, Float_t minRatioTPCclusters, Int_t  minNClustersITS, AliHFEextraCuts::ITSPixel_t pixel, Float_t EtaMin, Float_t EtaMax, Float_t DCAxy, Float_t DCAz, Int_t IsBcorr, Int_t IsDcorr)
{
    ///_______________________________________________________________________________________________________________
    ///Track selection: Cuts used to ensure a minimum quality level of the tracks selected to perform the analysis
    AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsMinBias","HFE Cuts");
    hfecuts->CreateStandardCuts();
    //cout<<minNClustersTPC<<endl;
    //cout<<minNClustersITS<<endl;
    
    //TPC Cuts
    hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetMinNClustersTPC(minNClustersTPC);							                //Minimum number of clusters on TPC
    hfecuts->SetMinNClustersTPCPID(minNClustersTPCPID);										//Minimum number of clusters for dE/dx
    hfecuts->SetMinRatioTPCclusters(minRatioTPCclusters);						                    //Number of clusters (Found/Findable)
    //ITS
    hfecuts->SetCutITSpixel(pixel);  //Require cluster in at least two layers of SPD: to reduce backg of conversion electrons.
    hfecuts->SetCheckITSLayerStatus(kFALSE);
    hfecuts->SetMinNClustersITS(minNClustersITS);								            //Minimum number of clusters on ITS
    //Additional Cuts
    hfecuts->SetPtRange(0.3, 20); //Transversal momentum range in GeV/c
    //testing this line for the DCA cut
    //hfecuts->SetRequireDCAToVertex();
    
    //cout<<DCAxy<<endl;
    //cout<<DCAz<<endl;
    
    hfecuts->SetMaxImpactParam(DCAxy,DCAz); //DCA to vertex
   
   
   cout<<isMCc<<"   "<<isAODc<<"    "<<isPPc<<"    "<<tpcPIDmincut<<"   "<<tpcPIDmaxcut<<"   "<<tofPID<<"   "<<minNClustersTPC<<"   "<<minNClustersTPCPID<<"   "<<minRatioTPCclusters<<"   "<<minNClustersITS<<"    "<<EtaMin<<"    "<<EtaMax<<"    "<<DCAxy<<"   "<<DCAz<<"   "<<IsBcorr<<"    "<<IsDcorr<<endl;
   
    ///Task config
    AliAnalysisHFEppTPCTOFBeauty *task = new AliAnalysisHFEppTPCTOFBeauty();
    printf("task ------------------------ %p\n ", task);
    task->SetHFECuts(hfecuts);
    task->SetAODanalysis(isAODc);
    task->SetPPanalysis(isPPc);
      
    //Setter for the PID cuts--------------
	task->SetPIDCuts(tpcPIDmincut, tpcPIDmaxcut, -tofPID, tofPID);
    //-----------------------------------------
    
    //Setter for the Eta cut--------------
	task->SetEtaCut(EtaMin, EtaMax);
    //-----------------------------------------
    
    //______________________________________
    //Particle identification
    AliHFEpid *pid = task->GetPID();
    
    //______________________________________
    //In the case of a simulation
    if(isMCc)
    {
        pid->SetHasMCData(kTRUE);
        task->SetMCanalysis();
    }
    //______________________________________
    
    if(isMCc){
    ///B correction
    ///Default FONLL
    if(IsBcorr == 0){
		TF1 *fBmesonShape = new TF1("fBmesonShape","(1.0/0.455737)*(5.78258e-03 / (TMath::Power(TMath::Exp( - 2.07792e+00*x[0] + 6.75483e-03 * x[0] * x[0] ) + x[0] / 1.66686e+01, 7.05243e-01))*(2.27740e-06 / (TMath::Power(TMath::Exp( - 1.72205e-01*x[0] + 3.69234e-03 * x[0] * x[0] ) + x[0] / 7.12479e+00, 4.34903e+00)))*TMath::Power(1.0+x[0]/1.79842e-05,6.74073e-01))/(1.87356e-02 / (TMath::Power(TMath::Exp( - 1.19647e+00*x[0] - 3.54687e-02 * x[0] * x[0] ) + x[0] / 4.67476e+00, 1.47150e+00))*(4.74056e-03 / (TMath::Power(TMath::Exp( - 5.31680e-01*x[0] + 1.02081e-02 * x[0] * x[0] ) + x[0] / 4.84483e+01, 9.27183e-01)))*TMath::Power(1.0+x[0]/1.74979e+03,-1.57702e+02))",0.0,50);
		task->SetBcorrFunction(fBmesonShape);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///FONLL Lower edge
    if(IsBcorr == 1){
		TF1 *fBmesonShape1 = new TF1("fBmesonShape1","(1.0/0.389241)*(3.77085e-03 / (TMath::Power(TMath::Exp( - 4.99937e+00*x[0] - 8.02258e-01 * x[0] * x[0] ) + x[0] / 2.60310e+04, 1.89019e-01))*(1.43903e-03 / (TMath::Power(TMath::Exp( - 1.82037e-01*x[0] + 3.88349e-03 * x[0] * x[0] ) + x[0] / 7.21515e+00, 4.49037e+00)))*TMath::Power(1.0+x[0]/2.91057e-01,4.05964e-01))/(1.87356e-02 / (TMath::Power(TMath::Exp( - 1.19647e+00*x[0] - 3.54687e-02 * x[0] * x[0] ) + x[0] / 4.67476e+00, 1.47150e+00))*(4.74056e-03 / (TMath::Power(TMath::Exp( - 5.31680e-01*x[0] + 1.02081e-02 * x[0] * x[0] ) + x[0] / 4.84483e+01, 9.27183e-01)))*TMath::Power(1.0+x[0]/1.74979e+03,-1.57702e+02))",0.0,50);
		task->SetBcorrFunction(fBmesonShape1);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///FONLL Upper edge
    if(IsBcorr == 2){
		TF1 *fBmesonShape2 = new TF1("fBmesonShape2","(1.0/0.532972)*(8.25352e-03 / (TMath::Power(TMath::Exp( - 6.11267e+00*x[0] - 1.05995e-01 * x[0] * x[0] ) + x[0] / 3.12399e+04, 2.13332e-01))*(1.70947e-03 / (TMath::Power(TMath::Exp( - 1.75188e-01*x[0] + 3.76833e-03 * x[0] * x[0] ) + x[0] / 6.82438e+00, 4.61154e+00)))*TMath::Power(1.0+x[0]/3.25350e-01,3.45204e-01))/(1.87356e-02 / (TMath::Power(TMath::Exp( - 1.19647e+00*x[0] - 3.54687e-02 * x[0] * x[0] ) + x[0] / 4.67476e+00, 1.47150e+00))*(4.74056e-03 / (TMath::Power(TMath::Exp( - 5.31680e-01*x[0] + 1.02081e-02 * x[0] * x[0] ) + x[0] / 4.84483e+01, 9.27183e-01)))*TMath::Power(1.0+x[0]/1.74979e+03,-1.57702e+02))",0.0,50);
		task->SetBcorrFunction(fBmesonShape2);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
 
    // Default
    if(IsDcorr == 0){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/0.296169)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.0,36);  
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/0.287189)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/0.27827)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/0.269485)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/0.260894)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/0.252547)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/0.232981)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/0.215513)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/0.200212)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/0.18697)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/0.175584)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/0.165817)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/0.150219)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/0.138568)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/0.129681)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/0.122738)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.117183)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.112639)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.108845)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.105620)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.100396)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.0962985)*(1.31661e-03/(TMath::Power(TMath::Exp( - 6.20935e-01*x[0] + 9.85658e-03 * x[0] * x[0] ) + x[0] / 3.59774e+00, 4.42801e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",9.0,36);  
		task->SetDcorrFunction15(fDmesonShape15);
	}
    
    // Data up tilted 
    
     if(IsDcorr == 1){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/0.318575)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.0,36);  
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/0.306756)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/0.295518)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/0.284868)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/0.274805)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/0.265321)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/0.244033)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/0.225915)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/0.210545)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/0.197494)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/0.186364)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/0.176812)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/0.16134)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/0.14936)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/0.139769)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/0.131863)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.125188)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.119445)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.114426)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.109988)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.102454)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.0962635)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",9.0,36);  
		task->SetDcorrFunction15(fDmesonShape15);
	}
   
    
    // Data Down tilted
    
     if(IsDcorr == 2){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/0.275971)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.0,36);  
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/0.267682)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/0.259608)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/0.251762)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/0.24416)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/0.236812)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/0.219616)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/0.204151)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/0.190411)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/0.178315)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/0.167736)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/0.15852)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/0.143546)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/0.132225)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/0.123608)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/0.116984)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.111835)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.107785)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.104563)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.101968)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.098111)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.0954295)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57474e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",9.0,36);  
		task->SetDcorrFunction15(fDmesonShape15);
	}
   
    
 	}
	
    
    ///Configure PID
    //_________________________
    //TPC PID
    pid->AddDetector("TPC", 0);				
    
    //_________________________
    ///Configure TPC cut
    //Defaul = -1 to 3 sigmas
    //Note that it is also possible to define a model instead of a constant
    //--------->For this change the "cut model"
    
    Double_t params[4];
    char *cutmodel;
    cutmodel = "pol0";
    params[0] = tpcPIDmincut;
    Double_t max= tpcPIDmaxcut;
    pid->ConfigureTPCdefaultCut(cutmodel,params,max);

    ///TOF PID
    pid->AddDetector("TOF", 1);
    pid->ConfigureTOF(tofPID);
    
    printf("*************************************\n");
    printf("Configuring standard Task:\n");
    pid->PrintStatus();
    printf("*************************************\n");
    
    return task;
}







