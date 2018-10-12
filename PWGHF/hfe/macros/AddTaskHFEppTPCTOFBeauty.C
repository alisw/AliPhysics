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
    ///Default
    if(IsBcorr == 0){
		TF1 *fBmesonShape = new TF1("fBmesonShape","(1.0/1.2853)*(5.27632e-03 / (TMath::Power(TMath::Exp( - 1.37634e+00*x[0] + 6.30089e-03 * x[0] * x[0] ) + x[0] / 8.00647e+00, 1.56482e+00))*(1.46817e-03 / (TMath::Power(TMath::Exp( - 2.30684e-01*x[0] + 4.99541e-03 * x[0] * x[0] ) + x[0] / 1.37772e+01, 2.92952e+00)))*TMath::Power(1.0+x[0]/2.20161e+02,2.26672e+00))/(9.25555e-03 / (pow(TMath::Exp( - 1.16122e+00*x[0] - 1.96195e-02 * x[0] * x[0] ) + x[0] / 5.21757e+00, 1.39344e+00))*(1.36994e-03 / (pow(TMath::Exp( - 2.12508e-01*x[0] + 9.07838e-04 * x[0] * x[0] ) + x[0] / 1.55251e+01, 3.04888e+00)))*TMath::Power(1.0+x[0]/1.12492e+02,1.17620e-01))",0.0,50);
		task->SetBcorrFunction(fBmesonShape);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///Var 1
  /*  if(IsBcorr == 1){
		//TF1 *fBmesonShape1 = new TF1("fBmesonShape1","(0.5/(1. + exp((x[0] - 7.) * 0.7)) + 0.5 + (x[0] - 15.)/300)+(1-(0.5/(1. + exp((x[0] - 7.) * 0.7)) + 0.5 + (x[0] - 15.)/300))/2", 0, 30);
		TF1 *fBmesonShape1 = new TF1("fBmesonShape1","(4.65644e-01) / (TMath::Power(TMath::Exp( - (2.63272e-01) * x[0] - (-7.24611e-03) * x[0] * x[0] ) + x[0] / (1.24435e+01), (2.09389e+00)))", 0, 30);
		task->SetBcorrFunction(fBmesonShape1);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///Var 2
    if(IsBcorr == 2){
		//TF1 *fBmesonShape2 = new TF1("fBmesonShape2","(0.5/(1. + exp((x[0] - 7.) * 0.7)) + 0.5 + (x[0] - 15.)/300)-(1-(0.5/(1. + exp((x[0] - 7.) * 0.7)) + 0.5 + (x[0] - 15.)/300))/2", 0, 30);
		TF1 *fBmesonShape2 = new TF1("fBmesonShape2","(4.23575e-01) / (TMath::Power(TMath::Exp( - (2.65700e-01) * x[0] - (4.35955e-02) * x[0] * x[0] ) + x[0] / (7.48510e+00), (2.13502e+00)))", 0, 30);
		task->SetBcorrFunction(fBmesonShape2);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///Weight from RAAxFONLL/MC
    if(IsBcorr == 3){
		TF1 *fBmesonShape3 = new TF1("fBmesonShape3","(1./0.925)*(3.74727e-01/(TMath::Power(TMath::Exp( - (3.65064e-01) * x[0] - (-9.29285e-03) * x[0] * x[0] ) + x[0] / (1.06879e+01), (2.01898e+00))))", 0, 30);
		task->SetBcorrFunction(fBmesonShape3);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    */
    
    if(IsDcorr == 0){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/31.3458)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.0,36);  
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/29.4832)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/27.7431)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/26.1229)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/24.6186)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/23.2251)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/20.1887)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/17.7144)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/15.707)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/14.0802)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/12.7599)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/11.6848)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/10.083)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/8.98936)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/8.22561)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/7.68037)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/7.28301)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/6.98789)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/6.76485)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/6.59356)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/6.3545)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/6.20141)*(3.05777e-03/(TMath::Power(TMath::Exp( - 4.61036e-01*x[0] + 1.19288e-02 * x[0] * x[0] ) + x[0] / 3.16263e+00, 4.60929e+00)))/(5.48342e-05 / (TMath::Power(TMath::Exp( - 4.48121e-01*x[0] + 1.11703e-02 * x[0] * x[0] ) + x[0] / 5.12305e+00, 4.57546e+00)))",9.0,36);  
		task->SetDcorrFunction15(fDmesonShape15);
	}
    
  /*  if(IsDcorr == 1){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/151.879)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.0,30);
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/115.442)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.1,30);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/89.7009)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.2,30);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/71.0357)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.3,30);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/57.1919)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.4,30);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("fDmesonShape17","(1/46.7188)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.5,30);
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/29.7105)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.75,30);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/20.079)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2,30);
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/14.2264)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2.25,30);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShape3","(1/10.4666)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2.5,30);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/7.94014)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2.75,30);   
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/6.17838)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",3,30);   
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/3.97542)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",3.5,30);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/2.7242)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",4,30);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/1.95827)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",4.5,30);   
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/1.46145)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",5,30);   
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/1.12401)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",5.5,30);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.88607)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",6,30);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.713003)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",6.5,30);   
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.583796)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",7,30);   
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.408415)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",8,30);   
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.298932)*(1.25915e+05*TMath::Power(1-(1-1.20176e+00)*x/6.15769e-02,1/(1-1.20176e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",9,30);   
		task->SetDcorrFunction15(fDmesonShape15);
	}
	
	if(IsDcorr == 2){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/86.5976)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.0,30);
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/65.8729)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.1,30);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/51.3564)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.2,30);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/40.8843)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.3,30);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/33.1368)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.4,30);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("fDmesonShape17","(1/27.2782)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.5,30);
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/17.7376)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",1.75,30);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/12.2852)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2,30);
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/8.92829)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2.25,30);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShape3","(1/6.73831)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2.5,30);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/5.24209)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",2.75,30);   
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/4.18062)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",3,30);   
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/2.81982)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",3.5,30);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/2.01957)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",4,30);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/1.51292)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",4.5,30);   
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/1.17351)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",5,30);   
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.935775)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",5.5,30);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.763158)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",6,30);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.634055)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",6.5,30);   
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.535073)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",7,30);   
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.395789)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",8,30);   
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.304682)*(2.99867e+05*TMath::Power(1-(1-1.22514e+00)*x/3.76340e-02,1/(1-1.22514e+00)))/(1.96376e+00*TMath::Power(1-(1-1.38770e+00)*x/6.84942e-01,1/(1-1.38770e+00)))",9,30);   
		task->SetDcorrFunction15(fDmesonShape15);
	}
	
	*/
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







