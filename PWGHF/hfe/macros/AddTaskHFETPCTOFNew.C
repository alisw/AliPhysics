AliAnalysisHFETPCTOFNew* AddTaskHFETPCTOFNew( ///-> to run locally
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
            Float_t Mass, 
            Float_t MinPt, 
            Float_t TpcNclus,
            Float_t EtaMin,
            Float_t EtaMax,
            AliVEvent::EOfflineTriggerTypes trigger,
            Int_t 	isCharPion 			= 0,
            Float_t DCAxy,
            Float_t DCAz,
            Int_t HadCont 			= 0
            
)           
{
    
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
	::Error("AddTaskHFETPCTOFRun2", "No analysis manager to connect to.");
	return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
	::Error("AddTaskHFETPCTOFRun2", "This task requires an input event handler");
	return NULL;
	}
	
	//_______________________
    AliAnalysisHFETPCTOFNew *task = ConfigHFETPCTOF(isMC,isAOD,isPP,tpcPIDmincut,tpcPIDmaxcut,tofPIDcut,MinNClustersTPC,MinNClustersTPCPID,MinRatioTPCclusters,MinNClustersITS,pixel,Mass,MinPt,TpcNclus,EtaMin,EtaMax,isCharPion,DCAxy,DCAz,HadCont);
    //_____________________________________________________
	//Trigger
		if(!isMC){
			task->SelectCollisionCandidates(trigger); //Selecting Minumum Bias events (selected randomlly)
		}
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


AliAnalysisHFETPCTOFNew* ConfigHFETPCTOF(Bool_t isMCc, Bool_t isAODc, Bool_t isPPc,Double_t tpcPIDmincut,Double_t tpcPIDmaxcut, Double_t tofPID, Int_t minNClustersTPC, Int_t minNClustersTPCPID, Float_t minRatioTPCclusters, Int_t  minNClustersITS, AliHFEextraCuts::ITSPixel_t pixel, Float_t Mass, Float_t MinPt, Float_t TpcNclus, Float_t EtaMin, Float_t EtaMax, Int_t isCharPion, Float_t DCAxy, Float_t DCAz, Int_t HadCont)
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
    hfecuts->SetPtRange(0.3, 6); //Transversal momentum range in GeV/c
    //testing this line for the DCA cut
    //hfecuts->SetRequireDCAToVertex();
    
    //cout<<DCAxy<<endl;
    //cout<<DCAz<<endl;
    
    hfecuts->SetMaxImpactParam(DCAxy,DCAz); 							                //DCA to vertex
    								//
    ///_______________________________________________________________________________________________________________
    
    //___________________________________________________________________________________________________________
    ///Task config
    AliAnalysisHFETPCTOFNew *task = new AliAnalysisHFETPCTOFNew();
    printf("task ------------------------ %p\n ", task);
    task->SetHFECuts(hfecuts);
    task->SetAODanalysis(isAODc);
    task->SetPPanalysis(isPPc);
    
    //Setter for the Partner cuts--------------
	task->SetPartnerCuts(Mass,MinPt,TpcNclus);
    //-----------------------------------------
    
    //Setter for the PID cuts--------------
	task->SetPIDCuts(tpcPIDmincut, tpcPIDmaxcut, -tofPID, tofPID);
    //-----------------------------------------
    
    //Setter for the Eta cut--------------
	task->SetEtaCut(EtaMin, EtaMax);
    //-----------------------------------------
    
   
    //______________________________________
    ///Hadron Contamination function
    
         
    ///Landau function
    if(HadCont == 0){
		TF1 *hBackground = new TF1("hadronicBackgroundFunction","[0]*TMath::Landau(x,[1],[2],0)", 0.5, 8);
		hBackground->FixParameter(0, 1.89864e-02);
		hBackground->FixParameter(1, 4.82510e+00);
		hBackground->FixParameter(2, 1.18288e+00);
		//hBackground->Eval(0,0,0);
		task->SetHadronFunction(hBackground);
	} 
     
    ///Error function
    if(HadCont == 1){
		TF1 *hBackground = new TF1("hadronicBackgroundFunction","[0]+[1]*TMath::Erf(([2]*x-[3]))", 0.5, 8);
		hBackground->FixParameter(0, 9.99999e-01);
		hBackground->FixParameter(1, 1.00000e+00);
		hBackground->FixParameter(2, 4.51097e-01);
		hBackground->FixParameter(3, 3.67968e+00);
		//hBackground->Eval(0,0,0);
		task->SetHadronFunction(hBackground);
	}
	
	///Landau function fit up
	if(HadCont == 2){
		TF1 *hBackground = new TF1("hadronicBackgroundFunction","[0]*TMath::Landau(x,[1],[2],0)", 0.5, 8);
		hBackground->FixParameter(0, 2.82262e-02);
		hBackground->FixParameter(1, 5.03114e+00);
		hBackground->FixParameter(2, 1.26631e+00);
		//hBackground->Eval(0,0,0);
		task->SetHadronFunction(hBackground);
	}
	
	///Landau function fit down
	if(HadCont == 3){
		TF1 *hBackground = new TF1("hadronicBackgroundFunction","[0]*TMath::Landau(x,[1],[2],0)", 0.5, 8);
		hBackground->FixParameter(0, 9.85918e-03);
		hBackground->FixParameter(1, 4.37298e+00);
		hBackground->FixParameter(2, 9.99751e-01);
		//hBackground->Eval(0,0,0);
		task->SetHadronFunction(hBackground);
		//cout<<"HadCont = "<<HadCont<<endl;
	}  
	
    //______________________________________
    
    
    ///Weight for pi0s and etas
    
    
    ///Weights tilted up and down calculated using the charged pion systematic uncertainty
    ///................................................................................................................................................
    /*
    if(isCharPion == 0){
		///Good weights for d20a2_extra (data over MB), where data = charged pions and eta from Mt-scaling (using |eta| < 1.2 in MC for the weight calculation)
		
		cout<<"-----------reference--------------------"<<endl;
		///Reference-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hpi0w1->FixParameter(0,7.96404e+00);
		hpi0w1->FixParameter(1,-9.14164e+00);
		hpi0w1->FixParameter(2,7.23065e+00);
		hpi0w1->FixParameter(3,2.73434e-02);
		hpi0w1->FixParameter(4,4.34588e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hpi0w2->FixParameter(0,1.15744e+00);
		hpi0w2->FixParameter(1,-8.46546e-01);
		hpi0w2->FixParameter(2,3.98888e-01);
		hpi0w2->FixParameter(3,7.19985e+00);
		hpi0w2->FixParameter(4,-1.50254e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,1.11960e+00);
		hpi0w3->FixParameter(1,9.98844e-01);
		hpi0w3->FixParameter(2,2.26967e-01);
		hpi0w3->FixParameter(3,7.49331e+00);
		hpi0w3->FixParameter(4,-6.41628e-03);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,3.02232e+00);
		hetaw1->FixParameter(1,-7.73934e+00);
		hetaw1->FixParameter(2,7.18446e+00);
		hetaw1->FixParameter(3,8.37699e-02);
		hetaw1->FixParameter(4,3.60050e-01);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,9.80086e-01);
		hetaw2->FixParameter(1,-4.79976e-01);
		hetaw2->FixParameter(2,3.50350e-01);
		hetaw2->FixParameter(3,9.47412e+00);
		hetaw2->FixParameter(4,-1.24664e-01);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0,1.04200e+00);
		hetaw3->FixParameter(1,-1.17124e+00);
		hetaw3->FixParameter(2,1.73411e-01);
		hetaw3->FixParameter(3,1.59364e+00);
		hetaw3->FixParameter(4,6.67547e-02);
		task->SetEtaWeight3(hetaw3);
		///--------------------------------------------------------------
		}
		*/
		
		
		if(isCharPion == 0){
		///Good weights for d20a2_extra (data over MB), where data = charged pions and eta from Mt-scaling (using |y| < 0.8 in MC for the weight calculation)
		
		//cout<<"-----------reference--------------------"<<endl;
		///Reference-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2.5);
		hpi0w1->FixParameter(0,1.00835e+01);
		hpi0w1->FixParameter(1,-8.79889e+00);
		hpi0w1->FixParameter(2,7.13270e+00);
		hpi0w1->FixParameter(3,3.30909e-02);
		hpi0w1->FixParameter(4,4.14692e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2.5,6);
		hpi0w2->FixParameter(0,1.07625e+00);
		hpi0w2->FixParameter(1,-3.99293e+00);
		hpi0w2->FixParameter(2,9.21251e-01);
		hpi0w2->FixParameter(3,8.76327e-02);
		hpi0w2->FixParameter(4,-1.00391e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,2.32912e+00);
		hpi0w3->FixParameter(1,-7.09878e-01);
		hpi0w3->FixParameter(2,8.37756e-02);
		hpi0w3->FixParameter(3,1.87128e+00);
		hpi0w3->FixParameter(4,1.88703e-01);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,1.46145e+00);
		hetaw1->FixParameter(1,-2.57763e+00);
		hetaw1->FixParameter(2,4.48690e+00);
		hetaw1->FixParameter(3,1.74529e+00);
		hetaw1->FixParameter(4,2.83520e-01);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,1.18725e+00);
		hetaw2->FixParameter(1,-1.87279e+00);
		hetaw2->FixParameter(2,6.36254e-01);
		hetaw2->FixParameter(3,1.60257e+00);
		hetaw2->FixParameter(4,-1.04365e-01);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0,1.36826e+00);
		hetaw3->FixParameter(1,1.00818e+00);
		hetaw3->FixParameter(2,3.23711e-01);
		hetaw3->FixParameter(3,8.42844e-08);
		hetaw3->FixParameter(4,-1.66152e-04);
		task->SetEtaWeight3(hetaw3);
		///--------------------------------------------------------------
		}
		
		if(isCharPion == 1){
			//cout<<"***********************tilt 1******************************"<<endl;
		///Tilt 1-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hpi0w1->FixParameter(0,8.70443e+00);
		hpi0w1->FixParameter(1,-7.49105e+00);
		hpi0w1->FixParameter(2,6.14059e+00);
		hpi0w1->FixParameter(3,4.76320e-02);
		hpi0w1->FixParameter(4,5.27575e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hpi0w2->FixParameter(0,1.49474e+00);
		hpi0w2->FixParameter(1,5.39796e-01);
		hpi0w2->FixParameter(2,1.69609e-01);
		hpi0w2->FixParameter(3,4.51710e+01);
		hpi0w2->FixParameter(4,-1.37151e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,1.12495e+00);
		hpi0w3->FixParameter(1,-9.96954e-02);
		hpi0w3->FixParameter(2,3.61900e-02);
		hpi0w3->FixParameter(3,1.06158e+01);
		hpi0w3->FixParameter(4,2.28435e-01);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,3.53156e+00);
		hetaw1->FixParameter(1,-7.42232e+00);
		hetaw1->FixParameter(2,6.96454e+00);
		hetaw1->FixParameter(3,8.35120e-02);
		hetaw1->FixParameter(4,4.10793e-01);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,9.38975e-01);
		hetaw2->FixParameter(1,-2.01011e+00);
		hetaw2->FixParameter(2,9.88095e-01);
		hetaw2->FixParameter(3,7.41818e+00);
		hetaw2->FixParameter(4,-6.80076e-02);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0,1.09401e+00);
		hetaw3->FixParameter(1,-6.60190e-01);
		hetaw3->FixParameter(2,1.15482e-01);
		hetaw3->FixParameter(3,2.74272e+00);
		hetaw3->FixParameter(4,1.44224e-01);
		task->SetEtaWeight3(hetaw3);
		}
		
		if(isCharPion == 2){
			//cout<<"-----------tilt 2--------------------"<<endl;
		///Tilt 2-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hpi0w1->FixParameter(0,5.13903e+00);
		hpi0w1->FixParameter(1,-7.25487e+00);
		hpi0w1->FixParameter(2,5.84647e+00);
		hpi0w1->FixParameter(3,6.02927e-02);
		hpi0w1->FixParameter(4,4.08393e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hpi0w2->FixParameter(0,1.10246e+00);
		hpi0w2->FixParameter(1,-1.09439e+00);
		hpi0w2->FixParameter(2,4.31339e-01);
		hpi0w2->FixParameter(3,5.61541e+00);
		hpi0w2->FixParameter(4,-1.59952e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,1.12386e+00);
		hpi0w3->FixParameter(1,1.10695e+00);
		hpi0w3->FixParameter(2,2.41931e-01);
		hpi0w3->FixParameter(3,7.41653e+00);
		hpi0w3->FixParameter(4,-5.66732e-02);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,2.78760e+00);
		hetaw1->FixParameter(1,-9.27114e+00);
		hetaw1->FixParameter(2,8.32542e+00);
		hetaw1->FixParameter(3,5.80300e-02);
		hetaw1->FixParameter(4,2.96995e-01);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,1.03427e+00);
		hetaw2->FixParameter(1,-1.03382e-01);
		hetaw2->FixParameter(2,1.99789e-01);
		hetaw2->FixParameter(3,9.97201e+00);
		hetaw2->FixParameter(4,-2.15465e-01);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0,9.34427e-01);
		hetaw3->FixParameter(1,8.84703e-01);
		hetaw3->FixParameter(2,2.48724e-01);
		hetaw3->FixParameter(3,8.12603e+00);
		hetaw3->FixParameter(4,-1.48239e-02);
		task->SetEtaWeight3(hetaw3);
		}
		///............................................................................................................................................
		
		
		///Weights tilted up and down calculated using the neutral pion systematic uncertainty
		///................................................................................................................................................
		if(isCharPion == 3){
		///Good weights for d20a2_extra (data over MB), where data = charged pions and eta from Mt-scaling (using |eta| < 1.2 in MC for the weight calculation)
		
		
		///Reference-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hpi0w1->FixParameter(0,6.28420e+00);
		hpi0w1->FixParameter(1,-4.73336e+00);
		hpi0w1->FixParameter(2,4.30371e+00);
		hpi0w1->FixParameter(3,2.05557e-01);
		hpi0w1->FixParameter(4,7.45169e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hpi0w2->FixParameter(0,1.06074e+00);
		hpi0w2->FixParameter(1,-1.65561e+00);
		hpi0w2->FixParameter(2,7.71429e-01);
		hpi0w2->FixParameter(3,3.93300e+00);
		hpi0w2->FixParameter(4,-1.86934e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,1.06619e+00);
		hpi0w3->FixParameter(1,-3.46033e-02);
		hpi0w3->FixParameter(2,3.15605e-02);
		hpi0w3->FixParameter(3,1.42625e+01);
		hpi0w3->FixParameter(4,2.94724e-01);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,2.65122e+00);
		hetaw1->FixParameter(1,-4.70401e+00);
		hetaw1->FixParameter(2,4.91259e+00);
		hetaw1->FixParameter(3,2.79371e-01);
		hetaw1->FixParameter(4,5.31804e-01);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,7.97299e-01);
		hetaw2->FixParameter(1,-1.82673e+00);
		hetaw2->FixParameter(2,8.46395e-01);
		hetaw2->FixParameter(3,2.51285e+00);
		hetaw2->FixParameter(4,-2.02592e-01);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0, 1.17724e+00);
		hetaw3->FixParameter(1, -5.18720e-01);
		hetaw3->FixParameter(2, 8.71248e-02);
		hetaw3->FixParameter(3, 2.80082e+00);
		hetaw3->FixParameter(4, 1.79024e-01);
		task->SetEtaWeight3(hetaw3);
		///--------------------------------------------------------------
		}
		
		
		if(isCharPion == 4){
			//cout<<"***********************tilt 1******************************"<<endl;
		///Tilt 1-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hpi0w1->FixParameter(0,9.11997e+00);
		hpi0w1->FixParameter(1,-5.55639e+00);
		hpi0w1->FixParameter(2,4.57881e+00);
		hpi0w1->FixParameter(3,7.49705e-02);
		hpi0w1->FixParameter(4,6.21982e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hpi0w2->FixParameter(0,1.11464e+00);
		hpi0w2->FixParameter(1,-1.79829e+00);
		hpi0w2->FixParameter(2,8.48301e-01);
		hpi0w2->FixParameter(3,4.52799e+00);
		hpi0w2->FixParameter(4,-1.32837e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,1.21393e+00);
		hpi0w3->FixParameter(1,-1.87943e-02);
		hpi0w3->FixParameter(2,2.66623e-02);
		hpi0w3->FixParameter(3,8.48609e+00);
		hpi0w3->FixParameter(4,5.39010e-01);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,4.44127e+00);
		hetaw1->FixParameter(1,-7.46162e+00);
		hetaw1->FixParameter(2, 6.89913e+00);
		hetaw1->FixParameter(3, 6.69980e-02);
		hetaw1->FixParameter(4, 4.57947e-01);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,8.15397e-01);
		hetaw2->FixParameter(1,-2.19968e+00);
		hetaw2->FixParameter(2,9.85861e-01);
		hetaw2->FixParameter(3,2.06653e+00);
		hetaw2->FixParameter(4,-1.53009e-01);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0,6.77599e-01);
		hetaw3->FixParameter(1,3.29809e-01);
		hetaw3->FixParameter(2,-2.49708e-03);
		hetaw3->FixParameter(3,1.83295e+01);
		hetaw3->FixParameter(4,4.76179e-01);
		task->SetEtaWeight3(hetaw3);
		}
		
		if(isCharPion == 5){
			//cout<<"-----------tilt 2--------------------"<<endl;
		///Tilt 2-----------------------------------------------------
		TF1 *hpi0w1 = new TF1("hpi0w1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hpi0w1->FixParameter(0,4.09359e+00);
		hpi0w1->FixParameter(1,-1.43022e+01);
		hpi0w1->FixParameter(2,1.05913e+01);
		hpi0w1->FixParameter(3,7.68920e-03);
		hpi0w1->FixParameter(4,2.06406e-01);
		task->SetPi0Weight(hpi0w1);
		
		TF1 *hpi0w2 = new TF1("hpi0w2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hpi0w2->FixParameter(0,1.00660e+00);
		hpi0w2->FixParameter(1,-2.97827e+00);
		hpi0w2->FixParameter(2,1.03280e+00);
		hpi0w2->FixParameter(3,3.68768e+00);
		hpi0w2->FixParameter(4,-1.46169e-01);
		task->SetPi0Weight2(hpi0w2);
		
		TF1 *hpi0w3 = new TF1("hpi0w3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hpi0w3->FixParameter(0,1.11703e+00);
		hpi0w3->FixParameter(1,1.28083e+00);
		hpi0w3->FixParameter(2,2.64689e-01);
		hpi0w3->FixParameter(3,7.14581e+00);
		hpi0w3->FixParameter(4,-1.95134e-01);
		task->SetPi0Weight3(hpi0w3);
		
		
		TF1 *hetaw1 = new TF1("hetaw1","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5,2);
		hetaw1->FixParameter(0,5.67745e-01);
		hetaw1->FixParameter(1, -1.95037e+01);
		hetaw1->FixParameter(2, 8.52115e+00);
		hetaw1->FixParameter(3, 3.79235e-05);
		hetaw1->FixParameter(4,-5.77702e-02);
		task->SetEtaWeight(hetaw1);
		
		
		TF1 *hetaw2 = new TF1("hetaw2","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 2,6);
		hetaw2->FixParameter(0,9.54625e-01);
		hetaw2->FixParameter(1,-1.62050e+00);
		hetaw2->FixParameter(2,7.18392e-01);
		hetaw2->FixParameter(3,7.91718e+00);
		hetaw2->FixParameter(4,-1.94481e-01);
		task->SetEtaWeight2(hetaw2);
		
		TF1 *hetaw3 = new TF1("hetaw3","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 6,20);
		hetaw3->FixParameter(0,7.58025e-01);
		hetaw3->FixParameter(1,6.31331e-02);
		hetaw3->FixParameter(2,-2.97603e-03);
		hetaw3->FixParameter(3,3.32068e+00);
		hetaw3->FixParameter(4,-1.94422e-01);
		task->SetEtaWeight3(hetaw3);
		}
		///............................................................................................................................................
		
		
		
		
		/*    
		if(!isCharPion){
		///Good weights for c3b2 (MB over enhanced+MB)
		TF1 *hpi0w = new TF1("hpi0w","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5, 5000);
		hpi0w->FixParameter(0, 1.01492e+00);
		hpi0w->FixParameter(1, 2.30839e-01);
		hpi0w->FixParameter(2, 4.07271e-02);
		hpi0w->FixParameter(3, 4.09513e+00);
		hpi0w->FixParameter(4, 4.50679e+00);
		task->SetPi0Weight(hpi0w);
    
		TF1 *hetaw = new TF1("hetaw","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5, 5000);
		hetaw->FixParameter(0, 9.92194e-01);
		hetaw->FixParameter(1, 2.64171e-01);
		hetaw->FixParameter(2, 4.52648e-02);
		hetaw->FixParameter(3, 3.68281e+00);
		hetaw->FixParameter(4, 4.43877e+00);
		task->SetEtaWeight(hetaw);
		}
		*/
		
		///Good weights for c3b2 (data over enhanced+MB), where data = charged pions and eta from Mt-scaling (using |eta| < 1.2 in MC for the weight calculation)
		/*
		TF1 *hpi0w = new TF1("hpi0w","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5, 100);
		hpi0w->FixParameter(0, 3.02495e+00);
		hpi0w->FixParameter(1, 1.54134e-01);
		hpi0w->FixParameter(2, 5.09308e-02);
		hpi0w->FixParameter(3, 3.35230e+00);
		hpi0w->FixParameter(4, 5.32343e+00);
		task->SetPi0Weight(hpi0w);
    
		TF1 *hetaw = new TF1("hetaw","[0] / (TMath::Power(TMath::Exp( - [1] * x[0] - [2] * x[0] * x[0] ) + x / [3], [4]))", 0.5, 100);
		hetaw->FixParameter(0, 4.75374e+00);
		hetaw->FixParameter(1, 1.99476e-01);
		hetaw->FixParameter(2, 4.24131e-02);
		hetaw->FixParameter(3, 3.08881e+00);
		hetaw->FixParameter(4, 5.34041e+00);
		task->SetEtaWeight(hetaw);
		*/
	
    
    
    
    
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
    
    //______________________________________________________
    //Configure PID
    
    //_________________________
    //TPC PID
    pid->AddDetector("TPC", 0);				//Add TPC PID
    
    //_________________________
    //Configure TPC cut
    //Defaul = -1 to 3 sigmas
    //Note that it is also possible to define a model instead of a constant
    //--------->For this change the "cut model"
    
    Double_t params[4];
    char *cutmodel;
    cutmodel = "pol0";
    params[0] = tpcPIDmincut;
    Double_t max= tpcPIDmaxcut;

    pid->ConfigureTPCdefaultCut(cutmodel,params,max);
    //________________________
    ///TOF PID
    pid->AddDetector("TOF", 1);
    //Automatizar cortes e colocar índices
    //___________________________
    //Configure	TOF cut: set number of sigmas of the TOF cut
    pid->ConfigureTOF(tofPID);
    ///_______________________________________________________________________________________________________________
    
    printf("*************************************\n");
    printf("Configuring standard Task:\n");
    pid->PrintStatus();
    printf("*************************************\n");
    
    return task;
}





