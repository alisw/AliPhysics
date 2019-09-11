#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH($ALICE_ROOT)

AliAnalysisHFEppTPCTOFBeauty* ConfigHFEppTPCTOFBeauty(Bool_t isMCc, Bool_t isAODc, Bool_t isPPc,Double_t tpcPIDmincut,Double_t tpcPIDmaxcut, Double_t tofPID, Int_t minNClustersTPC, Int_t minNClustersTPCPID, Float_t minRatioTPCclusters, Int_t  minNClustersITS, AliHFEextraCuts::ITSPixel_t pixel, Float_t EtaMin, Float_t EtaMax, Float_t DCAxy, Float_t DCAz, Int_t IsBcorr, Int_t IsDcorr);

#endif

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
            Int_t IsDcorr = 0,
            Bool_t fIsCalculateBDMesonpTWeights = kTRUE
            
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
    AliAnalysisHFEppTPCTOFBeauty *task = ConfigHFEppTPCTOFBeauty(isMC,isAOD,isPP,tpcPIDmincut,tpcPIDmaxcut,tofPIDcut,MinNClustersTPC,MinNClustersTPCPID,MinRatioTPCclusters,MinNClustersITS,pixel,EtaMin,EtaMax,DCAxy,DCAz,IsBcorr,IsDcorr,fIsCalculateBDMesonpTWeights);
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


AliAnalysisHFEppTPCTOFBeauty* ConfigHFEppTPCTOFBeauty(Bool_t isMCc, Bool_t isAODc, Bool_t isPPc,Double_t tpcPIDmincut,Double_t tpcPIDmaxcut, Double_t tofPID, Int_t minNClustersTPC, Int_t minNClustersTPCPID, Float_t minRatioTPCclusters, Int_t  minNClustersITS, AliHFEextraCuts::ITSPixel_t pixel, Float_t EtaMin, Float_t EtaMax, Float_t DCAxy, Float_t DCAz, Int_t IsBcorr, Int_t IsDcorr, Double_t fIsCalculateBDMesonpTWeights)
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
        task->SetBDMesonpTWeightCalc(fIsCalculateBDMesonpTWeights);
    }
    //______________________________________
   
   // Function for obtaining the Hadron Contamination 
    
   		//TF1 *fHadCont = new TF1("fHadCont","5.00674e+00*TMath::Landau(x[0],1.03540e+01,2.60371e+00) + 5.67850e-02*TMath::Gaus(x[0], 1.00000e+00, 9.25799e-02)",1.0,10.0);  
		//task->SetHCFunction(fHadCont);
    
    if(isMCc){
    //D correction
    
   // TFile *f = new TFile("DMesonWeight_pp5TeV_Aug30.root");
    
     TString filenameDweight = "DMesonWeight_pp5TeV_Aug30.root";
     //TFile *DweightFileEnh = TFile::Open(Form("%s", filenameDweight.Data()));
     TFile *DweightFileEnh = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", filenameDweight.Data()));
    
     TGraphErrors *fDmeson = (TGraphErrors *)DweightFileEnh->Get("DCentral");
   
   
     TString filenameStrweight = "KPi_ratio_5TeV_pp.root";
     //TFile *StrweightFileEnh = TFile::Open(Form("%s", filenameStrweight.Data()));
     TFile *StrweightFileEnh = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", filenameStrweight.Data()));
     //TFile *fSt = new TFile("KPi_ratio_5TeV_pp.root");
     TF1 *fStrangeWeight = (TF1 *)StrweightFileEnh->Get("f2");
     task->SetKtoPiFunction(fStrangeWeight);
 
    ///B correction
    //TFile *f = new TFile("FractionWeights_Graph_Aug01.root");
   // TFile *f = new TFile("FractionWeights_Graph_Sep02_Final.root");
    
     TString filenameFractweight = "FractionWeights_Graph_Sep02_Final.root";
   //TFile *FractweightFileEnh = TFile::Open(Form("%s", filenameFractweight.Data()));
   TFile *FractweightFileEnh = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", filenameFractweight.Data()));
    
    TGraphErrors *fLc = (TGraphErrors *)FractweightFileEnh->Get("hWeight_Lc");
    TGraphErrors *fDp = (TGraphErrors *)FractweightFileEnh->Get("hWeight_Dp");
    TGraphErrors *fDstar = (TGraphErrors *)FractweightFileEnh->Get("hWeight_Dstar");
    TGraphErrors *fDs = (TGraphErrors *)FractweightFileEnh->Get("hWeight_Ds");
    TGraphErrors *fD0 = (TGraphErrors *)FractweightFileEnh->Get("hWeight_D0");
    task->SetLcD0Function(fLc);
    task->SetDpD0Function(fDp);
    task->SetDstarD0Function(fDstar);
    task->SetDsD0Function(fDs);
    task->SetD0D0Function(fD0);
  // TString filenameBweight = "BMesonWeight_pp5TeV_ScaledToOne_July29.root";
    TString filenameBweight = "BMesonWeight_pp5TeV_Aug29.root";
   //TFile *BweightFileEnh = TFile::Open(Form("%s", filenameBweight.Data()));
   TFile *BweightFileEnh = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", filenameBweight.Data()));
    
    ///Default FONLL
    if(IsBcorr == 0){
    
    TGraphErrors *hWBEnhCentral = (TGraphErrors*)BweightFileEnh->Get("BFonllCentral");
		
		TF1 *fBmesonShape = new TF1("fBmesonShape","(1.0/0.530432)*(5.78258e-03 / (TMath::Power(TMath::Exp( - 2.07792e+00*x[0] + 6.75483e-03 * x[0] * x[0] ) + x[0] / 1.66686e+01, 7.05243e-01))*(2.27740e-06 / (TMath::Power(TMath::Exp( - 1.72205e-01*x[0] + 3.69234e-03 * x[0] * x[0] ) + x[0] / 7.12479e+00, 4.34903e+00)))*TMath::Power(1.0+x[0]/1.79842e-05,6.74073e-01))/(5.22379e-02 / (TMath::Power(TMath::Exp( - 1.27517e+00*x[0] - 3.54373e-02 * x[0] * x[0] ) + x[0] / 5.19647e+00, 1.31823e+00))*(1.77593e-03 / (TMath::Power(TMath::Exp( - 5.14768e-01*x[0] + 1.03322e-02 * x[0] * x[0] ) + x[0] / 4.18420e+01, 9.51497e-01)))*TMath::Power(1.0+x[0]/1.44862e+02,-1.59238e+01))",0.0,50);
		
		task->SetBcorrFunction(hWBEnhCentral);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///FONLL Lower edge
    if(IsBcorr == 1){
    
    TGraphErrors *hWBEnhLower = (TGraphErrors*)BweightFileEnh->Get("BFonllLower");
		
		TF1 *fBmesonShape1 = new TF1("fBmesonShape1","(1.0/0.441646)*(3.77085e-03 / (TMath::Power(TMath::Exp( - 4.99937e+00*x[0] - 8.02258e-01 * x[0] * x[0] ) + x[0] / 2.60310e+04, 1.89019e-01))*(1.43903e-03 / (TMath::Power(TMath::Exp( - 1.82037e-01*x[0] + 3.88349e-03 * x[0] * x[0] ) + x[0] / 7.21515e+00, 4.49037e+00)))*TMath::Power(1.0+x[0]/2.91057e-01,4.05964e-01))/(5.22379e-02 / (TMath::Power(TMath::Exp( - 1.27517e+00*x[0] - 3.54373e-02 * x[0] * x[0] ) + x[0] / 5.19647e+00, 1.31823e+00))*(1.77593e-03 / (TMath::Power(TMath::Exp( - 5.14768e-01*x[0] + 1.03322e-02 * x[0] * x[0] ) + x[0] / 4.18420e+01, 9.51497e-01)))*TMath::Power(1.0+x[0]/1.44862e+02,-1.59238e+01))",0.0,50);
		
		task->SetBcorrFunction(hWBEnhLower);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
    ///FONLL Upper edge
    if(IsBcorr == 2){
    
     TGraphErrors *hWBEnhUpper = (TGraphErrors*)BweightFileEnh->Get("BFonllUpper");
		
		TF1 *fBmesonShape2 = new TF1("fBmesonShape2","(1.0/0.58236)*(8.25352e-03 / (TMath::Power(TMath::Exp( - 6.11267e+00*x[0] - 1.05995e-01 * x[0] * x[0] ) + x[0] / 3.12399e+04, 2.13332e-01))*(1.70947e-03 / (TMath::Power(TMath::Exp( - 1.75188e-01*x[0] + 3.76833e-03 * x[0] * x[0] ) + x[0] / 6.82438e+00, 4.61154e+00)))*TMath::Power(1.0+x[0]/3.25350e-01,3.45204e-01))/(5.22379e-02 / (TMath::Power(TMath::Exp( - 1.27517e+00*x[0] - 3.54373e-02 * x[0] * x[0] ) + x[0] / 5.19647e+00, 1.31823e+00))*(1.77593e-03 / (TMath::Power(TMath::Exp( - 5.14768e-01*x[0] + 1.03322e-02 * x[0] * x[0] ) + x[0] / 4.18420e+01, 9.51497e-01)))*TMath::Power(1.0+x[0]/1.44862e+02,-1.59238e+01))",0.0,50);
		
		task->SetBcorrFunction(hWBEnhUpper);
		cout<<"-----------------------------------------------------IsBcorr"<<IsBcorr<<endl;
	}
 
    // Default
    if(IsDcorr == 0){
    		TF1 *fDmesonShape = new TF1("fDmesonShape","(1.44336e-03/(TMath::Power(TMath::Exp( - 6.08651e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50784e+00, 4.40766e+00)))/(1.60027e-02 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69999e+00, 4.51198e+00)))",0.0,36);  
		task->SetDcorrFunction(fDmeson);
    		/*TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/0.296109)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.0,36);
    
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/0.285340)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/0.274920)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/0.264896)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/0.255304)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/0.246166)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/0.2253789)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/0.2074933)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/0.19228)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/0.179479)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/0.1686951)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/0.159618)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/0.1454598)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/0.1351671)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/0.1274891)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/0.1216022)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.116966)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.113223)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.110133)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.1075296)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.103355)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.10011)*(1.42911e-03/(TMath::Power(TMath::Exp( - 6.08657e-01*x[0] + 9.8566e-03 * x[0] * x[0] ) + x[0] / 3.50791e+00, 4.40766e+00)))/(3.53468e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.20266e+00)))",9.0,36);  
		task->SetDcorrFunction15(fDmesonShape15);*/
	}
    
    // Data up tilted 
    
     if(IsDcorr == 1){
     
     		TF1 *fDmesonShape = new TF1("fDmesonShape","(1.92351e-03/(TMath::Power(TMath::Exp( - 5.46426e-01*x[0] + 9.65660e-03 * x[0] * x[0] ) + x[0] / 3.48352e+00, 4.65766e+00)))/(1.02198e-02 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69999e+00, 4.51197e+00)))",0.0,36);  
		task->SetDcorrFunction(fDmesonShape);
     
		/*TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/0.318842)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.0,36);  
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/0.307014)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/0.295766)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/0.285108)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/0.275036)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/0.265544)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/0.244238)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/0.226105)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/0.210722)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/0.19766)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/0.186521)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/0.176961)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/0.161475)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/0.149485)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/0.139886)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/0.131974)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.125294)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.119545)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.114522)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.11008)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.10254)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.0963444)*(1.63798e-03/(TMath::Power(TMath::Exp( - 5.59031e-01*x[0] + 9.65659e-03 * x[0] * x[0] ) + x[0] / 3.60055e+00, 4.65766e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",9.0,36);  
		task->SetDcorrFunction15(fDmesonShape15);*/
	}
   
    
    // Data Down tilted
    
     if(IsDcorr == 2){
		TF1 *fDmesonShape22 = new TF1("fDmesonShape22","(1/0.276203)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.0,36);  
		task->SetDcorrFunction22(fDmesonShape22);
    
		TF1 *fDmesonShape21 = new TF1("fDmesonShape21","(1/0.267907)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.1,36);
		task->SetDcorrFunction21(fDmesonShape21);
    
		TF1 *fDmesonShape20 = new TF1("fDmesonShape20","(1/0.259826)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.2,36);
		task->SetDcorrFunction20(fDmesonShape20);
    
		TF1 *fDmesonShape19 = new TF1("fDmesonShape19","(1/0.251974)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.3,36);
		task->SetDcorrFunction19(fDmesonShape19);
    
		TF1 *fDmesonShape18 = new TF1("fDmesonShape18","(1/0.244365)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.4,36);
		task->SetDcorrFunction18(fDmesonShape18);
       
		TF1 *fDmesonShape17 = new TF1("ffDmesonShape17","(1/0.237012)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.5,36); 
		task->SetDcorrFunction17(fDmesonShape17);
    
		TF1 *fDmesonShape16 = new TF1("fDmesonShape16","(1/0.2198)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",1.75,36);
		task->SetDcorrFunction16(fDmesonShape16);
    
		TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/0.204322)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.0,36); 
		task->SetDcorrFunction1(fDmesonShape1);
    
		TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/0.190571)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.25,36);   
		task->SetDcorrFunction2(fDmesonShape2);
    
		TF1 *fDmesonShape3 = new TF1("fDmesonShap3","(1/0.178465)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.5,36);   
		task->SetDcorrFunction3(fDmesonShape3);
    
		TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/0.167877)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",2.75,36);    
		task->SetDcorrFunction4(fDmesonShape4);
    
		TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/0.158653)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.0,36);  
		task->SetDcorrFunction5(fDmesonShape5);
    
		TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/0.143667)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",3.5,36);   
		task->SetDcorrFunction6(fDmesonShape6);
    
		TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/0.132336)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.0,36);   
		task->SetDcorrFunction7(fDmesonShape7);
    
		TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/0.123712)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",4.5,36);    
		task->SetDcorrFunction8(fDmesonShape8);
    
		TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/0.117083)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.0,36);  
		task->SetDcorrFunction9(fDmesonShape9);
    
		TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/0.111929)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",5.5,36);   
		task->SetDcorrFunction10(fDmesonShape10);
    
		TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/0.107876)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.0,36);   
		task->SetDcorrFunction11(fDmesonShape11);
    
		TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/0.10465)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",6.5,36);    
		task->SetDcorrFunction12(fDmesonShape12);
    
		TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/0.102053)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",7.0,36);    
		task->SetDcorrFunction13(fDmesonShape13);
    
		TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/0.0981935)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",8.0,36);  
		task->SetDcorrFunction14(fDmesonShape14);
    
		TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/0.0955098)*(1.31660e-03/(TMath::Power(TMath::Exp( - 6.08516e-01*x[0] - 1.38420e-02 * x[0] * x[0] ) + x[0] / 3.44002e+00, 4.24401e+00)))/(3.57173e-03 / (TMath::Power(TMath::Exp( - 6.05066e-01*x[0] + 1.62296e-02 * x[0] * x[0] ) + x[0] / 4.69154e+00, 4.14266e+00)))",9.0,36);  
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







