AliAnalysisTaskHFEmultTPCTOF *AddTaskHFEmultTPCTOF(
		TString HFEdir="",
		Bool_t isMC=kFALSE,
		TString estimatorFilename="",
		Double_t refMult=11.76,
		Char_t *periodName="16k_MB",
		
		AliVEvent::EOfflineTriggerTypes trigger=AliVEvent::kINT7 ,
		Int_t TPCNclus=100  ,
		Int_t ITSNclus= 3 ,
		Int_t TPCNclusPID= 80 ,
		Bool_t SPDBoth= kTRUE ,
		Bool_t SPDAny= kFALSE ,
		Bool_t SPDFirst= kFALSE ,
		Double_t DCAxyCut= 1 ,
		Double_t DCAzCut=2  ,
		Double_t Etarange= 0.8 ,
		Double_t TPCnsigmin= -1 ,
		Double_t TPCnsigmax= 3 ,
		Double_t TOFnsig= 3 ,
	
		Double_t InvmassCut= 0.14,		
		Int_t AssoTPCCluster= 60 ,
		Bool_t AssoITSRefit= kTRUE ,
		Double_t AssopTMin= 0.1  ,
		Double_t AssoEtarange= 0.9 ,
		Double_t AssoTPCnsig=  3.5,
		Int_t AssoTPCNclsForPID= 60 ,
		Int_t AssoITSNclus= 2 ,
		Float_t AssoDCAxy= 1 ,
		Float_t AssoDCAz= 2 
		)
{
  
  // AddTask for the AliAnalysisTaskHFEmultTPCTOF for HFE multiplicity Dependent study
  //====================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFEmultTPCTOF", "No analysis manager to connect to.");
    return NULL;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 
	TString outname1="Spectra";
	TString outname2="nEntries";

  TString filename="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":HFElectronOutput";
	filename +=HFEdir.Data();
	
  outname1 +=HFEdir.Data();
  outname2 +=HFEdir.Data();
 

       

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++
  TString taskname="ElecAnalysis";
  AliAnalysisTaskHFEmultTPCTOF *taskhfe = new AliAnalysisTaskHFEmultTPCTOF(taskname.Data());
  //taskhfe->SetDebugLevel(2);
  
	if(estimatorFilename.EqualTo("") )
 	{
    printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
  } 
  else
  {              
    TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator){
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
    
    taskhfe->SetReferenceMultiplicity(refMult);
    
    const Char_t* profilebasename="SPDtr";
    
    const Char_t* periodNames[9] = { "16k","18f1_extra","18f4b","18period","18_GPMC","18l5b","17period","17_GPMC","18l5a"};
   
    Int_t period=0;
    for(Int_t ip=0; ip< 9; ip++) {
		if(periodNames[ip]==periodName){ period =ip; break;}
     }
   
    TProfile* multEstimatorAvg[9];
    for(Int_t ip=0; ip<9; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvg[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }
    
    cout<<"==========================================================="<<endl;
    cout<<" period = "<<period<<" Name = "<<periodNames[period]<<endl;
    cout<<"==========================================================="<<endl;
    taskhfe->SetMultiplVsZProfile(multEstimatorAvg[period]);    
    taskhfe->SetEstimatorHistogram(period);
  }
  
	taskhfe->SetMCAnalysis(isMC);
	taskhfe->SetTrigger(trigger);
	taskhfe->SetEtaRange(Etarange);
	taskhfe->SetMinTPCCluster(TPCNclus);
	taskhfe->SetMinITSCluster(ITSNclus);
	taskhfe->SetMinTPCClusterPID(TPCNclusPID);
	taskhfe->SetHitsOnSPDLayers(SPDBoth,SPDAny,SPDFirst);
	taskhfe->SetDCACut(DCAxyCut,DCAzCut);
	taskhfe->SetTPCnsigma(TPCnsigmin,TPCnsigmax);
	taskhfe->SetTOFnsigma(TOFnsig);
	

	taskhfe->SetInvMassCut(InvmassCut);
	taskhfe->SetAssoTPCclus(AssoTPCCluster);
	taskhfe->SetAssoITSrefit(AssoITSRefit);
	taskhfe->SetAssopTMin(AssopTMin);
	taskhfe->SetAssoEtarange(AssoEtarange);
	taskhfe->SetAssoTPCnsig(AssoTPCnsig);
	
	taskhfe->SetAssoTPCNclsForPID(AssoTPCNclsForPID);
	taskhfe->SetAssoITSNclus(AssoITSNclus);
	taskhfe->SeAssoDCA(AssoDCAxy,AssoDCAz);

	
  mgr->AddTask(taskhfe);

  //_________Structure of Task O/P
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outname1,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(outname2,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); 
  

  //____________ Connect input/output
  mgr->ConnectInput(taskhfe,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskhfe,1,coutput1);
  mgr->ConnectOutput(taskhfe,2,coutput2);
  

  return taskhfe;
}
