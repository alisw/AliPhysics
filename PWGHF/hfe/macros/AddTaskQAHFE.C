AliAnalysisTaskQAHFE *AddTaskQAHFE(
		TString HFEdir="",
		Bool_t isMC=kFALSE,
		TString estimatorFilename="",
		Double_t refMult=13.0793,
		
		AliVEvent::EOfflineTriggerTypes trigger=AliVEvent::kINT7 ,
		//Int_t trigger=AliVEvent::kINT7 ,
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
		Double_t AssoTPCnsig=  3.5
		)
{
  
  // AddTask for the AliAnalysisTaskQAHFE for HFE multiplicity Dependent study
  //====================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskQAHFEID", "No analysis manager to connect to.");
    return NULL;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 
 	//TString HFEdir="";
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
  AliAnalysisTaskQAHFE *taskhfe = new AliAnalysisTaskQAHFE(taskname.Data());
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
    
    TProfile* multEstimatorAvg[2];                       
    //for(Int_t ip=0; ip<4; ip++)
    //{
    	multEstimatorAvg[0] = (TProfile*)(fileEstimator->Get("hProf1")->Clone("hProf1"));
			multEstimatorAvg[1] = (TProfile*)(fileEstimator->Get("hProf1")->Clone("hProf1"));
			if (!multEstimatorAvg[0] || !multEstimatorAvg[1])
		 	{
	  		AliFatal("Multiplicity estimator not found! Please check your estimator file");
	  		return;
			}
    //}
    taskhfe->SetMultiplVsZProfile_NoEtaCut(multEstimatorAvg[0]);
    taskhfe->SetMultiplVsZProfile(multEstimatorAvg[1]);
  }
	taskhfe->SetMCAnalysis(isMC);
	
	cout<<" trigger "<<trigger<<endl;
	//getchar();
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
