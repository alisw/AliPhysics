AliAnalysisTask_QA_EMCALElectrons *AddTask_QA_EMCALElectrons(
		TString HFEdir="",
		Bool_t isMC=kFALSE,
		AliVEvent::EOfflineTriggerTypes trigger=AliVEvent::kINT7 ,
		Bool_t isEG1=kFALSE,
		TString estimatorFilename="",
		Double_t refMult=11.76,
		
		Int_t TPCNclus=100  ,
		Int_t ITSNclus= 3 ,
		Int_t TPCNclusPID= 80 ,
		Bool_t SPDBoth= kTRUE ,
		Bool_t SPDAny= kFALSE ,
		Bool_t SPDFirst= kFALSE ,
		Double_t DCAxyCut= 1 ,
		Double_t DCAzCut=2  ,
		Double_t Etarange= 0.7 ,
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
  
  // AddTask for the AliAnalysisTask_QA_EMCALElectrons for HFE multiplicity Dependent study
  //====================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_QA_EMCALElectrons", "No analysis manager to connect to.");
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
  AliAnalysisTask_QA_EMCALElectrons *taskhfe = new AliAnalysisTask_QA_EMCALElectrons(taskname.Data());
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
    
    const Char_t* periodNames[4] = { "16l", "17d20a2_extra","16k", "17d20a1_extra"};
 
    TProfile* multEstimatorAvg[4];
    for(Int_t ip=0; ip<4; ip++) {
      cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
      multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
      if (!multEstimatorAvg[ip]) {
	AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	return;
      }
    }
    
    taskhfe->SetMultiplVsZProfile_16l(multEstimatorAvg[0]);
    taskhfe->SetMultiplVsZProfile_17d20a2_extra(multEstimatorAvg[1]);
    taskhfe->SetMultiplVsZProfile_16k(multEstimatorAvg[2]);
    taskhfe->SetMultiplVsZProfile_17d20a1_extra(multEstimatorAvg[3]);
  }
  
  	//taskhfe->SelectCollisionCandidates(AliVEvent::kINT5);
  	//taskhfe->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  	taskhfe->SelectCollisionCandidates(trigger);
  	taskhfe->SetMCAnalysis(isMC);
	taskhfe->SetTrigger(trigger);
	taskhfe->SetEtaRange(Etarange);
	taskhfe->SetMinTPCCluster(TPCNclus);
	taskhfe->SetMinITSCluster(ITSNclus);
	taskhfe->SetMinTPCClusterPID(TPCNclusPID);
	taskhfe->SetHitsOnSPDLayers(SPDBoth,SPDAny,SPDFirst);
	taskhfe->SetDCACut(DCAxyCut,DCAzCut);
	taskhfe->SetTPCnsigma(TPCnsigmin,TPCnsigmax);
	

	taskhfe->SetInvMassCut(InvmassCut);
	taskhfe->SetAssoTPCclus(AssoTPCCluster);
	taskhfe->SetAssoITSrefit(AssoITSRefit);
	taskhfe->SetAssopTMin(AssopTMin);
	taskhfe->SetAssoEtarange(AssoEtarange);
	taskhfe->SetAssoTPCnsig(AssoTPCnsig);
	
	
	if(trigger==AliVEvent::kINT7){
		
		isEG1=kFALSE;
		taskhfe->SetEMCalTriggerEG1(kFALSE);
		taskhfe->SetEMCalTriggerDG1(kFALSE);
		taskhfe->SetEMCalTriggerEG2(kFALSE);
		taskhfe->SetEMCalTriggerDG2(kFALSE);
		cout<<"1 trigger  "<<trigger<<"   "<< isEG1 <<endl;
		//getchar();
	}
	if(trigger==AliVEvent::kEMCEGA &&  isEG1==kTRUE){
		taskhfe->SetEMCalTriggerEG1(kTRUE);
		taskhfe->SetEMCalTriggerDG1(kTRUE);
		taskhfe->SetEMCalTriggerEG2(kFALSE);
		taskhfe->SetEMCalTriggerDG2(kFALSE);
		cout<<" 2 trigger  "<<trigger<<"   "<< isEG1 <<endl;
		//getchar();
	}
	if(trigger==AliVEvent::kEMCEGA && isEG1==kFALSE){
		taskhfe->SetEMCalTriggerEG1(kFALSE);
		taskhfe->SetEMCalTriggerDG1(kFALSE);
		taskhfe->SetEMCalTriggerEG2(kTRUE);
		taskhfe->SetEMCalTriggerDG2(kTRUE);
		cout<<"3 trigger  "<<trigger<<"   "<< isEG1 <<endl;
		//getchar();
	}
	
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
