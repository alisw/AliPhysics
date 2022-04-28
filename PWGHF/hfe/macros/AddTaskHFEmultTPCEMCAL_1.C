AliAnalysisTaskHFEmultTPCEMCAL_1 *AddTaskHFEmultTPCEMCAL_1(
		TString HFEdir="",
		Bool_t isMC=kFALSE,
		AliVEvent::EOfflineTriggerTypes trigger=AliVEvent::kINT7 ,
		Bool_t isEG1=kFALSE,
		TString estimatorFilename="",
		Double_t refMult=11.76,
		Char_t *periodName="16k_MB",
		
		Int_t TPCNclus=100  ,
		Int_t Ratioclus=0.8,
		Int_t ITSNclus= 3 ,
		Int_t TPCNclusPID= 80 ,
		Bool_t SPDBoth= kFALSE ,
		Bool_t SPDAny= kTRUE ,
		Bool_t SPDFirst= kFALSE ,
		Double_t DCAxyCut= 1 ,
		Double_t DCAzCut=2  ,
		Double_t Etarange= 0.7 ,
		Double_t TPCnsigmin= -1 ,
		Double_t TPCnsigmax= 3 ,
		Double_t TOFnsig= 3 ,
		Double_t EopEMin= 0.8 ,		
		Double_t EopEMax=  1.2,		
              Double_t  M20Min= 0.02 ,		
		Double_t M20Max1= 0.9,
              Double_t M20Max2= 0.7,
              Double_t M20Max3= 0.5,
		Double_t InvmassCut= 0.14,		
		Int_t AssoTPCCluster= 60 ,
		Bool_t AssoITSRefit= kTRUE ,
		Double_t AssopTMin= 0.1  ,
		Double_t AssoEtarange= 0.9 ,
		Double_t AssoTPCnsig=  3.5,
		Double_t Deltaeta = 0.01,
		Double_t Deltaphi = 0.01,
		Bool_t ClsTypeEMC=kTRUE, 
		Bool_t ClsTypeDCAL=kTRUE
		)
{
  
  // AddTask for the AliAnalysisTaskHFEmultTPCEMCAL_1 for HFE multiplicity Dependent study
  //====================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFEmultTPCEMCAL_1", "No analysis manager to connect to.");
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
  AliAnalysisTaskHFEmultTPCEMCAL_1 *taskhfe = new AliAnalysisTaskHFEmultTPCEMCAL_1(taskname.Data());
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
    
    const Char_t* periodNames[15] = { "16k_MB","16k_EG2","16k_EG1", "18f1_extra", "18f4b","18period_MB","18period_EG2","18period_EG1", "18_GPMC", "18l5b","17period_MB","17period_EG2","17period_EG1", "17_GPMC", "18l5a"};
    
    Int_t period=0;
	for(Int_t ip=0; ip<15; ip++) {
		if(periodNames[ip]==periodName){ period =ip; break;}
	}
	
    TProfile* multEstimatorAvg[15];
    for(Int_t ip=0; ip<15; ip++) {
      cout<<ip<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
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
  
  	//taskhfe->SelectCollisionCandidates(AliVEvent::kINT5);
  	//taskhfe->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  	taskhfe->SelectCollisionCandidates(trigger);
  	taskhfe->SetMCAnalysis(isMC);
	taskhfe->SetTrigger(trigger);
	taskhfe->SetEtaRange(Etarange);
	taskhfe->SetMinTPCCluster(TPCNclus);
	taskhfe->SetMinRatioCrossedRowOverFindable(Ratioclus);
	taskhfe->SetMinITSCluster(ITSNclus);
	taskhfe->SetMinTPCClusterPID(TPCNclusPID);
	taskhfe->SetHitsOnSPDLayers(SPDBoth,SPDAny,SPDFirst);
	taskhfe->SetDCACut(DCAxyCut,DCAzCut);
	taskhfe->SetTPCnsigma(TPCnsigmin,TPCnsigmax);
	taskhfe->SetEopE(EopEMin,EopEMax);
       taskhfe->SetShowerShapeEM20(M20Min,M20Max1,M20Max2,M20Max3);
       taskhfe->SetClusterTypeEMC(ClsTypeEMC);
       taskhfe->SetClusterTypeDCAL(ClsTypeDCAL);

	taskhfe->SetInvMassCut(InvmassCut);
	taskhfe->SetAssoTPCclus(AssoTPCCluster);
	taskhfe->SetAssoITSrefit(AssoITSRefit);
	taskhfe->SetAssopTMin(AssopTMin);
	taskhfe->SetAssoEtarange(AssoEtarange);
	taskhfe->SetAssoTPCnsig(AssoTPCnsig);
	taskhfe->SetDeltaEtaDeltaPhi(Deltaeta,Deltaphi);
	
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
	
	       if(ClsTypeEMC && ClsTypeDCAL){
		       taskhfe->SetEMCalTriggerEG1(kTRUE);
		       taskhfe->SetEMCalTriggerDG1(kTRUE);
		}
		if(ClsTypeEMC && !ClsTypeDCAL){
		       taskhfe->SetEMCalTriggerEG1(kTRUE);
		       taskhfe->SetEMCalTriggerDG1(kFALSE);
		}
		if(!ClsTypeEMC && ClsTypeDCAL){
		       taskhfe->SetEMCalTriggerEG1(kFALSE);
		       taskhfe->SetEMCalTriggerDG1(kTRUE);
		}
		
		taskhfe->SetEMCalTriggerEG2(kFALSE);
		taskhfe->SetEMCalTriggerDG2(kFALSE);
		
		cout<<" 2 trigger  "<<trigger<<"   "<< isEG1 <<endl;
		//getchar();
	}
	if(trigger==AliVEvent::kEMCEGA && isEG1==kFALSE){
	
		taskhfe->SetEMCalTriggerEG1(kFALSE);
		taskhfe->SetEMCalTriggerDG1(kFALSE);
		
		
		if(ClsTypeEMC && ClsTypeDCAL){
		       taskhfe->SetEMCalTriggerEG2(kTRUE);
		       taskhfe->SetEMCalTriggerDG2(kTRUE);
		}
		if(ClsTypeEMC && !ClsTypeDCAL){
		       taskhfe->SetEMCalTriggerEG2(kTRUE);
		       taskhfe->SetEMCalTriggerDG2(kFALSE);
		}
		if(!ClsTypeEMC && ClsTypeDCAL){
		       taskhfe->SetEMCalTriggerEG2(kFALSE);
		       taskhfe->SetEMCalTriggerDG2(kTRUE);
		}
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
