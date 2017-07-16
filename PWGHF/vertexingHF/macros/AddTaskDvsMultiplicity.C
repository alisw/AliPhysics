AliAnalysisTaskSEDvsMultiplicity *AddTaskDvsMultiplicity(Int_t system=0,
                                                         Bool_t readMC=kFALSE,
                                                         Int_t MCOption=0,
                                                         Int_t pdgMeson=411,
                                                         TString finDirname="Loose",
                                                         TString filename="",
                                                         TString finAnObjname="AnalysisCuts",
                                                         TString estimatorFilename="",
                                                         Double_t refMult=9.26,
                                                         Bool_t subtractDau=kFALSE,
                                                         Int_t NchWeight=0,
                                                         Int_t recoEstimator = AliAnalysisTaskSEDvsMultiplicity::kNtrk10,
                                                         Int_t MCEstimator = AliAnalysisTaskSEDvsMultiplicity::kEta10,
                                                         Bool_t isPPbData=kFALSE,
                                                         Int_t year = 16)
{
  //
  // Macro for the AliAnalysisTaskSE for D candidates vs Multiplicity
  // Invariant mass histogram in pt and multiplicity bins in a 3D histogram
  //   different estimators implemented
  //==============================================================================
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDvsMultiplicity", "No analysis manager to connect to.");
  }
    
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE;
  } else {
    filecuts=TFile::Open(filename.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      AliFatal("Input file not found : check your cut object");
    }
  }
    
    
  //Analysis Task
  AliRDHFCuts *analysiscuts=0x0;
    
  TString Name="";
  if(pdgMeson==411){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);
    Name="Dplus";
  }else if(pdgMeson==421){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(finAnObjname);
    Name="D0";
  }else if(pdgMeson==413){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(finAnObjname);
    Name="DStar";
  }else if(pdgMeson==431){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDstoKKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(finAnObjname);
    Name="Ds";
  }
    
  AliAnalysisTaskSEDvsMultiplicity *dMultTask = new AliAnalysisTaskSEDvsMultiplicity("dMultAnalysis",pdgMeson,analysiscuts,isPPbData);
  dMultTask->SetReadMC(readMC);
  dMultTask->SetDebugLevel(0);
  dMultTask->SetUseBit(kTRUE);
  dMultTask->SetDoImpactParameterHistos(kFALSE);
  dMultTask->SetSubtractTrackletsFromDaughters(subtractDau);
  dMultTask->SetMultiplicityEstimator(recoEstimator);
  dMultTask->SetMCPrimariesEstimator(MCEstimator);
  dMultTask->SetMCOption(MCOption);
  if(isPPbData) dMultTask->SetIsPPbData();
    
  if(NchWeight){
    TH1F *hNchPrimaries = NULL;
    TH1F *hMeasNchPrimaries = NULL;
    if(NchWeight==1){
      if(isPPbData) {
	hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithDWeight"); // MC distribution
      }
      else hNchPrimaries = (TH1F*)filecuts->Get("hGenPrimaryParticlesInelGt0");
      if(hNchPrimaries) {
	dMultTask->UseMCNchWeight(NchWeight);
	dMultTask->SetHistoNchWeight(hNchPrimaries);
      } else {
	AliFatal("Histogram for Nch multiplicity weights not found");
	return 0x0;
      }
      hMeasNchPrimaries = (TH1F*)filecuts->Get("hMeasNtrUnCorrEvWithD"); // data distribution
      if(hMeasNchPrimaries) {
	dMultTask->SetMeasuredNchHisto(hMeasNchPrimaries);
      }
    }
    else if(NchWeight==2){
      hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithDWeight"); // MC distribution
      hMeasNchPrimaries = (TH1F*)filecuts->Get("hMeasNtrUnCorrEvWithD"); // data distribution
      if(hNchPrimaries && hMeasNchPrimaries) {
	dMultTask->UseMCNchWeight(NchWeight);
	dMultTask->SetHistoNchWeight(hNchPrimaries);
	dMultTask->SetMeasuredNchHisto(hMeasNchPrimaries);
      } else {
	AliFatal("Histogram for Ntrk multiplicity weights not found");
	return 0x0;
      }
    }
  }
    
    
  if(pdgMeson==421) {
    dMultTask->SetMassLimits(1.5648,2.1648);
    dMultTask->SetNMassBins(200);
  }else if(pdgMeson==411)dMultTask->SetMassLimits(pdgMeson,0.2);
    
  if(estimatorFilename.EqualTo("") ) {
    printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
  } else{
        
    TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      AliFatal("File with multiplicity estimator not found\n");
      return;
    }
        
    dMultTask->SetReferenceMultiplcity(refMult);
        
    const Char_t* profilebasename="SPDmult10";
    if(recoEstimator==AliAnalysisTaskSEDvsMultiplicity::kVZEROA || recoEstimator==AliAnalysisTaskSEDvsMultiplicity::kVZEROAEq) profilebasename="VZEROAmult";
    else if(recoEstimator==AliAnalysisTaskSEDvsMultiplicity::kVZERO || recoEstimator==AliAnalysisTaskSEDvsMultiplicity::kVZEROEq) profilebasename="VZEROMmult";
    cout<<endl<<endl<<" profilebasename="<<profilebasename<<endl<<endl;
        
    if (isPPbData) {
      if(year == 16) {
	const Char_t* periodNames[4] = {"LHC16q_265499to265525_265309to265387", "LHC16q_265435","LHC16q_265388to265427","LHC16t_267163to267166"};
	TProfile* multEstimatorAvg[4];
	for(Int_t ip=0; ip<4; ip++) {
	  cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
	  multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	  if (!multEstimatorAvg[ip]) {
	    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	    return;
	  }
	}
	dMultTask->SetMultiplVsZProfileLHC16qt1stBunch(multEstimatorAvg[0]);
	dMultTask->SetMultiplVsZProfileLHC16qt2ndBunch(multEstimatorAvg[1]);
	dMultTask->SetMultiplVsZProfileLHC16qt3rdBunch(multEstimatorAvg[2]);
	dMultTask->SetMultiplVsZProfileLHC16qt4thBunch(multEstimatorAvg[3]);
      }
      else {
	//Only use two profiles if pPb 2013
	const Char_t* periodNames[2] = {"LHC13b", "LHC13c"};
	TProfile* multEstimatorAvg[2];
	for(Int_t ip=0; ip<2; ip++) {
	  cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
	  multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	  if (!multEstimatorAvg[ip]) {
	    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	    return;
	  }
	}
	dMultTask->SetMultiplVsZProfileLHC13b(multEstimatorAvg[0]);
	dMultTask->SetMultiplVsZProfileLHC13c(multEstimatorAvg[1]);
      }
    }
    else {
      const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
      TProfile* multEstimatorAvg[4];
      for(Int_t ip=0; ip<4; ip++) {
	multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	if (!multEstimatorAvg[ip]) {
	  AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
	  return;
	}
      }
      dMultTask->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
      dMultTask->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
      dMultTask->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
      dMultTask->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
    }
  }
  mgr->AddTask(dMultTask);
    
  // Create containers for input/output
    
  TString inname = "cinput";
  TString outname = "coutput";
  TString cutsname = "coutputCuts";
  TString normname = "coutputNorm";
  TString profname = "coutputProf";
    
  inname += Name.Data();
  outname += Name.Data();
  cutsname += Name.Data();
  normname += Name.Data();
  profname += Name.Data();
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  profname += finDirname.Data();
    
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DMult_";
  outputfile += Name.Data(); 
  outputfile += finDirname.Data(); 
    
  AliAnalysisDataContainer *coutputCuts = mgr->CreateContainer(cutsname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(normname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputProf = mgr->CreateContainer(profname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    
  mgr->ConnectInput(dMultTask,0,mgr->GetCommonInputContainer());
    
  mgr->ConnectOutput(dMultTask,1,coutput);
    
  mgr->ConnectOutput(dMultTask,2,coutputCuts);
    
  mgr->ConnectOutput(dMultTask,3,coutputNorm);  
    
  mgr->ConnectOutput(dMultTask,4,coutputProf);
    
  return dMultTask;
}
