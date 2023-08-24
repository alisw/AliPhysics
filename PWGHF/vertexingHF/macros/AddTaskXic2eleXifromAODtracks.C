AliAnalysisTaskSEXic2eleXifromAODtracks *AddTaskXic2eleXifromAODtracks(TString finname="",
								   Bool_t theMCon=kFALSE,
								   Int_t iscoltype= 0,
								   Bool_t writeVariableTree=kTRUE,
								   Bool_t domixing=kFALSE,
								   Bool_t reconstructPrimVert=kFALSE,
								   Bool_t writeEachVariableTree=kFALSE,
								   Bool_t writeMCVariableTree=kFALSE,
				                   TString estimatorFilename="",
								   Int_t nTour=0
								   )

{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLc2V0YW", "No analysis manager to connect to.");
    return NULL;
  }  

  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( finname.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
    filecuts=TFile::Open(finname.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      Printf("FATAL: Input file not found : check your cut object");
      return NULL;
    }
  }

  AliRDHFCutsXictoeleXifromAODtracks* RDHFCutsXic2eleXianal = new AliRDHFCutsXictoeleXifromAODtracks();
  if (stdcuts) RDHFCutsXic2eleXianal->SetStandardCutsPP2010();
  else RDHFCutsXic2eleXianal = (AliRDHFCutsXictoeleXifromAODtracks*)filecuts->Get("eleXiAnalysisCuts");
  RDHFCutsXic2eleXianal->SetName("eleXiAnalysisCuts");
  RDHFCutsXic2eleXianal->SetMinPtCandidate(-1.);
  RDHFCutsXic2eleXianal->SetMaxPtCandidate(10000.);
  if (!RDHFCutsXic2eleXianal) {
    cout << "Specific AliRDHFCutsXic2eleXianal not found\n";
    return NULL;
  }


  //CREATE THE TASK

  printf("CREATE TASK\n");
  AliAnalysisTaskSEXic2eleXifromAODtracks *task = new AliAnalysisTaskSEXic2eleXifromAODtracks("AliAnalysisTaskSEXic2eleXifromAODtracks",RDHFCutsXic2eleXianal,writeVariableTree);
  task->SetMC(theMCon);
	if(iscoltype==0){
		task->SetUseCentralityV0M(kFALSE);
		task->SetUseCentralitySPDTracklet(kFALSE);
	}else{
		task->SetUseCentralityV0M(kTRUE);
		task->SetUseEventPlane(4);
	}
  task->SetDebugLevel(1);
  task->SetReconstructPrimVert(reconstructPrimVert);
  task->SetWriteEachVariableTree(writeEachVariableTree);
  task->SetWriteMCVariableTree(writeMCVariableTree);
  if(domixing)
    task->SetEventMixingWithPools();
  else
    task->SetEventMixingOff();

	//PVz Binning for pool PP or PbPb	
	Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
	Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;

	//Cent Binning for pool	pPb 
	Double_t cent_mult_binlimitspPb[] = { 0,10,20,30,40,50,60,70,80,90,100};
	Int_t cent_mult_bin_numbpPb = sizeof(cent_mult_binlimitspPb)/sizeof(Double_t) - 1;
	Double_t cent_mult_binlimitspp[] = { 0,100};
	Int_t cent_mult_bin_numbpp = sizeof(cent_mult_binlimitspp)/sizeof(Double_t) - 1;
	Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
	Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;

	task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
	if(iscoltype==0){
		task->SetPoolCentBinLimits(cent_mult_bin_numbpp,cent_mult_binlimitspp);
		task->SetNumberOfEventsForMixing(10);//pp
  }else if(iscoltype==0){
		task->SetPoolCentBinLimits(cent_mult_bin_numbpPb,cent_mult_binlimitspPb);
		task->SetNumberOfEventsForMixing(10);//pPb
	}else{
		task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);
		task->SetNumberOfEventsForMixing(10);//PbPb
	}

  if(iscoltype==0 || iscoltype==1){
    Int_t nrpbin = 1.;
    Double_t rpbinlimits[2] = {-3.2,3.2};
    task->SetPoolRPBinLimits(nrpbin,rpbinlimits);
  }else{
    Int_t nrpbin = 8;
    Double_t rpbinlimits[9];
    Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
    for(Int_t ir=0;ir<9;ir++){
      rpbinlimits[ir] = -1.*TMath::Pi()/2. + steprp * (Double_t) ir;
    }
    task->SetPoolRPBinLimits(nrpbin,rpbinlimits);
  }

  //multiplicity study
  if(iscoltype==0){
    if(estimatorFilename.EqualTo("") ) {
      printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
    } else{
      TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
      if(!fileEstimator)  {
        Printf("FATAL: File with multiplicity estimator not found\n");
        return NULL;
      }
      task->SetReferenceMultiplcity(9.26);
      const Char_t* profilebasename="SPDmult10";
      const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
      TProfile* multEstimatorAvg[4];                       
      for(Int_t ip=0; ip<4; ip++) {
        multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
          Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
          return NULL;
        }
      }
      task->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
      task->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
      task->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
      task->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
    }
  }

	//=================== for 5 TeV reweighting==========
	TF1 * weightfit = new TF1("weightfit","expo");
	weightfit -> SetParameter(0,8.66109e-01);
	weightfit -> SetParameter(1,-2.59812e-01);
    task -> SetFunction(weightfit);
	//=================== for 13 TeV reweighting===========
	TF1 * weightfit13 = new TF1("weightfit13","expo");
	weightfit13 -> SetParameter(0,1.04603e+00);
	weightfit13 -> SetParameter(1,-2.46920e-01);
    task -> SetFunction13(weightfit13);


	//============== weight process for the acceptance phi distribution
	// =========== electron =================================
	  TF1 * AccWeight = new TF1("AccWeight","expo+pol4(2)");
      AccWeight -> SetParameter(0,7.25182e+00);
      AccWeight -> SetParameter(1,-1.59753e+01);
      AccWeight -> SetParameter(2,1.48939e-01);
      AccWeight -> SetParameter(3,-2.29063e-01);
      AccWeight -> SetParameter(4,1.78563e-01);
      AccWeight -> SetParameter(5,-5.70464e-02);
      AccWeight -> SetParameter(6,6.25614e-03);
      task -> SetFunctionElectron(AccWeight);
   // ============== positron ================================
      TF1 * AccWeightPositron = new TF1("AccWeightPositron","expo+pol4(2)");
      AccWeightPositron -> SetParameter(0,7.23466e+00);
      AccWeightPositron -> SetParameter(1,-1.59978e+01);
      AccWeightPositron -> SetParameter(2,1.25472e-01);
      AccWeightPositron -> SetParameter(3,-1.80599e-01);
      AccWeightPositron -> SetParameter(4,1.43474e-01);
      AccWeightPositron -> SetParameter(5,-4.66011e-02);
      AccWeightPositron -> SetParameter(6,5.17285e-03);
 
      task -> SetFunctionPositron(AccWeightPositron);


	mgr->AddTask(task);

  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_Xic2eleXi_";
  outputfile += nTour;

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("Xichist%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutputLc2 = mgr->CreateContainer(Form("Xic2eleXiCuts%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  mgr->ConnectOutput(task,2,coutputLc2);

  AliAnalysisDataContainer *coutputLc3 = mgr->CreateContainer(Form("eleXiHisto%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,3,coutputLc3);

  AliAnalysisDataContainer *coutputLc4 = mgr->CreateContainer(Form("eleXivariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,4,coutputLc4);
  AliAnalysisDataContainer *coutputLc5 = mgr->CreateContainer(Form("eleXi_elevariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,5,coutputLc5);
  AliAnalysisDataContainer *coutputLc6 = mgr->CreateContainer(Form("eleXi_cascvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,6,coutputLc6);
  AliAnalysisDataContainer *coutputLc7 = mgr->CreateContainer(Form("eleXi_mcvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,7,coutputLc7);
  AliAnalysisDataContainer *coutputLc8 = mgr->CreateContainer(Form("eleXiCounter%1d",nTour),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //counter
  mgr->ConnectOutput(task,8,coutputLc8);
  AliAnalysisDataContainer *coutputLc9 = mgr->CreateContainer(Form("eleXi_mcelevariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,9,coutputLc9);
  AliAnalysisDataContainer *coutputLc10 = mgr->CreateContainer(Form("eleXi_mccascvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,10,coutputLc10);
  //AliAnalysisDataContainer *coutputLc11 = mgr->CreateContainer(Form("eleXi_mcgenpairvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  //mgr->ConnectOutput(task,11,coutputLc11);
  AliAnalysisDataContainer *coutputLc11 = mgr->CreateContainer(Form("eleXi_singlevariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,11,coutputLc11);
  AliAnalysisDataContainer *coutputLc12 = mgr->CreateContainer(Form("eleXi_correlationvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,12,coutputLc12);

  return task;

}
