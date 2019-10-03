AliAnalysisTaskSELc2eleLambdafromAODtracks *AddTaskLc2eleLambdafromAODtracks(TString finname="",
								   Bool_t theMCon=kFALSE,
									 Int_t iscoltype= 0,
								   Bool_t writeVariableTree=kFALSE,
									 Int_t domixing=0,
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
      AliFatal("Input file not found : check your cut object");
    }
  }

  AliRDHFCutsLctoeleLambdafromAODtracks* RDHFCutsLc2eleLambdaanal = new AliRDHFCutsLctoeleLambdafromAODtracks();
  if (stdcuts) RDHFCutsLc2eleLambdaanal->SetStandardCutsPP2010();
  else RDHFCutsLc2eleLambdaanal = (AliRDHFCutsLctoeleLambdafromAODtracks*)filecuts->Get("eleLambdaAnalysisCuts");
  RDHFCutsLc2eleLambdaanal->SetName("eleLambdaAnalysisCuts");
  RDHFCutsLc2eleLambdaanal->SetMinPtCandidate(-1.);
  RDHFCutsLc2eleLambdaanal->SetMaxPtCandidate(10000.);
  if (!RDHFCutsLc2eleLambdaanal) {
    cout << "Specific AliRDHFCutsLc2eleLambdaanal not found\n";
    return;
  }

  //CREATE THE TASK

  printf("CREATE TASK\n");
  AliAnalysisTaskSELc2eleLambdafromAODtracks *task = new AliAnalysisTaskSELc2eleLambdafromAODtracks("AliAnalysisTaskSELc2eleLambdafromAODtracks",RDHFCutsLc2eleLambdaanal,writeVariableTree);
  task->SetMC(theMCon);
	if(iscoltype==0){
		task->SetUseCentralityV0M(kFALSE);
		task->SetUseCentralitySPDTracklet(kFALSE);
	}else{
		task->SetUseCentralityV0M(kTRUE);
		task->SetUseEventPlane(14);
	}
  task->SetDebugLevel(1);
  task->SetReconstructPrimVert(reconstructPrimVert);
  task->SetWriteEachVariableTree(writeEachVariableTree);
  task->SetWriteMCVariableTree(writeMCVariableTree);

	if(iscoltype==0){
		if(domixing==0){
			task->SetEventMixingOff();
		}else{
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitspp[] = { 0,100};
			Int_t cent_mult_bin_numbpp = sizeof(cent_mult_binlimitspp)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbpp,cent_mult_binlimitspp);

			Int_t nrpbin = 1.;
			Double_t rpbinlimits[2] = {-3.2,3.2};
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//pp
		}
	}else if(iscoltype==1){
		if(domixing==0){
			task->SetEventMixingOff();
		}else{
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitspPb[] = { 0,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbpPb = sizeof(cent_mult_binlimitspPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbpPb,cent_mult_binlimitspPb);

			Int_t nrpbin = 1.;
			Double_t rpbinlimits[2] = {-3.2,3.2};
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//pPb
		}
	}else if(iscoltype==2){
		if(domixing==0){
			task->SetEventMixingOff();
		}else if(domixing==1){
			//Standard
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else if(domixing==2){
			//Depth x 1/2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(5);//PbPb
		}else if(domixing==3){
			//Depth x 2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(20);//PbPb
		}else if(domixing==4){
			//z binning x1/2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-6,-2,2,6,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else if(domixing==5){
			//z binning  x2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else if(domixing==6){
			//Cent x 2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,15.,20,25.,30,35.,40,45.,50,55.,60,65.,70,75.,80,85.,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else if(domixing==7){
			//Cent x 1/2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,5,10,30,50,70,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 8;
			Double_t rpbinlimits[9];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<9;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else if(domixing==8){
			//RP x 2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 16;
			Double_t rpbinlimits[17];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<17;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else if(domixing==9){
			//RP x 2
			task->SetEventMixingWithPools();
			Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
			Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
			task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

			Double_t cent_mult_binlimitsPbPb[] = { 0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
			Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
			task->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);

			Int_t nrpbin = 4;
			Double_t rpbinlimits[5];
			Double_t steprp = TMath::Pi()/(Double_t)nrpbin;
			for(Int_t ir=0;ir<5;ir++){
				rpbinlimits[ir] = steprp * (Double_t) ir;
			}
			task->SetPoolRPBinLimits(nrpbin,rpbinlimits);

			task->SetNumberOfEventsForMixing(10);//PbPb
		}else{
			cout<<"Such pool binning option is not supported"<<endl;
			exit(1);
		}
	}

  //multiplicity study
  if(iscoltype==0){
    if(estimatorFilename.EqualTo("") ) {
      printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
    } else{
      TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
      if(!fileEstimator)  {
        AliFatal("File with multiplicity estimator not found\n");
        return;
      }
      task->SetReferenceMultiplcity(9.26);
      const Char_t* profilebasename="SPDmult10";
      const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
      TProfile* multEstimatorAvg[4];                       
      for(Int_t ip=0; ip<4; ip++) {
        multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
          AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
          return;
        }
      }
      task->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
      task->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
      task->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
      task->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
    }
  }

  mgr->AddTask(task);

  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_Lc2eleLambda_";
  outputfile += nTour;

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("Lchist%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutputLc2 = mgr->CreateContainer(Form("Lc2eleLambdaCuts%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  mgr->ConnectOutput(task,2,coutputLc2);

  AliAnalysisDataContainer *coutputLc3 = mgr->CreateContainer(Form("eleLambdaHisto%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,3,coutputLc3);

  AliAnalysisDataContainer *coutputLc4 = mgr->CreateContainer(Form("eleLambdavariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,4,coutputLc4);
  AliAnalysisDataContainer *coutputLc5 = mgr->CreateContainer(Form("eleLambda_elevariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,5,coutputLc5);
  AliAnalysisDataContainer *coutputLc6 = mgr->CreateContainer(Form("eleLambda_v0variables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,6,coutputLc6);
  AliAnalysisDataContainer *coutputLc7 = mgr->CreateContainer(Form("eleLambda_mcvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,7,coutputLc7);
  AliAnalysisDataContainer *coutputLc8 = mgr->CreateContainer(Form("eleLambdaCounter%1d",nTour),AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //counter
  mgr->ConnectOutput(task,8,coutputLc8);
  AliAnalysisDataContainer *coutputLc9 = mgr->CreateContainer(Form("eleLambda_mcelevariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,9,coutputLc9);
  AliAnalysisDataContainer *coutputLc10 = mgr->CreateContainer(Form("eleLambda_mcv0variables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,10,coutputLc10);
  //AliAnalysisDataContainer *coutputLc11 = mgr->CreateContainer(Form("eleLambda_mcgenpairvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  //mgr->ConnectOutput(task,11,coutputLc11);
  AliAnalysisDataContainer *coutputLc11 = mgr->CreateContainer(Form("eleLambda_singlevariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,11,coutputLc11);
  AliAnalysisDataContainer *coutputLc12 = mgr->CreateContainer(Form("eleLambda_correlationvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,12,coutputLc12);

  return task;

}
