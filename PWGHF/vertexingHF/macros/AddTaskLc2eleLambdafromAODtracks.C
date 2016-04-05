AliAnalysisTaskSELc2eleLambdafromAODtracks *AddTaskLc2eleLambdafromAODtracks(TString finname="",
								   Bool_t theMCon=kFALSE,
									 Bool_t ispp= kFALSE,
								   Bool_t writeVariableTree=kFALSE,
									 Bool_t domixing=kFALSE,
									 Bool_t reconstructPrimVert=kFALSE,
								   Bool_t writeEachVariableTree=kFALSE,
								   Bool_t writeMCVariableTree=kFALSE,
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
	if(ispp)
		task->SetUseCentralityV0M(kFALSE);
	else
		task->SetUseCentralityV0M(kTRUE);
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

	task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
	if(ispp){
		task->SetPoolCentBinLimits(cent_mult_bin_numbpp,cent_mult_binlimitspp);
		task->SetNumberOfEventsForMixing(1000);//pp
	}else{
		task->SetPoolCentBinLimits(cent_mult_bin_numbpPb,cent_mult_binlimitspPb);
		task->SetNumberOfEventsForMixing(30);//pPb
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
  AliAnalysisDataContainer *coutputLc11 = mgr->CreateContainer(Form("eleLambda_mcgenpairvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,11,coutputLc11);
  AliAnalysisDataContainer *coutputLc12 = mgr->CreateContainer(Form("eleLambda_correlationvariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
  mgr->ConnectOutput(task,12,coutputLc12);

  return task;

}
