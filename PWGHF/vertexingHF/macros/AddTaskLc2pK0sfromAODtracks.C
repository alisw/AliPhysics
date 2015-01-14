AliAnalysisTaskSELc2pK0sfromAODtracks *AddTaskLc2pK0sfromAODtracks(TString finname="",
								   Bool_t theMCon=kFALSE,
								   Bool_t writeVariableTree=kFALSE,
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

  AliRDHFCutsLctopK0sfromAODtracks* RDHFCutsLc2pK0sprod = new AliRDHFCutsLctopK0sfromAODtracks();
  if (stdcuts) RDHFCutsLc2pK0sprod->SetStandardCutsPP2010();
  else RDHFCutsLc2pK0sprod = (AliRDHFCutsLctopK0sfromAODtracks*)filecuts->Get("LcProductionCuts");
  RDHFCutsLc2pK0sprod->SetName("LcProductionCuts");
  RDHFCutsLc2pK0sprod->SetMinPtCandidate(-1.);
  RDHFCutsLc2pK0sprod->SetMaxPtCandidate(10000.);
  if (!RDHFCutsLc2pK0sprod) {
    cout << "Specific AliRDHFCutsLc2pK0sprod not found\n";
    return;
  }
  
  AliRDHFCutsLctopK0sfromAODtracks* RDHFCutsLc2pK0sanal = new AliRDHFCutsLctopK0sfromAODtracks();
  if (stdcuts) RDHFCutsLc2pK0sanal->SetStandardCutsPP2010();
  else RDHFCutsLc2pK0sanal = (AliRDHFCutsLctopK0sfromAODtracks*)filecuts->Get("LcAnalysisCuts");
  RDHFCutsLc2pK0sanal->SetName("LcAnalysisCuts");
  RDHFCutsLc2pK0sanal->SetMinPtCandidate(-1.);
  RDHFCutsLc2pK0sanal->SetMaxPtCandidate(10000.);
  if (!RDHFCutsLc2pK0sanal) {
    cout << "Specific AliRDHFCutsLc2pK0sanal not found\n";
    return;
  }


  //CREATE THE TASK

  printf("CREATE TASK\n");
  AliAnalysisTaskSELc2pK0sfromAODtracks *task = new AliAnalysisTaskSELc2pK0sfromAODtracks("AliAnalysisTaskSELc2pK0sfromAODtracks",RDHFCutsLc2pK0sprod,RDHFCutsLc2pK0sanal,writeVariableTree);
  task->SetIspp(kFALSE);
  task->SetIspA(kFALSE);
  task->SetIsAA(kTRUE);
  task->SetMC(theMCon);
  task->SetDebugLevel(1);


  mgr->AddTask(task);

  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_Lc2pK0s_";
  outputfile += nTour;

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("Lchist%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutputLc2 = mgr->CreateContainer(Form("Lc2pK0sCuts%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // cuts
  mgr->ConnectOutput(task,2,coutputLc2);

  if (writeVariableTree) {
    AliAnalysisDataContainer *coutputLc3 = mgr->CreateContainer(Form("Lc2pK0svariables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
    mgr->ConnectOutput(task,3,coutputLc3);
  } else {
    AliAnalysisDataContainer *coutputLc3 = mgr->CreateContainer(Form("LcAll%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // variables tree
    mgr->ConnectOutput(task,3,coutputLc3);
  }

  return task;

}
