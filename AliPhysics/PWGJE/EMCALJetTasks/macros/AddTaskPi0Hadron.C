// $Id$

AliAnalysisTaskPi0Hadron* AddTaskPi0Hadron(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *ncells             = "EMCALCells",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  const char *cutType            = "TPC",
  const char *taskname           = "AliAnalysisTaskPi0Hadron"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskPi0Hadron", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskPi0Hadron", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  
  TString name(taskname);
  if (strcmp(ntracks,"")) {
    name += "_";
    name += ntracks;
  }
  if (strcmp(nclusters,"")) {
    name += "_";
    name += nclusters;
  }
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  
  name += "_";
  name += cutType;


  //The next XX lines are implemented to get the reader and
  //AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();



//  AliAnalysisTaskPi0Hadron* qaTask = new AliAnalysisTaskPi0Hadron(name);
  AliAnalysisTaskPi0Hadron* qaTask = new AliAnalysisTaskPi0Hadron();
  qaTask->SetCaloCellsName(ncells);
  qaTask->SetRhoName(nrho,-1);
  qaTask->SetVzRange(-10,10);

  AliParticleContainer *partCont = qaTask->AddParticleContainer(ntracks);
  if (partCont) partCont->SetParticlePtCut(trackptcut);
  AliClusterContainer *clusCont = qaTask->AddClusterContainer(nclusters);
  if (clusCont) {
    clusCont->SetClusPtCut(clusptcut);
    qaTask->SetNeedEmcalGeom(kTRUE);
  }
  else {
    qaTask->SetNeedEmcalGeom(kFALSE);
  }


  AliJetContainer *jetCont = qaTask->AddJetContainer(njets,cutType,jetradius);
  if (jetCont) {
    jetCont->SetJetPtCut(jetptcut);
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetRhoName(nrho);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  cout<<"ELI: test1 in task"<<endl;

  mgr->AddTask(qaTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  cout<<"ELI: test1 in task"<<endl;

  TString contName(name);
  contName += "_histosgrams";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );

  return qaTask;
}
