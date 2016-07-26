// $Id$

AliAnalysisTaskGammaHadron* AddTaskGammaHadron(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *ncells             = "EMCALCells",
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  const char *taskname           = "AliAnalysisTaskGammaHadron"
)
{  
	cout<<"in AddTaskGammaHadron.C(...)"<<endl;
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskGammaHadron", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskGammaHadron", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  
  TString name(taskname);
  if(strcmp(ntracks,""))
  {
    name += "_";
    name += ntracks;
  }
  if(strcmp(nclusters,""))
  {
    name += "_";
    name += nclusters;
  }

  AliAnalysisTaskGammaHadron* AnalysisTask = new AliAnalysisTaskGammaHadron();

  //perform some settings
  AnalysisTask->SetCaloCellsName(ncells);
  AnalysisTask->SetVzRange(-10,10);

 //for later AnalysisTask->SetEffHistGamma(THnF *h);
 //for later AnalysisTask->SetEffHistHadron(THnF *h);


  //add particle and cluster container
  AliParticleContainer *partCont = AnalysisTask->AddParticleContainer(ntracks);
  AliClusterContainer  *clusCont = AnalysisTask->AddClusterContainer(nclusters);
  if(partCont)
  {
	  partCont->SetParticlePtCut(trackptcut);
  }
  if(clusCont)
  {
    clusCont->SetClusPtCut(clusptcut);
    AnalysisTask->SetNeedEmcalGeom(kTRUE);
  }
  else
  {
	  AnalysisTask->SetNeedEmcalGeom(kFALSE);
  }


  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(AnalysisTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  //cout<<"ELI: test1 in task"<<endl;

  TString contName(name);
  contName += "_histosgrams";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (AnalysisTask, 0,  cinput1 );
  mgr->ConnectOutput (AnalysisTask, 1, coutput1 );

  return AnalysisTask;
}
