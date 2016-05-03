// $Id$

AliAnalysisTaskGammaHadron* AddTaskGammaHadron(
  Bool_t      InputGammaOrPi0        = 0,      //gamma analysis=0, pi0 analyis=1
  Bool_t      InputSameEventAnalysis = 1,      //same event=1 mixed event =0 (currently only used to throw out event in  Run() function)
  const char *tracksName             = "Tracks",
  const char *clustersName           = "CaloClusters",
  const char *ncells                 = "EMCALCells",   //probably delete this this is nowhere used
  Double_t    trackptcut             = 0.15,
  Double_t    clusptcut              = 0.30,
  const char *taskname               = "AliAnalysisTask",
  const char *suffix                 = ""
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
  // Built the name of the Task together
  //-------------------------------------------------------
  TString GammaPi0Name;
  if(InputGammaOrPi0 == 0)
  {
	  GammaPi0Name += "GH";
  }
  else
  {
	  GammaPi0Name += "Pi0H";
  }
  TString SameMixName;
  if(InputSameEventAnalysis == 1)
  {
	  SameMixName += "SE";
  }
  else
  {
	  SameMixName += "ME";
  }

  TString combinedName;
  combinedName.Form("%s_%s_%s_%s_%s",taskname,(const char*)GammaPi0Name,(const char*)SameMixName,tracksName,clustersName);
  if(suffix!="")
  {
	  combinedName += "_";
	  combinedName += suffix;
  }

  cout<<"combinedName: "<<combinedName<<endl;
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  AliAnalysisTaskGammaHadron* AnalysisTask = new AliAnalysisTaskGammaHadron(InputGammaOrPi0,InputSameEventAnalysis);

  //perform some settings
  AnalysisTask->SetCaloCellsName(ncells);
  AnalysisTask->SetVzRange(-10,10);

  /*
  // Used for physics selection
  task->SetUseAliAnaUtils(kTRUE);
  task->DoVertexRCut(doVertexRCut);
  */

 //for later AnalysisTask->SetEffHistGamma(THnF *h);
 //for later AnalysisTask->SetEffHistHadron(THnF *h);


  //add particle and cluster container
  AliParticleContainer *partCont = AnalysisTask->AddParticleContainer(tracksName);
  AliClusterContainer  *clusCont = AnalysisTask->AddClusterContainer(clustersName);
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

  TString contName(combinedName);
  contName += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (AnalysisTask, 0,  cinput1 );
  mgr->ConnectOutput (AnalysisTask, 1, coutput1 );

  return AnalysisTask;
}
