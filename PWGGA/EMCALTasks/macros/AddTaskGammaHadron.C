// $Id$

AliAnalysisTaskGammaHadron* AddTaskGammaHadron(
  Bool_t      InputGammaOrPi0        = 0,      //gamma analysis=0, pi0 analyis=1
  Bool_t      InputSameEventAnalysis = 1,      //same event=1 mixed event =0 (currently only used to throw out event in  Run() function)
  const char *trackName              = "usedefault",
  const char *clusName               = "usedefault",
  const char *cellName               = "usedefault",   //probably delete this this is nowhere used
  Double_t    trackptcut             = 0.15,
  Double_t    clusptcut              = 0.30,
  const char *taskname               = "AliAnalysisTask",
  const char *suffix                 = ""
)
{  
  //cout<<"in AddTaskGammaHadron.C(...)"<<endl;
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskGammaHadron", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskGammaHadron", "This task requires an input event handler");
    return 0;
  }
  if (handler->InheritsFrom("AliESDInputHandler"))
  {
	  ::Error("AddTaskGammaHadron", "We have never taken care if this works for ESDs");
	  return 0;
   }

  //in case of AOD the default names are:
  if(trackName=="usedefault")trackName = "tracks";   //could be hybrid tracks
  if(clusName =="usedefault")clusName  = "caloClusters";
  if(cellName =="usedefault")cellName  = "emcalCells";

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
  combinedName.Form("%s_%s_%s_%s_%s",taskname,(const char*)GammaPi0Name,(const char*)SameMixName,trackName,clusName);
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

  //..Add the containers and set the names
  AnalysisTask->SetCaloCellsName(cellName);
  AnalysisTask->AddClusterContainer(clusName);

  if (trackName == "mcparticles")
  {
	  AliMCParticleContainer* mcpartCont = AnalysisTask->AddMCParticleContainer(trackName);
	  mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks")
  {
	  AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer(trackName);
	  trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
  }
  //..check that condition!! maybe for mixed events its different!!!!!!
  if(!AnalysisTask->GetTrackContainer(trackName) || !AnalysisTask->GetClusterContainer(clusName))
  {
	 cout<<"Task can not run like this!"<<endl;
	 return 0;
  }

  //..Add some selection criteria
  AnalysisTask->SetVzRange(-10,10);
  if(AnalysisTask->GetTrackContainer(trackName))
  {
	  AnalysisTask->GetTrackContainer(trackName)->SetParticlePtCut(trackptcut);
  }
  if(AnalysisTask->GetClusterContainer(clusName))
  {
	  AnalysisTask->GetClusterContainer(clusName)->SetClusPtCut(clusptcut);
  }

  /*
  // Used for physics selection
  task->SetUseAliAnaUtils(kTRUE);
  task->DoVertexRCut(doVertexRCut);
  */



  //..some additional input for the analysis
  AnalysisTask->SetNeedEmcalGeom(kTRUE);
  //for later AnalysisTask->SetEffHistGamma(THnF *h);
  //for later AnalysisTask->SetEffHistHadron(THnF *h);

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
