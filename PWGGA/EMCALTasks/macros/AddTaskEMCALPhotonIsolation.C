AliAnalysisTaskEMCALPhotonIsolation* AddTaskEMCALPhotonIsolation(
 const char *ntracks            = "EmcalTracks",
  const char *nclusters          = "EmcalClusters",
Bool_t		bHisto		= kTRUE,
Int_t		iOutput		= 1,
Bool_t		bIsMC		= kFALSE
)
{

  printf("Preparing neutral cluster analysis\n");
/*  // #### Detect the demanded trigger with its readable name
  TString triggerName(Form("Trigger_%i", trigger));
  if (trigger == AliVEvent::kAnyINT)
    triggerName = "kAnyINT";
  else if (trigger == AliVEvent::kAny)
    triggerName = "kAny";
  else if(trigger == AliVEvent::kINT7)
    triggerName = "kINT7";
  else if(trigger == AliVEvent::kMB)
    triggerName = "kMB";
  else if(trigger == AliVEvent::kEMC7)
    triggerName = "kEMC7";
  else if(trigger == AliVEvent::kEMCEJE)
    triggerName = "kEMCEJE";
  else if(trigger == AliVEvent::kEMCEGA)
    triggerName = "kEMCEGA";
*/
  // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskNeutralCluster", "No analysis manager to connect to.");
    return NULL;
  }






  printf("Creating container names for cluster analysis\n");
    TString myContName("");
  if(bIsMC)
    myContName = Form("Analysis_Neutrals_MC");
  else
    myContName = Form("Analysis_Neutrals");
    //myContName = Form("AnalysisR0%2.0f_%s%s%s", jetRadius*100, triggerName.Data(), stringPtHard.Data(), containerNameSuffix.Data());

    // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);

  // #### Task preferences
  task->SetOutputFormat(iOutput);
 task->SetLCAnalysis(kFALSE);
  task->SetQA(kTRUE);
  task->SetIsoMethod(1);
  task->SetUEMethod(1);
  task->SetUSEofTPC(kFALSE);


  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliParticleContainer *clusterCont = task->AddParticleContainer(nclusters);

  printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
  // #### Add analysis task
  manager->AddTask(task);

  // AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:NeutralCluster", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralCluster",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *cinput  = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);

  /*if(isEMCalTrain)
    RequestMemory(task,200*1024);

  // #### Do some nasty piggybacking on demand
  if (externalMacro)
    gROOT->Macro(externalMacro);
*/
  return task;
}
