///\file AddTaskEMCALPhotonIsolation.C
///\brief Configuration of the PhotonIsolation Task
///
/// \author Lucile Ronflette <lucile.ronflette@cern.ch>, SUBATECH, Nantes
/// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
/// \author Marco Marquard <marco.marquard@cern.ch>, University Frankfurt am Main



AliAnalysisTaskEMCALPhotonIsolation* AddTaskEMCALPhotonIsolation(
 const char *ntracks            = "EmcalTracks",
 const char *nclusters          = "EmcalClusters",
 const char *nhybtracks         = "EmcalTracksAna",
Bool_t		bHisto		= kTRUE,
 Int_t		iOutput		= 0,
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
    ::Error("AddTaskEMCALPhotonIsolation", "No analysis manager to connect to.");
    return NULL;
  }






  printf("Creating container names for cluster analysis\n");
    TString myContName("");
  if(bIsMC)
    myContName = Form("Analysis_Neutrals_MC");
  else
    myContName = Form("Analysis_Neutrals");

    // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);

  // #### Task preferences
  task->SetOutputFormat(iOutput);
 task->SetLCAnalysis(kFALSE);
 task->SetIsoConeRadius(0.4);
 task->SetEtIsoThreshold(2.);
task->SetCTMdeltaEta(0.02);
task->SetCTMdeltaPhi(0.03);
  task->SetQA(kTRUE);
  task->SetIsoMethod(1);
  task->SetEtIsoMethod(0);
  task->SetUEMethod(1);
  task->SetUSEofTPC(kFALSE);
   task->SetMC(bIsMC);


  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliParticleContainer *clusterCont = task->AddParticleContainer(nclusters);
  //  AliParticleContainer *hybTrackCont = task->AddParticleContainer(nhybtracks);

  printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
  // #### Add analysis task
  manager->AddTask(task);


  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralCluster",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *cinput  = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);

  //if(isEMCalTrain)
    //    RequestMemory(task,200*1024);


  return task;
}
