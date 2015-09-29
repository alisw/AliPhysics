AliAnalysisTaskNeutralMesons *AddTaskEMCALNeutralMesons(TString trigger = "MC",
							TString geoname = "EMCAL_COMPLETEV1",
							Int_t mctype = 1,
							Int_t mcpart = 111,
							Double_t lowEnrgCut = 0.3,
							Bool_t ncellsInCl = "kFALSE",
							TString clusterBranch = "V1Unfold_Ecell50_Eseed280_DT250_WT1000"
							)
{

/// analysis manager

// Get the pointer to the existing analysis manager via the static access method.

  //trigger: MC, EMC, INT7, MB

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALNeutralMesons", "No analysis manager to connect to.");
    return NULL;
  }  
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALNeutralMesons", "This task requires an input event handler");
    return NULL;
  }
  
   AliAnalysisTaskNeutralMesons *task = new AliAnalysisTaskNeutralMesons("task");
   task->SetMCtype(mctype);
   task->SetMCprt(mcpart);
   task->SetGeoName(geoname);
   task->SetLowEnergyCut(lowEnrgCut);
   task->SetNcellsPerClusterCut(ncellsInCl);
   task->SetNewClBranchName(clusterBranch);
   
   if(trigger=="EMC7") task->SelectCollisionCandidates(AliVEvent::kEMC7);
   if(trigger=="INT7") task->SelectCollisionCandidates(AliVEvent::kINT7);
   if(trigger=="MB") task->SelectCollisionCandidates(AliVEvent::kMB);

   mgr->AddTask(task);

   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   mgr->ConnectInput(task,0,cinput);

   AliAnalysisDataContainer *co1 = mgr->CreateContainer(Form("out_%s_%d",trigger.Data(),mcpart),
                                                        TList::Class(),
                                                        AliAnalysisManager::kOutputContainer,
                                                        Form("%s:NeutralMesons",AliAnalysisManager::GetCommonFileName()));
   mgr->ConnectOutput(task,1,co1);

   return task; 
}
