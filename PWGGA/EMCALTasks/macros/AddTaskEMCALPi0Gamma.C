// Author: B. Sahlmueller
// Analysis task for neutral pions (into two gammas), and for direct photons by subtraction method

AliAnalysisTask *AddTaskEMCALPi0Gamma(const UInt_t triggermask = AliVEvent::kMB,
                                      Bool_t mcmode = 0, 
                                      Bool_t addsig = 0, 
                                      const char geoname[] = "EMCAL_COMPLETEV1",
                                      Bool_t qf = 0, 
				      Double_t asym1 = 0.3, 
                                      Double_t asym2 = 0.7, 
                                      Double_t minE = 0.3, 
                                      Double_t minecc = -10, 
                                      Int_t ncells = 2, 
                                      Bool_t bunfold = 0)
{

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPi0Gamma", "No analysis manager to connect to.");
    return NULL;
  }  

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPi0Gamma", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  if(bunfold){
    AliAnalysisTaskEMCALClusterize * clusterize = new AliAnalysisTaskEMCALClusterize("EMCALClusterize");
    clusterize->SetOCDBPath("local://$ALICE_ROOT/OCDB");
    //  clusterize->SetAODBranchName("newEMCALClusters"); //Make sure you use this same name to recover the list during the analysis.
    //  clusterize->FillAODFile(kTRUE); //if true aod.root file filled, in any case list of clusters is accessible after in the event analysis.
    clusterize->JustUnfold(kTRUE); // if TRUE, do just unfolding, do not recluster cells
    //  clusterize->SetImportGeometryFromFile(kTRUE,"$ALICE_ROOT/OADB/EMCAL/geometry_2011.root"); // change only in case 2010 to geometry_2010.root
    //  clusterize->SwitchOnLoadOwnGeometryMatrices();
    //////Set the parameters used in standard reconstruction
    
    AliEMCALRecParam * params = clusterize->GetRecParam();
    params->SetClusterizerFlag( AliEMCALRecParam ::kClusterizerv1); //Select the clusterization algorithm this or kClusterizerv1, kClusterizerNxN
    params->SetClusteringThreshold(0.3); // 100 MeV
    params->SetMinECut(0.05); //10 MeV , no effect
    params->SetUnfold(kTRUE); //Do the unfolding or not like in the standard reconstruction
    params->SetW0(4.5);
    mgr->AddTask(clusterize);
    
    mgr->ConnectInput (clusterize, 0, cinput);  //connect input
  }
  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskEMCALPi0Gamma* task = new  AliAnalysisTaskEMCALPi0Gamma("Pi0GammaTask");
  task->SelectCollisionCandidates(triggermask);
  task->SetUseQualFlag(qf);
  task->SetAsymMax1(asym1);
  task->SetAsymMax2(asym2);
  task->SetMinClusEnergy(minE);
  task->SetMinEcc(minecc);
  task->SetNminCells(ncells);
  task->SetGeoName(geoname);
  task->SetMcMode(mcmode);
  task->SetAddedSignal(addsig);
  task->SetFillNtuple(0);
  mgr->AddTask(task);
  
  //RequestMemory(task,320*1024); // request 0.5GB memory for task 

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosPi0Gamma",
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
							    //							    "_pi0gamma.root");
  
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, coutput1 );
  
  return task;
  
}
