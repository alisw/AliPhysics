AliAnalysisTaskCaloFilter * AddTaskCaloFilter(const Bool_t  bias      = kTRUE, 
                                              const Float_t minE      = 6, 
                                              const Float_t minN      = 2, 
                                              const Float_t vz        = 10.,
                                              const Int_t   opt       = AliAnalysisTaskCaloFilter::kEMCAL, //kPHOS, kEMCAL or kBoth
                                              const Bool_t  correct   = kFALSE,
                                              const Bool_t  fillTrack = kTRUE,
                                              const Bool_t  fillAOD   = kTRUE)
{

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTaskCaloFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTaskCaloFilter", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("CaloFilter");
  
  filter->SetCaloFilter(opt); //kPHOS, kEMCAL or kBoth
  
  filter->SetVzCut(vz);
  
  if(bias) // select events with significant signal in EMCAL
  {
    filter->SetEnergyCut(minE);
    filter->SetNcellsCut(minN);
    
    filter->SelectCollisionCandidates(AliVEvent::kAny) ;
    
    printf("--- Select events with 1 cluster with at least %2.2f GeV and N = %d ---\n",minE,minN);
  }
  else // Do not bias the signal in EMCAL, select MB events 
  {
    filter->SetEnergyCut(0);
    filter->SetNcellsCut(0);    

    filter->SelectCollisionCandidates(AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral ) ;
    
    printf("--- Select Min Bias events ---\n");
  }

  if(correct)   filter->SwitchOnClusterCorrection();
  else          filter->SwitchOffClusterCorrection();  
  
  if(fillTrack) filter->SwitchOnFillTracks();
  else          filter->SwitchOffFillTracks();
  
  if(fillAOD)   filter->SwitchOnFillAODFile();
  else          filter->SwitchOffFillAODFile();
    
  filter->PrintInfo();
  
  mgr->AddTask(filter);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
  mgr->ConnectInput  (filter, 0, cinput1);
  mgr->ConnectOutput (filter, 0, coutput1 );
  
  return filter;

}