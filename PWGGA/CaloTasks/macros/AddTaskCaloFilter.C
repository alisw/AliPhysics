AliAnalysisTaskCaloFilter * AddTaskCaloFilter(const Bool_t  bias      = kTRUE, 
                                              const Float_t minE      = 5, 
                                              const Float_t minN      = 3, 
                                              const Float_t vz        = 10.,
                                              const Int_t   opt       = AliAnalysisTaskCaloFilter::kBoth, //kPHOS, kEMCAL or kBoth
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
  
  //filter->SetDebugLevel(2);
  
  filter->SetCaloFilter(opt); //kPHOS, kEMCAL or kBoth
  
  filter->SetVzCut(vz);
  
  if(bias) // select events with significant signal in EMCAL or TPC or PHOS
  {
    filter->SetEventSelection(1,1,1); // Select events depending on EMCAL, PHOS and tracks criteria
    filter->SwitchOnAcceptAllMBEvent();
    
    filter->SetEMCALEnergyCut(minE);
    filter->SetEMCALNcellsCut(minN);
 
    filter->SetPHOSEnergyCut(minE);
    filter->SetPHOSNcellsCut(minN);
    
    filter->SetTrackPtCut(minE);

    filter->SelectCollisionCandidates(AliVEvent::kAny) ;
    
    printf("--- Select events with 1 cluster with at least %2.2f GeV and N = %d ---\n",minE,minN);
  }
  else // Do not bias the signal in EMCAL, select MB events 
  {
    
    filter->SetEventSelection(0,0,0);
    filter->SwitchOnAcceptAllMBEvent();
    
    filter->SetEMCALEnergyCut(-1);
    filter->SetEMCALNcellsCut(0);  
    
    filter->SetPHOSEnergyCut(-1);
    filter->SetPHOSNcellsCut(0); 
    
    filter->SetTrackPtCut(-1);
    
    filter->SelectCollisionCandidates(AliVEvent::kMB);// | AliVEvent::kCentral | AliVEvent::kSemiCentral ) ;
    
    printf("--- Select Min Bias events ---\n");
  }

  if(correct)   filter->SwitchOnClusterCorrection();
  else          filter->SwitchOffClusterCorrection();  
  
  AliEMCALRecoUtils * reco = filter->GetEMCALRecoUtils();
  reco->SwitchOnRejectExoticCluster() ;
  reco->SetExoticCellFractionCut(0.97);
  reco->SetExoticCellMinAmplitudeCut(2.);

  if(fillTrack) { filter->SwitchOnFillTracks()  ; filter->SwitchOnFillHybridTracks()  ; }
  else          { filter->SwitchOffFillTracks() ; filter->SwitchOffFillHybridTracks() ; }
  
  filter->SwitchOffFillv0s() ; // Put ON if you know what you do.
  
  filter->SwitchOnFillVZERO(); // Be able to recalculate centrality and event plane afterwards even it is stored in header
  
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