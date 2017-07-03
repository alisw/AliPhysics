/// \file AddTaskCaloFilter.C
/// \ingroup EMCALPerfAddTaskMacros
/// \brief Configuration analysis task filtering events and calorimeter data into AOD format.
///
/// The parameters for the analysis are:
/// \param bias: bias or not depending on EMCal signal.
/// \param mc: simulation or data.
/// \param minE: at least one cluster in EMCal with this energy.
/// \param minN: at least one cluster in EMCal with this number of cells.
/// \param vz: z vertex cut.
/// \param opt: filter EMCal only, PHOS only or both.
/// \param correct: calibrate the EMCal, remove bad cells ...
/// \param fillTrack: Fill event with hybrid tracks
/// \param fillAOD: Output AOD is filled, not used only for a secondary task analysis at the same time.
///
/// The configuration is EMCal oriented, if PHOS oriented modifications are needed.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

AliAnalysisTaskCaloFilter * AddTaskCaloFilter(const Bool_t  bias      = kTRUE,
                                              const Bool_t  mc        = kFALSE,
                                              const Float_t minE      = 6, 
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
  
  if(mc)
  {
    filter->SetEventSelection(1,0,0); // Select events depending on EMCAL, PHOS and tracks criteria
    filter->SwitchOnAcceptAllMBEvent();
    
    filter->SwitchOnFillMCParticles();
    
    filter->SetEMCALEnergyCut(minE);
    filter->SetEMCALNcellsCut(minN);
    
    filter->SetPHOSEnergyCut(minE);
    filter->SetPHOSNcellsCut(minN);
    
    filter->SetTrackPtCut(minE);
  }
  else if(bias) // select events with significant signal in EMCAL or TPC or PHOS
  {
    filter->SetEventSelection(1,0,0); // Select events depending on EMCAL, PHOS and tracks criteria
    filter->SwitchOnAcceptAllMBEvent();
    
    filter->SetEMCALEnergyCut(minE);
    filter->SetEMCALNcellsCut(minN);
 
    filter->SetPHOSEnergyCut(minE);
    filter->SetPHOSNcellsCut(minN);
    
    filter->SetTrackPtCut(minE);

    //filter->SetMBTriggerMask(AliVEvent::kAnyINT);
    filter->SetMBTriggerMask(AliVEvent::kMB); // not working for all productions
    
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

  
  //    filter->SelectCollisionCandidates(AliVEvent::kAny) ;
      
  if(correct)   filter->SwitchOnClusterCorrection();
  else          filter->SwitchOffClusterCorrection();  
  
  AliEMCALRecoUtils * reco = filter->GetEMCALRecoUtils();
  reco->SwitchOnRejectExoticCluster() ;
  reco->SetExoticCellFractionCut(0.97);
  reco->SetExoticCellMinAmplitudeCut(4.);

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