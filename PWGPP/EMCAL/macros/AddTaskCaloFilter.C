/// \file AddTaskCaloFilter.C
/// \ingroup EMCALPerfAddTaskMacros
/// \brief Configuration analysis task filtering events and calorimeter data into AOD format.
///
/// The configuration is EMCal oriented, if PHOS oriented modifications are needed.
/// Events filtered with at least EMCal or PHOS clusters/cells are stored in a new AOD file.
/// Optionally tracks and other event info can be stored.
/// Event selection depending on EMCal/PHOS/TPC info can be set.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>
#include <TROOT.h>

#include "AliAnalysisTaskCaloFilter.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliEMCALRecoUtils.h"

#endif

///
/// Main method calling all the configuration
///
/// The parameters passed to the analysis are:
/// \param bias      : bias or not depending on EMCal signal.
/// \param mc        : simulation or data.
/// \param minE      : at least one cluster in EMCal with this energy.
/// \param minN      : at least one cluster in EMCal with this number of cells.
/// \param vz        : z vertex cut.
/// \param opt       : filter EMCal only, PHOS only or both.
/// \param correct   : calibrate the EMCal, remove bad cells ...
/// \param fillTrack : Fill event with hybrid tracks
/// \param fillAOD   : Output AOD is filled, not used only for a secondary task analysis at the same time.
///
AliAnalysisTaskCaloFilter * AddTaskCaloFilter
(  const Bool_t  bias      = kTRUE
 , const Bool_t  mc        = kFALSE
 , const Float_t minE      = 6.
 , const Int_t   minN      = 3 
 , const Float_t vz        = 10.
 , const Int_t   opt       = 0
 //AliAnalysisTaskCaloFilter::kBoth,kPHOS, kEMCAL or kBoth
 , const Bool_t  correct   = kFALSE
 , const Bool_t  fillTrack = kFALSE
 , const Bool_t  fillAOD   = kTRUE)
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
  
  // Set the output AOD handler
  //
  printf("AddTaskCaloFilter --- Init output handler ---\n");
  
  AliAODHandler* aodoutHandler   = new AliAODHandler();
  aodoutHandler->SetOutputFileName("AliAOD.EMCAL.root");
  //aodoutHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodoutHandler);

  // Configure the task
  //
  printf("AddTaskCaloFilter --- Init task ---\n");

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
    printf("AddTaskCaloFilter --- Select MC events with bias in EMCal ---\n");
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
    filter->SetMBTriggerMask(AliVEvent::kINT7); // not working for all productions
    
    filter->SelectCollisionCandidates(AliVEvent::kAny) ;
    
    printf("AddTaskCaloFilter --- Select events with bias in EMCal ---\n");
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
    
    filter->SelectCollisionCandidates(AliVEvent::kINT7);// | AliVEvent::kCentral | AliVEvent::kSemiCentral ) ;
    
    printf("AddTaskCaloFilter --- Select INT7 events ---\n");
  }

  // Activate the cluster corrections (calibration, bad map...)     
  //
  if(correct)   filter->SwitchOnClusterCorrection();
  else          filter->SwitchOffClusterCorrection();  
  
  // Exoticity cut settings
  //
  AliEMCALRecoUtils * reco = filter->GetEMCALRecoUtils();
  reco->SwitchOnRejectExoticCluster() ;
  reco->SetExoticCellFractionCut(0.97);
  reco->SetExoticCellMinAmplitudeCut(4.);

  // Track storing
  //
  if(fillTrack) { filter->SwitchOnFillTracks()  ; filter->SwitchOnFillHybridTracks()  ; }
  else          { filter->SwitchOffFillTracks() ; filter->SwitchOffFillHybridTracks() ; }
  
  // Other options to store in event
  //
  filter->SwitchOffFillv0s() ; // Put ON if you know what you do.
  
  filter->SwitchOnFillVZERO(); // Be able to recalculate centrality and event plane 
                               // afterwards even it is stored in header
  
  // AOD output storing
  //
  if(fillAOD)   filter->SwitchOnFillAODFile();
  else          filter->SwitchOffFillAODFile();
  
  // Pass the task to the manager, print first set parameters
  filter->PrintInfo();
  
  mgr->AddTask(filter);
  
  // Create containers for input/output
  //
  printf("AddTaskCaloFilter --- Created input/output containers ---\n");

  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
  printf("AddTaskCaloFilter --- Created containers, add them ---\n");
  
  mgr->ConnectInput  (filter, 0, cinput1  );
  mgr->ConnectOutput (filter, 0, coutput1 );
  
  printf("AddTaskCaloFilter --- End ---\n");

  return filter;

}

