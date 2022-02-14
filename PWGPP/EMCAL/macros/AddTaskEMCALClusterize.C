/// \file AddTaskEMCALClusterize.C
/// \ingroup EMCALPerfAddTaskMacros
/// \brief Configuration of EMCal re-clusterization analysis task.
///
/// This task reclusterizes on the fly EMCal clusters, creates a new
/// branch with those clusters so that it can be used by another analysis
/// task accessing this cluster branch.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>

#include "AliAnalysisManager.h"

#include "AliAnalysisTaskEMCALClusterize.h"

//#include "ConfigureEMCALRecoUtils.C"

#endif // CINT

/// Main method
///
/// The parameters for the analysis are:
/// \param clusArrTit: char string with default name of new clusters container, not needed but leave for backward compatibility 
/// \param bFillAOD: Bool, keep the new clusters in output file.
/// \param bMC: Bool, simulation or data.
/// \param exotic: Bool, remove exotic clusters.
/// \param name: TString, name of clusterizer: V1, V2, V3 (faster V2), V1Unfold, NxN.
/// \param trigger: TString, name of triggered events to be analyzed.
/// \param tm: Bool, perform track matching recalculation.
/// \param minEcell: float, minimum cell energy entering the cluster.
/// \param minEseed: float, minimum cell energy seed of the cluster
/// \param maxDeltaT: float, maximum difference in time of cells in cluster, keep it rather open.
/// \param timeWindow: float, maximum/minimum time of the clusters/cells, after time recalibration.
/// \param minEUnf: minimum energy cut for unfolding (check what it does)
/// \param minFrac: minimum fraction of energy cut for unfolding (check what it does)
/// \param bRecalE: Bool, recalibrate EMCal energy
/// \param bBad: Bool, remove bad channels
/// \param bRecalT: Bool, recalibrate EMCal time
/// \param iNonLine: Int, correct cluster non linearity
/// \param minCen: Integer, minimum centrality, -1 no selection
/// \param maxCen: Integer, maximum centrality, -1 no selection
/// \param clusterEnergyCutEvent: Float, in case of event filtering, select events with at least one EMCal cluster with this energy
/// \param nRowDiff: Integer, number of rows for NxM clusterizer
/// \param nColDiff: Integer, number of collumns for NxM clusterizer
/// \param skipOrReject: Bool, for unfolding (check)
/// \param tCardMimic: option for TCard correlation emulation, MC only. Options: 0 - no emulation; 1 - just add energy to adjacent cells; >1 - add energy to adjacent cells and subtract added energy to reference cell
/// \param cellUpd: update cells list with cuts
///
AliAnalysisTaskEMCALClusterize* AddTaskEMCALClusterize(const char  * clusArrTit = "EMCAL_Clusters_New",
                                                       const Bool_t  bFillAOD   = kFALSE,                                                
                                                       const Int_t   bMC        = kFALSE,
                                                       const Bool_t  exotic     = kTRUE,
                                                       const TString name       = "V1Unfold", 
                                                       const TString trigger    = "", 
                                                       const Int_t   tm         = 1, 
                                                       const Int_t   minEcell   = 50,
                                                       const Int_t   minEseed   = 100,
                                                       const Int_t   maxDeltaT  = 250,
                                                       const Int_t   timeWindow = 1000,
                                                       const Int_t   minEUnf    = 15, 
                                                       const Int_t   minFrac    = 1,
                                                       const Bool_t  bRecalE    = kTRUE,
                                                       const Bool_t  bBad       = kTRUE,
                                                       const Bool_t  bRecalT    = kTRUE,
                                                       const Int_t   iNonLine   = 0,
                                                       const Int_t   minCen     = -1,
                                                       const Int_t   maxCen     = -1,
                                                       const Float_t clusterEnergyCutEvent = -1,
                                                       const Int_t   nRowDiff   = 1,
                                                       const Int_t   nColDiff   = 1,
                                                       const Bool_t  skipOrReject = kFALSE,
                                                       const Int_t   tCardMimic = 0,
                                                       const Bool_t  cellUpd    = kTRUE
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEMCALClusterize", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEMCALClusterize", "This clusterize requires an input event handler");
    return NULL;
  }
  
  printf("AddTaskEMCALClusterize() - Passed Settings :\n");
  printf("\t mc %d, exo %d, name %s, trigger %s, tm %d\n",
         bMC,exotic,name.Data(),trigger.Data(),tm);
  printf("\t Ecell %d, Eseed %d, dT %d, wT %d, minUnf %d, minFrac %d \n",
         minEcell, minEseed,maxDeltaT,timeWindow,minEUnf,minFrac);
  printf("\t recalE %d, bad %d, recalT %d, nonlin %d, minCen %d, maxCen %d, rowDiff %d, colDiff %d, t-card %d, cell update %d \n",
         bRecalE,bBad,bRecalT,iNonLine,minCen,maxCen,nRowDiff,nColDiff,tCardMimic,cellUpd);

  // Create name of task and AOD branch depending on different settings
  TString arrayName = clusArrTit;
  if(name.Contains("NxN")) arrayName = Form("%dx%d_Ecell%d_Eseed%d_DT%d_WT%d",2*nRowDiff+1,2*nColDiff+1,minEcell,minEseed,maxDeltaT,timeWindow);
  else                     arrayName = Form(   "%s_Ecell%d_Eseed%d_DT%d_WT%d",              name.Data(),minEcell,minEseed,maxDeltaT,timeWindow);
  
  if(minCen != -1 && maxCen != -1)
    arrayName+=Form("_Cen%d_%d",minCen,maxCen);

  printf("AddTaskEMCALClusterize() - Created Branch Name: %s \n",arrayName.Data());
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskEMCALClusterize* clusterize = new AliAnalysisTaskEMCALClusterize(Form("EMCALClusterize%s",arrayName.Data()));

  clusterize->SetAODBranchName(arrayName);
  
  //clusterize->SetOCDBPath("local://$ALICE_ROOT/OCDB");

  // Centrality range
  clusterize->SetCentralityBin(minCen, maxCen); 
  
  // Some general settings to create AOD file in case we want to keep it
  clusterize->SwitchOffFillAODCaloCells();
  clusterize->SwitchOffFillAODHeader();
  clusterize->FillAODFile(bFillAOD); // fill aod.root with clusters?, not really needed for analysis.

  // Update cells list after cuts
  if ( cellUpd ) clusterize->SwitchOnUpdateCell();
  
  //-------------------------------------------------------
  // Set clusterization parameters via rec param
  //-------------------------------------------------------

  AliEMCALRecParam * params = clusterize->GetRecParam();

  //
  // Position and SS weight parameter
  //
  params->SetW0(4.5);

  //
  // Time cuts
  // Be careful using time cuts, best thing is to leave them open.
  //
  if(maxDeltaT > 1) params->SetTimeCut(maxDeltaT*1.e-9);
  else            { params->SetTimeCut(250*1.e-9); printf("AddTaskEMCALClusterize() - default maxDeltaT = 250 ns\n"); }// Same as in reco
  
  if(timeWindow > 1)
  {
    params->SetTimeMin(-1*timeWindow*1.e-9);
    params->SetTimeMax(timeWindow*1.e-9);
  }
  else
  {
    if(!bMC)
    {
      if(bRecalT)
      {
        params->SetTimeMin(-250*1.e-9);
        params->SetTimeMax( 250*1.e-9);
        printf("AddTaskEMCALClusterize() - default time window for calibrated time -250 ns < T < 250 ns\n");
      }
      else
      {
        // same as in reco, USE IF NO TIME RECALIBRATION
        params->SetTimeMin(425*1.e-9);
        params->SetTimeMax(825*1.e-9);
        printf("AddTaskEMCALClusterize() - default time window 425 ns < T < 825 ns\n");
      }
    }
    else
    {
      params->SetTimeMin(-100000000);
      params->SetTimeMax( 100000000);
      printf("AddTaskEMCALClusterize() - open time cut\n");
    }
  }

  //
  // Energy cuts
  //
  params->SetClusteringThreshold(minEseed/1.e3);
  params->SetMinECut            (minEcell/1.e3); 

  //
  // Clusterizer type
  //
  if(name.Contains("V3")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv3); // faster V2
  if(name.Contains("V2")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
  if(name.Contains("V1")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
  if(name.Contains("NxN"))
  {
    params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerNxN);
    printf("AddTaskEMCALClusterize() - Set NxN cluster size to %dx%d (row diff %d, col diff %d)\n",
           2*nRowDiff+1,2*nColDiff+1,nRowDiff,nColDiff);
    params->SetNxM(nRowDiff, nColDiff);
  }

  //-------------------------------------------------------
  // Unfolding, 2 options :
  //-------------------------------------------------------

  //    1) Just unfold existing clusters
  if(name.Contains("JustUnfold"))
    clusterize->JustUnfold(kTRUE); // if TRUE, do just unfolding, do not recluster cells
  else  
    clusterize->JustUnfold(kFALSE); 

  //   2) Unfold clusters created in the clusterize (revise settings)
  if (name.Contains("Unfold"))
  {
    clusterize->SwitchOnCellEnergySelection();
    clusterize->SetCellCuts(minEUnf/1000.,minFrac/10000.);
    clusterize->SetRejectBelowThreshold(skipOrReject);
    printf("AliAnalysisTaskEMCALClusterize() - Unfolding Cuts: min E %f, frac %f\n",
           minEUnf/1000.,minFrac/10000.);
    //clusterize->SwitchOffCellEnergySelection(); 
    
    if(!name.Contains("Just"))
      params->SetUnfold(kTRUE);
    else 
      params->SetUnfold(kFALSE);
    
  } // unfold
  
  //-------------------------------------------------------
  // Configure AliEMCALRecoUtils
  //-------------------------------------------------------

//  AliEMCALRecoUtils * reco = clusterize->GetRecoUtils();
//  
//  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
//  
  clusterize->ConfigureEMCALRecoUtils(bMC,exotic,iNonLine,bRecalE,bBad,bRecalT);
  
  //-------------------------------------------------------
  // Do track matching after clusterization
  //-------------------------------------------------------
  if ( tm > 0 ) 
  {
    clusterize->SwitchOnTrackMatching();
    if      ( tm == 2 ) clusterize->GetRecoUtils()->SwitchOnAODHybridTracksMatch();
    else if ( tm == 1 ) clusterize->GetRecoUtils()->SwitchOnAODTPCOnlyTracksMatch();
    else                clusterize->GetRecoUtils()->SetAODTrackFilterMask(tm);
  }
  else   clusterize->SwitchOffTrackMatching();

  //-------------------------------------------------------
  // Alignment matrices
  //-------------------------------------------------------

  //clusterize->SetImportGeometryFromFile(kTRUE,"$ALICE_PHYSICS/OADB/EMCAL/geometry_2011.root"); // done automatically, set here to force

  clusterize->SwitchOnLoadOwnGeometryMatrices();
  
  //-------------------------------------------------------
  // Clusterize events with some significant signal
  //-------------------------------------------------------
  
  if(clusterEnergyCutEvent > 0)
  {
    clusterize->SwitchOnSelectEMCALEvent();
    clusterize->SetEMCALEnergyCut(clusterEnergyCutEvent);
    clusterize->SetEMCALNcellsCut(3);
  }
  else clusterize->SwitchOffSelectEMCALEvent();
  
  //-------------------------------------------------------
  // Cluster MC labels recalculation
  //-------------------------------------------------------
  
  if(bMC)
  {
    //printf("Recalculate MC labels\n");
    clusterize->SwitchOnUseClusterMCLabelForCell(0) ; // Take the cell MC label as basis (only possible in recent productions, from 2012?)
    clusterize->SwitchOnRemapMCLabelForAODs()  ;      // Only in case 0, and for productions where the re-mapping of cell label in AODs was not done (productions before March 2013?)

    //clusterize->SwitchOnUseClusterMCLabelForCell(1) ; // Assign to each cell the same MC label as the original cluster to which it belonged
    //clusterize->SwitchOnUseClusterMCLabelForCell(2) ; // Find the original clusters that have the same cells as the new cluster,
                                                        // assign the labels of the original clusters to the new cluster.
                                                        // only interesting if output is V1
    
    // clusterize->SwitchOnUseMCEdepFracLabelForCell(); // For Run2 MC, switch all the above off
  
  }
  
  //-------------------------------------------------------
  // T-Card cell correlation
  // See https://alice-notes.web.cern.ch/node/837
  //-------------------------------------------------------
  
  clusterize->SwitchOffTCardCorrelation();

  if(bMC && tCardMimic > 0)
  {
    if ( tCardMimic == 1 ) clusterize->SwitchOnTCardCorrelation(kFALSE);
    else                   clusterize->SwitchOnTCardCorrelation(kTRUE);
        
    // Parameters setting    
    // Optional emulation, all SM have cross talk all the time
    for(Int_t ism = 0; ism < 22; ism++)
      clusterize->SetInducedEnergyLossProbabilityPerSM(1.0, ism);
    
    clusterize->SetInducedTCardMinimumCellEnergy(0.01) ;
    clusterize->SetInducedTCardMaximum(100) ;
    clusterize->SetInducedTCardMaximumLowE(-1);
    clusterize->SetInducedTCardMinimum(0.1);
    
    clusterize->SwitchOnRandomizeTCardInducedEnergy() ;
    
    clusterize->SetInducedEnergyLossMinimumFraction(0.35/100.);
    clusterize->SetInducedEnergyLossMaximumFraction(1.6/100.);
    
    // SM0,4,5,6,8,9,12,13,14,15,16,17,18,19 (set it first for all SM equal)
    Float_t mu1 = 0.80/100.;
    Float_t mu2 =-0.11/100.;
    Float_t wid = 0.50/100.;
    clusterize->SetInducedEnergyLossFraction  (mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1(mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossFractionWidth(wid,wid,wid,0.00000);
    
    // SM3,7
    mu1 = 1.20/100.;
    mu2 =-0.11/100.;
    clusterize->SetInducedEnergyLossFractionPerSM  (3,mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1PerSM(3,mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossMinimumFractionPerSM(0.6/100.,3);
    clusterize->SetInducedEnergyLossMaximumFractionPerSM(1.8/100.,3);
    
    mu1 = 1.20/100.;
    mu2 =-0.11/100.;
    clusterize->SetInducedEnergyLossFractionPerSM  (7,mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1PerSM(7,mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossMinimumFractionPerSM(0.6/100.,7);
    clusterize->SetInducedEnergyLossMaximumFractionPerSM(1.8/100.,7);
    
    // SM1,2,10,11
    mu1 = 1.20/100.;
    mu2 =-0.11/100.;
    clusterize->SetInducedEnergyLossFractionPerSM  (1,mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1PerSM(1,mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossMinimumFractionPerSM(0.5/100.,1);
 
    clusterize->SetInducedEnergyLossFractionPerSM  (10,mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1PerSM(10,mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossMinimumFractionPerSM(0.5/100.,10);
    
    clusterize->SetInducedEnergyLossFractionPerSM  (11,mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1PerSM(11,mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossMinimumFractionPerSM(0.5/100.,11);
    
    mu1 = 1.15/100.;
    mu2 =-0.11/100.;
    clusterize->SetInducedEnergyLossFractionPerSM  (2,mu1, mu1, mu1, 0.00000);
    clusterize->SetInducedEnergyLossFractionP1PerSM(2,mu2, mu2, mu2, 0.00000); 
    clusterize->SetInducedEnergyLossMinimumFractionPerSM(0.45/100.,2);
  }
  
  //-------------------------------------------------------
  // Trigger options
  //-------------------------------------------------------

  if(trigger=="EMC7")
  {
    printf("AddTaskEMCALClusterize() - trigger EMC7\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMC7);
  }
  else if (trigger=="INT7")
  {
    printf("AddTaskEMCALClusterize() - trigger INT7\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else if(trigger=="EMC1")
  {
    printf("AddTaskEMCALClusterize() - trigger EMC1\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(trigger=="MB")
  {
    printf("AddTaskEMCALClusterize() - trigger MB\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kMB);
  }  
  else if(trigger=="PHOS")
  {
    printf("AddTaskEMCALClusterize() - trigger PHOS\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kPHI7);
  }  
  else if(trigger=="PHOSPb")
  {
    printf("AddTaskEMCALClusterize() - trigger PHOSPb\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kPHOSPb);
  }
  else if(trigger=="AnyINT")
  {
    printf("AddTaskEMCALClusterize() - trigger AnyINT\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kAnyINT);
  }  
  else if(trigger=="INT")
  {
    printf("AddTaskEMCALClusterize() - trigger AnyINT\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kAny);
  }
  else if(trigger=="EMCEGA")
  {
    printf("AddTaskEMCALClusterize() - EMC Gamma\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  } 
  else if(trigger=="EMCEJE")
  {
    printf("AddTaskEMCALClusterize() - trigger EMC Jet\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  }
  else if(trigger=="Central")
  {
    printf("AddTaskEMCALClusterize() - trigger Central\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kCentral);
  } 
  else if(trigger=="SemiCentral")
  {
    printf("AddTaskEMCALClusterize() - trigger SemiCentral\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  }
  
  // Set clusters branch name, make sure the analysis after this one reads this name
    
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(clusterize);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (clusterize, 0,  cinput1 );
  
  if(bFillAOD)  
  {
    printf("AddTaskEMCALClusterize() - Fill output AOD\n");
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer() ;
    mgr->ConnectOutput (clusterize, 0, coutput1 );
  }
  return clusterize;
}
