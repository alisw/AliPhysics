/// \file AddTaskEMCALPi0Calibration.C
/// \brief Configuration of task AliAnalysisTaskEMCALPi0CalibSelection.
///
/// Configuration of task AliAnalysisTaskEMCALPi0CalibSelection, which fills invariant mass
/// histograms for each of the EMCal channels. It has to be executed in several iterations.
///
/// The parameters for the analysis are:
/// \param arrayName: TString name of new cluster branch.
/// \param bFillAOD: Bool, keep the new clusters in output file.
/// \param bMC: Bool, simulation or data.
/// \param exotic: Bool, remove exotic clusters.
/// \param name: TString, name of clusterizer: V1, V2, V1Unfold, NxN.
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
/// \param bNonLine: Bool, correct cluster non linearity
/// \param minCen: Integer, minimum centrality, -1 no selection
/// \param maxCen: Integer, maximum centrality, -1 no selection
/// \param clusterEnergyCutEvent: Float, in case of event filtering, select events with at least one EMCal cluster with this energy
/// \param nRowDiff: Integer, number of rows for NxM clusterizer
/// \param nColDiff: Integer, number of collumns for NxM clusterizer
/// \param skipOrReject: Bool, for unfolding (check)
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

AliAnalysisTaskEMCALPi0CalibSelection * AddTaskEMCALPi0Calibration(TString outputFile = "", // AnalysisResults.root
                                                                   TString trigger ="CEMC7",
                                                                   Bool_t recalE = kFALSE, 
                                                                   Bool_t recalT = kFALSE,
                                                                   Bool_t rmBad  = kFALSE,
                                                                   Bool_t nonlin = kTRUE,
                                                                   Bool_t simu   = kFALSE)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTaskEMCALTriggerQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTaskEMCALPi0Calibration", "This task requires an input event handler");
    return NULL;
  }
    
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  
  AliAnalysisTaskEMCALPi0CalibSelection * pi0calib = new AliAnalysisTaskEMCALPi0CalibSelection ("EMCALPi0Calibration");
  //pi0calib->SetDebugLevel(10); 
  //pi0calib->UseFilteredEventAsInput();
  pi0calib->SetClusterMinEnergy(0.3);
  pi0calib->SetClusterMaxEnergy(10.);
  pi0calib->SetClusterLambda0Cuts(0.1,0.5);
  pi0calib->SetAsymmetryCut(1.);
  pi0calib->SetClusterMinNCells(1);
  pi0calib->SetNCellsGroup(0);
  pi0calib->SwitchOnSameSM();
  pi0calib->SetPairDTimeCut(20);
  pi0calib->SetClusterMinTime(560);
  pi0calib->SetClusterMaxTime(610);

  pi0calib->SetTriggerName(trigger);
  
  // Cluster recalculation, Reco Utils configuration
  
  
  AliEMCALRecoUtils * reco = pi0calib->GetEMCALRecoUtils();
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  
  ConfigureEMCALRecoUtils(reco,
                          simu,                             
                          kTRUE, // exotic
                          nonlin,
                          recalE, 
                          rmBad,
                          recalT); 
  
  reco->SetNumberOfCellsFromEMCALBorder(0); // Do not remove clusters in borders!
  
  // recalibrate energy and do corrections because of Temperature corrections 
  pi0calib->SwitchOnClusterCorrection();
  reco->SwitchOnRecalibration();
  reco->SwitchOnRunDepCorrection();
  
  //reco->Print("");
  
  //---------------------
  // Geometry alignment
  //---------------------
  
  pi0calib->SetGeometryName("EMCAL_COMPLETE12SMV1_DCAL_8SM");

  pi0calib->SwitchOnLoadOwnGeometryMatrices();
  
  
  //---------------------
  // Recalibration
  //---------------------
  
  if(recalE)
  {
    TFile * f = new TFile("RecalibrationFactors.root","read");
    for(Int_t ism = 0; ism < 12; ism++)
    {
      TH2F * h = (TH2F*)f->Get("EMCALRecalFactors_SM0");
      reco->SetEMCALChannelRecalibrationFactors(0,h);
    }
  }
  
  pi0calib->PrintInfo();
  
  mgr->AddTask(pi0calib);
  
  if(outputFile.Length()==0) outputFile = AliAnalysisManager::GetCommonFileName(); 

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("Pi0Calibration_Trig%s",trigger.Data()), 
                                                           TList::Class(), AliAnalysisManager::kOutputContainer,  
                                                           outputFile.Data());
  
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("ParamsPi0Calibration_Trig%s",trigger.Data()), 
                                                             TList::Class(), AliAnalysisManager::kOutputContainer, 
                                                             "AnalysisParameters.root");
  
  mgr->AddTask(pi0calib);
                                                             
  mgr->ConnectInput  (pi0calib, 0, cinput1);
  mgr->ConnectOutput (pi0calib, 1, coutput);
  mgr->ConnectOutput (pi0calib, 2, cout_cuts);

  return pi0calib;
}
