// $Id: AliAnalysisTaskEMCALPi0CalibSelection 56196 2012-05-03 20:53:20Z gconesab $

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
  pi0calib->SelectCollisionCandidates(); 
  //pi0calib->SetDebugLevel(10); 
  //pi0calib->UseFilteredEventAsInput();
  pi0calib->SetClusterMinEnergy(0.3);
  pi0calib->SetClusterMaxEnergy(10.);
  pi0calib->SetClusterLambda0Cuts(0.1,0.5);
  pi0calib->SetAsymmetryCut(1.);
  pi0calib->SetClusterMinNCells(1);
  pi0calib->SetNCellsGroup(0);
  pi0calib->SwitchOnSameSM();
  pi0calib->SetPairDTimeCut(40);
  pi0calib->SetClusterMinTime(560);
  pi0calib->SetClusterMaxTime(610);

  pi0calib->SetTriggerName(trigger);
  
  // Cluster recalculation, Reco Utils configuration
  
  pi0calib->SwitchOnClusterCorrection();

  
  AliEMCALRecoUtils * reco = pi0calib->GetEMCALRecoUtils();
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigureEMCALRecoUtils.C");
  
  ConfigureEMCALRecoUtils(reco,
                          simu,                             
                          kTRUE, // exotic
                          nonlin,
                          recalE, 
                          rmBad,
                          recalT);   
  
  reco->SetNumberOfCellsFromEMCALBorder(1);

  //reco->Print("");
  
  //cu->SwitchOnEMCALOADB();

  
  //---------------------
  // Geometry alignment
  //---------------------
  
  //pi0calib->SetGeometryName("EMCAL_COMPLETE12SMV1");
  pi0calib->SetGeometryName("EMCAL_COMPLETEV1");

  pi0calib->SwitchOffLoadOwnGeometryMatrices();
  
  
  //---------------------
  // Recalibration
  //---------------------
  
  if(recalE)
  {
    reco->SwitchOnRecalibration();
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
