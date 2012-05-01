// $Id$

AliTender *AddTaskEMCALTender(
  const char *geoname="EMCAL_COMPLETEV1", 
  const char* datatype="pp")
{
  // Parameters: geoname = "EMCAL_FIRSTYEARV1" or "EMCAL_COMPLETEV1" or ""

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTrgContam", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliTender* ana = new  AliTender("AliTender");

  ana->SelectCollisionCandidates( AliVEvent::kEMC1 | AliVEvent::kMB | AliVEvent::kEMC7 | AliVEvent::kINT7);

  
  Bool_t ismc = (mgr->GetMCtruthEventHandler() != NULL);

  
  mgr->AddTask(ana);
  // Adding EMCAL supply
  AliEMCALTenderSupply *EMCALSupply=new AliEMCALTenderSupply("EMCALtender");  
  EMCALSupply->SetDebugLevel(2);

  AliEMCALRecParam *params = new AliEMCALRecParam();
  params->SetClusteringThreshold(0.1); // 100 MeV
  params->SetMinECut(0.05); //50 MeV  
  params->SetW0(4.5);
  params->SetTimeCut(1e6);//Open this cut for AODs
  params->SetTimeMin(-1);//Open this cut for AODs
  params->SetTimeMax(1e6);//Open this cut for AODs
  EMCALSupply->SetRecParam(params);

  EMCALSupply->SetEMCALGeometryName(geoname);  

  EMCALSupply->SwitchOffCellFiducialRegion(); //Does NOT remove edge clusters
  if (!ismc) {
    if (1){//datatype == "pp") {
      //::Info("AddTaskEMCALTender", "USING pp data configuration...");
      //params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2); //Std for pbpb
      EMCALSupply->SwitchOnRecalDistBadChannel();
      EMCALSupply->SwitchOnReCalibrateCluster();
      EMCALSupply->SwitchOnRecalculateClusPos();
      //EMCALSupply->SetNonLinearityFunction(AliEMCALTenderSupply::kBeamTestCorrected);
      //EMCALSupply->SwitchOnUpdateCell(); // will update cells and reclusterize
      //EMCALSupply->SwitchOnReclustering(); //SwitchOnReclustering if needed      
    } else {
      //::Info("AddTaskEMCALTender", "USING PbPb data configuration...");
      //params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2); //Std for pbpb
      EMCALSupply->SwitchOnCalibrateTime(); // 
      EMCALSupply->SwitchOnUpdateCell(); // will update cells and reclusterize
      //EMCALSupply->SwitchOnReclustering(); //SwitchOnReclustering if needed      
    }
  } else {
    ::Info("AddTaskEMCALTender", "USING MC configuration...");
  }
  EMCALSupply->SetMass(0.139);
  //EMCALSupply->SetStep(5);
  //EMCALSupply->SwitchOnCutEtaPhiSum(); 
  //EMCALSupply->SetRCut(0.0025);

  EMCALSupply->SwitchOnCutEtaPhiSeparate();
  EMCALSupply->SetEtaCut(0.025);
  EMCALSupply->SetPhiCut(0.05);

  ana->AddSupply(EMCALSupply);
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("tender_event", 
                         AliESDEvent::Class(), 
                         AliAnalysisManager::kExchangeContainer,
                         "default_tender");
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
