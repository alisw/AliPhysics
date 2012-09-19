// $Id$

AliEMCALTenderSupply* ConfigEmcalTenderSupply(
  Bool_t timeCut       = kTRUE,
  Bool_t distBC        = kTRUE, 
  Bool_t recalibClus   = kTRUE, 
  Bool_t recalcClusPos = kTRUE, 
  Bool_t nonLinearCorr = kTRUE, 
  Bool_t remExotic     = kTRUE,
  Bool_t fidRegion     = kFALSE,
  Bool_t calibEnergy   = kTRUE,
  Bool_t calibTime     = kTRUE,
  UInt_t nonLinFunct   = AliEMCALRecoUtils::kBeamTestCorrected)
{
  AliEMCALTenderSupply *EMCALSupply = new AliEMCALTenderSupply("EMCALtender");  
  EMCALSupply->SetDebugLevel(2);

  AliEMCALRecParam *params = new AliEMCALRecParam();
  params->SetClusteringThreshold(0.1); // 100 MeV
  params->SetMinECut(0.05);            //50 MeV  
  params->SetW0(4.5);
  if (timeCut) {
    params->SetTimeCut(1e6); //Open this cut for AODs
    params->SetTimeMin(-1);  //Open this cut for AODs
    params->SetTimeMax(1e6); //Open this cut for AODs
  }
  EMCALSupply->SetRecParam(params);

  if (distBC)
    EMCALSupply->SwitchOnRecalDistBadChannel();
  else
    EMCALSupply->SwitchOffRecalDistBadChannel();
  
  if (recalibClus) {
    EMCALSupply->SwitchOnReCalibrateCluster();
    EMCALSupply->SwitchOnUpdateCell();
  }
  else {
    EMCALSupply->SwitchOffReCalibrateCluster();
    EMCALSupply->SwitchOffUpdateCell();
  }
  
  if (recalcClusPos)
    EMCALSupply->SwitchOnRecalculateClusPos();
  else
    EMCALSupply->SwitchOffRecalculateClusPos();
  
  if (nonLinearCorr) {
    EMCALSupply->SetNonLinearityFunction(nonLinFunct);
    EMCALSupply->SwitchOnNonLinearityCorrection();
  }
  else {
    EMCALSupply->SwitchOffNonLinearityCorrection();
  }
  
  if (remExotic) 
      EMCALSupply->SwitchOnClusterExoticChannelCheck();
  else
    EMCALSupply->SwitchOffClusterExoticChannelCheck();
  
  if (fidRegion)
    EMCALSupply->SwitchOnCellFiducialRegion();
  else
    EMCALSupply->SwitchOffCellFiducialRegion();

  if (calibTime)
    EMCALSupply->SwitchOnCalibrateTime();
  else
    EMCALSupply->SwitchOffCalibrateTime();

  if (calibEnergy)
    EMCALSupply->SwitchOnCalibrateEnergy();
  else
    EMCALSupply->SwitchOffCalibrateEnergy();

  EMCALSupply->SetMass(0.139);
  EMCALSupply->SwitchOnCutEtaPhiSeparate();
  EMCALSupply->SetEtaCut(0.025);
  EMCALSupply->SetPhiCut(0.05);
  
  return EMCALSupply;
}
