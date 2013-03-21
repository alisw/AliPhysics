// $Id$

AliEMCALTenderSupply* ConfigEmcalTenderSupply(
  Bool_t timeCut       = kFALSE,
  Bool_t distBC        = kTRUE, 
  Bool_t recalibClus   = kTRUE, 
  Bool_t recalcClusPos = kTRUE, 
  Bool_t nonLinearCorr = kTRUE, 
  Bool_t remExotic     = kTRUE,
  Bool_t fidRegion     = kFALSE,
  Bool_t calibEnergy   = kTRUE,
  Bool_t calibTime     = kTRUE,
  Bool_t remBC         = kTRUE,
  UInt_t nonLinFunct   = AliEMCALRecoUtils::kBeamTestCorrected,
  Bool_t reclusterize  = kFALSE,
  UInt_t clusterizer   = AliEMCALRecParam::kClusterizerNxN,
  Bool_t trackMatch    = kFALSE,
  Bool_t updateCellOnly= kFALSE)
{
  AliEMCALTenderSupply *EMCALSupply = new AliEMCALTenderSupply("EMCALtender");  
  EMCALSupply->SetDebugLevel(2);

  AliEMCALRecParam *params = new AliEMCALRecParam();
  params->SetClusteringThreshold(0.1); // 100 MeV
  params->SetMinECut(0.05);            // 50 MeV  
  params->SetW0(4.5);
  if (reclusterize) {
    params->SetClusterizerFlag(clusterizer);
    if (clusterizer == AliEMCALRecParam::kClusterizerNxN)
      params->SetNxM(3,3);
  }
  //if (timeCut) {
  // No time cut
  params->SetTimeCut(1e6);
  params->SetTimeMin(-1);
  params->SetTimeMax(1e6);
  //}
  EMCALSupply->SetRecParam(params);

  if (reclusterize) {
    EMCALSupply->SwitchOnReclustering();
  }
  else {
    EMCALSupply->SwitchOffReclustering();
  }

  if (remBC) {
    EMCALSupply->SwitchOnClusterBadChannelCheck();
    EMCALSupply->SwitchOnBadCellRemove();
  }
  else {
    EMCALSupply->SwitchOffClusterBadChannelCheck();
    EMCALSupply->SwitchOffBadCellRemove();
  }

  if (distBC) {
    EMCALSupply->SwitchOnRecalDistBadChannel();
  }
  else {
    EMCALSupply->SwitchOffRecalDistBadChannel();
  }
  
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
  
  if (remExotic) {
    EMCALSupply->SwitchOnExoticCellRemove();
    EMCALSupply->SwitchOnClusterExoticChannelCheck();
  }
  else {
    EMCALSupply->SwitchOffExoticCellRemove();
    EMCALSupply->SwitchOffClusterExoticChannelCheck();
  }
  
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

  if (trackMatch) {
    EMCALSupply->SetMass(0.139);
    EMCALSupply->SwitchOnCutEtaPhiSeparate();
    EMCALSupply->SetEtaCut(0.025);
    EMCALSupply->SetPhiCut(0.05);
    EMCALSupply->SwitchOnTrackMatch();
  }
  else {
    EMCALSupply->SwitchOffTrackMatch();
  }
  
  if (updateCellOnly) 
    EMCALSupply->SwitchOnUpdateCellOnly();
  else
    EMCALSupply->SwitchOffUpdateCellOnly();

  return EMCALSupply;
}
