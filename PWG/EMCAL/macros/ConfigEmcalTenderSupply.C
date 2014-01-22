// $Id$

AliEMCALTenderSupply* ConfigEmcalTenderSupply(
  Bool_t distBC         = kTRUE,   //distance to bad channel
  Bool_t recalibClus    = kTRUE,   //recalibrate cluster energy
  Bool_t recalcClusPos  = kTRUE,   //recalculate cluster position
  Bool_t nonLinearCorr  = kTRUE,   //apply non-linearity
  Bool_t remExotic      = kTRUE,   //remove exotic cells
  Bool_t fidRegion      = kFALSE,  //apply fiducial cuts
  Bool_t calibEnergy    = kTRUE,   //calibrate energy
  Bool_t calibTime      = kTRUE,   //calibrate timing
  Bool_t remBC          = kTRUE,   //remove bad channels
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected,
  Bool_t reclusterize   = kTRUE,   //reclusterize
  Float_t seedthresh    = 0.100,   //seed threshold
  Float_t cellthresh    = 0.050,   //cell threshold
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Bool_t trackMatch     = kTRUE,   //track matching
  Bool_t updateCellOnly = kFALSE,  //only change if you run your own clusterizer task
  Float_t timeMin       = 100e-9,  //minimum time of physical signal in a cell/digit (s)
  Float_t timeMax       = 900e-9,  //maximum time of physical signal in a cell/digit (s)
  Float_t timeCut       = 900e-9   //maximum time difference between the digits inside EMC cluster (s)
)
{
  AliEMCALTenderSupply *EMCALSupply = new AliEMCALTenderSupply("EMCALtender");  
  EMCALSupply->SetDebugLevel(2);

  AliEMCALRecParam *params = new AliEMCALRecParam();
  params->SetClusteringThreshold(seedthresh); 
  params->SetMinECut(cellthresh);             
  params->SetW0(4.5);
  if (reclusterize) {
    params->SetClusterizerFlag(clusterizer);
    if (clusterizer == AliEMCALRecParam::kClusterizerNxN)
      params->SetNxM(3,3);
  }
  params->SetTimeMin(timeMin);
  params->SetTimeMax(timeMax);
  params->SetTimeCut(timeCut);

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
