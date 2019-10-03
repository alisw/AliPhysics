// $Id$

AliEMCALTenderSupply* ConfigEmcalTenderSupply(
  const Bool_t distBC             = kFALSE,     // distance to bad channel
  const Bool_t recalibClus        = kFALSE,     // recalibrate cluster energy
  const Bool_t recalcClusPos      = kFALSE,     // recalculate cluster position
  const Bool_t nonLinearCorr      = kFALSE,     // apply non-linearity
  const Bool_t remExoticCell      = kFALSE,     // remove exotic cells
  const Bool_t remExoticClus      = kTRUE,      // remove exotic clusters
  const Bool_t fidRegion          = kFALSE,     // apply fiducial cuts
  const Bool_t calibEnergy        = kFALSE,     // calibrate energy
  const Bool_t calibTime          = kFALSE,     // calibrate timing
  const Bool_t remBC              = kFALSE,     // remove bad channels
  const UInt_t nonLinFunct        = AliEMCALRecoUtils::kNoCorrection,
  const Bool_t reclusterize       = kTRUE,      // reclusterize
  const Float_t seedthresh        = 0.100,      // seed threshold
  const Float_t cellthresh        = 0.050,      // cell threshold
  const UInt_t clusterizer        = AliEMCALRecParam::kClusterizerv1,
  const Bool_t trackMatch         = kTRUE,      // track matching
  const Bool_t updateCellOnly     = kFALSE,     // only change if you run your own clusterizer task
  const Float_t timeMin           = -1,         // minimum time of physical signal in a cell/digit (s)
  const Float_t timeMax           = +1,         // maximum time of physical signal in a cell/digit (s)
  const Float_t timeCut           =  1,         // maximum time difference between the digits inside EMC cluster (s)
  const Float_t diffEAggregation  =  0.         // difference E in aggregation of cells (i.e. stop aggregation if E_{new} > E_{prev} + diffEAggregation)
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
      params->SetNxM(1,1); // -> (1,1) means 3x3!
  }
  params->SetTimeMin(timeMin);
  params->SetTimeMax(timeMax);
  params->SetTimeCut(timeCut);
  params->SetLocMaxCut(diffEAggregation); // Set minimum energy difference to start new cluster

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
    EMCALSupply->SwitchOnUpdateCell();
  } else {
    EMCALSupply->SwitchOffClusterBadChannelCheck();
    EMCALSupply->SwitchOffBadCellRemove();
  }

  if (distBC) {
    EMCALSupply->SwitchOnRecalDistBadChannel();
  } else {
    EMCALSupply->SwitchOffRecalDistBadChannel();
  }
  
  if (recalibClus) {
    EMCALSupply->SwitchOnReCalibrateCluster();
    EMCALSupply->SwitchOnUpdateCell();
  } else {
    EMCALSupply->SwitchOffReCalibrateCluster();
    //EMCALSupply->SwitchOffUpdateCell();
  }
  
  if (recalcClusPos)
    EMCALSupply->SwitchOnRecalculateClusPos();
  else
    EMCALSupply->SwitchOffRecalculateClusPos();
  
  if (nonLinearCorr) {
    EMCALSupply->SetNonLinearityFunction(nonLinFunct);
    EMCALSupply->SwitchOnNonLinearityCorrection();
  } else {
    EMCALSupply->SwitchOffNonLinearityCorrection();
  }
  
  if (remExoticCell) {
    EMCALSupply->SwitchOnExoticCellRemove();
  } else {
    EMCALSupply->SwitchOffExoticCellRemove();
  }

  if (remExoticClus) {
    EMCALSupply->SwitchOnClusterExoticChannelCheck();
  } else {
    EMCALSupply->SwitchOffClusterExoticChannelCheck();
  }
  
  if (fidRegion)
    EMCALSupply->SwitchOnCellFiducialRegion();
  else
    EMCALSupply->SwitchOffCellFiducialRegion();

  if (calibTime) {
    EMCALSupply->SwitchOnCalibrateTime();
    EMCALSupply->SwitchOnUpdateCell();
  }
  else
    EMCALSupply->SwitchOffCalibrateTime();

  if (calibEnergy) {
    EMCALSupply->SwitchOnCalibrateEnergy();
    EMCALSupply->SwitchOnUpdateCell();
  }
  else
    EMCALSupply->SwitchOffCalibrateEnergy();

  if (trackMatch) {
    EMCALSupply->SetMass(0.139);
    EMCALSupply->SwitchOnCutEtaPhiSeparate();
    EMCALSupply->SetEtaCut(0.025);
    EMCALSupply->SetPhiCut(0.05);
    EMCALSupply->SwitchOnTrackMatch();
  } else {
    EMCALSupply->SwitchOffTrackMatch();
  }
  
  if (updateCellOnly) 
    EMCALSupply->SwitchOnUpdateCellOnly();
  else
    EMCALSupply->SwitchOffUpdateCellOnly();

  return EMCALSupply;
}
