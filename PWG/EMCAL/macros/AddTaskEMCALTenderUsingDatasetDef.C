// $Id$

AliAnalysisTaskSE *AddTaskEMCALTenderUsingDatasetDef(
  const char *dstr    = "LHC11h",
  const char *pass    = 0 /*should not be needed*/
) 
{
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEMCALTender.C");

  Bool_t distBC         = kTRUE;   //distance to bad channel
  Bool_t recalibClus    = kTRUE;   //recalibrate cluster energy
  Bool_t recalcClusPos  = kTRUE;   //recalculate cluster position
  Bool_t nonLinearCorr  = kTRUE;   //apply non-linearity
  Bool_t remExotic      = kTRUE;   //remove exotic cells
  Bool_t fidRegion      = kFALSE;  //apply fiducial cuts
  Bool_t calibEnergy    = kTRUE;   //calibrate energy
  Bool_t calibTime      = kTRUE,   //calibrate timing
  Bool_t remBC          = kTRUE,   //remove bad channels
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected;
  Bool_t reclusterize   = kFALSE,  //reclusterize
  Float_t seedthresh    = 0.3,     //seed threshold
  Float_t cellthresh    = 0.05,    //cell threshold
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Bool_t trackMatch     = kFALSE,  //track matching
  Bool_t updateCellOnly = kFALSE,  //only change if you run your own clusterizer task
  Float_t timeMin       = 100e-9,  //minimum time of physical signal in a cell/digit (s)
  Float_t timeMax       = 900e-9,  //maximum time of physical signal in a cell/digit (s)
  Float_t timeCut       = 50e-9,   //maximum time difference between the digits inside EMC cluster (s)


AliAnalysisTaskSE *task = AddTaskEMCALTender(
  "tothink",distBC,recalibClus,recalcClusPos,nonLinearCorr,remExotic,fidRegion,calibEnergy,calibTime,
  remBC, nonLinFunct, reclusterize, seedthresh, cellthresh, clusterizer, trackMatch, updateCellOnly,
  timeMin, timeMax, timeCut, 0);

  return et;
}
