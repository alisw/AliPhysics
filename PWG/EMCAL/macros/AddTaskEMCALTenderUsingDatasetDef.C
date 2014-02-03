// $Id$

AliAnalysisTaskSE *AddTaskEMCALTenderUsingDatasetDef(
  const char *perstr  = "LHC11h",
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
  Bool_t calibTime      = kTRUE;   //calibrate timing
  Bool_t remBC          = kTRUE;   //remove bad channels
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected;
  Bool_t reclusterize   = kFALSE;  //reclusterize
  Float_t seedthresh    = 0.100;   //seed threshold
  Float_t cellthresh    = 0.050;   //cell threshold
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2;
  Bool_t trackMatch     = kTRUE;   //track matching
  Bool_t updateCellOnly = kFALSE;  //only change if you run your own clusterizer task
  Float_t timeMin       = -50e6;   //minimum time of physical signal in a cell/digit (s)
  Float_t timeMax       =  50e6;   //maximum time of physical signal in a cell/digit (s)
  Float_t timeCut       =   1e6;   //maximum time difference between the digits inside EMC cluster (s)

  Bool_t isMC = kFALSE;

  TString period(perstr);
  period.ToLower();

  if (period == "lhc11h") {
    trackMatch = kFALSE;
    timeMin = 100e-9;
    timeMax = 900e-9;
    timeCut = 900e-9;
  }
  else if (period == "lhc12a15e" || period == "lhc12a15e_fix") {
    nonLinFunct = AliEMCALRecoUtils::kPi0MCv3;
    isMC = kTRUE;
  }
  else if (period == "lhc12a15a") {
    nonLinFunct = AliEMCALRecoUtils::kPi0MCv2;
    isMC = kTRUE;
  }
  else if(period == "lhc12a" || period == "lhc12b" || period == "lhc12c" || period == "lhc12d" || period == "lhc12e" || period == "lhc12f" || period == "lhc12g" || period == "lhc12h" || period == "lhc12i") {
    reclusterize = kTRUE;
    seedthresh = 0.3;
    cellthresh = 0.05;
    timeMin = 485e-9; //no timing calib available
    timeMax = 685e-9;
    timeCut =  1e6;
  }
  else if(period == "lhc13b" || period == "lhc13c" || period == "lhc13d" || period == "lhc13e" || period == "lhc13f" || period == "lhc13g" || 
	  period == "lhc13b4" || period == "lhc13b4_fix" || period == "lhc13b4_plus") {
    reclusterize = kTRUE;
    seedthresh = 0.3;
    cellthresh = 0.05;
    if(period == "lhc13b" || period == "lhc13c" || period == "lhc13d" || period == "lhc13e" || period == "lhc13f" || period == "lhc13g") {
      timeMin = -50e-9;
      timeMax =  50e-9;
      timeCut =  1e6;
    }
    else if(period == "lhc13b4" || period == "lhc13b4_fix" || period == "lhc13b4_plus") {
      nonLinFunct = AliEMCALRecoUtils::kMCPi0v3;
      isMC = kTRUE;
    }
  }

  if(isMC) { //no timing cuts when running on MC
    timeMin = -1;
    timeMax =  1e6;
    timeCut =  1e6;
  }

  AliAnalysisTaskSE *task = AddTaskEMCALTender(
    distBC, recalibClus, recalcClusPos, nonLinearCorr,
    remExotic, fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, 
    reclusterize, seedthresh, cellthresh, clusterizer, trackMatch, updateCellOnly,
    timeMin, timeMax, timeCut, pass);

  return task;
}
