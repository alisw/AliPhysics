#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

const Int_t fknTasks = 7;
enum AliTRDrecoTasks{
  kInfoGen = 0
  ,kTrackingEfficiency = 1
  ,kTrackingCombinedEfficiency = 2
  ,kTrackingResolution = 3
  ,kCalibration = 4
  ,kPIDChecker = 5
  ,kCheckDetector = 6
};


#endif

