#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

#define BIT(n)        (1 << (n))
#define SETBIT(n,i)   ((n) |= BIT(i))
#define TESTBIT(n,i)  ((Bool_t)(((n) & BIT(i)) != 0))
#define CLEARBIT(n,i) ((n) &= ~BIT(i))

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

Char_t *fTaskClass[fknTasks] = {
  "AliTRDtrackInfoGen"
  ,"AliTRDtrackingEfficiency"
  ,"AliTRDtrackingEfficiencyCombined"
  ,"AliTRDtrackingResolution"
  ,"AliTRDcalibration"
  ,"AliTRDpidChecker"
  ,"AliTRDcheckDetector"
};

#endif

