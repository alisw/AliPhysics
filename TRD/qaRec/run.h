#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

#define BIT(n)        (1 << (n))
#define SETBIT(n,i)   ((n) |= BIT(i))
#define TESTBIT(n,i)  ((Bool_t)(((n) & BIT(i)) != 0))
#define CLEARBIT(n,i) ((n) &= ~BIT(i))

const Int_t fknTasks = 8;
enum AliTRDrecoTasks{
   kInfoGen = 0
  ,kCheckDetector = 1
  ,kTrackingEfficiency = 2
  ,kTrackingCombinedEfficiency = 3
  ,kTrackingResolution = 4
  ,kCalibration = 5
  ,kPIDChecker = 6
  ,kPIDRefMaker = 7
};

Char_t *fTaskClass[fknTasks] = {
  "AliTRDtrackInfoGen"
  ,"AliTRDcheckDetector"
  ,"AliTRDtrackingEfficiency"
  ,"AliTRDtrackingEfficiencyCombined"
  ,"AliTRDtrackingResolution"
  ,"AliTRDcalibration"
  ,"AliTRDpidChecker"
  ,"AliTRDpidRefMaker"
};

Char_t *fTaskOpt[fknTasks+2] = {
  "ALL"
  ,"DET"
  ,"EFF"
  ,"EFFC"
  ,"RES"
  ,"CAL"
  ,"PID"
  ,"PIDR"
  ,"NOFR"
  ,"NOMC"
};

#endif

