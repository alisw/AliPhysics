#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

#define BIT(n)      (1 << (n))
#define SETBIT(n,i) ((n) |= BIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
#define CLRBIT(n,i) ((n) &= ~BIT(i))

#define NTRDTASKS 7

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

Char_t* fgkTRDtaskClassName[NTRDTASKS] = {
  "AliTRDcheckDetector"
  ,"AliTRDtrackingEfficiency"
  ,"AliTRDtrackingEfficiencyCombined"
  ,"AliTRDtrackingResolution"
  ,"AliTRDpidChecker"
  ,"AliTRDpidRefMaker"
  ,"AliTRDcalibration"
};

Char_t *fgkTRDtaskOpt[NTRDTASKS+3] = {
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

