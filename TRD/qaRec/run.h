#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

#define BIT(n)      (1 << (n))
#define SETBIT(n,i) ((n) |= BIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
#define CLRBIT(n,i) ((n) &= ~BIT(i))

#define NTRDTASKS 8

enum AliTRDrecoTasks{
   kInfoGen = 0
  ,kCheckDetector = 1
  ,kTrackingEfficiency = 2
  ,kTrackingCombinedEfficiency = 3
  ,kTrackingResolution = 4
  ,kCalibration = 5
  ,kPIDChecker = 6
  ,kPIDRefMaker = 7
  ,kClusterErrorParam = 8
};

const Char_t* fgkTRDtaskClassName[NTRDTASKS] = {
  "AliTRDcheckDetector"
  ,"AliTRDtrackingEfficiency"
  ,"AliTRDtrackingEfficiencyCombined"
  ,"AliTRDtrackingResolution"
  ,"AliTRDcalibration"
  ,"AliTRDpidChecker"
  ,"AliTRDpidRefMaker"
  ,"AliTRDclusterResolution"
};

const Char_t *fgkTRDtaskOpt[NTRDTASKS+3] = {
  "ALL"
  ,"DET"
  ,"EFF"
  ,"EFFC"
  ,"RES"
  ,"CAL"
  ,"PID"
  ,"PIDR"
  ,"CLRES"
  ,"NOFR"
  ,"NOMC"
};

#endif

