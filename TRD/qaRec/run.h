#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

#define BIT(n)      (1 << (n))
#define SETBIT(n,i) ((n) |= BIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
#define CLRBIT(n,i) ((n) &= ~BIT(i))

#define NQATASKS 6
#define NCALIBTASKS 5
const Int_t NTRDTASKS = NQATASKS+NCALIBTASKS;

enum AliTRDrecoTasks{
   kInfoGen = 0
  ,kCheckDetector = 1
  ,kTrackingEff = 2
  ,kTrackingEffMC = 3
  ,kResolution = 4
  ,kPIDChecker = 5
  ,kCalibration = 6
  ,kAlignment = 7
  ,kPIDRefMaker = 8
  ,kClErrParam = 9
  ,kMultiplicity = 10
};

const Char_t* fgkTRDtaskClassName[NQATASKS+NCALIBTASKS] = {
  "AliTRDtrackInfoGen"
  ,"AliTRDcheckDetector"
  ,"AliTRDtrackingEfficiency"
  ,"AliTRDtrackingEfficiencyCombined"
  ,"AliTRDresolution"
  ,"AliTRDpidChecker"
  ,"AliTRDcalibration"
  ,"AliTRDalignmentTask"
  ,"AliTRDpidRefMaker"
  ,"AliTRDclusterResolution"
  ,"AliTRDmultiplicity"
};

const Char_t *fgkTRDtaskOpt[NQATASKS+NCALIBTASKS+1] = {
  "GEN"
  ,"DET"
  ,"EFF"
  ,"EFFC"
  ,"RES"
  ,"PID"
  ,"CAL"
  ,"ALGN"
  ,"PIDR"
  ,"CLRES"
  ,"MULT"
  ,"ALL"
};

#define NTPCPERFORMANCE 6
#define NTPCCALIBRATION 0
const Int_t NTPCTASKS = NTPCPERFORMANCE+NTPCCALIBRATION;

Char_t *fgkTPCtaskClassName[NTPCTASKS] = {
  "AliPerformanceTask"
  ,"AliPerformanceEff"
  ,"AliPerformanceRes"
  ,"AliPerformanceTPC"
  ,"AliPerformanceDEdx"
  ,"AliPerformanceDCA"
};

const Char_t *fgkTPCtaskOpt[NTPCTASKS+1] = {
  "GEN"
  ,"EFF"
  ,"RES"
  ,"TPC"
  ,"DEDX"
  ,"DCA"
  ,"ALL"
};
#endif

