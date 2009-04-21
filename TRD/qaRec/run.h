#ifndef TRDRECONSTRUCTIONTRAIN_H
#define TRDRECONSTRUCTIONTRAIN_H

#define BIT(n)      (1 << (n))
#define SETBIT(n,i) ((n) |= BIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
#define CLRBIT(n,i) ((n) &= ~BIT(i))

#define NQATASKS 6
#define NCALIBTASKS 6
const Int_t NTRDTASKS = NQATASKS+NCALIBTASKS;

enum AliTRDrecoTasks{
   kCheckESD      = 0
  ,kInfoGen       = 1
  ,kCheckDetector = 2
  ,kEfficiency    = 3
  ,kResolution    = 4
  ,kPID           = 5
  ,kCalibration   = 6
  ,kEfficiencyMC  = 7
  ,kAlignment     = 8
  ,kPIDRefMaker   = 9
  ,kClErrParam    =10
  ,kMultiplicity  =11
};

const Char_t* fgkTRDtaskClassName[NQATASKS+NCALIBTASKS] = {
  "AliTRDcheckESD"
  ,"AliTRDinfoGen"
  ,"AliTRDcheckDetector"
  ,"AliTRDefficiency"
  ,"AliTRDresolution"
  ,"AliTRDpidChecker"
  ,"AliTRDcalibration"
  ,"AliTRDefficiencyMC"
  ,"AliTRDalignmentTask"
  ,"AliTRDpidRefMaker"
  ,"AliTRDclusterResolution"
  ,"AliTRDmultiplicity"
};

const Char_t *fgkTRDtaskOpt[NQATASKS+NCALIBTASKS+1] = {
  ""
  ,"GEN"
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

