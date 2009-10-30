#ifndef ALITRDPERFORMANCETRAIN_H
#define ALITRDPERFORMANCETRAIN_H

#define BITBIT(n)      (1 << (n))
#define SETBITT(n,i) ((n) |= BITBIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BITBIT(i)) != 0))
#define CLRBITT(n,i) ((n) &= ~BITBIT(i))

#define NTRDQATASKS 6
#define NTRDCALIBTASKS 6
/* #define NTRDCALIBTASKS 7 */
const Int_t NTRDTASKS = NTRDQATASKS+NTRDCALIBTASKS;

enum AliTRDrecoTasks{
   kCheckESD      = 0
  ,kInfoGen       = 1
  ,kCheckDET      = 2
  ,kEfficiency    = 3
  ,kResolution    = 4
  ,kCheckPID      = 5
  ,kCalibration   = 6
  ,kEfficiencyMC  = 7
  ,kAlignment     = 8
  ,kPIDRefMaker = 9
/*   ,kPIDRefMakerLQ = 9 */
/*   ,kPIDRefMakerNN =10 */
/*   ,kClErrParam    =11 */
/*   ,kMultiplicity  =12 */
  ,kClErrParam    =10
  ,kMultiplicity  =11
};

const Char_t* fgkTRDtaskClassName[NTRDTASKS] = {
  "AliTRDcheckESD"
  ,"AliTRDinfoGen"
  ,"AliTRDcheckDET"
  ,"AliTRDefficiency"
  ,"AliTRDresolution"
  ,"AliTRDcheckPID"
  ,"AliTRDcalibration"
  ,"AliTRDefficiencyMC"
  ,"AliTRDalignmentTask"
  ,"AliTRDpidRefMaker"
/*   ,"AliTRDpidRefMakerLQ" */
/*   ,"AliTRDpidRefMakerNN" */
  ,"AliTRDclusterResolution"
  ,"AliTRDmultiplicity"
};

const Char_t *fgkTRDtaskOpt[NTRDTASKS+1] = {
  ""
  ,"GEN"
  ,"DET"
  ,"EFF"
  ,"RES"
  ,"PID"
  ,"CAL"
  ,"EFFC"
  ,"ALGN"
  ,"PIDR"
/*   ,"LQR" */
/*   ,"NNR" */
  ,"CLRES"
  ,"MULT"
  ,"ALL"
};

#include <cstring>

//____________________________________________
Bool_t HasReadMCData(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOMC");
}

//____________________________________________
Bool_t HasReadFriendData(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOFR");
}

#endif

