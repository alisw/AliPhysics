#ifndef ALITRDPERFORMANCETRAIN_H
#define ALITRDPERFORMANCETRAIN_H

#define BITBIT(n)      (1 << (n))
#define SETBITT(n,i) ((n) |= BITBIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BITBIT(i)) != 0))
#define CLRBITT(n,i) ((n) &= ~BITBIT(i))

#define NTRDQATASKS 6
#define NTRDCALIBTASKS 6
const Int_t NTRDTASKS = NTRDQATASKS+NTRDCALIBTASKS;

enum ETRDinfoGenOutSlots {
   kEventInfo     = 1
  ,kTracksBarrel  = 2
  ,kTracksSA      = 3
  ,kTracksKink    = 4
  ,kV0List        = 5
  ,kMonitor       = 6
  ,kNOutSlots     = 7
};

enum ETRDrecoTasks{
   kCheckESD      = 0
  ,kInfoGen       = 1
  ,kCheckDET      = 2
  ,kEfficiency    = 3
  ,kResolution    = 4
  ,kCheckPID      = 5
  ,kCalibration   = 6
  ,kEfficiencyMC  = 7
  ,kAlignment     = 8
  ,kPIDRefMaker   = 9
  ,kClErrParam    =10
  ,kMultiplicity  =11
};

Char_t const* fgkTRDtaskClassName[NTRDTASKS] = {
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
  ,"AliTRDclusterResolution"
  ,"AliTRDmultiplicity"
};

Char_t const* fgkTRDtaskOpt[NTRDTASKS+1] = {
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

