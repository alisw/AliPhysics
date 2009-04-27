#ifndef ALITRDPERFORMANCETRAIN_H
#define ALITRDPERFORMANCETRAIN_H

#define BIT(n)      (1 << (n))
#define SETBIT(n,i) ((n) |= BIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
#define CLRBIT(n,i) ((n) &= ~BIT(i))

#define NTRDQATASKS 6
#define NTRDCALIBTASKS 6
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
  ,kPIDRefMaker   = 9
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
  ,"CLRES"
  ,"MULT"
  ,"ALL"
};

#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TError.h"
#endif

//____________________________________________
Int_t ParseOptions(Char_t *trd)
{
  Int_t fSteerTask = 1;
  TObjArray *tasksArray = TString(trd).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 0; itask < NTRDQATASKS; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 2; itask < NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask); SETBIT(fSteerTask, kInfoGen);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("run.C", Form("TRD task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TSTBIT(fSteerTask, kCalibration)) SETBIT(fSteerTask, kCheckDET);
  if(TSTBIT(fSteerTask, kMultiplicity)) SETBIT(fSteerTask, kEfficiency);
  if(TSTBIT(fSteerTask, kEfficiencyMC)) SETBIT(fSteerTask, kEfficiency);
  if(TSTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kPIDRefMaker)) SETBIT(fSteerTask, kCheckPID);

  return fSteerTask;
}


#endif

