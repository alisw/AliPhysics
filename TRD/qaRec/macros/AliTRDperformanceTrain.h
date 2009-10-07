#ifndef ALITRDPERFORMANCETRAIN_H
#define ALITRDPERFORMANCETRAIN_H

#define BIT(n)      (1 << (n))
#define SETBIT(n,i) ((n) |= BIT(i))
#define TSTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))
#define CLRBIT(n,i) ((n) &= ~BIT(i))

#define NTRDQATASKS 6
#define NTRDCALIBTASKS 7
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
  ,kPIDRefMakerLQ = 9
  ,kPIDRefMakerNN =10
  ,kClErrParam    =11
  ,kMultiplicity  =12
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
  ,"AliTRDpidRefMakerLQ"
  ,"AliTRDpidRefMakerNN"
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
  ,"LQR"
  ,"NNR"
  ,"CLRES"
  ,"MULT"
  ,"ALL"
};

#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TFileMerger.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TError.h"
#endif

#include <cstring>

//____________________________________________
Bool_t HasReadMCData(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOMC");
}

//____________________________________________
Bool_t HasReadFriendData(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOFR");
}

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
//   if(TSTBIT(fSteerTask, kCalibration)) SETBIT(fSteerTask, kCheckDET);
//   if(TSTBIT(fSteerTask, kMultiplicity)) SETBIT(fSteerTask, kEfficiency);
//   if(TSTBIT(fSteerTask, kEfficiencyMC)) SETBIT(fSteerTask, kEfficiency);
//   if(TSTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
//   if(TSTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);
//   if(TSTBIT(fSteerTask, kPIDRefMaker)) SETBIT(fSteerTask, kCheckPID);

  return fSteerTask;
}

//______________________________________________________
void mergeProd(const Char_t *mark, const Char_t *files)
{
  const Int_t kBatch = 20;

  TFileMerger *fFM = new TFileMerger(1);
  fFM->OutputFile(Form("%s/0_%s",  gSystem->ExpandPathName("$PWD"), mark));

  Int_t jbatch = 0, nbatch = 0;
  string filename;
  ifstream file(files);
  while(getline(file, filename)){
    if(Int_t(filename.find(mark)) < 0) continue;
    fFM->AddFile(filename.c_str()); nbatch++;
    if(nbatch==kBatch){
      //printf("MERGE BATCH %d [%d]\n", jbatch, nbatch);
      fFM->Merge(); jbatch++;
      new(fFM) TFileMerger(kTRUE);
      fFM->OutputFile(Form("%s/%d_%s",  gSystem->ExpandPathName("$PWD"), jbatch, mark));
      nbatch=0;
    }
  }
  if(nbatch){
    //printf("MERGE BATCH %d[%d]\n", jbatch, nbatch);
    fFM->Merge();
    jbatch++;
  }
  if(!jbatch){
    delete fFM;
    return;
  }

  new(fFM) TFileMerger(kTRUE);
  fFM->OutputFile(Form("%s/%s",  gSystem->ExpandPathName("$PWD"), mark));
  for(Int_t ib=jbatch; ib--;){
    fFM->AddFile(Form("%s/%d_%s",  gSystem->ExpandPathName("$PWD"), ib, mark));
    gSystem->Exec(Form("rm -f %s/%d_%s", gSystem->ExpandPathName("$PWD"), ib, mark));
  }
  fFM->Merge();
  delete fFM;
}

#endif

