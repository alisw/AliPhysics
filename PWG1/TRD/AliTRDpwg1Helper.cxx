#include <Rtypes.h>
#include <TError.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRandom.h>
#include <TString.h>
#include <TSystem.h>

#include <cstring>
#include <fstream>

#include "AliTRDpwg1Helper.h"

Char_t const* AliTRDpwg1Helper::fgkTRDtaskClassName[AliTRDpwg1Helper::kNTRDTASKS] = {
  "AliTRDcheckESD"
  ,"AliTRDinfoGen"
  ,"AliTRDcheckDET"
  ,"AliTRDefficiency"
  ,"AliTRDresolution"
  ,"AliTRDcheckPID"
  ,"AliTRDv0Monitor"
  ,"AliTRDcalibration"
  ,"AliTRDefficiencyMC"
  ,"AliTRDalignmentTask"
  ,"AliTRDpidRefMaker"
  ,"AliTRDclusterResolution"
  ,"AliTRDmultiplicity"
};

Char_t const* AliTRDpwg1Helper::fgkTRDtaskOpt[AliTRDpwg1Helper::kNTRDTASKS+1] = {
  "ESD"
  ,"GEN"
  ,"DET"
  ,"EFF"
  ,"RES"
  ,"PID"
  ,"V0"
  ,"CAL"
  ,"EFFC"
  ,"ALGN"
  ,"PIDR"
  ,"CLRES"
  ,"MULT"
  ,"ALL"
};

//______________________________________________________
Int_t AliTRDpwg1Helper::ParseOptions(Char_t *trd){
  Int_t fSteerTask = 0;
  TObjArray *tasksArray = TString(trd).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 0; itask < kNTRDQATASKS; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOMC") == 0 || s.CompareTo("NOFR") == 0){
      continue; // taken care by special functions
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 0; itask < kNTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask); 
        if(itask>1) SETBIT(fSteerTask, kInfoGen);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("ParseOptions()", Form("TRD task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TESTBIT(fSteerTask, kCalibration)) SETBIT(fSteerTask, kCheckDET);
  if(TESTBIT(fSteerTask, kMultiplicity)) SETBIT(fSteerTask, kEfficiency);
  if(TESTBIT(fSteerTask, kEfficiencyMC)) SETBIT(fSteerTask, kEfficiency);
  if(TESTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
  if(TESTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);
  if(TESTBIT(fSteerTask, kPIDRefMaker)) SETBIT(fSteerTask, kCheckPID);
  if(TESTBIT(fSteerTask, kV0Monitor)) SETBIT(fSteerTask, kCheckPID);

  return fSteerTask;

}

//______________________________________________________
void AliTRDpwg1Helper::MergeProd(const Char_t *mark, const Char_t *files, const Int_t nBatch = 20)
{


  // Clear first predefines
  Char_t MERGE[8]; sprintf(MERGE, "%d.lst", (Int_t)gRandom->Uniform(9999.));
  Char_t PURGE[8]; sprintf(PURGE, "%d.lst", (Int_t)gRandom->Uniform(9999.));
  gSystem->Exec("mkdir -p merge; rm -rf merge/*");

  // purge file list
  std::string filename;
  std::ifstream file(files);
  Int_t iline(0);
  while(getline(file, filename)){
    if(Int_t(filename.find(mark)) < 0) continue;
    gSystem->Exec(Form("echo %s >> %s", filename.c_str(), PURGE));
    iline++;
  }
  Int_t nBatches(iline/nBatch);

  for(Int_t ibatch(0); ibatch<nBatches; ibatch++){
    Int_t first(ibatch*nBatch);
    if(!gSystem->Exec(Form("root.exe -b -q \'$ALICE_ROOT/PWG1/TRD/macros/mergeBatch.C(\"%s\", \"%s\", %d, %d)\'", mark, PURGE, nBatch, first))) continue;
    gSystem->Exec(Form("mv %d_%s merge/", first, mark));
    gSystem->Exec(Form("echo %s/merge/%d_%s >> %s", gSystem->ExpandPathName("$PWD"), first, mark, MERGE));
  }
  gSystem->Exec(Form("root.exe -b -q \'$ALICE_ROOT/PWG1/TRD/macros/mergeBatch.C(\"%s\", \"%s\", %d, 0, kFALSE, kTRUE)\'", mark, MERGE, nBatches));
  gSystem->Exec(Form("mv 0_%s %s", mark, mark));
  
  gSystem->Exec(Form("rm -rfv %s %s merge", MERGE, PURGE));
}

//______________________________________________________
Int_t AliTRDpwg1Helper::GetTaskIndex(const Char_t *name){
  for(Int_t it(0); it<kNTRDTASKS; it++){
    if(strcmp(fgkTRDtaskClassName[it], name)==0) return it;
  }
  return -1;
}

//______________________________________________________
Bool_t AliTRDpwg1Helper::HasReadMCData(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOMC");
}

//____________________________________________
Bool_t AliTRDpwg1Helper::HasReadFriendData(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOFR");
}

