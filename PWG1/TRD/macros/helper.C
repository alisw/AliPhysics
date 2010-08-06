#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TFileMerger.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TError.h"
#include "TRandom.h"
#include <fstream>
#endif

#ifndef HELPER_C
#define HELPER_C

//____________________________________________
Int_t ParseOptions(Char_t *trd)
{
  Int_t fSteerTask = 1;
  TObjArray *tasksArray = TString(trd).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 0; itask < NTRDQATASKS; itask++) SETBITT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOMC") == 0 || s.CompareTo("NOFR") == 0){
      continue; // taken care by special functions
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 0; itask < NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBITT(fSteerTask, itask); 
        if(itask>1) SETBITT(fSteerTask, kInfoGen);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("ParseOptions()", Form("TRD task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TSTBIT(fSteerTask, kCalibration)) SETBITT(fSteerTask, kCheckDET);
  if(TSTBIT(fSteerTask, kMultiplicity)) SETBITT(fSteerTask, kEfficiency);
  if(TSTBIT(fSteerTask, kEfficiencyMC)) SETBITT(fSteerTask, kEfficiency);
  if(TSTBIT(fSteerTask, kClErrParam)) SETBITT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kAlignment)) SETBITT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kPIDRefMaker)) SETBITT(fSteerTask, kCheckPID);
  if(TSTBIT(fSteerTask, kV0Monitor)) SETBITT(fSteerTask, kCheckPID);

  return fSteerTask;
}

//______________________________________________________
void mergeProd(const Char_t *mark, const Char_t *files, const Int_t nBatch = 20)
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

#endif
