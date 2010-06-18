#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TFileMerger.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TError.h"
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
      for(Int_t itask = 2; itask < NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBITT(fSteerTask, itask); SETBITT(fSteerTask, kInfoGen);
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

  return fSteerTask;
}

//______________________________________________________
void mergeProd(const Char_t *mark, const Char_t *files, const Int_t kBatch = 20)
{
  TFileMerger *fFM = 0x0;
  Bool_t kSVN = kFALSE;

  Int_t jbatch = 0, nbatch = 0;
  std::string filename;
  std::ifstream file(files);
  while(getline(file, filename)){
    if(Int_t(filename.find(mark)) < 0) continue;
    if(!nbatch){
      if(fFM){ 
        delete fFM;
        fFM = new TFileMerger(kTRUE);
      } else fFM = new TFileMerger(kTRUE);
      fFM->OutputFile(Form("%s/%d_%s",  gSystem->ExpandPathName("$PWD"), jbatch, mark));
    }
    if(!kSVN){ // download SVN info for trending
      string base=filename.substr(0, filename.find_last_of('/'));
      if(gSystem->Exec(Form("cp -v %s/svnInfo.log .", base.c_str())) == 0) kSVN=kTRUE;
    }
    fFM->AddFile(filename.c_str()); nbatch++;
    if(nbatch==kBatch){
      printf("MERGING BATCH %d [%d] ... \n", jbatch, nbatch);
      fFM->Merge(); jbatch++;
      nbatch=0;
    }
  }
  if(nbatch){
    printf("MERGING INCOMPLETE BATCH %d [%d] ... \n", jbatch, nbatch);
    fFM->Merge(); jbatch++;
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
