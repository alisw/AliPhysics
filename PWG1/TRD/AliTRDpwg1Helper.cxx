////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Helper class for PWG1 TRD train                                       //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include <Rtypes.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TFileMerger.h>
#include <TRandom.h>
#include <TString.h>
#include <TSystem.h>

#include <string>
#include <cstring>
#include <fstream>

#include "AliTRDpwg1Helper.h"

const Char_t * AliTRDpwg1Helper::fgkTRDtaskClassName[AliTRDpwg1Helper::kNTRDTASKS] = {
  "AliTRDcheckESD"
  ,"AliTRDinfoGen"
  ,"AliTRDcheckDET"
  ,"AliTRDefficiency"
  ,"AliTRDresolution"
  ,"AliTRDcheckPID"
  ,"AliTRDv0Monitor"
  ,"AliTRDcheckTRK"
  ,"AliTRDcalibration"
  ,"AliTRDefficiencyMC"
  ,"AliTRDalignmentTask"
  ,"AliTRDpidRefMaker"
  ,"AliTRDclusterResolution"
  ,"AliTRDmultiplicity"
};

const Char_t * AliTRDpwg1Helper::fgkTRDtaskOpt[AliTRDpwg1Helper::kNTRDTASKS+1] = {
  "ESD"
  ,"GEN"
  ,"DET"
  ,"EFF"
  ,"RES"
  ,"PID"
  ,"V0"
  ,"TRK"
  ,"CAL"
  ,"EFFC"
  ,"ALGN"
  ,"PIDR"
  ,"CLRES"
  ,"MULT"
  ,"ALL"
};

//______________________________________________________
Bool_t AliTRDpwg1Helper::DoTask(Int_t idx, Int_t map)
{
  return TESTBIT(map, idx);
}

//______________________________________________________
Int_t AliTRDpwg1Helper::ParseOptions(Char_t *trd)
{
// Parse space separated options.
// Possible options are:
//      "ALL" : [default] all performance (no calibration) tasks
// ------- Performance tasks ----------
//     "ESD"  : Basic TRD Detector checks on ESD only (no TRD tracks analysed)
//     "DET"  : Basic TRD Detector checks
//     "RES"  : TRD tracking Resolution
//     "EFF"  : TRD Tracking Efficiency
//     "PID"  : TRD PID - pion efficiency
//     "V0"   : monitor V0 performance for use in TRD PID calibration
// ------- Calibration tasks ----------
//     "TRK"  : multidimensional tracking performance resolution
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "MULT"  : TRD single track selection
//     "CLRES": clusters Resolution
//     "CAL"  : TRD calibration
//     "ALGN" : TRD alignment
//     "PIDR" : TRD PID - reference data
// ------- SPECIAL OPTIONS -----------
//     "NOFR" : Data set does not have AliESDfriends.root
//     "NOMC" : Data set does not have Monte Carlo Informations (default have MC),

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
      if(!foundOpt) Warning("AliTRDpwg1Helper::ParseOptions()", Form("TRD task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
//  if(TESTBIT(fSteerTask, kCheckTRK)) SETBIT(fSteerTask, kResolution);
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
void AliTRDpwg1Helper::MergeProd(const Char_t *mark, const Char_t *files, const Int_t nBatch, Int_t level)
{
// Recursively merge files named "mark" from list in "files" in groups of "nBatch" files.
// parameter "level" is used to index recurent calls of this function.

  Char_t lMERGE[8]; snprintf(lMERGE, 8, "%04d.lst", (Int_t)gRandom->Uniform(9999.));
  Char_t lPURGE[8]; snprintf(lPURGE, 8, "%04d.lst", (Int_t)gRandom->Uniform(9999.));

  // purge file list
  std::string filename;
  std::ifstream file(files);
  Int_t iline(0);
  while(getline(file, filename)){
    if(Int_t(filename.find(mark)) < 0) continue;
    gSystem->Exec(Form("echo %s >> %s", filename.c_str(), lPURGE));
    iline++;
  }
  Int_t nBatches=Int_t(TMath::Ceil(Double_t(iline)/nBatch));
  Info("MergeProd()", Form("Merge %d files in %d batches.", iline, nBatches));

  Int_t first(0);
  for(Int_t ibatch(0); ibatch<nBatches; ibatch++){
    first = ibatch*nBatch;
    if(gSystem->Exec(Form("aliroot -b -q \'$ALICE_ROOT/PWG1/TRD/macros/mergeBatch.C(\"%s\", \"%s\", %d, %d)\'", mark, lPURGE, nBatch, first))) continue;
    gSystem->Exec(Form("mv %d_%s merge/%d_%d_%s", first, mark, level, first, mark));
    gSystem->Exec(Form("echo %s/merge/%d_%d_%s >> %s", gSystem->ExpandPathName("$PWD"), level, first, mark, lMERGE));
  }

  if(nBatches==1){
    Info("AliTRDpwg1Helper::MergeProd()", "Rename 1 merged file.");
    gSystem->Exec(Form("mv merge/%d_%d_%s %s", level, first, mark, mark));
  } else if(nBatches<=nBatch){
    Info("AliTRDpwg1Helper::MergeProd()", Form("Merge %d files in 1 batch.", nBatches));
    if(!gSystem->Exec(Form("aliroot -b -q \'$ALICE_ROOT/PWG1/TRD/macros/mergeBatch.C(\"%s\", \"%s\", %d, 0, kFALSE)\'", mark, lMERGE, nBatches))) return;
    gSystem->Exec(Form("mv 0_%s %s", mark, mark));
  } else {
    level++;
    Info("AliTRDpwg1Helper::MergeProd()", Form("Merge level %d.", level));
    MergeProd(mark, lMERGE, nBatch, level);
  }
  gSystem->Exec(Form("rm -fv %s %s", lMERGE, lPURGE));
}


//______________________________________________________
Int_t AliTRDpwg1Helper::MergeBatch(const Char_t *mark, const Char_t *files, const Int_t nfiles, const Int_t first, Bool_t kSVN, Bool_t kCLEAR)
{
// Merge files specified in the file list "files" by the token "mark".
// The script will merge "nfiles" files starting from the "first" file.
// If the file "svnInfo.log" is found together with the files to be merged it is copied locally
// if option "kSVN". The input files are removed from disk if option "kCLEAR".
//
// On return the name of the merged file is return or NULL in case of failure.
//
  TObjArray arr(nfiles); arr.SetOwner(kTRUE);
  TFileMerger fFM(kTRUE);
  fFM.OutputFile(Form("%s/%d_%s",  gSystem->ExpandPathName("$PWD"), first, mark));
  Int_t iline(0), nbatch(0);
  std::string filename;
  std::ifstream file(files);
  while(getline(file, filename)){
    if(Int_t(filename.find(mark)) < 0) continue;
    if(iline<first){
      iline++;
      continue;
    }
    if(kSVN){ // download SVN info for trending
      if(gSystem->Exec(Form("if [ ! -f svnInfo.log ]; then cp -v %s/svnInfo.log %s; fi", Dirname(filename.c_str()), gSystem->ExpandPathName("$PWD"))) == 0) kSVN=kFALSE;
    }
    Info("AliTRDpwg1Helper::MergeBatch()", filename.c_str());
    if(!fFM.AddFile(filename.c_str())) return NULL;
    arr.Add(new TObjString(filename.c_str()));
    nbatch++;
    if(nbatch==nfiles) break;
  }
  if(!nbatch){
    Warning("AliTRDpwg1Helper::MergeBatch()", "NOTHING TO MERGE"); return NULL;
  } else {
    Info("AliTRDpwg1Helper::MergeBatch()", "MERGING FILES[%d] START[%d] %s ... ", nbatch, first, ((nbatch<nfiles)?"INCOMPLETE":""));
  }
  if(!fFM.Merge()){
    Info("AliTRDpwg1Helper::MergeBatch()", "Failed [%s]", fFM.GetOutputFileName());
    return 1;
  }
  Info("AliTRDpwg1Helper::MergeBatch()", "Done [%s]", fFM.GetOutputFileName());

  if(kCLEAR){
    for(Int_t ifile(0); ifile<arr.GetEntries(); ifile++){
      gSystem->Exec(Form("rm -fv %s", ((TObjString*)arr.At(ifile))->GetString().Data()));
    }
  }
  return 0;
}

//______________________________________________________
const Char_t* AliTRDpwg1Helper::Basename(const char* filepath)
{
// Implementation of shell "basename" builtin
  TString s(filepath);
  Int_t idx(s.Last('/')+1);
  s=s(idx, idx+100);
  return s;
}

//______________________________________________________
const Char_t* AliTRDpwg1Helper::Dirname(const char* filepath)
{
// Implementation of shell "dirname" builtin
  TString s(filepath);
  Int_t idx(s.Last('/'));
  s=s(0, idx);
  return s;
}

//______________________________________________________
Int_t AliTRDpwg1Helper::GetTaskIndex(const Char_t *name)
{
//  Give index in TRD train of task class "name"
  for(Int_t it(0); it<kNTRDTASKS; it++){
    if(strcmp(fgkTRDtaskClassName[it], name)==0) return it;
  }
  return -1;
}

//______________________________________________________
Bool_t AliTRDpwg1Helper::HasReadMCData(Char_t *opt)
{
// Use MC data option
  return !(Bool_t)strstr(opt, "NOMC");
}

//____________________________________________
Bool_t AliTRDpwg1Helper::HasReadFriendData(Char_t *opt)
{
// Use friends data option
  return !(Bool_t)strstr(opt, "NOFR");
}

