////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Helper class for PWGPP TRD train                                       //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TError.h>
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

#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliESDtrackCuts.h>
#include <AliTRDtrackerV1.h>

#include "AliTRDpwgppHelper.h"

#include "AliTRDtrackInfo.h"
#include "AliTRDv0Info.h"
#include "AliTRDchmbInfo.h"
#include "AliTRDtriggerInfo.h"
#include "AliTRDeventInfo.h"
#include "AliTRDeventCuts.h"
#include "AliTRDinfoGen.h"
#include "AliTRDcheckESD.h"
#include "AliTRDcheckDET.h"
#include "AliTRDcheckPID.h"
#include "AliTRDcheckTRK.h"
#include "AliTRDcalibration.h"
#include "AliTRDalignmentTask.h"
#include "AliTRDresolution.h"
#include "AliTRDclusterResolution.h"
#include "AliTRDpidRefMaker.h"
#include "AliTRDpidRefMakerNN.h"
#include "AliTRDpidRefMakerLQ.h"
#include "AliTRDefficiency.h"
#include "AliTRDefficiencyMC.h"
#include "AliTRDmultiplicity.h"
#include "AliTRDv0Monitor.h"
#include "AliTRDpwgppHelper.h"

Int_t AliTRDpwgppHelper::fgYear = 0;
const Char_t * AliTRDpwgppHelper::fgkTRDtaskClassName[AliTRDpwgppHelper::kNTRDTASKS] = {
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

const Char_t * AliTRDpwgppHelper::fgkTRDtaskOpt[AliTRDpwgppHelper::kNTRDTASKS+1] = {
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
Bool_t AliTRDpwgppHelper::AddTrainPerformanceTRD(const Char_t *trd, const Char_t *addMacroPath)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTrainPerformanceTRD", "AliAnalysisManager not set!");
    return kFALSE;
  }

  // TRD data containers
  AliAnalysisDataContainer *ci[kNOutSlots];
  AliAnalysisDataContainer *ce[6];
  Char_t macro[100], pars[100]; Int_t npars(0);
  Info("AddTrainPerformanceTRD", "Add Macros taken from %s", addMacroPath);
  Info("AddTrainPerformanceTRD", "TRD wagons \"%s\"", trd);
  Int_t bitmap = ParseOptions(trd);
  for(Int_t it=0; it<kNTRDQATASKS; it++){
    if(!DoTask(it, bitmap)) continue;
    snprintf(macro, 100, "%s/Add%s.C", addMacroPath, TString(TaskClassName(it))(3,20).Data());
    npars=0;
    switch(it){
    case kCheckESD:
      break;
    case kInfoGen:
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p", (void*)ci); npars=1; 
      break;
    case kCheckDET:
      // map slots
      ce[0]=ci[kTracksBarrel];
      ce[1]=ci[kTracksSA];
      ce[2]=ci[kTracksKink];
      ce[3]=ci[kEventInfo];
      ce[4]=ci[kTracklets];
      ce[5]=ci[kClusters];
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p, %d", (void*)ce, TESTBIT(bitmap, AliTRDpwgppHelper::kCalibration)); npars=2; 
      break;
    case kEfficiency:
      // map slots
      ce[0]=ci[kTracksBarrel];
      ce[1]=ci[kTracksITS];
      ce[2]=ci[kTracksKink];
      ce[3]=ci[kEventInfo];
      ce[4]=ci[kTracklets];
      ce[5]=ci[kClusters];
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p, %d, %d", (void*)ce, TESTBIT(bitmap, AliTRDpwgppHelper::kEfficiencyMC), TESTBIT(bitmap, AliTRDpwgppHelper::kMultiplicity)); npars=3; 
      break;
    case kResolution:
      // map slots
      ce[0]=ci[kTracksBarrel];
      ce[1]=ci[kTracksITS];
      ce[2]=ci[kTracksKink];
      ce[3]=ci[kEventInfo];
      ce[4]=ci[kTracklets];
      ce[5]=ci[kClusters];
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p, %d, %d", (void*)ce, TESTBIT(bitmap, AliTRDpwgppHelper::kClErrParam), TESTBIT(bitmap, AliTRDpwgppHelper::kAlignment)); npars=3; 
      break;
    case kCheckPID:
      // map slots
      ce[0]=ci[kTracksBarrel];
      ce[1]=ci[kEventInfo];
      ce[2]=ci[kTracklets];
      ce[3]=ci[kV0List];
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p, (AliAnalysisDataContainer**)%p, %d", (void*)ce, (void*)&ce[4], TESTBIT(bitmap, AliTRDpwgppHelper::kPIDRefMaker)); npars=3; 
      break;
    case kCheckTRK:
      // map slots
      ce[0]=ci[kTracksBarrel];
      ce[1]=ci[kEventInfo];
      ce[2]=ci[kTracklets];
      ce[3]=ci[kClusters];
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p", (void*)ce); npars=1; 
      break;
    case kV0Monitor:
      // slots already mapped by checkPID
      ce[2] = ce[3];
      ce[3] = ce[4];
      snprintf(pars, 100, "(AliAnalysisDataContainer**)%p", (void*)ce); npars=1; 
      break;
    default:
      Warning("AddTrainPerformanceTRD()", "No performance task registered at slot %d.", it);
      npars=-1;
    }
    if(npars>0){
      if(gROOT->Macro(Form("%s(%s)", macro, pars))!=0) return kFALSE;
    } else if(npars==0){
      if(gROOT->Macro(macro)!=0) return kFALSE;
    } else continue;     
 
  }
  return kTRUE;
}

//______________________________________________________
const Char_t* AliTRDpwgppHelper::Translate(Bool_t doCheckESD, Bool_t doCheckDET, Bool_t doEffic, Bool_t doResolution, Bool_t doCheckPID, Bool_t doCheckV0)
{
  TString *opt = new TString("");
  if( doCheckESD==kTRUE &&
      doCheckDET==kTRUE &&
      doEffic==kTRUE &&
      doResolution==kTRUE &&
      doCheckPID==kTRUE &&
      doCheckV0==kTRUE
  ){
    (*opt)="ESD DET EFF RES PID V0";
  } else {
    Bool_t kINDENT(kFALSE);
    if(doCheckESD){ 
      opt->Append("ESD");
      kINDENT=kTRUE;
    }
    if(doCheckDET){ 
      if(kINDENT) opt->Append(" ");
      opt->Append("DET"); 
      kINDENT = kTRUE;
    }
    if(doEffic){ 
      if(kINDENT) opt->Append(" ");
      opt->Append("EFF");
      kINDENT=kTRUE;
    }
    if(doResolution){ 
      if(kINDENT) opt->Append(" ");
      opt->Append("RES");
      kINDENT=kTRUE;
    }
    if(doCheckPID){ 
      if(kINDENT) opt->Append(" ");
      opt->Append("PID");
      kINDENT=kTRUE;
    }
    if(doCheckV0){ 
      if(kINDENT) opt->Append(" ");
      opt->Append("V0");
      kINDENT=kTRUE;
    }
  }

  return (const Char_t*)opt->Data();
}

//______________________________________________________
Bool_t AliTRDpwgppHelper::DoTask(Int_t idx, Int_t map)
{
  return TESTBIT(map, idx);
}

//______________________________________________________
Int_t AliTRDpwgppHelper::ParseOptions(const Char_t *trd)
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
      if(!foundOpt) Warning("AliTRDpwgppHelper::ParseOptions()", "TRD task %s not implemented (yet).", s.Data());
    }
  }
  tasksArray->Delete(); delete tasksArray;

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
void AliTRDpwgppHelper::MergeProd(const Char_t *mark, const Char_t *files, const Int_t nBatch, Int_t level)
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
  Info("MergeProd()", "Merge %d files in %d batches.", iline, nBatches);

  Int_t first(0);
  for(Int_t ibatch(0); ibatch<nBatches; ibatch++){
    first = ibatch*nBatch;
    if(gSystem->Exec(Form("aliroot -b -q \'$ALICE_PHYSICS/PWGPP/TRD/macros/mergeBatch.C(\"%s\", \"%s\", %d, %d)\'", mark, lPURGE, nBatch, first))) continue;
    gSystem->Exec(Form("mv %d_%s merge/%d_%d_%s", first, mark, level, first, mark));
    gSystem->Exec(Form("echo %s/merge/%d_%d_%s >> %s", gSystem->ExpandPathName("$PWD"), level, first, mark, lMERGE));
  }

  if(nBatches==1){
    Info("AliTRDpwgppHelper::MergeProd()", "Rename 1 merged file.");
    gSystem->Exec(Form("mv merge/%d_%d_%s %s", level, first, mark, mark));
  } else if(nBatches<=nBatch){
    Info("AliTRDpwgppHelper::MergeProd()", "Merge %d files in 1 batch.", nBatches);
    if(!gSystem->Exec(Form("aliroot -b -q \'$ALICE_PHYSICS/PWGPP/TRD/macros/mergeBatch.C(\"%s\", \"%s\", %d, 0, kFALSE)\'", mark, lMERGE, nBatches))) return;
    gSystem->Exec(Form("mv 0_%s %s", mark, mark));
  } else {
    level++;
    Info("AliTRDpwgppHelper::MergeProd()", "Merge level %d.", level);
    MergeProd(mark, lMERGE, nBatch, level);
  }
  gSystem->Exec(Form("rm -fv %s %s", lMERGE, lPURGE));
}


//______________________________________________________
Int_t AliTRDpwgppHelper::MergeBatch(const Char_t *mark, const Char_t *files, const Int_t nfiles, const Int_t first, Bool_t kSVN, Bool_t kCLEAR)
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
    Info("AliTRDpwgppHelper::MergeBatch()", "%s", filename.c_str());
    if(!fFM.AddFile(filename.c_str())) return 1;
    arr.Add(new TObjString(filename.c_str()));
    nbatch++;
    if(nbatch==nfiles) break;
  }
  if(!nbatch){
    Warning("AliTRDpwgppHelper::MergeBatch()", "NOTHING TO MERGE"); return 0;
  } else {
    Info("AliTRDpwgppHelper::MergeBatch()", "MERGING FILES[%d] START[%d] %s ... ", nbatch, first, ((nbatch<nfiles)?"INCOMPLETE":""));
  }
  if(!fFM.Merge()){
    Info("AliTRDpwgppHelper::MergeBatch()", "Failed [%s]", fFM.GetOutputFileName());
    return 1;
  }
  Info("AliTRDpwgppHelper::MergeBatch()", "Done [%s]", fFM.GetOutputFileName());

  if(kCLEAR){
    for(Int_t ifile(0); ifile<arr.GetEntries(); ifile++){
      gSystem->Exec(Form("rm -fv %s", ((TObjString*)arr.At(ifile))->GetString().Data()));
    }
  }
  return 0;
}

//______________________________________________________
const Char_t* AliTRDpwgppHelper::Basename(const char* filepath)
{
// Implementation of shell "basename" builtin
  TString s(filepath);
  Int_t idx(s.Last('/')+1);
  s=s(idx, idx+100);
  return s;
}

//______________________________________________________
const Char_t* AliTRDpwgppHelper::Dirname(const char* filepath)
{
// Implementation of shell "dirname" builtin
  TString s(filepath);
  Int_t idx(s.Last('/'));
  s=s(0, idx);
  return s;
}

//______________________________________________________
Int_t AliTRDpwgppHelper::GetTaskIndex(const Char_t *name)
{
//  Give index in TRD train of task class "name"
  for(Int_t it(0); it<kNTRDTASKS; it++){
    if(strcmp(fgkTRDtaskClassName[it], name)==0) return it;
  }
  return -1;
}

//______________________________________________________
Bool_t AliTRDpwgppHelper::HasReadMCData(const Char_t *opt)
{
// Use MC data option
  return !(Bool_t)strstr(opt, "NOMC");
}

//____________________________________________
Bool_t AliTRDpwgppHelper::HasReadFriendData(const Char_t *opt)
{
// Use friends data option
  return !(Bool_t)strstr(opt, "NOFR");
}

