// Usage:
//   makeResults.C("tasks", "file_list", kGRID)
//   tasks : "ALL" or one/more of the following separated by space:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "RES"  : TRD tracking Resolution
//     "CAL"  : TRD calibration
//     "PID"  : TRD PID - pion efficiency 
//     "PIDR" : TRD PID - reference data
//     "DET"  : Basic TRD Detector checks
//     "NOFR" : Data set does not have AliESDfriends.root 
//     "NOMC" : Data set does not have Monte Carlo Informations (real data), so all tasks which rely
//              on MC information are switched off
//   file_list : is the list of the files to be processed. 
//     They should contain the full path. Here is an example:
// /lustre_alpha/alice/TRDdata/HEAD/1.0GeV/RUN0/TRD.TaskResolution.root
// or for GRID alien:///alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/1/TRD.TaskDetChecker.root
//   kGRID : specify if files are comming from a GRID collection [default kFALSE]
//
// HOW TO BUILD THE FILE LIST
//   1. locally
// ls -1 BaseDir/TRD.Task*.root > files.lst
// 
//   2. on Grid
// char *BaseDir="/alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/";
// TGrid::Connect("alien://");
// TGridResult *res = gGrid->Query(BaseDir, "%/TRD.Task%.root");
// TGridCollection *col = gGrid->OpenCollectionQuery(res);
// col->Reset();
// TMap *map = 0x0;
// while(map = (TMap*)col->Next()){
//   TIter it((TCollection*)map);
//   TObjString *info = 0x0;
//   while(info=(TObjString*)it()){
//     printf("alien://%s\n", col->GetLFN(info->GetString().Data()));
//   }
// }
//
// The files which will be processed are the intersection between the
// condition on the tasks and the files in the file list.
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 
//

#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <fstream>
#include "TError.h"
#include <TClass.h>
#include <TFileMerger.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>

#include "qaRec/AliTRDrecoTask.h"

#endif

#include "run.h"

Char_t *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libTRDqaRec.so"};

void makeResults(Char_t *opt = "ALL", const Char_t *files=0x0, Bool_t kGRID=kFALSE)
{
  if(kGRID){
    if(!gSystem->Getenv("GSHELL_ROOT")){
      printf("alien not initialized.\n");
      return;
    }
    TGrid::Connect("alien://");
  }

	// Load Libraries in interactive mode
  Int_t nlibs = static_cast<Int_t>(sizeof(libs)/sizeof(Char_t *));
  for(Int_t ilib=0; ilib<nlibs; ilib++){
    if(gSystem->Load(libs[ilib]) >= 0) continue;
    printf("Failed to load %s.\n", libs[ilib]);
    return;
  }


  gStyle->SetOptStat(0);
  Bool_t mc      = kTRUE;
  Bool_t friends = kTRUE;

  // select tasks to process; we should move it to an 
  // individual function and move the task identifiers 
  // outside the const space
  TString tasks(opt);
  TObjArray *tasksArray = tasks.Tokenize(" ");
  Int_t fSteerTask = 0;
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 1; itask < NQATASKS; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOFR") == 0){ 
      friends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      mc = kFALSE;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 1; itask < NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("makeResults.C", Form("Task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TSTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kPIDRefMaker)) SETBIT(fSteerTask, kPIDChecker);
  if(TSTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);

  // file merger object
  TFileMerger *fFM = 0x0;
  TClass *ctask = new TClass;
  AliTRDrecoTask *task = 0x0;
  Int_t nFiles;

  if(gSystem->AccessPathName(Form("%s/merge",  gSystem->ExpandPathName("$PWD")))) gSystem->Exec(Form("mkdir -v %s/merge",  gSystem->ExpandPathName("$PWD")));

  for(Int_t itask = 1; itask<NTRDTASKS; itask++){
    if(!TSTBIT(fSteerTask, itask)) continue;

    new(ctask) TClass(fgkTRDtaskClassName[itask]);
    task = (AliTRDrecoTask*)ctask->New();
    task->SetDebugLevel(0);
    task->SetMCdata(mc);
    task->SetFriends(friends);

     // setup filelist
    nFiles = 0;
    TString mark("TRD.Performance.root");
    string filename;
    if(files){
      ifstream filestream(files);
      while(getline(filestream, filename)){
        if(Int_t(filename.find(mark.Data())) < 0) continue;
        nFiles++;
      }
    } else {
      nFiles = !gSystem->AccessPathName("TRD.Performance.root");
    }

    if(!nFiles){
      Info("makeResults.C", Form("No Files found for Task %s", task->GetName()));
      delete task;
      continue;
    }
    Info("makeResults.C", Form("  Processing %d files for task %s ...", nFiles, task->GetName()));

    if(files){
      fFM = new TFileMerger(kTRUE);
      fFM->OutputFile(Form("%s/merge/TRD.Performance.root",  gSystem->ExpandPathName("$PWD")));

      ifstream file(files);
      while(getline(file, filename)){
        if(Int_t(filename.find(mark.Data())) < 0) continue;
        fFM->AddFile(filename.c_str());
      }
      fFM->Merge();
      delete fFM;
      if(!task->Load(Form("%s/merge/TRD.Performance.root", gSystem->ExpandPathName("$PWD")))){
        delete task;
        break;
      }
    } else{
      if(!task->Load(Form("%s/TRD.Performance.root", gSystem->ExpandPathName("$PWD")))){
        delete task;
        break;
      }
    }

    printf("Processing ...\n");
    task->PostProcess();
    TCanvas *c=new TCanvas("c", "TRD Performance", 10, 10, 800, 500);
    for(Int_t ipic=0; ipic<task->GetNRefFigures(); ipic++){
      if(!task->GetRefFigure(ipic)) continue;
      c->SaveAs(Form("%s_fig%d.gif", task->GetName(), ipic));
      c->Clear();
    }
    delete c;
    delete task;
  }
  delete ctask;
}

