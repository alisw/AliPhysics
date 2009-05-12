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
// /lustre/alice/local/TRDdata/SIM/P-Flat/TRUNK/RUN0/TRD.Performance.root
// or for GRID alien:///alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/1/TRD.Performance.root
//   kGRID : specify if files are comming from a GRID collection [default kFALSE]
//
// HOW TO BUILD THE FILE LIST
//   1. locally
// ls -1 BaseDir/RUN*/TRD.Performance.root > files.lst
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

#include "macros/AliTRDperformanceTrain.h"
//#include "../../PWG1/macros/AddPerformanceTask.h"

Char_t *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libTRDqaRec.so"};

void mergeProd(const Char_t *mark="TRD.Performance.root", const Char_t *files=0);
void makeResults(Char_t *opt = "ALL", const Char_t *files=0x0, Bool_t kGRID=kFALSE)
{
  if(kGRID){
    if(!gSystem->Getenv("GSHELL_ROOT")){
      Error("makeResults.C", "AliEn not initialized.");
      return;
    }
    TGrid::Connect("alien://");
  }

	// Load Libraries in interactive mode
  Int_t nlibs = static_cast<Int_t>(sizeof(libs)/sizeof(Char_t *));
  for(Int_t ilib=0; ilib<nlibs; ilib++){
    if(gSystem->Load(libs[ilib]) >= 0) continue;
    Error("makeResults.C", Form("Failed to load %s.", libs[ilib]));
    return;
  }

  // define setup 
  gStyle->SetOptStat(0);
  Bool_t mc      = kTRUE;
  Bool_t friends = kTRUE;


  if(files) mergeProd("TRD.Performance.root", files);
  Int_t fSteerTask = ParseOptions(opt);
  TCanvas *c=new TCanvas("c", "TRD Performance", 10, 10, 800, 500);

  TClass *ctask = new TClass;
  AliTRDrecoTask *task = 0x0;
  for(Int_t itask = NTRDQATASKS; itask--;){
    if(!TSTBIT(fSteerTask, itask)) continue;

    new(ctask) TClass(fgkTRDtaskClassName[itask]);
    task = (AliTRDrecoTask*)ctask->New();
    task->SetDebugLevel(0);
    task->SetMCdata(mc);
    task->SetFriends(friends);

    if(!task->Load(Form("%s/TRD.Performance.root", gSystem->ExpandPathName("$PWD")))){
      Error("makeResults.C", Form("Load data container for task %s failed.", task->GetName()));
      delete task;
      break;
    }

    if(!task->PostProcess()){
      Error("makeResults.C", Form("Processing data container for task %s failed.", task->GetName()));
      delete task;
      break;
    }
    for(Int_t ipic=0; ipic<task->GetNRefFigures(); ipic++){
      c->Clear();
      if(!task->GetRefFigure(ipic)) continue;
      c->SaveAs(Form("%s_Fig%02d.gif", task->GetName(), ipic));
    }
    delete task;
  }
  delete ctask;
  delete c;
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

