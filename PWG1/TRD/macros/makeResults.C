// Usage:
//   makeResults.C("tasks", "file_list", ""task_id, kGRID)
//   tasks : "ALL" or one/more of the following separated by space:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "RES"  : TRD tracking Resolution
//     "PID"  : TRD PID - pion efficiency 
//     "DET"  : Basic TRD Detector checks
//     "NOFR" : Data set does not have AliESDfriends.root 
//     "NOMC" : Data set does not have Monte Carlo Informations (real data), so all tasks which rely
//              on MC information are switched off
//   file_list : is the list of the files to be processed. 
//     They should contain the full path. Here is an example:
// /lustre/alice/local/TRDdata/SIM/P-Flat/TRUNK/RUN0/TRD.Performance.root
// or for GRID alien:///alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/1/TRD.Performance.root
//   task_id : identifier of task speciality as defined by the AddMacro.C. 
//             (e.g. AddTRDresolution.C defines "" for barrel tracks, "K" for kink tracks and "SA" for stand alone tracks)
//   kGRID : specify if files are comming from a GRID collection [default kFALSE]
//
// HOW TO BUILD THE FILE LIST
//   1. locally
// ls -1 BaseDir/RUN*/TRD.Performance.root > files.lst
// 
//   2. on Grid
// char *BaseDir="/alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/";
// TGrid::Connect("alien://");
// TGridResult *res = gGrid->Query(BaseDir, "%/TRD.Performance.root");
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
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TGridCollection.h>

#include "AliLog.h"

#include "PWG1/TRD/AliTRDrecoTask.h"
#include "PWG1/TRD/AliTRDcheckESD.h"

#endif

#include "AliTRDperformanceTrain.h"
#include "helper.C"
//#include "../../PWG1/macros/AddPerformanceTask.h"

Char_t *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libANALYSISalice.so", "libTENDER.so", "libPWG1.so"};
// define setup
TCanvas *c = 0x0;
Bool_t mc(kFALSE), friends(kFALSE);

void processTRD(TNamed* task, const Char_t *filename);
void processESD(TNamed* task, const Char_t *filename);
void makeResults(Char_t *opt = "ALL", const Char_t *files="QAResults.root", Char_t *cid = "", Bool_t kGRID=kFALSE)
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

  mc = HasReadMCData(opt);
  friends = HasReadFriendData(opt);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString outputFile;
  if(!TString(files).EndsWith(".root")){ 
    outputFile = Form("%s/QAResults.root", gSystem->ExpandPathName("$PWD"));
    mergeProd("QAResults.root", files);
  } else {
    outputFile = files;
  }
  Int_t fSteerTask = ParseOptions(opt);

  if(!c) c=new TCanvas("c", "Performance", 10, 10, 800, 500);

  TClass *ctask = new TClass;
  AliAnalysisTask *task = 0x0;
  for(Int_t itask = NTRDQATASKS; itask--;){
    if(!TSTBIT(fSteerTask, itask)) continue;
    new(ctask) TClass(fgkTRDtaskClassName[itask]);
    task = (AliAnalysisTask*)ctask->New();
    task->SetName(Form("%s%s", task->GetName(), cid));
    printf("task %s, output file %s\n", task->GetName(), outputFile.Data());
    if(task->IsA()->InheritsFrom("AliTRDrecoTask")) processTRD(task, outputFile.Data());
    else processESD(task, outputFile.Data());
  }
  delete ctask;
  delete c;
}


//______________________________________________________
void processTRD(TNamed *otask, const Char_t *filename)
{
  printf("process[%s] : %s\n", otask->GetName(), otask->GetTitle());
  Int_t debug(0);
  AliTRDrecoTask *task = dynamic_cast<AliTRDrecoTask*>(otask);
  task->SetDebugLevel(debug);
  AliLog::SetClassDebugLevel(otask->IsA()->GetName(), debug);
  task->SetMCdata(mc);
  task->SetFriends(friends);

  //if(!task->Load(Form("%s/AnalysisResults.root", gSystem->ExpandPathName("$PWD")))){
  if(!task->Load(filename)){
    Error("makeResults.C", Form("Load data container for task %s failed.", task->GetName()));
    delete task;
    return;
  }

  if(!task->PostProcess()){
    Error("makeResults.C", Form("Processing data container for task %s failed.", task->GetName()));
    delete task;
    return;
  }
  for(Int_t ipic=0; ipic<task->GetNRefFigures(); ipic++){
    c->Clear();
    if(!task->GetRefFigure(ipic)) continue;
    c->SaveAs(Form("%s_Fig%02d.gif", task->GetName(), ipic), "gif");
  }
  delete task;
}

//______________________________________________________
void processESD(TNamed *otask, const Char_t *filename)
{
  printf("process[%s] : %s\n", otask->GetName(), otask->GetTitle());

  AliTRDcheckESD *esd = dynamic_cast<AliTRDcheckESD*>(otask);
  if(!esd){
    Info("makeResults.C", Form("Processing of task %s failed.", otask->GetName()));
    delete otask;
    return;
  }
  //if(!esd->Load(Form("%s/AnalysisResults.root", gSystem->ExpandPathName("$PWD")), "TRD_Performance")){
  if(!esd->Load(filename, "TRD_Performance")){
    Error("makeResults.C", Form("Load data container for task %s failed.", esd->GetName()));
    delete esd;
    return;
  }
  esd->Terminate(NULL);

  for(Int_t ipic(0); ipic<esd->GetNRefFigures(); ipic++){
    c->Clear(); 
    if(!esd->GetRefFigure(ipic)) continue;
    c->SaveAs(Form("%s_Fig%02d.gif", esd->GetName(), ipic));
  }
  delete esd;
}
