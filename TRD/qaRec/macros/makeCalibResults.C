// Usage:
//   makeCalibResults.C("task", "file_list", kGRID)
//   tasks : "one/more of the following separated by space:
//     "CAL"  : TRD calibration
//     "ALGN" : TRD alignment
//     "PIDR" : TRD PID - reference data
//     "CLRES": Cluster position and error parameterization
//     "NOFR" : Data set does not have AliESDfriends.root 
//     "NOMC" : Data set does not have Monte Carlo Informations (real data), so all tasks which rely
//              on MC information are switched off
//   file_list : is the list of the files to be processed. 
//     They should contain the full path. Here is an example:
// /lustre/alice/local/TRDdata/SIM/P-Flat/TRUNK/RUN0/TRD.CalibName.root
// or for GRID alien:///alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/1/TRD.Calib%.root
//   kGRID : specify if files are comming from a GRID collection [default kFALSE]
//
// HOW TO BUILD THE FILE LIST
//   1. locally
// ls -1 BaseDir/RUN*/TRD.Calib*.root > files.lst
// 
//   2. on Grid
// char *BaseDir="/alice/cern.ch/user/m/mfasel/MinBiasProd/results/ppMinBias80000/";
// TGrid::Connect("alien://");
// TGridResult *res = gGrid->Query(BaseDir, "%/TRD.Calib%.root");
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
// In compiled mode : 
// Don't forget to load first the libraries
// gSystem->Load("libMemStat.so")
// gSystem->Load("libMemStatGui.so")
// gSystem->Load("libANALYSIS.so")
// gSystem->Load("libANALYSISalice.so")
// gSystem->Load("libTRDqaRec.so")
// gSystem->Load("libSTAT.so")
// gSystem->Load("libPWG1.so");
// gSystem->Load("libNetx.so") ;
// gSystem->Load("libRAliEn.so");
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 
//

#if ! defined (__CINT__) || defined (__MAKECINT__)

#include "qaRec/AliTRDrecoTask.h"
#include <fstream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TGridCollection.h>

#endif

#include "AliTRDperformanceTrain.h"
//#include "../../PWG1/macros/AddPerformanceTask.h"

Char_t *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libTRDqaRec.so", "libSTAT.so"};
// define setup
TCanvas *c = 0x0;
Bool_t mc(kFALSE), friends(kFALSE);

void calibrateTRD(TNamed* task);
void makeCalibResults(Char_t *opt, const Char_t *files=0x0, Bool_t kGRID=kFALSE)
{
  if(kGRID){
    if(!gSystem->Getenv("GSHELL_ROOT")){
      Error("makeCalibResults.C", "AliEn not initialized.");
      return;
    }
    TGrid::Connect("alien://");
  }

	// Load Libraries in interactive mode
  Int_t nlibs = static_cast<Int_t>(sizeof(libs)/sizeof(Char_t *));
  for(Int_t ilib=0; ilib<nlibs; ilib++){
    if(gSystem->Load(libs[ilib]) >= 0) continue;
    Error("makeCalibResults.C", Form("Failed to load %s.", libs[ilib]));
    return;
  }

  mc = HasReadMCData(opt);
  friends = HasReadFriendData(opt);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  Int_t fSteerTask = ParseOptions(opt);

  if(!c) c=new TCanvas("c", "Calibration", 10, 10, 800, 500);

  TClass *ctask = new TClass;
  AliAnalysisTask *task = 0x0;
  for(Int_t itask = NTRDQATASKS; itask<NTRDTASKS; itask++){
    if(!TSTBIT(fSteerTask, itask)) continue;
    new(ctask) TClass(fgkTRDtaskClassName[itask]);
    task = (AliAnalysisTask*)ctask->New();
    if(files) mergeProd(Form("TRD.Calib%s.root", task->GetName()), files);

    if(task->IsA()->InheritsFrom("AliTRDrecoTask")) calibrateTRD(task);
  }
  delete ctask;
  delete c;
}


//______________________________________________________
void calibrateTRD(TNamed *otask)
{
  AliTRDrecoTask *task = dynamic_cast<AliTRDrecoTask*>(otask);
  task->SetDebugLevel(0);
 
  AliLog::SetClassDebugLevel(Form("AliTRD%s", task->GetName()), 3); 
  task->SetMCdata(mc);
  task->SetFriends(friends);

  if(!task->Load(Form("%s/TRD.Calib%s.root", gSystem->ExpandPathName("$PWD"), task->GetName()))){
    Error("makeCalibResults.C", Form("Load data container for task %s failed.", task->GetName()));
    delete task;
    return;
  }

  if(!task->PostProcess()){
    Error("makeCalibResults.C", Form("Processing data container for task %s failed.", task->GetName()));
    delete task;
    return;
  }
  for(Int_t ipic=0; ipic<task->GetNRefFigures(); ipic++){
    c->Clear();
    if(!task->GetRefFigure(ipic)) continue;
    c->SaveAs(Form("%s_Fig%02d.gif", task->GetName(), ipic));
  }
  delete task;
}

