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
// gSystem->Load("libTENDER.so")
// gSystem->Load("libSTAT.so")
// gSystem->Load("libPWGPP.so");
// gSystem->Load("libNetx.so") ;
// gSystem->Load("libRAliEn.so");
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 
//

#if ! defined (__CINT__) || defined (__MAKECINT__)

#include "AliLog.h"
#include "PWGPP/TRD/AliTRDrecoTask.h"
#include <fstream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TGridCollection.h>

#endif

Char_t const *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libANALYSISalice.so", "libTENDER.so", "libSTAT.so", "libPWGPP.so", "libPWGmuon.so"};

// define setup
TClass *ctask = new TClass;
TCanvas *c(NULL);
Bool_t fMC(kFALSE), fFriends(kFALSE);
Char_t const *fFiles(NULL);

void calibrateTRD(Int_t itask, Char_t const* ntask=NULL, Char_t const* nfile=NULL);
void makeCalibResults(Char_t *opt, Char_t const *files=NULL, Bool_t kGRID=kFALSE)
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

  fMC = HasReadMCData(opt);
  fFriends = HasReadFriendData(opt);
  fFiles = files;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  Int_t fSteerTask = ParseOptions(opt);

  if(!c) c=new TCanvas("c", "Calibration", 10, 10, 800, 500);

  for(Int_t itask = AliTRDpwgppHelper::kNTRDQATASKS; itask<AliTRDpwgppHelper::NTRDTASKS; itask++){
    if(!TESTBIT(fSteerTask, itask)) continue;
    switch(itask){
    case AliTRDpwgppHelper::kPIDRefMaker:
      calibrateTRD(itask, "AliTRDpidRefMakerLQ", "PIDrefMaker");
      //calibrateTRD(itask, "AliTRDpidRefMakerNN", "PIDrefMaker");
      break;
    default:
      calibrateTRD(itask);
      break;
    }
  }
  delete ctask;
  delete c;
}


//______________________________________________________
void calibrateTRD(Int_t itask, Char_t const* ntask, Char_t const* nfile)
{
  if(!ntask) ntask=AliTRDpwgppHelper::fgkTRDtaskClassName[itask];
  new(ctask) TClass(ntask);
  if(!ctask){
    Error("makeCalibResults.C", Form("Asking for wrong class name [%s].", ntask));
    return;
  }
  AliTRDrecoTask *task = static_cast<AliTRDrecoTask*>(ctask->New());
  if(!task){
    Error("makeCalibResults.C", Form("Class name [%s] does not inherit from AliTRDrecoTask.", ntask));
    return;
  }
  if(!nfile) nfile=Form("TRD.Calib%s.root", task->GetName());
  if(fFiles) mergeProd(Form("TRD.Calib%s.root", nfile), fFiles);
  task->SetDebugLevel(0);
  task->SetMCdata(fMC);
  task->SetFriends(fFriends);
  AliLog::SetClassDebugLevel(ntask, 3); 

  if(!task->Load(nfile)){
    Error("makeCalibResults.C", Form("Loading data container for task[%s] failed.", ntask));
    delete task;
    return;
  }
  if(!task->PostProcess()){
    Error("makeCalibResults.C", Form("Processing data container for task[%s] failed.", ntask));
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

