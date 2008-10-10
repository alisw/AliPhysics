#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TError.h"
#include <TClass.h>
#include <TFileMerger.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPython.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>

#include "qaRec/AliTRDrecoTask.h"

#endif

#include "run.h"

Char_t *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libTRDqaRec.so", "libPyROOT"};

void makeResults(Char_t *tasks = "ALL", Char_t* dir=0x0)
{
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
  TObjArray *tasksArray = TString(tasks).Tokenize(" ");
  Int_t fSteerTask = 0;
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 1; itask < fknTasks; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOFR") == 0){ 
      friends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      mc = kFALSE;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 1; itask < fknTasks; itask++){
        if(s.CompareTo(fTaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("makeResults.C", Form("Task %s not implemented (yet).", s.Data()));
    }
  }


  // catch the list of files using the ROOT Python Interface
  TPython *pyshell = new TPython();
  pyshell->Exec("import commands");

  // file merger object
  TFileMerger *fFM = new TFileMerger();
  TClass *ctask = 0x0;
  TObject *o = 0x0;
  TObjArray *fContainer = 0x0;
  AliTRDrecoTask *task = 0x0;

  if(gSystem->AccessPathName(Form("%s/merge",  gSystem->ExpandPathName("$PWD")))) gSystem->Exec(Form("mkdir -v %s/merge",  gSystem->ExpandPathName("$PWD")));

  printf("\n\tPROCESSING DATA FOR TASKS [%b]:\n", fSteerTask);
  for(Int_t itask = 1; itask <fknTasks; itask++){
    if(!TESTBIT(fSteerTask, itask)) continue;

    ctask = new TClass(fTaskClass[itask]);
    task = (AliTRDrecoTask*)ctask->New();
    task->SetDebugLevel(0);
    task->SetMCdata(mc);
    task->SetFriends(friends);
    printf("\t%s [%s]\n", task->GetTitle(), task->GetName());

     // setup filelist
    TString pathname = gSystem->ExpandPathName( dir ? dir : "$PWD");
    TString filestring((const Char_t*) pyshell->Eval(Form("commands.getstatusoutput(\"find %s | grep TRD.Task%s.root\")[1]", pathname.Data(), task->GetName())));
    TObjArray *filenames = filestring.Tokenize("\n");
    Int_t nFiles = filenames->GetEntriesFast();
    if(!nFiles){
      printf("No Files found for Task %s\n", task->GetName());
      delete task;
      delete ctask;
      continue;
    }

    if(nFiles>1){
      fFM = new(fFM) TFileMerger(kTRUE);
      fFM->OutputFile(Form("%s/merge/TRD.Task%s.root",  gSystem->ExpandPathName("$PWD"), task->GetName()));
      for(Int_t ifile = 0; ifile < nFiles; ifile++){
        TString filename = (dynamic_cast<TObjString *>(filenames->UncheckedAt(ifile)))->String();
        if(filename.Contains("merge")) continue;
        //printf("\tProcessing %s ...\n", filename.Data());
        fFM->AddFile(filename.Data());
      }
      fFM->Merge();
      fFM->~TFileMerger();
      task->Load(Form("%s/merge/TRD.Task%s.root", gSystem->ExpandPathName("$PWD"), task->GetName()));
    } else task->Load((dynamic_cast<TObjString *>(filenames->UncheckedAt(0)))->String().Data());

    if(!(fContainer = task->Container())) {
      delete task;
      delete ctask;
      continue;
    } 
    
    task->PostProcess();
    for(Int_t ipic=0; ipic<task->GetNRefFigures(); ipic++){
      TCanvas *c = new TCanvas("c", "", 500, 500);
      Int_t ifirst, ilast; Option_t *opt;
      TH1 *h = 0x0; TGraph *g = 0x0;
      task->GetRefFigure(ipic, ifirst, ilast, opt);
      if(!(o = fContainer->At(ifirst))) continue;
      
      if(o->InheritsFrom("TH1")){ 
        h = dynamic_cast<TH1*>(o);
        h->Draw(opt);
      } else if(o->InheritsFrom("TGraph")){ 
        g = dynamic_cast<TGraph*>(o);
        g->Draw(Form("a%s", opt));
      } else{
        printf("No idea how to plot object of type %s.\n", o->IsA()->GetName());
        printf("Please teach me.\n");
        continue;
      }

      for(Int_t ig=ifirst+1; ig<ilast; ig++){
        if(!(o = fContainer->At(ig))) continue;
        if(o->InheritsFrom("TH1")){
          h = dynamic_cast<TH1*>(o);
          h->Draw(Form("%ssame", opt));
        } else if(o->InheritsFrom("TGraph")){
          g = dynamic_cast<TGraph*>(o);
          g->Draw(opt);
        }
      }
      c->SaveAs(Form("%s_fig%d.gif", task->GetName(), ipic));
      delete c;
    }
    delete task;
    delete ctask;
  }
  delete pyshell;
}

