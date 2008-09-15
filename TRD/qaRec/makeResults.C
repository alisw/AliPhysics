#include "run.h"
#define BIT(n)        (1 << (n))
#define TESTBIT(n,i)  ((Bool_t)(((n) & BIT(i)) != 0))


const Int_t fSteerTask = 0xffffffff;
Char_t *fTaskClass[fknTasks] = {
  "AliTRDtrackInfoGen"
  ,"AliTRDtrackingEfficiency"
  ,"AliTRDtrackingEfficiencyCombined"
  ,"AliTRDtrackingResolution"
  ,"AliTRDcalibration"
  ,"AliTRDpidChecker"
  ,"AliTRDcheckDetector"
};

const Int_t nlibs = 3;
Char_t *libs[] = {"libProofPlayer.so", "libANALYSIS.so", "libTRDqaRec.so"};
void makeResults(Char_t *dir=".")
{
  for(Int_t ilib=0; ilib<nlibs; ilib++){
    if(!gSystem->Load(libs[ilib])) continue;
    printf("Failed to load %s.\n", libs[ilib]);
    return;
  }
  gStyle->SetOptStat(0);

  Int_t ndir = 1;

  // file merger object
  TFileMerger *fFM = 0x0;
  TClass *ctask = 0x0;
  TObject *o = 0x0;
  TObjArray *fContainer = 0x0;
  AliTRDrecoTask *task = 0x0;

  printf("\n\tPROCESSING DATA FOR TASKS:\n");
  for(Int_t itask = 3; itask <4/*fknTasks*/ ; itask++){
    if(!TESTBIT(fSteerTask, itask)) continue;

    ctask = new TClass(fTaskClass[itask]);
    task = (AliTRDrecoTask*)ctask->New();
    task->SetDebugLevel(2);
    printf("\t%s [%s]\n", task->GetTitle(), task->GetName());

    fFM = new TFileMerger(kTRUE);
    fFM->OutputFile(Form("merge/TRD.Task%s.root",  task->GetName()));
    
    Int_t idir = 0;
    while(idir<ndir){
      fFM->AddFile(Form("./TRD.Task%s.root",  task->GetName()));
      idir++;
    }
    fFM->Merge();
    delete fFM;
    task->Load(Form("merge/TRD.Task%s.root", task->GetName()));
    task->PostProcess();
    if(!(fContainer = task->Container())) {
      delete task;
      delete ctask;
      continue;
    } 
    for(Int_t ipic=0; ipic<task->GetNRefFigures(); ipic++){
      Int_t ifirst, ilast;
      task->GetRefFigure(ipic, ifirst, ilast);
      if(!(o = fContainer->At(ifirst))) continue;
      
      if(o->InheritsFrom("TH1")){ 
        h = dynamic_cast<TH1*>(o);
        h->Draw("pl");
      } else if(o->InheritsFrom("TGraph")){ 
        g = dynamic_cast<TGraph*>(o);
        g->Draw("apl");
      } else{
        printf("No idea how to plot object of type %s.\n", o->IsA()->GetName());
        printf("Please teach me.\n");
        continue;
      }

      for(Int_t ig=ifirst+1; ig<ilast; ig++){
        if(!(o = fContainer->At(ig))) continue;
        if(o->InheritsFrom("TH1")){
          h = dynamic_cast<TH1*>(o);
          h->Draw("plsame");
        } else if(o->InheritsFrom("TGraph")){
          g = dynamic_cast<TGraph*>(o);
          g->Draw("pl");
        }
      }
      if(gPad) gPad->SaveAs(Form("%s_fig%d.gif", task->GetName(), ipic));
    }
    delete task;
    delete ctask;
  }
}