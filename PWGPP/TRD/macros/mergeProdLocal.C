#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TError.h"
#include <THnSparse.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2.h>
#endif

const Int_t ntasks(3);
const Char_t *taskConfig[]={
  "TRDefficiency_EFF"
  ,"TRDresolution_Cluster2Track:Tracklet2Track:Tracklet2TRDin:Cluster2MC:Tracklet2MC:TRDin2MC:TRD2MC"
  ,"TRDresolutionK_Cluster2Track:Tracklet2Track:Tracklet2TRDin:Cluster2MC:Tracklet2MC:TRDin2MC:TRD2MC"
};
Char_t taskName[ntasks][100]={{""}};
TObjArray *plots[ntasks] = {NULL};
THnSparse **H[ntasks]={NULL};
TObjArray *arr(NULL);
TFile *file(NULL);
void addPerformance(const Int_t it);
void mergeProdLocal(const Char_t *list)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTender.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGPP.so");

  FILE *fp(NULL);
  if(!(fp = fopen(list, "rt"))){
    Error("mergeProd.C", "Missing input list \"%s\".", list);
    return;
  }
  // read task descriptors
  TString s, s1; 
  for(Int_t it(0); it<ntasks; it++){
    s=taskConfig[it]; 
    s1=s(0, s.Index('_')); sprintf(taskName[it], "%s", s1.Data());
    s = s(s.Index('_')+1, 200); plots[it] = s.Tokenize(":");
    H[it] = new THnSparse*[plots[it]->GetEntries()];
    memset(H[it], 0, plots[it]->GetEntries()*sizeof(THnSparse*));
  }
  TObjArray *infoGenArr(NULL);

  printf("Merging ");
  Char_t fn[200]; TString slist;
  while(slist.Gets(fp)){
    if(slist.EndsWith("zip")) snprintf(fn, 200, "%s#AnalysisResults.root", slist.Data());
    else if(slist.EndsWith("root")) snprintf(fn, 200, "%s", slist.Data());
    else{
      Warning("mergeProd", "Don't know what to do with file \"%s\". Skip file.", slist.Data());
      continue;
    }
    if(!(file = TFile::Open(fn))) continue;
    //Info("mergeProd", "Adding file %s ...", fn);
    printf(".");fflush(stdout);
    if(!file->cd("TRD_Performance")){ 
      Warning("mergeProd", "Missing TRD_Performance. Skip file.");
      continue;
    }
      
    if(!(arr = (TObjArray*)gDirectory->Get("TRDinfoGen"))){
      Warning("mergeProd", "Missing TRDinfoGen. Skip file.");
      continue;
    }
    if(!infoGenArr) infoGenArr=(TObjArray*)arr->Clone();
    else{ 
      ((TH1*)infoGenArr->At(0))->Add((TH1*)arr->At(0)); 
      ((TH1*)infoGenArr->At(1))->Add((TH1*)arr->At(1)); 
      ((TH2*)infoGenArr->At(2))->Add((TH2*)arr->At(2)); 
      //((AliTRDtriggerInfo*)infoGenArr->At(3))->Merge((AliTRDtriggerInfo*)arr->At(3)); 
    }
    for(Int_t itask(0); itask<ntasks; itask++){
      if(!(arr = (TObjArray*)gDirectory->Get(taskName[itask]))){
        Warning("mergeEntry", "Missing %s. Skip task.", taskName[itask]);
        continue;
      }
      addPerformance(itask);
    }
    file->Close();
    delete file;
  }
  
  slist=list;
  for(Int_t i(slist.Length()); --i;) {
    if(slist(i)!='.') continue;
    s1 = slist(0, i);
    snprintf(fn, 200, "Merged_%s.root", s1.Data());
    break;
  }
  printf("\nWritting to \"%s\"...\n", fn);
  TFile::Open(fn, "RECREATE");
  gFile->mkdir("TRD_Performance")->Write();
  gFile->cd("TRD_Performance");
  infoGenArr->Write("TRDinfoGen", 1);
  for(Int_t it(0); it<ntasks; it++){
    Int_t nplots(plots[it]->GetEntries());
    arr = new TObjArray(nplots); arr->SetOwner();
    for(Int_t ih(0); ih<nplots; ih++) arr->AddAt(H[it][ih], ih);
    arr->Write(taskName[it], 1);
  }  
  gFile->Close();
}

void addPerformance(const Int_t it)
{ 
  THnSparse *h(NULL);
  Int_t nplots(plots[it]->GetEntries());
  char sname[50];
  for(Int_t ih(0); ih<nplots; ih++){
    snprintf(sname, 50, "h%s", ((TObjString*)(*plots[it])[ih])->GetName());
    if(!(h = (THnSparse*)arr->FindObject(sname))) {
      Warning("addPerformance", "Missing %s.", sname);
      continue;
    }
    gROOT->cd();
    if(H[it][ih]) H[it][ih]->Add(h);
    else H[it][ih] = (THnSparse*)h->Clone(sname);
    file->cd("TRD_Performance");
  }
}
