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
TObjArray *arr(NULL), *narr[ntasks] = {NULL};
TFile *file(NULL), *fOut(NULL);
Bool_t mc(kFALSE), 
       kWrite(kFALSE),
       kVerbose(kFALSE);
void addPerformance(const Int_t it);
void mergeProdLocal(const Char_t *list, Bool_t doMC=kFALSE, const Int_t nwrite=10000)
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
  printf("Usage : mergeProdLocal(list, doMC, n)\n"
    "  list : the list of AnalysisResults file to be merged. Files can be root or zip archived.\n"
    "       : using \"%s\".\n"
    "  doMC : true if production is MC. Default false.\n"
    "       : using \"%s\".\n"
    "  n    : Optional parameter to specify the number of files to be added until a flush of the\n"
    "         merged file is performed. Use only for large productions.\n"
    "       : using %d.\n", list, (doMC?"true":"false"), nwrite
  );
  
  // build out filename
  TString s, s1, slist=list;
  Char_t fn[200];
  for(Int_t i(slist.Length()); --i;) {
    if(slist(i)!='.') continue;
    s1 = slist(0, i);
    snprintf(fn, 200, "Merged_%s.root", s1.Data());
    break;
  }
  printf("Opening merged file \"%s\"...\n", fn);
  fOut = TFile::Open(fn, "RECREATE");
  fOut->mkdir("TRD_Performance")->Write();

  // read task descriptors a
  for(Int_t it(0); it<ntasks; it++){
    s=taskConfig[it]; 
    s1=s(0, s.Index('_')); sprintf(taskName[it], "%s", s1.Data());
    s = s(s.Index('_')+1, 200); plots[it] = s.Tokenize(":");
  }
  TObjArray *infoGenArr(NULL);
  mc = doMC;
  Int_t ifile(0);
  if(!kVerbose) printf("Merging ");
  while(slist.Gets(fp)){
    if(slist.EndsWith("zip")) snprintf(fn, 200, "%s#AnalysisResults.root", slist.Data());
    else if(slist.EndsWith("root")) snprintf(fn, 200, "%s", slist.Data());
    else{
      Warning("mergeProd", "Don't know what to do with file \"%s\". Skip file.", slist.Data());
      continue;
    }
    if(!(file = TFile::Open(fn))) continue;
    if(kVerbose) Info("mergeProd", "Adding file %s ...", fn);
    else printf(".");fflush(stdout);
    if(!file->cd("TRD_Performance")){ 
      Warning("mergeProd", "Missing TRD_Performance. Skip file.");
      continue;
    }
      
    if(!(arr = (TObjArray*)gDirectory->Get("TRDinfoGen"))){
      Warning("mergeProd", "Missing TRDinfoGen. Skip file.");
      continue;
    }

    if((ifile%nwrite)==(nwrite-1)){ 
      kWrite = kTRUE;
      if(!kVerbose) printf("\nMerging ");
      else printf(" ... Flushing merged file ...");
    }
    fOut->cd("TRD_Performance");
    if(!infoGenArr) infoGenArr=(TObjArray*)arr->Clone();
    else { 
      ((TH1*)infoGenArr->At(0))->Add((TH1*)arr->At(0)); 
      ((TH1*)infoGenArr->At(1))->Add((TH1*)arr->At(1)); 
      ((TH2*)infoGenArr->At(2))->Add((TH2*)arr->At(2)); 
    }
    arr->Delete(); delete arr;
    
    for(Int_t itask(0); itask<ntasks; itask++){
      file->cd("TRD_Performance");
      if(!(arr = (TObjArray*)gDirectory->Get(taskName[itask]))){
        Warning("mergeEntry", "Missing %s. Skip task.", taskName[itask]);
        continue;
      }
      fOut->cd("TRD_Performance");
      addPerformance(itask);
      arr->Delete(); delete arr;
    }
    file->Close(); delete file;
    ifile++; kWrite = kFALSE;
    //if(ifile>=35) break; 
  }
  
  if(!kVerbose) printf("\n");
  Info("mergeProd", "Flushing merged file ...");
  fOut->cd("TRD_Performance");
  if(infoGenArr) infoGenArr->Write("TRDinfoGen", TObject::kSingleKey);
  for(Int_t it(0); it<ntasks; it++){
    if(!(arr = (TObjArray*)gDirectory->Get(taskName[it]))) continue;
    arr->Write(taskName[it], TObject::kSingleKey);
  }  
  fOut->Close();
}

void addPerformance(const Int_t it)
{ 
  if(!narr[it]){
    narr[it] = (TObjArray*)arr->Clone(); narr[it]->SetOwner();
    return;
  }
  
  THnSparse *h(NULL), *H(NULL);
  Int_t nplots(plots[it]->GetEntries());
  TString sname;
  for(Int_t ih(0); ih<nplots; ih++){
    sname= ((TObjString*)(*plots[it])[ih])->String();
    if(!mc && sname.EndsWith("MC")) continue;
    sname.Prepend("h");
    if(!(h = (THnSparse*)arr->FindObject(sname.Data()))) {
      Warning("addPerformance", "Missing %s. from \"%s\"", sname.Data(), file->GetName());
      continue;
    }
    if(!(H = (THnSparse*)narr[it]->FindObject(sname.Data()))) {
      Error("addPerformance", "Missing %s. from Merging.", sname.Data());
      return;
    }
    //printf("  -> %s[%d]\n", sname.Data(), (Int_t)H->GetEntries());
    H->Add(h);
  }
  if(kWrite){
    narr[it]->Write(taskName[it], /*TObject::kOverwrite|*/TObject::kSingleKey);
    narr[it]->Delete(); delete narr[it]; narr[it] = NULL;
  }
}
