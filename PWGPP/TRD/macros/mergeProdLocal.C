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

const Char_t *nameTask="TRD_Performance";
const Int_t ntasks(4);
const Char_t *taskConfig[]={
  "TRDinfoGen_hStat:hEv:hBCtrack"
  ,"TRDefficiency_hEFF"
  ,"TRDresolution_hCluster2Track:hTracklet2Track:hTracklet2TRDin:hCluster2MC:hTracklet2MC:hTRDin2MC:hTRD2MC"
  ,"TRDresolutionK_hCluster2Track:hTracklet2Track:hTracklet2TRDin:hCluster2MC:hTracklet2MC:hTRDin2MC:hTRD2MC"
};
Char_t taskName[ntasks][100]={{""}};
TObjArray *plots[ntasks] = {NULL};
TObjArray *arr(NULL), *narr[ntasks] = {NULL};
TFile *file(NULL), *fOut(NULL);
Bool_t mc(kFALSE), 
       kWrite(kFALSE),
       kVerbose(kTRUE);
void addPerformance(const Int_t it);
void mergeProdLocal(const Char_t *list, Bool_t doMC=kFALSE, const Int_t nwrite=10000)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTender.so");
  //gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGLFspectra.so");

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
  fOut->mkdir(nameTask)->Write();

  // read task descriptors a
  for(Int_t it(0); it<ntasks; it++){
    s=taskConfig[it]; 
    s1=s(0, s.Index('_')); sprintf(taskName[it], "%s", s1.Data());
    s = s(s.Index('_')+1, 200); plots[it] = s.Tokenize(":");
  }
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
    if(!file->cd(nameTask)){ 
      Warning("mergeProd", "Missing %s. Skip file.", nameTask);
      continue;
    }
      
    if((ifile%nwrite)==(nwrite-1)){ 
      kWrite = kTRUE;
      if(!kVerbose) printf("\nMerging ");
      else printf(" ... Flushing merged file ...");
    }
    
    for(Int_t itask(0); itask<ntasks; itask++){
      file->cd(nameTask);
      if(!(arr = (TObjArray*)gDirectory->Get(taskName[itask]))){
        Warning("mergeEntry", "Missing %s. Skip task.", taskName[itask]);
        continue;
      }
      fOut->cd(nameTask);
      addPerformance(itask);
      arr->Delete(); delete arr;
    }
    file->Close(); delete file;
    ifile++; kWrite = kFALSE;
    //if(ifile>=35) break; 
  }
  
  if(!kVerbose) printf("\n");
  Info("mergeProd", "Flushing merged file ...");
  fOut->cd(nameTask);
  for(Int_t it(0); it<ntasks; it++){
    //if(!(arr = (TObjArray*)gDirectory->Get(taskName[it]))) continue;
    //arr->Write(taskName[it], TObject::kSingleKey);
    narr[it]->Write(taskName[it], TObject::kSingleKey);
  }  
  fOut->Close();
}

void addPerformance(const Int_t it)
{ 
  if(!narr[it]){
    narr[it] = (TObjArray*)arr->Clone(); narr[it]->SetOwner();
    gDirectory->ls();
    return;
  }

  TObject *oo(NULL);
  TH1 *h1(NULL), *H1(NULL);  
  THnSparse *h(NULL), *H(NULL);
  Int_t nplots(plots[it]->GetEntries());
  TString sname;
  for(Int_t ih(0); ih<nplots; ih++){
    sname= ((TObjString*)(*plots[it])[ih])->String();
    if(!mc && sname.EndsWith("MC")) continue;

    if(!(oo = (TObject*)arr->FindObject(sname.Data()))) {
      Warning("addPerformance", "Missing %s from \"%s\"", sname.Data(), file->GetName());
      continue;
    }
    if(strstr(oo->IsA()->GetName(), "THnSparse")){ 
      h = dynamic_cast<THnSparse*>(oo);
      if(!(H = (THnSparse*)narr[it]->FindObject(sname.Data()))) {
        Error("addPerformance", "Missing %s from Merging.", sname.Data());
        return;
      }
      H->Add(h);
      if(kWrite) printf("  -> %s[%d]\n", sname.Data(), (Int_t)H->GetEntries());
    } else if(strstr(oo->IsA()->GetName(), "TH")){ 
      h1 = dynamic_cast<TH1*>(oo);
      if(!(H1 = (TH1*)narr[it]->FindObject(sname.Data()))) {
        Error("addPerformance", "Missing %s. from Merging.", sname.Data());
        return;
      }
      H1->Add(h1);
      if(kWrite) printf("  -> %s[%d]\n", sname.Data(), (Int_t)H1->GetEntries());
    } else{
      Warning("addPerformance", "Object %s. from \"%s\" of unknown class %s", sname.Data(), file->GetName(), oo->IsA()->GetName());
      continue;
    }
  }
  if(kWrite){
    narr[it]->Write(taskName[it], /*TObject::kOverwrite|*/TObject::kSingleKey);
    narr[it]->Delete(); delete narr[it]; narr[it] = NULL;
  }
}

void finalMerge(const Char_t *fn, Int_t ih=0)
{
  TFile *fIn(NULL);
  if(!(fIn=TFile::Open(fn))) return;
  if(!gFile->cd(nameTask)){
    Error("finalMerge", "message");
    return;
  }
  TIterator *i=gDirectory->GetListOfKeys()->MakeIterator();
  THnSparse *H1(NULL);

  TFile *fOut=TFile::Open(Form("Merged_%d.root", ih), "RECREATE");

  TKey *k(NULL); 
  TObjArray *arr(NULL);
  while(k=(TKey*)i->Next()){
    //printf("%s;%d\n", k->GetName(), k->GetCycle());
    fIn->cd(nameTask);
    arr=(TObjArray*)gDirectory->Get(Form("%s;%d", k->GetName(), k->GetCycle()));
    printf("Processing %s;%d ...\n", k->GetName(), k->GetCycle());
    fOut->cd();
    if(H1) H1->Add((THnSparse*)arr->At(ih));
    else H1=(THnSparse*)((THnSparse*)arr->At(ih))->Clone();
    arr->Delete(); delete arr;
  }
  fOut->cd();
  if(H1) H1->Write();
  fOut->Close();
}