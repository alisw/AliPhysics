#include "UtilFiles.h"
#include "UtilMessages.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TSystem.h"
#include <iostream>

using namespace std;

//_________________________________________________________________________________________________
TFile * GetFile(TString fname, TString opt){
  TFile * f = 0x0;
  if(fname.IsNull() || fname.IsWhitespace()) Fatalmsg("GetFile", "Provided empty string");
  f = TFile::Open(fname, opt);
  if(!f || !f->IsOpen()){
    TString path = "./";
    path += fname;
    const Int_t last = path.Last('/');
    if(last == 1) path = "./";
    else while(path.Last('/') == last && path.Sizeof() > 1) path.Chop();
    cout<<"In path:"<<endl;
    gSystem->Exec(Form("ls %s", path.Data()));
    Fatalmsg("GetFile", Form("Cannot find File %s", fname.Data()));
  }

  return f;
}

//_________________________________________________________________________________________________
void GetListFromFile(TFile *fin, TString name, TList *& lin){
  fin->GetObject(name, lin);
  if(!lin){
    fin->ls();
    Fatalmsg("GetListFromFile", Form("Cannot find TList %s from file %s", name.Data(), fin->GetName()));
    
  }
}

//_________________________________________________________________________________________________
void GetListFromDirectory(TDirectory *dir, TString name, TList *& lin){
  dir->GetObject(name, lin);
  if(!lin){
    dir->ls();
    Fatalmsg("GetListFromDirectory", Form("Cannot find TList %s from file %s", name.Data(), dir->GetName()));
    
  }
}

//_________________________________________________________________________________________________
void GetDirectoryFromFile(TFile *fin, TString name, TDirectory *& dir){
  fin->GetObject(name, dir);
  if(!dir){
    fin->ls();
    Fatalmsg("GetDirectoryFromFile", Form("Cannot find TDirectory %s from file %s", name.Data(), fin->GetName()));
    
  }
}

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH1F *& histo){//Gets the histogram from the list
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  histo = (TH1F *) (lin->FindObject(hname));
  if(!histo){
    lin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH1 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH1D *& histo){//Gets the histogram from the list
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  histo = (TH1D *) (lin->FindObject(hname));
  if(!histo){
    lin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH1 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH1F *& histo){//Gets the histogram from the file
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  fin->GetObject(hname, histo);
  if(!histo){
    fin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH1 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH1D *& histo){//Gets the histogram from the file
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  fin->GetObject(hname, histo);
  if(!histo){
    fin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH1 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH2I *& histo){//Gets the histogram from the file
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  fin->GetObject(hname, histo);
  if(!histo){
    histo = (TH2I*) fin->Get(hname);
    if(histo){
      Warningmsg("GetHistogram", Form("Requested histogram class is not found! Casting %s to TH2I !!!", fin->Get(hname)->ClassName()));
      return;
    }
    fin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH2 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH2F *& histo){//Gets the histogram from the file
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  fin->GetObject(hname, histo);
  if(!histo){
    fin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH2 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH2F *& histo){//Gets the histogram from the list
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  histo = (TH2F *) (lin->FindObject(hname));
  if(!histo){
    lin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH2 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH2I *& histo){//Gets the histogram from the list
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  histo = (TH2I *) (lin->FindObject(hname));
  if(!histo){
    lin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH2 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH3I *& histo){//Gets the histogram from the list
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  histo = (TH3I *) (lin->FindObject(hname));
  if(!histo){
    lin->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH2 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
void GetHistogram(TDirectory *dir, const TString hname, TH1F *& histo){//Gets the histogram from the directory
  //   if(histo) Fatalmsg("GetHistogram", Form("Histogram %s already exists!", histo->GetName()));
  dir->GetObject(hname, histo);
  if(!histo){
    dir->ls();
    Fatalmsg("GetHistogram", Form("Cannot find TH1 %s !", hname.Data()));
  }
  
}

//_________________________________________________________________________________________________
TList *ReduceList(TList *lin, const TString criteria){//Macro to produce multiple lists from one -> Useful for writing to file
  TList *result = new TList();
  result->SetOwner();
  
  TIter next(lin);
  TObject *nextobj;
  TString objname, objtitle, objclass;
  
  while ((nextobj = next())) {
    objname = nextobj->GetName();
    objtitle = nextobj->GetTitle();
    objclass = nextobj->ClassName();
    if(!objname.Contains(criteria)) continue;
    
    result->Add(nextobj);
    lin->Remove(nextobj);
    
  }
  
  return result;
  
}

//_________________________________________________________________________________________________
TList *FormListFromFile(TFile *fin, const TString criteria, const TString checklists){
  TList *result = new TList();
  result->SetOwner();
  
  TList *lkeys = fin->GetListOfKeys();
  const Int_t max = lkeys->GetEntries();
  Int_t counter = 0;
  TIter next(lkeys);
  TKey *key;
  while ((key = (TKey*)next())) {
    TString objclass  = key->GetClassName();
    TString objname  = key->GetName();
    counter++;
    if(counter > max) break;
    if(!checklists.IsNull() && objclass.Contains("TList") && objname.Contains(checklists)){
      TList * lin = (TList*)key->ReadObj();
      TIter nextinlist(lin);
      TObject *nextobjinlist;
      TString objnameinlist, objtitleinlist, objclassinlist;
      
      while ((nextobjinlist = nextinlist())) {
        objnameinlist = nextobjinlist->GetName();
        objtitleinlist = nextobjinlist->GetTitle();
        objclassinlist = nextobjinlist->ClassName();
        if(!objnameinlist.Contains(criteria)) continue;
        
        if(!criteria.IsNull() && !objnameinlist.Contains(criteria)) continue;
        if(!objclassinlist.Contains("TH1")) continue;
        TH1F *h = (TH1F*)nextobjinlist->Clone();
        h->SetDirectory(0);
        result->Add(static_cast<TH1F*>(h));
        //     PrintProgress(counter, max);
        if(!criteria.IsNull()) Infomsg("FormListFromFile", Form("%i/%i %s", counter, max, objnameinlist.Data()));
        
      }
      
    }
    else{
      if(!criteria.IsNull() && !objname.Contains(criteria)) continue;
      if(!objclass.Contains("TH1")) continue;
      TH1F *h = (TH1F*)key->ReadObj();
      h->SetDirectory(0);
      result->Add(static_cast<TH1F*>(h));
      //     PrintProgress(counter, max);
      if(!criteria.IsNull()) Infomsg("FormListFromFile", Form("%i/%i %s", counter, max, objname.Data()));
    }
  }
  return result;
}
  
