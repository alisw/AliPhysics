#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TClass.h>
#include <AliCounterCollection.h>

//read the file and take list and stat

Bool_t ReadFile(TList* &list,TH1F* &hstat,TString listname,TString partname,TString path="./",TString filename="AnalysisResults.root");

Bool_t ReadFile(TList* &list,TH1F* &hstat,TString listname,TString partname,TString path,TString filename){

  TString hstatname="nEntriesQA",dirname="PWG3_D2H_QA";
  filename.Prepend(path);
  listname+=partname;
  hstatname+=partname;

  TFile* f=new TFile(filename.Data());
  if(!f){
    cout<<filename.Data()<<" not found"<<endl;
    return kFALSE;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname);
  if(!f){
    cout<<dirname.Data()<<" not found  in "<<filename.Data()<<endl;
    return kFALSE;
  }

  list=(TList*)dir->Get(listname);
  if(!list){
    cout<<"List "<<listname.Data()<<" not found"<<endl;
    dir->ls();
    return kFALSE;
  }

  hstat=(TH1F*)dir->Get(hstatname);
  if(!hstat){
    cout<<hstatname.Data()<<" not found"<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//draw "track related" histograms (list "outputTrack")
void DrawOutputTrack(TString partname="D0",TString textleg="",TString path="./"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputTrack";

  TList* list;
  TH1F * hstat;

  Bool_t isRead=ReadFile(list,hstat,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }

  for(Int_t i=0;i<list->GetEntries();i++){
    TH1F* h=(TH1F*)list->At(i);
    if(!h){
      cout<<"Histogram "<<i<<" not found"<<endl;
      continue;
    }
    TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
    c->cd();
    c->SetGrid();
    TString hname=h->GetName();
    if(!hname.Contains("nCls")){
      c->SetLogy();
      h->Draw();
    } else h->Draw("htext0");
    c->SaveAs(Form("%s%s.png",c->GetName(),textleg.Data()));
  }
  
  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  hstat->Draw("htext0");
  cst->SaveAs(Form("%s%s.png",hstat->GetName(),textleg.Data()));


}

//draw "pid related" histograms (list "outputPID")
void DrawOutputPID(TString partname="D0",TString textleg="",TString path="./"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputPid";

  TList* list;
  TH1F * hstat;

  Bool_t isRead=ReadFile(list,hstat,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }

  for(Int_t i=0;i<list->GetEntries();i++){
    TClass* objtype=list->At(i)->IsA();
    TString tpname=objtype->GetName();

    if(tpname=="TH1F"){
      TH1F* h=(TH1F*)list->At(i);

      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      //h->Scale(1./h->Integral("width"));
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      c->SetLogz();
      c->cd();
      h->Draw();
    
      //write
      c->SaveAs(Form("%s.png",h->GetName()));
      TFile* fout=new TFile(Form("%s.root",h->GetName()),"recreate");
      fout->cd();
      c->Write();
    }
  
    if(tpname=="TH2F"){
      TH2F* h=(TH2F*)list->At(i);
      
      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      h->Scale(1./h->Integral("width"));
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      c->SetLogz();
      c->cd();
      
      h->Draw("colz");
      //write
      c->SaveAs(Form("%s.png",h->GetName()));
      TFile* fout=new TFile(Form("%s.root",h->GetName()),"recreate");
      fout->cd();
      c->Write();
    }
  }
}

void DrawOutputCentrality(TString partname="D0",TString textleg="",TString path="./"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);

  TString listname="outputCentrCheck";

  TList* list;
  TH1F * hstat;

  Bool_t isRead=ReadFile(list,hstat,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  
  for(Int_t i=0;i<list->GetEntries();i++){
    TH1F* h=(TH1F*)list->At(i);
    if(!h){
      cout<<"Histogram "<<i<<" not found"<<endl;
      continue;
    }
    TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
    c->cd();
    c->SetGrid();
    h->Draw();
    c->SaveAs(Form("%s%s.png",c->GetName(),textleg.Data()));
  }
  
  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  hstat->Draw("htext0");
  cst->SaveAs(Form("%s%s.png",hstat->GetName(),textleg.Data()));
  
  listname="countersCentrality";

  isRead=ReadFile(list,hstat,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  for(Int_t i=0;i<list->GetEntries();i++){
    AliCounterCollection* coll=(AliCounterCollection*)list->At(i);
    
    coll->SortRubric("run");//sort by run number
    //coll->PrintKeyWords();
    Int_t ncentr=10;//check this
    TH1F* h020=0x0;
    TH1F* h2080=0x0;
    TCanvas *ccent=new TCanvas(Form("ccent%s",coll->GetName()),Form("Centrality vs Run (%s)",coll->GetName()),1400,800);
    ccent->Divide(5,3);

    for(Int_t ic=0;ic<ncentr;ic++){

      TH1F* h=(TH1F*)coll->Get("run",Form("centralityclass:%d_%d",ic*10,ic*10+10));
      h->SetName(Form("h%d%d",i,ic));
      if(ic>=0 && ic<=1){ //0-20
	if(!h020) {
	  h020=(TH1F*)h->Clone(Form("h020%s",coll->GetName()));
	  h020->SetTitle(Form("Centrality 0-20 %s",coll->GetName()));
	}
	else h020->Add(h);
      }
      if(ic>=2 && ic<=7){ //20-80
	if(!h2080) {
	  h2080=(TH1F*)h->Clone(Form("h2080%s",coll->GetName()));
	  h2080->SetTitle(Form("Centrality 20-80 %s",coll->GetName()));
	}
	else h2080->Add(h);
      }

      ccent->cd(ic+1);
      h->DrawClone();
      // ccent->cd(1);
      // h->SetLineColor(ic+1);
      // if(ic==0)h->DrawClone();
      // else h->DrawClone("sames");
    }

    ccent->cd(ncentr+1);
    h020->DrawClone();

    ccent->cd(ncentr+2);
    h2080->DrawClone();

    ccent->SaveAs(Form("%s%s.png",ccent->GetName(),textleg.Data()));
  }

}
