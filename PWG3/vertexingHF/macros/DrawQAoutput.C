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

  /*
  TCanvas* citscls=new TCanvas("citscls","ITS clusters");
  citscls->SetGridy();
  TCanvas* cspdcls=new TCanvas("cspdcls","SPD clusters");
  cspdcls->SetGridy();
  TCanvas* cngoodtr=new TCanvas("cngoodtr","Good tracks per event");
  TCanvas* cptgoodtr=new TCanvas("cptgoodtr","Pt distribution good tracks");
  TCanvas* cd0=new TCanvas("cd0","Impact parameter distribution (track selected)");



  TString hname="hnClsITS";
  TH1F* hnClsITS=(TH1F*)list->FindObject(hname);
  if(!hnClsITS){
    cout<<hname<<" not found"<<endl;
  }else{
    citscls->cd();
    hnClsITS->Draw();
    citscls->SaveAs(Form("%s%s.png",hname.Data(),textleg.Data()));
  }

  hname="hnClsSPD";
  TH1F* hnClsSPD=(TH1F*)list->FindObject(hname);
  if(!hnClsSPD){
    cout<<hname<<" not found"<<endl;
  }else{
    cspdcls->cd();
    hnClsSPD->Draw();
    cspdcls->SaveAs(Form("%s%s.png",hname.Data(),textleg.Data()));
  }

  hname="hdistrGoodTr";
  TH1F* hngoodtr=(TH1F*)list->FindObject(hname);
  if(!hngoodtr){
    cout<<hname<<" not found"<<endl;
  }else{
    cngoodtr->cd();
    hngoodtr->Draw();
    cngoodtr->SaveAs(Form("%s%s.png",hname.Data(),textleg.Data()));
  }

  hname="hptGoodTr";
  TH1F* hptgoodtr=(TH1F*)list->FindObject(hname);
  if(!hptgoodtr){
    cout<<hname<<" not found"<<endl;
  }else{
    cptgoodtr->cd();
    hptgoodtr->Draw();
    cptgoodtr->SaveAs(Form("%s%s.png",hname.Data(),textleg.Data()));
  }

  hname="hd0";
  TH1F* hd0=(TH1F*)list->FindObject(hname);
  if(!hd0){
    cout<<hname<<" not found"<<endl;
  }else{
    cd0->cd();
    hd0->Draw();
    cd0->SaveAs(Form("%s%s.png",hname.Data(),textleg.Data()));
  }
  */
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
    TH1F* h=(TH1F*)list->At(i);
    if(!h){
      cout<<"Histogram "<<i<<" not found"<<endl;
      continue;
    }
    TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
    c->cd();
    if(i<3)h->Draw();
    else h->Draw("colz");
  }

}
