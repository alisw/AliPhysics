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
      c->SaveAs(Form("%s%s.png",h->GetName(),textleg.Data()));
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
      h->Sumw2();
      h->Scale(1./h->Integral("width"));
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      c->SetLogz();
      //c->SetLogx();
      c->cd();
      
      h->Draw("colz");

      //mean and pull, code from Jens Wiechula
      TF1 fg("fg","gaus",-2.,2.); // fit range +- 2 sigma
      TLine l;
      TObjArray arr;

      h->Draw("colz");
      fg.SetParameters(1,0,1);
      h->FitSlicesY(&fg,0,-1,0,"NQR",&arr);

      TH1 *hM=(TH1*)arr.At(1);
      hM->SetMarkerStyle(20);
      hM->SetMarkerSize(.5);
      hM->Draw("same");

      TH1 *hS=(TH1*)arr.At(2);
      hS->SetMarkerStyle(20);
      hS->SetMarkerSize(.5);
      hS->SetMarkerColor(kRed);
      hS->SetLineColor(kRed);
      hS->Draw("same");

      l.SetLineColor(kBlack);
      l.DrawLine(.2,0,20,0);
      l.SetLineColor(kRed);
      l.DrawLine(.2,1,20,1);

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

  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  Int_t nevents=hstat->Integral(1,2);
  hstat->Draw("htext0");
  cst->SaveAs(Form("%s%s.png",hstat->GetName(),textleg.Data()));
  Int_t nevents080=1;

  //TCanvas *spare=new TCanvas("sparecv","Spare");

  for(Int_t i=0;i<list->GetEntries();i++){

    TClass* objtype=list->At(i)->IsA();
    TString tpname=objtype->GetName();

    if(tpname=="TH1F"){

      TH1F* h=(TH1F*)list->At(i);
      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
      pvtxt->SetBorderSize(0);
      pvtxt->SetFillStyle(0);

      c->cd();
      c->SetGrid();
      Int_t entries=h->Integral();
      pvtxt->AddText(Form("%.1f per cent of the events",(Double_t)entries/(Double_t)nevents));
      h->Draw();
      pvtxt->Draw();
      c->SaveAs(Form("%s%s.png",c->GetName(),textleg.Data()));
    }
    if(tpname=="TH2F"){
      TH2F* h=(TH2F*)list->At(i);
      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      TCanvas* c=new TCanvas(Form("c%s",h->GetName()),h->GetName());
      TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
      pvtxt->SetBorderSize(0);
      pvtxt->SetFillStyle(0);

      c->cd();
      c->SetGrid();
      Int_t entries=h->Integral();
      pvtxt->AddText(Form("%.1f per cent of the events",(Double_t)entries/(Double_t)nevents));
      h->Draw("colz");
      c->SetLogz();
      pvtxt->Draw();
      c->SaveAs(Form("%s%s.png",c->GetName(),textleg.Data()));
    }
  }
  
  
  listname="countersCentrality";

  isRead=ReadFile(list,hstat,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }

  TH1F* hallcntr=0x0;
  cout<<"normalizing to 0-80% as a check"<<endl;

  TCanvas *cvnocnt=new TCanvas("cvnocnt","No Centrality estimation",800,400);
  cvnocnt->Divide(2,1);

  for(Int_t i=0;i<list->GetEntries();i++){
    AliCounterCollection* coll=(AliCounterCollection*)list->At(i);
    
    coll->SortRubric("run");//sort by run number
    //coll->PrintKeyWords();
    Int_t ncentr=10;//check this
    TH1F* h020=0x0;
    TH1F* h2080=0x0;
    hallcntr=0x0; 

    TH1F* hbad=(TH1F*)coll->Get("run",Form("centralityclass:-990_-980"));
    cvnocnt->cd(i+1);
    if(hbad) hbad->Draw();

    TCanvas *ccent=new TCanvas(Form("ccent%s",coll->GetName()),Form("Centrality vs Run (%s)",coll->GetName()),1400,800);
    ccent->Divide(5,3);

    for(Int_t ic=0;ic<8/*ncentr*/;ic++){ //normalizing to 0-80% as a check

      TH1F* h=(TH1F*)coll->Get("run",Form("centralityclass:%d_%d",ic*10,ic*10+10));
      h->SetName(Form("h%d%d",i,ic));
      if(!hallcntr) {
	hallcntr=(TH1F*)h->Clone("hallcntr");
	hallcntr->Sumw2();
      }
      else hallcntr->Add(h);
      nevents080+=h->Integral();
    }

    for(Int_t ic=0;ic<ncentr;ic++){

      TH1F* h=(TH1F*)coll->Get("run",Form("centralityclass:%d_%d",ic*10,ic*10+10));
      h->SetName(Form("h%d%d",i,ic));
      h->Sumw2();

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

      h->Divide(hallcntr);

      ccent->cd(ic+1);
      h->GetYaxis()->SetRangeUser(0.,0.15);
      h->DrawClone();
      /*
      if(ic==0&&i==0){
	spare->cd();
	h->Draw();
      }
      */
      // ccent->cd(1);
      // h->SetLineColor(ic+1);
      // if(ic==0)h->DrawClone();
      // else h->DrawClone("sames");
    }

    ccent->cd(ncentr+1);
    h020->Divide(hallcntr);
    h020->DrawClone();
    TCanvas* cv020=new TCanvas(Form("cv020-%d",i),"0-20% vs run number");
    cv020->cd();
    h020->GetYaxis()->SetRangeUser(0.,1.);
    h020->DrawClone();
    cv020->SaveAs(Form("cv020-%d.png",i));

    ccent->cd(ncentr+2);
    h2080->Divide(hallcntr);
    h2080->DrawClone();

    TCanvas* cv2080=new TCanvas(Form("cv2080-%d",i),"20-80% vs run number");
    cv2080->cd();
    h2080->GetYaxis()->SetRangeUser(0.,1.);
    h2080->DrawClone();
    cv2080->SaveAs(Form("cv2080-%d.png",i));

    ccent->SaveAs(Form("%s%s.png",ccent->GetName(),textleg.Data()));
  }

}

void DrawProjections(TString partname="D0",TString h2dname="hMultvsPercentile",Int_t nsteps=0,TString direction="X",TString path="./"){
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  TString listname="outputCentrCheck";

  TList* list;
  TH1F * hstat;

  Bool_t isRead=ReadFile(list,hstat,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat){
    cout<<":-( null pointers..."<<endl;
    return;
  }

  TH2F* h2=(TH2F*)list->FindObject(h2dname);
  if(!h2){
    cout<<h2dname.Data()<<" not found"<<endl;
    return;
  }
  Int_t kbins=1;
  if(nsteps==0){
    if(direction=="X") nsteps=h2->GetNbinsY();
    if(direction=="Y") nsteps=h2->GetNbinsX();
  }
  else{
    if(direction=="X") kbins=h2->GetNbinsY()/nsteps;
    if(direction=="Y") kbins=h2->GetNbinsX()/nsteps;
  }

  TCanvas *cvpj=new TCanvas(Form("cvpj%s%s",direction.Data(),h2dname.Data()),Form("cvpj%s",direction.Data()),1000,800);
  cvpj->Divide((Int_t)(nsteps/3)+1,3);
  TFile* fout=new TFile(Form("proj%s%s.root",direction.Data(),h2dname.Data()), "recreate");
  for(Int_t i=0;i<nsteps;i++){
    TH1F* h=0x0;
    if(direction=="X")h=(TH1F*)h2->ProjectionX(Form("px%d",i),i+kbins,i+2*kbins);
    if(direction=="Y")h=(TH1F*)h2->ProjectionY(Form("py%d",i),i+kbins,i+2*kbins);

    TPaveText *pvtxt=new TPaveText(0.6,0.6,0.9,0.9,"NDC");
    pvtxt->SetBorderSize(0);
    pvtxt->SetFillStyle(0);
    pvtxt->AddText(Form("%d - %d",((i+kbins-1)*10),(i+2*kbins-1)*10));

    cvpj->cd(i+1);
    h->Draw();
    pvtxt->Draw();
    fout->cd();
    h->Write();
  }
  cvpj->SaveAs(Form("cvpj%s%s.png",direction.Data(),h2dname.Data()));

}
