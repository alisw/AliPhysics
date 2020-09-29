#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include "AliPID.h"
#endif

TH1F* CreateHisto(TString nam, Int_t tote);
void SetupPadStyle();
void DoMakeUp(TH1* h, const TString det, const TString part, const TString partL);//Function to perform style the drawn TH1

Bool_t DrawTrendingPIDQA(TString mergedTrendFile = "trending.root"){
  TFile *fin = TFile::Open(mergedTrendFile.Data());
  if(!fin){
    printf("Cannot open file with PID QA trending: %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  TTree * ttree = (TTree*) fin->Get("trending");
  if (!ttree){
    printf("Trending tree not found in file %s\n",mergedTrendFile.Data());
    return kFALSE;
  }
  Int_t nrun;
  ttree->SetBranchAddress("nrun",&nrun);
  TObjArray* lb=(TObjArray*)ttree->GetListOfBranches();
  Int_t nVars=lb->GetEntries();
  Float_t* vect=new Float_t[nVars-1];
  TH1F** htr=new TH1F*[nVars-1];
  Int_t kv=0;
  Int_t totEnt=ttree->GetEntries();
  for(Int_t j=0; j<nVars; j++){
    TBranch* br=(TBranch*)lb->At(j);
    printf("Branch %d  %s\n",j,br->GetName());
    TString bnam=br->GetName();
    if(!bnam.Contains("nrun")){
      ttree->SetBranchAddress(bnam,&vect[kv]);
      htr[kv]=CreateHisto(bnam,totEnt);
      ++kv;
    }
  }
  for(Int_t je=0; je<totEnt; je++){
    ttree->GetEvent(je);
    printf(" Run %d\n",nrun);
    for(Int_t k=0; k<nVars-1; k++){
      htr[k]->SetBinContent(je+1,vect[k]);
      htr[k]->GetXaxis()->SetBinLabel(je+1,Form("%d",nrun));
    }
  }


  TCanvas** c=new TCanvas*[20];
  TString detStrings[4]={"ITS","TPC_Basic","TPC_TOF","TOF"};
  TString partStrings[5]={"electron","pion","kaon","proton","deuteron"};
  const TString partStringslatex[5]={"e^{#pm}","#pi^{#pm}","K^{#pm}","p, #bar{p}","d, #bar{d}"};
  TString canname[20];
  for(Int_t k=0; k<20; k++){ 
    c[k]=new TCanvas(Form("c%d",k),Form("c%d",k),1500,600);
    c[k]->Divide(2,2);
    canname[k]="PIDqaTrend";
  }
  for(Int_t jdet=0; jdet<4; jdet++){
    for(Int_t jpar=0; jpar<5; jpar++){
      Int_t k=jdet*5+jpar;
      canname[k].Append(detStrings[jdet].Data());
      canname[k].Append(partStrings[jpar].Data());
    }
  }

  for(Int_t jvar=0; jvar<nVars-1; jvar++){
    TString hname=htr[jvar]->GetName();
    Int_t theDet=-1;
    for(Int_t jdet=0; jdet<4; jdet++){
      if(hname.Contains(Form("nSigma%s",detStrings[jdet].Data()))) theDet=jdet;
    }
    Int_t thePar=-1;
    for(Int_t jpar=0; jpar<5; jpar++){     
      if(hname.Contains(partStrings[jpar].Data())) thePar=jpar;
    }
    if(theDet>=0 && thePar>=0){
      Int_t theCanv=theDet*5+thePar;
      Int_t thePad=-1;
      if(hname.Contains("meannSigma")) thePad=1;
      else if(hname.Contains("signSigma")) thePad=2;
      if(thePad>0){
	TObjArray* arr=hname.Tokenize("_");
	TObjString* lasts=(TObjString*)arr->At(arr->GetEntries()-1);
	TString stmom=lasts->GetString();
	if(stmom.Contains("MeV")){
	  stmom.ReplaceAll("p","");
	  stmom.ReplaceAll("MeV","");
	  Int_t mom=stmom.Atoi();
	  if(theDet<=2 && mom>800) thePad+=2;
	  if(theDet==3 && mom>1100) thePad+=2;
	  c[theCanv]->cd(thePad);
          DoMakeUp(htr[jvar]->DrawCopy(""), detStrings[theDet], partStrings[thePar], partStringslatex[thePar]);
          htr[jvar]->Draw("psame");
        }
      }
    }
  }
  for(Int_t k=0; k<20; k++) c[k]->SaveAs(Form("%s.png",canname[k].Data()));
  delete [] vect;

  TDirectoryFile* pidtr=(TDirectoryFile*)fin->Get("PIDinTr");
  if(pidtr && pidtr->GetNkeys()>0){
    TCanvas* ctrpid=new TCanvas("ctrpid","PID-in-track",1000,700);
    TCanvas* ctrpidall=new TCanvas("ctrpidall","PID-in-track",1500,700);
    TLegend* leg=new TLegend(0.7,0.35,0.89,0.87);
    leg->SetHeader("PID in tracking");
    ctrpidall->Divide(3,3);
    Int_t cols[9]={kGreen+2,kGray,1,2,4,kMagenta,kOrange+1,kYellow,kCyan};
    Bool_t drawLeg=kFALSE;
    for(Int_t j=0; j<9; j++){
      TString histoname=Form("hSigP_TPC_TrackedAs_%sforMerge",AliPID::ParticleName(j));
      TH2* histo = (TH2*)pidtr->Get(histoname.Data());
      if(histo){
	drawLeg=kTRUE;
	TH2* histo2=(TH2*)histo->Clone(Form("%scolor",histoname.Data()));
	histo2->SetTitle(" ");
	histo2->SetStats(0);
	histo2->SetMarkerColor(cols[j]);
	leg->AddEntry(histo2,Form("%s",AliPID::ParticleName(j)),"")->SetTextColor(histo2->GetMarkerColor());
	ctrpid->cd();
	SetupPadStyle();
	if(j==0) histo2->Draw();
	else histo2->Draw("same");
	ctrpidall->cd(j+1);
	SetupPadStyle();
	histo->Draw("colz");
	ctrpidall->Update();
      }
    }
    ctrpid->cd();
    if(drawLeg) leg->Draw();
    ctrpidall->SaveAs("TPCdEdx-PIDinTracking-9pads-AllRuns.png");
    ctrpid->SaveAs("TPCdEdx-PIDinTracking-AllRuns.png");
  }
  return kTRUE;
}


TH1F* CreateHisto(TString nam, Int_t tote){
  TH1F* h=new TH1F(nam.Data(),Form(" ; run ; %s",nam.Data()),tote,0.5,tote+0.5);
  h->SetStats(0);
  h->SetMarkerStyle(20);
  return h;
}

void SetupPadStyle()
{
  gPad->SetLogx();
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();
}

TGraph* ShowLimit(const TH1* h) //Function to show if a datapoint out of range is up or down
{
  TGraph* g = 0x0;
  const Double_t lim[2] = { h->GetMinimum(), h->GetMaximum() };
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    const Double_t y = h->GetBinContent(i);
    if (y < lim[1] && y > lim[0])
      continue;
    const Double_t x = h->GetXaxis()->GetBinCenter(i);
    if (!g) {
      g = new TGraph();
      g->SetMarkerColor(kRed);
      g->SetMarkerStyle(27);
    }
    if (y > lim[1])
      g->SetPoint(g->GetN(), x, lim[1]);
    else
      g->SetPoint(g->GetN(), x, lim[0]);
  }
  if (g)
    g->Draw("psame");
  return g;
}

void DoMakeUp(TH1* h, const TString det, const TString part, const TString partL)
{
  const TString hname = h->GetName();
  if (det.EqualTo("TOF")) { //Setting histogram range for TOF plots + styling
    if (hname.Contains("meannSigma"))
      h->GetYaxis()->SetRangeUser(-.5, .5);
    else if (hname.Contains("signSigma"))
      h->GetYaxis()->SetRangeUser(.5, 1.5);
    h->GetYaxis()->SetTitleSize(0.09);
    h->GetYaxis()->SetTitleOffset(0.6);
    h->GetYaxis()->SetLabelSize(0.07);
    TString title = h->GetYaxis()->GetTitle();
    title.ReplaceAll("_", " ");
    title.ReplaceAll("meannSigma" + det, "#mu n#sigma" + det);
    title.ReplaceAll("signSigma" + det, "#sigma n#sigma" + det);
    title.ReplaceAll(" " + part, "_{" + partL + "}");
    title.ReplaceAll("p1000MeV", "#it{p} 1.0 GeV/#it{c}");
    title.ReplaceAll("p2000MeV", "#it{p} 2.0 GeV/#it{c}");
    h->GetYaxis()->SetTitle(title);
    ShowLimit(h);
  }
  if (det.EqualTo("ITS")) { //Setting histogram range for TOF plots + styling
    if (hname.Contains("meannSigma"))
      h->GetYaxis()->SetRangeUser(-1.2, 1.2);
    else if (hname.Contains("signSigma"))
      h->GetYaxis()->SetRangeUser(.5, 1.5);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetYaxis()->SetTitleOffset(0.5);
    h->GetYaxis()->SetLabelSize(0.07);
    TString title = h->GetYaxis()->GetTitle();
    title.ReplaceAll("_", " ");
    title.ReplaceAll("meannSigma" + det, "#mu n#sigma" + det);
    title.ReplaceAll("signSigma" + det, "#sigma n#sigma" + det);
    title.ReplaceAll(" " + part, "_{" + partL + "}");
    title.ReplaceAll("p300MeV", "#it{p} 0.3 GeV/#it{c}");
    title.ReplaceAll("p1200MeV", "#it{p} 1.2 GeV/#it{c}");
    h->GetYaxis()->SetTitle(title);
    ShowLimit(h);
  }
}
