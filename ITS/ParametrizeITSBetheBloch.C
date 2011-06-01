#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include "AliITSPIDResponse.h"
#endif

/*  $Id$    */


// Macro to extract the parameterizations of dE/dx in ITS vs. p
// using Phobos BB parameterization


Double_t BetheBlochPhobos(const Double_t *x, const Double_t *par);
Double_t BetheBlochPhobosBG(const Double_t *x, const Double_t *par);
Double_t BetheBlochAleph(const Double_t *x, const Double_t *par);


void ParametrizeITSBetheBloch(Bool_t optFromITSsa=kFALSE){

  TFile *fil=new TFile("AnalysisResults.root");
  TH2F *hsigvsp=0x0;
  TString dfName="PWG2SpectraITSsa";
  TString listName="clistITSsaMult-1to-1";
  TString hisName="fHistDEDX";

  if(optFromITSsa){
    dfName="ITSsaTracks";
    listName="clistITSsaTracks";
    hisName="hdedxvsPt4clsITSpureSA";
  }
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dfName.Data());
  if(!df){
    printf("ERROR: directory file %s not found\n",dfName.Data());
    return;
  }
  TList* l=(TList*)df->Get(listName.Data());
  if(!l){
    printf("ERROR: TList %s not found\n",listName.Data());
    return;
  }
  hsigvsp=(TH2F*)l->FindObject(hisName.Data());
  if(!hsigvsp){ 
    printf("ERROR: historgam %s not found\n",hisName.Data());
    return;
  }
  hsigvsp->GetXaxis()->SetTitle("momentum (GeV/c)");
  hsigvsp->GetYaxis()->SetTitle("dE/dx (keV/300 #mum)");

  Double_t minpPionFit=0.1;
  Double_t maxpPionFit=0.8;
  Double_t minbgGlobFit=0.3;
  Double_t maxbgGlobFit=2.5;
  Double_t minbgPions=2.;  
  Double_t maxbgPions=10.;  
  Double_t minpKaons=0.2;
  Double_t maxpKaons=0.4;
  Double_t maxbgKaons=0.7;
  Double_t minpProtons=0.35;
  Double_t maxpProtons=0.57;
  Double_t maxbgProtons=0.7;

  Double_t massPion=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t massKaon=TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t massProton=TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t massDeuteron=1.8756;

  Int_t iPad=1;
  Int_t iCanvas=0;
  Int_t nbinsX=hsigvsp->GetNbinsX();
  printf("Total Number of bins =%d\n",nbinsX);
  Int_t group=2;
 
  Int_t totPads=(Int_t)((Float_t)nbinsX/(Float_t)group+0.5);
  Int_t totCanvas=totPads/16;
  TCanvas** c0=new TCanvas*[totCanvas];
  c0[iCanvas]=new TCanvas(Form("c0_%d",iCanvas),Form("c0_%d",iCanvas));
  c0[iCanvas]->Divide(4,4);
  TCanvas** c0k=new TCanvas*[totCanvas];
  c0k[iCanvas]=new TCanvas(Form("c0k_%d",iCanvas),Form("c0k_%d",iCanvas));
  c0k[iCanvas]->Divide(4,4);
  TCanvas** c0p=new TCanvas*[totCanvas];
  c0p[iCanvas]=new TCanvas(Form("c0p_%d",iCanvas),Form("c0p_%d",iCanvas));
  c0p[iCanvas]->Divide(4,4);

  TGraphErrors* gdedxpi=new TGraphErrors(0);
  TGraphErrors* gdedxk=new TGraphErrors(0);
  TGraphErrors* gdedxp=new TGraphErrors(0);
  TGraphErrors* gdedxbg=new TGraphErrors(0);
  gdedxbg->SetTitle("");

  Int_t nPpi=0;
  Int_t nPk=0;
  Int_t nPp=0;
  Int_t nPbg=0;
  for(Int_t iBin=1; iBin<=nbinsX; iBin+=group){
    if(iPad>16){
      iCanvas++;
      c0[iCanvas]=new TCanvas(Form("c0_%d",iCanvas),Form("c0_%d",iCanvas));
      c0[iCanvas]->Divide(4,4);
      c0k[iCanvas]=new TCanvas(Form("c0k_%d",iCanvas),Form("c0k_%d",iCanvas));
      c0k[iCanvas]->Divide(4,4);
      c0p[iCanvas]=new TCanvas(Form("c0p_%d",iCanvas),Form("c0p_%d",iCanvas));
      c0p[iCanvas]->Divide(4,4);
      iPad=1;
    }

    Float_t sum=0;
    Float_t sumwei=0;
    Double_t pmed=0.;
    for(Int_t ibx=iBin; ibx<iBin+group;ibx++){
      for(Int_t iby=1; iby<=hsigvsp->GetNbinsY();iby++){
	sum+=hsigvsp->GetXaxis()->GetBinCenter(ibx)*hsigvsp->GetBinContent(ibx,iby);
	sumwei+=hsigvsp->GetBinContent(ibx,iby);
      }
    }
    if(sumwei>0.) pmed=sum/sumwei;
    if(pmed<0.1) continue;
    Int_t iPmed=(Int_t)(pmed*1000);
    TH1D* hProj=hsigvsp->TH2F::ProjectionY(Form("proj%d",iPmed),iBin,iBin+group-1);
    TH1D* hKaon=(TH1D*)hProj->Clone(Form("projk%d",iPmed));
    c0[iCanvas]->cd(iPad);
    hProj->SetLineColor(2);
    hProj->Draw();    

    Float_t peakpos=hProj->GetBinCenter(hProj->GetMaximumBin());
    Float_t peaksig=0.1*peakpos;
    hProj->Fit("gaus","LL","",peakpos-3*peaksig,peakpos+3.*peaksig);
    TF1* fgaus=(TF1*)hProj->GetListOfFunctions()->FindObject("gaus");
    gdedxpi->SetPoint(nPpi,sum/sumwei,fgaus->GetParameter(1));
    gdedxpi->SetPointError(nPpi,0.,fgaus->GetParError(1));
    nPpi++;
    if(pmed/massPion>minbgPions && pmed/massPion<maxbgPions){
      gdedxbg->SetPoint(nPbg,pmed/massPion,fgaus->GetParameter(1));
      gdedxbg->SetPointError(nPbg,0.,fgaus->GetParError(1));
      nPbg++;
    }
    for(Int_t iBinK=1; iBinK<=hKaon->GetNbinsX(); iBinK++){
      Double_t cont=hKaon->GetBinContent(iBinK)-fgaus->Eval(hKaon->GetBinCenter(iBinK));
      if(cont<0) cont=0.;
      if(hKaon->GetBinCenter(iBinK)<fgaus->GetParameter(1)) cont=0;
      hKaon->SetBinContent(iBinK,cont);
    }
    c0k[iCanvas]->cd(iPad);
    hKaon->SetLineColor(4);
    hKaon->Draw();
    TH1D* hProton=(TH1D*)hKaon->Clone(Form("projp%d",iPmed));
    TF1* fgausk=0x0;
    Double_t mink=250-600*(pmed-0.2);
    Double_t maxk=mink+mink;
    hKaon->Fit("gaus","LL","",mink,maxk);
    fgausk=(TF1*)hKaon->GetListOfFunctions()->FindObject("gaus");
    if(pmed>minpKaons && pmed<maxpKaons){
      gdedxk->SetPoint(nPk,pmed,fgausk->GetParameter(1));
      gdedxk->SetPointError(nPk,0.,fgausk->GetParError(1));
      nPk++;
      if(pmed/massKaon<maxbgKaons){
	gdedxbg->SetPoint(nPbg,pmed/massKaon,fgausk->GetParameter(1));
	gdedxbg->SetPointError(nPbg,0.,fgausk->GetParError(1));
	nPbg++;
      }    
    }
    if(fgausk){
      for(Int_t iBinP=1; iBinP<=hProton->GetNbinsX(); iBinP++){
	Double_t cont=hKaon->GetBinContent(iBinP)-fgausk->Eval(hKaon->GetBinCenter(iBinP));
	if(cont<0) cont=0.;
	if(hProton->GetBinCenter(iBinP)<fgausk->GetParameter(1)) cont=0;
	hProton->SetBinContent(iBinP,cont);
      }
    }
    c0p[iCanvas]->cd(iPad);
    hProton->SetLineColor(kGreen+1);
    hProton->Draw();
    if(pmed>minpProtons && pmed<maxpProtons){
      Double_t minP=300-600*(pmed-0.35);
      Double_t maxP=minP+minP;
      hProton->Fit("gaus","LL","",minP,maxP);
      TF1* fgausp=(TF1*)hProton->GetListOfFunctions()->FindObject("gaus");
      gdedxp->SetPoint(nPp,pmed,fgausp->GetParameter(1));
      gdedxp->SetPointError(nPp,0.,fgausp->GetParError(1));
      nPp++;
      if(pmed/massProton<maxbgProtons){
	gdedxbg->SetPoint(nPbg,pmed/massProton,fgausp->GetParameter(1));
	gdedxbg->SetPointError(nPbg,0.,fgausp->GetParError(1));
	nPbg++;
      }    
    }
    ++iPad;
  }

  gStyle->SetPalette(1);
  hsigvsp->GetXaxis()->SetRangeUser(0.07,1.8);
  hsigvsp->SetStats(0);
  TCanvas* c1=new TCanvas("c1","Pion Fit");
  c1->SetLogz();
  hsigvsp->Draw("col");
  gdedxpi->Draw("Psame");
  gdedxk->Draw("Psame");
  gdedxp->Draw("Psame");

  TF1 *fPhobos=new TF1("fPhobos",BetheBlochPhobos,minpPionFit,maxpPionFit,5);
  fPhobos->SetParLimits(0,0.,1.E9);
  fPhobos->SetParLimits(1,0.,100.);
  fPhobos->SetParLimits(2,0.,1.);
  fPhobos->SetParLimits(3,0.,1.);
  fPhobos->SetParLimits(4,0.,1.);
  gdedxpi->Fit("fPhobos","R");

  TCanvas* c1g=new TCanvas("c1g","Global Fit");
  TF1 *fPhobosbg=new TF1("fPhobosbg",BetheBlochPhobosBG,minbgGlobFit,maxbgGlobFit,5);
  fPhobosbg->SetParLimits(0,0.,1.E9);
  fPhobosbg->SetParLimits(1,0.,100.);
  fPhobosbg->SetParLimits(2,0.,1.);
  fPhobosbg->SetParLimits(3,0.,1.);
  fPhobosbg->SetParLimits(4,0.,1.);
  gdedxbg->Draw("AP");
  gdedxbg->GetXaxis()->SetTitle("#beta#gamma");
  gdedxbg->GetYaxis()->SetTitle("dE/dx (keV/300 #mum)");
  gdedxbg->Fit("fPhobosbg","R");
  c1g->Update();
  

  AliITSPIDResponse* itsresGlobFit=new AliITSPIDResponse(kFALSE);
  AliITSPIDResponse* itsres=new AliITSPIDResponse(kFALSE);
  Double_t parameters[5],parametersGlob[5];

  printf("-------- Pion Fit Parameters -----------\n");
  for(Int_t iPar=0; iPar<5; iPar++){
    parameters[iPar]=fPhobos->GetParameter(iPar);
    printf("parameters[%d]=%g;\n",iPar,parameters[iPar]);
  }
  printf("-------- Global Fit Parameters -----------\n");
  for(Int_t iPar=0; iPar<5; iPar++){
    parametersGlob[iPar]=fPhobosbg->GetParameter(iPar);
    printf("parametersGlob[%d]=%g;\n",iPar,parametersGlob[iPar]);
  }

  itsresGlobFit->SetBetheBlochParamsITSsa(parametersGlob);
  itsres->SetBetheBlochParamsITSsa(parameters);
  
  TGraph* gPion=new TGraph(0);
  TGraph* gKaon=new TGraph(0);
  TGraph* gProton=new TGraph(0);
  TGraph* gPionGlobFit=new TGraph(0);
  TGraph* gKaonGlobFit=new TGraph(0);
  TGraph* gProtonGlobFit=new TGraph(0);

  TGraphAsymmErrors* gBBpions=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gBBkaons=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gBBprotons=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gBBpionsGlobFit=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gBBkaonsGlobFit=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gBBprotonsGlobFit=new TGraphAsymmErrors(0);

  Double_t numberofsigmapions=2.;
  Float_t resodedx[4];
  resodedx[0]=0.0001; //default value
  resodedx[1]=0.0001; //default value
  resodedx[2]=0.116;
  resodedx[3]=0.103;
  

  for(Int_t ip=0; ip<400; ip++){
    Double_t mom=0.1+0.02*ip;
    Double_t dedxpions=itsres->Bethe(mom,massPion,kTRUE);
    Double_t dedxkaons=itsres->Bethe(mom,massKaon,kTRUE);
    Double_t dedxprotons=itsres->Bethe(mom,massProton,kTRUE);
    Double_t dedxdeutons=itsres->Bethe(mom,massDeuteron,kTRUE);
    Double_t dedxpionsGlobFit=itsresGlobFit->Bethe(mom,massPion,kTRUE);
    Double_t dedxkaonsGlobFit=itsresGlobFit->Bethe(mom,massKaon,kTRUE);
    Double_t dedxprotonsGlobFit=itsresGlobFit->Bethe(mom,massProton,kTRUE);
    Double_t dedxdeutonsGlobFit=itsresGlobFit->Bethe(mom,massDeuteron,kTRUE);

    gPion->SetPoint(ip,mom,dedxpions);
    gKaon->SetPoint(ip,mom,dedxkaons);
    gProton->SetPoint(ip,mom,dedxprotons);
    gPionGlobFit->SetPoint(ip,mom,dedxpionsGlobFit);
    gKaonGlobFit->SetPoint(ip,mom,dedxkaonsGlobFit);
    gProtonGlobFit->SetPoint(ip,mom,dedxprotonsGlobFit);

    gBBpions->SetPoint(ip,mom,dedxpions);
    gBBpions->SetPointError(ip,0,0,(dedxpions*resodedx[3]*numberofsigmapions),TMath::Abs((dedxpions-dedxkaons)/2));
    gBBkaons->SetPoint(ip,mom,dedxkaons);
    gBBkaons->SetPointError(ip,0,0,TMath::Abs((dedxpions-dedxkaons)/2),TMath::Abs((dedxprotons-dedxkaons)/2));
    gBBprotons->SetPoint(ip,mom,dedxprotons);
    gBBprotons->SetPointError(ip,0,0,TMath::Abs((dedxprotons-dedxkaons)/2),TMath::Abs((dedxprotons-dedxdeutons)/2));

    gBBpionsGlobFit->SetPoint(ip,mom,dedxpionsGlobFit);
    gBBpionsGlobFit->SetPointError(ip,0,0,(dedxpionsGlobFit*resodedx[3]*numberofsigmapions),TMath::Abs((dedxpions-dedxkaonsGlobFit)/2));
    gBBkaonsGlobFit->SetPoint(ip,mom,dedxkaonsGlobFit);
    gBBkaonsGlobFit->SetPointError(ip,0,0,TMath::Abs((dedxpionsGlobFit-dedxkaonsGlobFit)/2),TMath::Abs((dedxprotonsGlobFit-dedxkaonsGlobFit)/2));
    gBBprotonsGlobFit->SetPoint(ip,mom,dedxprotonsGlobFit);
    gBBprotonsGlobFit->SetPointError(ip,0,0,TMath::Abs((dedxprotonsGlobFit-dedxkaonsGlobFit)/2),TMath::Abs((dedxprotonsGlobFit-dedxdeutonsGlobFit)/2));
  }


  TCanvas* c2=new TCanvas("c2","Parameterizations");
  c2->SetLogz();
  hsigvsp->SetMinimum(20);
  hsigvsp->Draw("col");
  gPion->Draw("LSAME");
  gKaon->Draw("LSAME");
  gProton->Draw("LSAME");
  gPionGlobFit->SetLineColor(kRed+1);
  gKaonGlobFit->SetLineColor(kRed+1);
  gProtonGlobFit->SetLineColor(kRed+1);
  gPionGlobFit->Draw("LSAME");
  gKaonGlobFit->Draw("LSAME");
  gProtonGlobFit->Draw("LSAME");
  TLatex* t1=new TLatex(0.7,0.85,"Pion Fit");
  t1->SetNDC();
  t1->Draw();
  TLatex* tG=new TLatex(0.7,0.75,"Global Fit");
  tG->SetNDC();
  tG->SetTextColor(kRed+1);
  tG->Draw();


  gBBpions->SetLineColor(2);
  gBBpions->SetLineWidth(2);
  gBBpions->SetFillColor(2);
  gBBpions->SetFillStyle(3001);
  gBBkaons->SetLineColor(3);
  gBBkaons->SetLineWidth(2);
  gBBkaons->SetFillColor(3);
  gBBkaons->SetFillStyle(3001);
  gBBprotons->SetLineColor(4);
  gBBprotons->SetLineWidth(2);
  gBBprotons->SetFillColor(4);
  gBBprotons->SetFillStyle(3001);

  TCanvas* c3=new TCanvas("c3","Asymm Bands, PionFit");
  c3->SetLogz();
  c3->SetLogx();
  hsigvsp->SetMinimum(20);
  hsigvsp->Draw("col");
  gBBpions->Draw("L4same");
  gBBkaons->Draw("L4same");
  gBBprotons->Draw("L4same");

  gBBpionsGlobFit->SetLineColor(2);
  gBBpionsGlobFit->SetLineWidth(2);
  gBBpionsGlobFit->SetFillColor(2);
  gBBpionsGlobFit->SetFillStyle(3001);
  gBBkaonsGlobFit->SetLineColor(3);
  gBBkaonsGlobFit->SetLineWidth(2);
  gBBkaonsGlobFit->SetFillColor(3);
  gBBkaonsGlobFit->SetFillStyle(3001);
  gBBprotonsGlobFit->SetLineColor(4);
  gBBprotonsGlobFit->SetLineWidth(2);
  gBBprotonsGlobFit->SetFillColor(4);
  gBBprotonsGlobFit->SetFillStyle(3001);

  TCanvas* c3g=new TCanvas("c3g","Asymm Bands, GlobFit");
  c3g->SetLogz();
  c3g->SetLogx();
  hsigvsp->SetMinimum(20);
  hsigvsp->Draw("col");
  gBBpionsGlobFit->Draw("L4same");
  gBBkaonsGlobFit->Draw("L4same");
  gBBprotonsGlobFit->Draw("L4same");

}






Double_t BetheBlochPhobos(const Double_t *x, const Double_t *par){

  //
  // returns AliExternalTrackParam::BetheBloch normalized to 
  // fgMIP at the minimum
  //

  Double_t massPion=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t bg=x[0]/massPion;
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t gamma=bg/beta;
  Double_t eff=1.0;
  if(bg<par[2])
    eff=(bg-par[3])*(bg-par[3])+par[4];
  else
    eff=(par[2]-par[3])*(par[2]-par[3])+par[4];

  Double_t bb=0.;
  if(gamma>=0. && beta>0.){
    bb=(par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(beta*beta))*eff;
  }
  return bb;
}

Double_t BetheBlochPhobosBG(const Double_t *x, const Double_t *par){

  //
  // returns AliExternalTrackParam::BetheBloch normalized to 
  // fgMIP at the minimum
  //

  Double_t bg=x[0];
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t gamma=bg/beta;
  Double_t eff=1.0;
  if(bg<par[2])
    eff=(bg-par[3])*(bg-par[3])+par[4];
  else
    eff=(par[2]-par[3])*(par[2]-par[3])+par[4];

  Double_t bb=0.;
  if(gamma>=0. && beta>0.){
    bb=(par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(beta*beta))*eff;
  }
  return bb;
}

Double_t BetheBlochAleph(const Double_t *x, const Double_t *par){

  //
  // This is the empirical ALEPH parameterization of the Bethe-Bloch formula.
  // It is normalized to 1 at the minimum.
  //
  // bg - beta*gamma
  // 
  // The default values for the kp* parameters are for ALICE TPC.
  // The returned value is in MIP units
  //

  Double_t massPion=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t bg=x[0]/massPion;
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);

  Double_t aa = TMath::Power(beta,par[3]);
  Double_t bb = TMath::Power(1./bg,par[4]);

  bb=TMath::Log(par[2]+bb);
  
  return (par[1]-aa-bb)*par[0]/aa;
}
