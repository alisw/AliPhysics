#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TGrid.h>
#include <TFile.h>
#include <TList.h>
#include <TPaveStats.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#endif

/*  $Id$    */

//-------------------------------------------------------
//
// Macro do plot the output histograms of the QA task for ITS standalone tracks
// General Plots: ratios between ITSsa, ITSpureSA and TPC+ITS tracks 
//                eta phi distributions of tracks
//                number of clusters per track
// Pt resolution (matching ITSèureSA with TPC+ITS tracks)
// d0 resolution and bias 
//
// Authors: Leonardo Milano, Francesco Prino
//
//-------------------------------------------------------

enum{kDoGeneral,kDoPt,kDoImpPar,kDoAll};

void PlotITSsa(TList* l);
void PlotPtResol(TList* l, Bool_t optFromMC);
void PlotImpPar(TList* l);
void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1);

void PlotOutputQATaskITSsa(TString filename, Int_t analysisType=kDoAll, Bool_t isMC=kFALSE){

  gROOT->SetStyle("Plain");

  if(filename.Contains("alien")) TGrid::Connect("alien:");
  TFile* fil=TFile::Open(filename.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("TracksITSsa");
  if(!df) df=(TDirectoryFile*)fil->Get("ITSsaTracks");
  TList* l=(TList*)df->Get("clistITSsaTracks");
  if(analysisType==kDoAll || analysisType==kDoGeneral) PlotITSsa(l);
  if(analysisType==kDoAll || analysisType==kDoPt) PlotPtResol(l,isMC);
  if(analysisType==kDoAll || analysisType==kDoImpPar) PlotImpPar(l);

}

//-----------------------------------------------------

void PlotITSsa(TList* l){

  TH1F* hPtTPCITS=(TH1F*)l->FindObject("hPtTPCITS");
  TH1F* hPtITSsa=(TH1F*)l->FindObject("hPtITSsa");
  TH1F* hPtITSpureSA=(TH1F*)l->FindObject("hPtITSpureSA");

  TH2F* hEtaPhiTPCITS=(TH2F*)l->FindObject("hEtaPhiTPCITS");
  TH2F* hEtaPhiITSsa=(TH2F*)l->FindObject("hEtaPhiITSsa");
  TH2F* hEtaPhiITSpureSA=(TH2F*)l->FindObject("hEtaPhiITSpureSA");

  TH1F* hChi2TPCITS=(TH1F*)l->FindObject("hChi2TPCITS");
  TH1F* hChi2ITSsa=(TH1F*)l->FindObject("hChi2ITSsa");
  TH1F* hChi2ITSpureSA=(TH1F*)l->FindObject("hChi2ITSpureSA");

  TH1F* hNcluTPCITS=(TH1F*)l->FindObject("hNcluTPCITS");
  TH1F* hNcluITSsa=(TH1F*)l->FindObject("hNcluITSsa");
  TH1F* hNcluITSpureSA=(TH1F*)l->FindObject("hNcluITSpureSA");


  TH1F* hRatio=(TH1F*)hPtTPCITS->Clone("hRatio");
  hRatio->Add(hPtITSsa);
  hRatio->Divide(hPtITSpureSA);
  hRatio->SetStats(0);

  TCanvas* c1=new TCanvas("c1","Pt",800,1000);
  c1->Divide(1,2);
  c1->cd(1);
  hPtITSpureSA->Draw();
  hPtITSpureSA->GetXaxis()->SetTitle("Pt (GeV/c)");
  gPad->Update();
  TPaveStats *st1=(TPaveStats*)hPtITSpureSA->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.71);
  st1->SetY2NDC(0.9);
  hPtTPCITS->SetLineColor(2);
  hPtTPCITS->Draw("sames");
  gPad->Update();
  TPaveStats *st2=(TPaveStats*)hPtTPCITS->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.51);
  st2->SetY2NDC(0.7);
  st2->SetTextColor(2);

  hPtITSsa->SetLineColor(4);
  hPtITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *st3=(TPaveStats*)hPtITSsa->GetListOfFunctions()->FindObject("stats");
  st3->SetY1NDC(0.31);
  st3->SetY2NDC(0.5);
  st3->SetTextColor(4);
  TLegend* leg=new TLegend(0.5,0.5,0.69,0.69);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hPtITSpureSA,"ITS pureSA","L");
  ent->SetTextColor(hPtITSpureSA->GetLineColor());
  ent=leg->AddEntry(hPtTPCITS,"TPC+ITS","L");
  ent->SetTextColor(hPtTPCITS->GetLineColor());
  ent=leg->AddEntry(hPtITSsa,"ITSsa","L");
  ent->SetTextColor(hPtITSsa->GetLineColor());
  leg->Draw();
  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hRatio->Draw();
  hRatio->GetXaxis()->SetTitle("Pt (GeV/c)");
  hRatio->GetYaxis()->SetTitle("(TPCITS+ITSsa)/ITSpureSA");

  hChi2ITSpureSA->Scale(1./hChi2ITSpureSA->GetEntries());
  hChi2ITSsa->Scale(1./hChi2ITSsa->GetEntries());
  hChi2TPCITS->Scale(1./hChi2TPCITS->GetEntries());

  TCanvas* c2=new TCanvas("c2","Chi2");
  hChi2ITSpureSA->Draw();
  hChi2ITSpureSA->GetXaxis()->SetTitle("Chi2");
  gPad->Update();
  TPaveStats *stc1=(TPaveStats*)hChi2ITSpureSA->GetListOfFunctions()->FindObject("stats");
  stc1->SetY1NDC(0.71);
  stc1->SetY2NDC(0.9);
  hChi2TPCITS->SetLineColor(2);
  hChi2TPCITS->Draw("sames");
  gPad->Update();
  TPaveStats *stc2=(TPaveStats*)hChi2TPCITS->GetListOfFunctions()->FindObject("stats");
  stc2->SetY1NDC(0.51);
  stc2->SetY2NDC(0.7);
  stc2->SetTextColor(2);
  c2->Update();
  hChi2ITSsa->SetLineColor(4);
  hChi2ITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *stc3=(TPaveStats*)hChi2ITSsa->GetListOfFunctions()->FindObject("stats");
  stc3->SetY1NDC(0.31);
  stc3->SetY2NDC(0.5);
  stc3->SetTextColor(4);
  leg->Draw();

  hNcluITSpureSA->Scale(1./hNcluITSpureSA->GetEntries());
  hNcluITSsa->Scale(1./hNcluITSsa->GetEntries());
  hNcluTPCITS->Scale(1./hNcluTPCITS->GetEntries());

  TCanvas* c3=new TCanvas("c3","Nclu");
  c3->SetRightMargin(0.22);
  hNcluITSpureSA->Draw();
  hNcluITSpureSA->GetXaxis()->SetTitle("n. ITS clusters");
  gPad->Update();
  TPaveStats *stn1=(TPaveStats*)hNcluITSpureSA->GetListOfFunctions()->FindObject("stats");
  stn1->SetY1NDC(0.71);
  stn1->SetY2NDC(0.9);
  hNcluTPCITS->SetLineColor(2);
  hNcluTPCITS->Draw("sames");
  gPad->Update();
  TPaveStats *stn2=(TPaveStats*)hNcluTPCITS->GetListOfFunctions()->FindObject("stats");
  stn2->SetY1NDC(0.51);
  stn2->SetY2NDC(0.7);
  stn2->SetTextColor(2);

  hNcluITSsa->SetLineColor(4);
  hNcluITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *stn3=(TPaveStats*)hNcluITSsa->GetListOfFunctions()->FindObject("stats");
  stn3->SetY1NDC(0.31);
  stn3->SetY2NDC(0.5);
  stn3->SetTextColor(4);
  leg->Draw();

  gStyle->SetPalette(1);
  hEtaPhiITSpureSA->SetStats(0);
  hEtaPhiITSpureSA->SetTitle("ITS pureSA");
  hEtaPhiITSsa->SetStats(0);
  hEtaPhiITSsa->SetTitle("ITSsa");
  hEtaPhiTPCITS->SetStats(0);
  hEtaPhiTPCITS->SetTitle("TPC+ITS");
  TCanvas* c4=new TCanvas("c4","EtaPhi",1000,700);
  c4->Divide(3,1);
  c4->cd(1);
  hEtaPhiITSpureSA->Draw("colz");
  hEtaPhiITSpureSA->GetXaxis()->SetTitle("Eta");
  hEtaPhiITSpureSA->GetYaxis()->SetTitle("Phi");
  c4->cd(2);
  hEtaPhiITSsa->Draw("colz");
  hEtaPhiITSsa->GetXaxis()->SetTitle("Eta");
  hEtaPhiITSsa->GetYaxis()->SetTitle("Phi");
  c4->cd(3);
  hEtaPhiTPCITS->Draw("colz");
  hEtaPhiTPCITS->GetXaxis()->SetTitle("Eta");
  hEtaPhiTPCITS->GetYaxis()->SetTitle("Phi");
}

//-----------------------------------------------------

void PlotPtResol(TList* l, Bool_t optFromMC){
  TString hNameR,hNameA;
  TString partName[3]={"Pion","Kaon","Proton"};
  TString prefix;
  if(optFromMC) prefix="hMC";
  else prefix="h";
  
  TCanvas* c2d[3];
  TCanvas* c1dA[3];
  TCanvas* c1dR[3];

  TH2F* h2DA[3];
  TH2F* h2DR[3];
  TH1F* hptres[3][40];
  TH1F* h1ptrelres[3][40];
  TH1F* hptreco[3][40];

  TGraphErrors* gbias[3];
  TGraphErrors* grelresol[3];

  gStyle->SetPalette(1);

  for(Int_t iSpec=0; iSpec<3; iSpec++){
    hNameA=Form("%sPtResid%s",prefix.Data(),partName[iSpec].Data());
    hNameR=Form("%sInvPtRelResid%s",prefix.Data(),partName[iSpec].Data());
    printf("%s %s\n",hNameA.Data(),hNameR.Data());
    h2DA[iSpec]=(TH2F*)l->FindObject(hNameA.Data());
    h2DR[iSpec]=(TH2F*)l->FindObject(hNameR.Data());
    c2d[iSpec]=new TCanvas(Form("c2d%s",partName[iSpec].Data()),Form("c2d%s",partName[iSpec].Data()));
    c2d[iSpec]->Divide(2,1);
    c2d[iSpec]->cd(1);
    h2DA[iSpec]->Draw("colz"); 
    c2d[iSpec]->cd(2);
    h2DR[iSpec]->Draw("colz");

    Int_t nptbins=h2DR[iSpec]->GetNbinsX();

    Int_t nybinsA=h2DA[iSpec]->GetNbinsY();
    Float_t minyA=h2DA[iSpec]->GetYaxis()->GetBinLowEdge(1);
    Float_t maxyA=h2DA[iSpec]->GetYaxis()->GetBinUpEdge(nybinsA);

    Int_t nybinsR=h2DR[iSpec]->GetNbinsY();
    Float_t minyR=h2DR[iSpec]->GetYaxis()->GetBinLowEdge(1);
    Float_t maxyR=h2DR[iSpec]->GetYaxis()->GetBinUpEdge(nybinsR);
    printf("%d   %d %f %f   %d %f %f\n",nptbins,nybinsA,minyA,maxyA,nybinsR,minyR,maxyR);

    c1dA[iSpec]=new TCanvas(Form("c1dA%s",partName[iSpec].Data()),Form("c1dA%s",partName[iSpec].Data()));
    c1dA[iSpec]->Divide(6,5);
    c1dR[iSpec]=new TCanvas(Form("c1dR%s",partName[iSpec].Data()),Form("c1dR%s",partName[iSpec].Data()));
    c1dR[iSpec]->Divide(6,5);


    gbias[iSpec]=new TGraphErrors(0);
    grelresol[iSpec]=new TGraphErrors(0);
    gbias[iSpec]->SetTitle("");
    grelresol[iSpec]->SetTitle("");

    for(Int_t iptbin=0; iptbin<nptbins;iptbin++){
      Float_t avept=h2DA[iSpec]->GetXaxis()->GetBinCenter(iptbin+1);
      Float_t widpt=0.5*h2DA[iSpec]->GetXaxis()->GetBinWidth(iptbin+1);
      Int_t minptbinmev=(Int_t)(h2DA[iSpec]->GetXaxis()->GetBinLowEdge(iptbin+1)*1000.+0.5);

      hptres[iSpec][iptbin]=new TH1F(Form("hptres%s_%d",partName[iSpec].Data(),minptbinmev),
				     Form("hptres%s_%d",partName[iSpec].Data(),minptbinmev),
				     nybinsA,minyA,maxyA);
      h1ptrelres[iSpec][iptbin]=new TH1F(Form("h1ptrelres%s_%d",partName[iSpec].Data(),minptbinmev),
					 Form("h1ptrelres%s_%d",partName[iSpec].Data(),minptbinmev),
					 nybinsR,minyR,maxyR);
      hptreco[iSpec][iptbin]=new TH1F(Form("hptreco%s_%d",partName[iSpec].Data(),minptbinmev),
				      Form("hptreco%s_%d",partName[iSpec].Data(),minptbinmev),
				      400,0.,2.);
      for(Int_t iBin=1; iBin<=nybinsA; iBin++){
	hptres[iSpec][iptbin]->SetBinContent(iBin,h2DA[iSpec]->GetBinContent(iptbin+1,iBin));
	hptres[iSpec][iptbin]->SetBinError(iBin,h2DA[iSpec]->GetBinError(iptbin+1,iBin));
      }
      for(Int_t iBin=1; iBin<=nybinsR; iBin++){
	h1ptrelres[iSpec][iptbin]->SetBinContent(iBin,h2DR[iSpec]->GetBinContent(iptbin+1,iBin));
	h1ptrelres[iSpec][iptbin]->SetBinError(iBin,h2DR[iSpec]->GetBinError(iptbin+1,iBin));
      }

      c1dA[iSpec]->cd(iptbin+1);
      hptres[iSpec][iptbin]->Draw();
      if(hptres[iSpec][iptbin]->Integral()>50){
	hptres[iSpec][iptbin]->Fit("gaus");
	hptres[iSpec][iptbin]->GetXaxis()->SetTitle("Pt residuals (GeV/c)");
	hptres[iSpec][iptbin]->GetXaxis()->CenterTitle();
	TF1* fgaus= (TF1*)hptres[iSpec][iptbin]->GetListOfFunctions()->FindObject("gaus");
	Int_t nPoint=gbias[iSpec]->GetN();
	gbias[iSpec]->SetPoint(nPoint, avept, fgaus->GetParameter(1));//hptres[iSpec][iptbin]->GetMean());
	gbias[iSpec]->SetPointError(nPoint, widpt, fgaus->GetParError(1)); //hptres[iSpec][iptbin]->GetMeanError());
      }
      c1dR[iSpec]->cd(iptbin+1);
      h1ptrelres[iSpec][iptbin]->Draw();
      if(h1ptrelres[iSpec][iptbin]->Integral()>50){
	h1ptrelres[iSpec][iptbin]->Fit("gaus");//,"L");
	h1ptrelres[iSpec][iptbin]->GetXaxis()->SetTitle("1/Pt relative residuals");
	h1ptrelres[iSpec][iptbin]->GetXaxis()->CenterTitle();
	TF1* fgaus= (TF1*)h1ptrelres[iSpec][iptbin]->GetListOfFunctions()->FindObject("gaus");    
	Int_t nPoint=grelresol[iSpec]->GetN();
	grelresol[iSpec]->SetPoint(nPoint, avept, fgaus->GetParameter(2));
	grelresol[iSpec]->SetPointError(nPoint, widpt, fgaus->GetParError(2));    
      }

    }
  }
  

  TCanvas* cb=new TCanvas("cb","Bias");
  gbias[2]->SetMarkerStyle(22);
  gbias[2]->SetMarkerColor(4);
  gbias[2]->SetLineColor(4);
  gbias[2]->Draw("PA");
  gbias[0]->SetMarkerStyle(20);
  gbias[0]->SetMarkerColor(1);
  gbias[0]->SetLineColor(1);
  gbias[0]->Draw("PSAME");
  gbias[1]->SetMarkerStyle(25);
  gbias[1]->SetMarkerColor(2);
  gbias[1]->SetLineColor(2);
  gbias[1]->Draw("PSAME");
  gbias[2]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  if(optFromMC) gbias[2]->GetYaxis()->SetTitle("<p_{T}(ITSsa)-p_{T}(MC)> (GeV/c)");
  else gbias[2]->GetYaxis()->SetTitle("<p_{T}(ITSsa)-p_{T}(TPCITS)> (GeV/c)");
  gbias[2]->GetYaxis()->SetTitleOffset(1.2);
  cb->Update();

  TCanvas* cr=new TCanvas("cr","Resol");
  grelresol[2]->SetMinimum(0.);
  grelresol[2]->SetMaximum(0.2);
  grelresol[2]->SetMarkerStyle(22);
  grelresol[2]->SetMarkerColor(4);
  grelresol[2]->SetLineColor(4);
  grelresol[2]->Draw("PA");
  grelresol[0]->SetMarkerStyle(20);
  grelresol[0]->SetMarkerColor(1);
  grelresol[0]->SetLineColor(1);
  grelresol[0]->Draw("PSAME");
  grelresol[1]->SetMarkerStyle(25);
  grelresol[1]->SetMarkerColor(2);
  grelresol[1]->SetLineColor(2);
  grelresol[1]->Draw("PSAME");
  grelresol[2]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  grelresol[2]->GetYaxis()->SetTitle("1/Pt relative resolution (%)");
  grelresol[2]->GetYaxis()->SetTitleOffset(1.2);
  cr->Update();

}

//-----------------------------------------------------

void PlotImpPar(TList* l){

  TString particle[3]={"Pion","Kaon","Proton"};
  Int_t   colors[3]={1,2,4};
  TCanvas *cImpPar=new TCanvas("ImpParRes","ImpParRes",1000,800);
  cImpPar->SetGridx();
  cImpPar->SetLogx();
  cImpPar->SetLeftMargin(0.14);
  TCanvas *cImpParMean=new TCanvas("ImpParMean","ImpParMean",1000,800);
  cImpParMean->SetGridx();
  cImpParMean->SetLogx();
  cImpParMean->SetLeftMargin(0.14);

  TLegend* leg1=new TLegend(0.6,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegendEntry* ent;
  
  TH2F *hd0rphiITSpureSA[3];
  //binning
  const Int_t nbins = 29;
  Double_t xbins[nbins+1]={0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,
			   0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,
			   0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,1.90,2.00};
  
  TH1F *fHistDCA[nbins];  
  
  for(Int_t iparticle=0;iparticle<3;iparticle++){
    hd0rphiITSpureSA[iparticle]=(TH2F*)l->FindObject(Form("hd0rphiITSpureSA%s",particle[iparticle].Data()));
    
    
    TH1F *fHistImpParRes = new TH1F(Form("fHistImpParRes%s",particle[iparticle].Data()),"",nbins,xbins);
    fHistImpParRes->SetStats(0);
    fHistImpParRes->GetXaxis()->SetTitle("Pt [GeV/c]");
    fHistImpParRes->GetYaxis()->SetTitle("d0 r#phi resolution [#mum]");
    fHistImpParRes->GetYaxis()->SetTitleOffset(1.4);
    TH1F *fHistImpParMean = new TH1F(Form("fHistImpParMean%s",particle[iparticle].Data()),"",nbins,xbins);
    fHistImpParMean->SetStats(0);
    fHistImpParMean->GetXaxis()->SetTitle("Pt [GeV/c]");
    fHistImpParMean->GetYaxis()->SetTitle("d0 r#phi mean [#mum]");
    fHistImpParMean->GetYaxis()->SetTitleOffset(1.4);
    TF1 *fPar = new TF1("fPar","gaus",-1,1);
    for(Int_t m=0;m<nbins;m++){
      
      fHistDCA[m]= (TH1F*)hd0rphiITSpureSA[iparticle]->ProjectionY(Form("%s%i",particle[iparticle].Data(),m),hd0rphiITSpureSA[iparticle]->GetXaxis()->FindBin(xbins[m]+0.000001),hd0rphiITSpureSA[iparticle]->GetXaxis()->FindBin(xbins[m+1]-0.000001));
      fHistDCA[m]->Rebin();
    }
    
    TCanvas *cgaus=new TCanvas(Form("cgaus%s",particle[iparticle].Data()),Form("DCA dist %s ",particle[iparticle].Data()),1000,800);
    cgaus->Divide(8,4,0.001,0.001);
    for(Int_t i=0; i<nbins; i++){
      cgaus->cd(i+1);
      fHistDCA[i]->SetLineColor(colors[iparticle]);
      fHistDCA[i]->SetMarkerColor(colors[iparticle]);
      fPar->SetLineColor(colors[iparticle]);
      
      fHistDCA[i]->Draw();
      fHistDCA[i]->SetFillColor(16);
      printf("\n\n\n\n\n\n first fit step\n\n\n\n\n");
      fHistDCA[i]->Fit(fPar,"NM","",-1,1);
      printf("\n\n\n\n\n\n second fit step\n\n\n\n\n");
      Float_t nsigmas=1.;
      Float_t sigma=fPar->GetParameter(2);
      fHistDCA[i]->Fit(fPar,"NM","",fPar->GetParameter(1)-nsigmas*sigma,fPar->GetParameter(1)+nsigmas*sigma);
      fHistDCA[i]->GetXaxis()->SetRangeUser(fPar->GetParameter(1)-nsigmas*sigma,fPar->GetParameter(1)+nsigmas*sigma);
      fPar->DrawClone("same");
      fHistImpParRes->Fill((xbins[i]+xbins[i+1])/2,10000*fPar->GetParameter(2));
      fHistImpParRes->SetBinError(fHistImpParRes->FindBin((xbins[i]+xbins[i+1])/2),10000*fPar->GetParError(2));
      fHistImpParMean->Fill((xbins[i]+xbins[i+1])/2,10000*fPar->GetParameter(1));
      fHistImpParMean->SetBinError(fHistImpParMean->FindBin((xbins[i]+xbins[i+1])/2),10000*fPar->GetParError(1));
    }
    fHistImpParRes->SetMaximum(1000);
    fHistImpParRes->SetMinimum(0);
    fHistImpParMean->SetMaximum(80);
    fHistImpParMean->SetMinimum(-80);
    
    SetDrawAtt(iparticle+20,colors[iparticle],1,colors[iparticle],1,fHistImpParRes);
    SetDrawAtt(iparticle+20,colors[iparticle],1,colors[iparticle],1,fHistImpParMean);
        
    cImpPar->cd();
    if(iparticle==0)fHistImpParRes->DrawCopy("p");
    else fHistImpParRes->DrawCopy("psame");

    cImpParMean->cd();
    if(iparticle==0)fHistImpParMean->DrawCopy("p");
    else fHistImpParMean->DrawCopy("psame");
    
    ent=leg1->AddEntry(fHistImpParMean,Form("%s",particle[iparticle].Data()),"PL");
    ent->SetTextColor(fHistImpParMean->GetLineColor());
  }  
  cImpPar->cd();
  leg1->Draw();
  cImpParMean->cd();
  leg1->Draw();
  
} 


void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1){ 


  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
}


