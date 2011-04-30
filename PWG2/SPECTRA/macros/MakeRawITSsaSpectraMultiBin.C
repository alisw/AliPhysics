#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TLatex.h>
#include <TImage.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TList.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <Rtypes.h>
#include "AliITSsadEdxFitter.h"
#endif


Double_t BetheBloch(Double_t *mom, Double_t *mass);
void Logo();
//

//______________________________________________________________________
void MakeRawITSsaSpectraMultiBin(Bool_t optMC=kTRUE, Int_t multibin=0){


  AliITSsadEdxFitter *ITSsa=new AliITSsadEdxFitter(optMC);

  TString filename, dirName, ps0, ps1, ps2;
  Char_t listname[50];
  if(optMC) dirName=(Form("../gridmultiplicitybins/LHC10d1_1.5sigma_7DCA_negmag/Spectra_MC_negmag_MultBin%d",multibin));
  else dirName=(Form("../gridmultiplicitybins/data_1.5sigma_7DCA_negmag/Spectra_data_negmag_MultBin%d",multibin));
  switch(multibin){
  case 0:
    sprintf(listname,"clistITSsaMult0to9999");
    break;
  case 1:
    sprintf(listname,"clistITSsaMult0to5");
    break;
  case 2:
    sprintf(listname,"clistITSsaMult6to9");
    break;
  case 3:
    sprintf(listname,"clistITSsaMult10to14");
    break;
  case 4:
    sprintf(listname,"clistITSsaMult15to22");
    break;
  case 5:
    sprintf(listname,"clistITSsaMult23to9999");
    break;
  }
  if(optMC){
    filename=Form("%s/AnalysisResults.root",dirName.Data());
    ps0 = Form("%s/outSIM.ps[",dirName.Data());
    ps1 = Form("%s/outSIM.ps",dirName.Data());
    ps2 = Form("%s/outSIM.ps]",dirName.Data());
  }
  else{
    filename=Form("%s/AnalysisResults.root",dirName.Data());
    ps0 = Form("%s/outDATA.ps[",dirName.Data());
    ps1 = Form("%s/outDATA.ps",dirName.Data());
    ps2 = Form("%s/outDATA.ps]",dirName.Data());
  }
  TString openfilename="./";
  openfilename+=filename;
  cout<<openfilename<<endl;
  TFile *fi=new TFile(openfilename.Data());
  if(!fi){
    cout<<"TFile loading failed"<<endl;
    return;
  }
  TDirectoryFile *di=(TDirectoryFile*) fi->Get("PWG2SpectraITSsa");
  if(!di){
    cout<<"TDirectory loading failed!"<<endl;
    return;
  }
  TList *li=(TList*)di->Get(listname);
  if(!li){
    cout<<"TList loading failed"<<endl;
    return;
  }
  cout<<"File loaded"<<endl;
  
  TCanvas *cdummy=new TCanvas("dummy","IntegralMethod",1000,800);
  cdummy->Print(ps0.Data());
  
  //binning
  const Int_t nbins = 22;
  Double_t xbins[nbins+1]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  
  //histograms
  TH1F *fHistMCPosPi[nbins]; 
  TH1F *fHistMCPosK[nbins]; 
  TH1F *fHistMCPosP[nbins]; 
  TH1F *fHistMCNegPi[nbins]; 
  TH1F *fHistMCNegK[nbins]; 
  TH1F *fHistMCNegP[nbins]; 
  TH1F *fHistPosPi[nbins]; 
  TH1F *fHistPosK[nbins]; 
  TH1F *fHistPosP[nbins]; 
  TH1F *fHistNegPi[nbins]; 
  TH1F *fHistNegK[nbins]; 
  TH1F *fHistNegP[nbins]; 
  for(Int_t m=0;m<nbins;m++){
    fHistMCNegPi[m] = (TH1F*)li->FindObject(Form("fHistMCNegPi%d",m));
    fHistMCNegK[m] = (TH1F*)li->FindObject(Form("fHistMCNegK%d",m));
    fHistMCNegP[m] = (TH1F*)li->FindObject(Form("fHistMCNegP%d",m));
    fHistMCPosPi[m] = (TH1F*)li->FindObject(Form("fHistMCPosPi%d",m));
    fHistMCPosK[m] = (TH1F*)li->FindObject(Form("fHistMCPosK%d",m));
    fHistMCPosP[m] = (TH1F*)li->FindObject(Form("fHistMCPosP%d",m));
    fHistPosPi[m] = (TH1F*)li->FindObject(Form("fHistPosPi%d",m));
    fHistPosK[m] = (TH1F*)li->FindObject(Form("fHistPosK%d",m));
    fHistPosP[m] = (TH1F*)li->FindObject(Form("fHistPosP%d",m));
    fHistNegPi[m] = (TH1F*)li->FindObject(Form("fHistNegPi%d",m));
    fHistNegK[m] = (TH1F*)li->FindObject(Form("fHistNegK%d",m));
    fHistNegP[m] = (TH1F*)li->FindObject(Form("fHistNegP%d",m));
  }
  
  TH1F *fHistNEvents = (TH1F*)li->FindObject("fHistNEvents");
  TH2F *fHistDEDX = (TH2F*)li->FindObject("fHistDEDX");
  TH2F *fHistDEDXdouble = (TH2F*)li->FindObject("fHistDEDXdouble");
  TH1F *fHistCharge[4];
  for(Int_t j=0;j<4;j++) fHistCharge[j] = (TH1F*)li->FindObject((Form("fHistChargeLay%d",j)));

  TH1F *hEffPos[3];
  TH1F *hEffNeg[3];
  TH1F *hCorrFacPos[3];
  TH1F *hCorrFacNeg[3];
  TH1F *hEffMCPIDPos[3];
  TH1F *hEffMCPIDNeg[3];
  TH1F *hCorrFacMCPIDNeg[3];
  TH1F *hCorrFacMCPIDPos[3];
  TH1F *hSpectraPrimPosMC[3];
  TH1F *hSpectraPrimNegMC[3];
  TH1F *hSpectraPrimPosMCBefEvSel[3];
  TH1F *hSpectraPrimNegMCBefEvSel[3];
  TH1F *hSpectraPos[3];
  TH1F *hSpectraNeg[3];
  TH1F *hSpectraMCPIDPos[3];
  TH1F *hSpectraMCPIDNeg[3];
  TH1F* hMeanPos[3];
  TH1F* hMeanNeg[3];
  TH1F* hSigmaPos[3];
  TH1F* hSigmaNeg[3];
  TGraph *gres[6][22];

  for(Int_t i=0; i<3; i++){
    hSpectraPrimPosMC[i]=(TH1F*)li->FindObject(Form("fHistPrimMCpos%d",i));
    hSpectraPrimNegMC[i]=(TH1F*)li->FindObject(Form("fHistPrimMCneg%d",i));
    hSpectraPrimPosMCBefEvSel[i]=(TH1F*)li->FindObject(Form("fHistPrimMCposBefEvSel%d",i));
    hSpectraPrimNegMCBefEvSel[i]=(TH1F*)li->FindObject(Form("fHistPrimMCnegBefEvSel%d",i));
    hSpectraMCPIDPos[i]=new TH1F(Form("hSpectraMCPIDPos%d",i),Form("hSpectraMCPIDPos%d",i),nbins,xbins);
    hSpectraMCPIDNeg[i]=new TH1F(Form("hSpectraMCPIDNeg%d",i),Form("hSpectraMCPIDNeg%d",i),nbins,xbins);
    hSpectraPos[i]=new TH1F(Form("hSpectraPos%d",i),Form("hSpectraPos%d",i),nbins,xbins);
    hSpectraNeg[i]=new TH1F(Form("hSpectraNeg%d",i),Form("hSpectraNeg%d",i),nbins,xbins);
    hMeanPos[i]=new TH1F(Form("hMeanPos%d",i),Form("hMeanPos%d",i),nbins,xbins);
    hMeanNeg[i]=new TH1F(Form("hMeanNeg%d",i),Form("hMeanNeg%d",i),nbins,xbins);
    hSigmaPos[i]=new TH1F(Form("hSigmaPos%d",i),Form("hSigmaPos%d",i),nbins,xbins);
    hSigmaNeg[i]=new TH1F(Form("hSigmaNeg%d",i),Form("hSigmaNeg%d",i),nbins,xbins);
  }

  //division for DeltaPt (for MC spectra)
  for(Int_t ipt=0;ipt<3;ipt++){
    for(Int_t bin=1; bin <= hSpectraPrimPosMC[ipt]->GetNbinsX(); bin++){ 
      Float_t binSize=hSpectraPrimPosMC[ipt]->GetBinLowEdge(bin+1) - hSpectraPrimPosMC[ipt]->GetBinLowEdge(bin);
      hSpectraPrimPosMC[ipt]->SetBinContent(bin, hSpectraPrimPosMC[ipt]->GetBinContent(bin) / binSize);
      hSpectraPrimPosMC[ipt]->SetBinError(bin, 0);
      hSpectraPrimNegMC[ipt]->SetBinContent(bin, hSpectraPrimNegMC[ipt]->GetBinContent(bin) / binSize);
      hSpectraPrimNegMC[ipt]->SetBinError(bin, 0);
      if(hSpectraPrimPosMCBefEvSel[ipt]){
	hSpectraPrimPosMCBefEvSel[ipt]->SetBinContent(bin, hSpectraPrimPosMCBefEvSel[ipt]->GetBinContent(bin) / binSize);
	hSpectraPrimPosMCBefEvSel[ipt]->SetBinError(bin, 0);
      }
      if(hSpectraPrimNegMCBefEvSel[ipt]){
	hSpectraPrimNegMCBefEvSel[ipt]->SetBinContent(bin, hSpectraPrimNegMCBefEvSel[ipt]->GetBinContent(bin) / binSize);
	hSpectraPrimNegMCBefEvSel[ipt]->SetBinError(bin, 0);
      }
    }
  }
  cout<<"All plots loaded"<<endl;

  //open output file to store the dedx distribution and the spectra
  TString savename=filename.Data();  
  savename.ReplaceAll("AnalysisResults","SpectraReco");
  TFile *fout=new TFile(savename.Data(),"recreate");
  fout->cd();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetPalette(1);

  //propaganda plot
  Float_t pdgmass[5]={0.13957,0.493677,0.938272,1.875612762,0.00596}; //mass for pi, K, P, d (Gev/c^2)
  TF1 *funPos[5];
  TF1 *funNeg[5];
  for(Int_t m=0;m<5;m++){
    funPos[m] = new TF1(Form("funPos%d",m),BetheBloch,0.02,5,1);
    funPos[m]->SetParameter(0,pdgmass[m]);
    funPos[m]->SetLineWidth(2);
    funNeg[m] = new TF1(Form("funNeg%d",m),BetheBloch,-5,0.02,1);
    funNeg[m]->SetParameter(0,-pdgmass[m]);
    funNeg[m]->SetLineWidth(2);
  }

  TCanvas *cEvents=new TCanvas("cEvents","cEvents",500,400);
  cEvents->cd();
  fHistNEvents->Draw("text");
  fHistNEvents->Write();

  TCanvas *cDEDX=new TCanvas("cDEDX","DEDX",1000,800);
  cDEDX->cd();
  gPad->SetLogx();
  gPad->SetLogz();
  fHistDEDX->GetXaxis()->SetRangeUser(0.08,5);
  fHistDEDX->GetYaxis()->SetRangeUser(0.,700);
  fHistDEDX->GetXaxis()->SetTitle("momentum [GeV/c]");
  fHistDEDX->GetYaxis()->SetTitle("dE [keV/300#mum]");
  fHistDEDX->Draw("colz");
  fHistDEDX->Write();
  for(Int_t m=0;m<3;m++) funPos[m]->Draw("same");
  Logo();

  TCanvas *cDEDXdouble=new TCanvas("cDEDXdouble","DEDXdouble",1000,800);
  cDEDXdouble->cd();
  gPad->SetLogz();
  fHistDEDXdouble->GetXaxis()->SetRangeUser(-3,3);
  fHistDEDXdouble->GetXaxis()->SetTitle("momentum * sign [GeV/c]");
  fHistDEDXdouble->GetYaxis()->SetTitle("dE [keV/300#mum]");
  fHistDEDXdouble->Draw("colz");
  fHistDEDXdouble->Write();
  for(Int_t m=0;m<3;m++) {
    funPos[m]->Draw("same");
    funNeg[m]->Draw("same");
  }	
  Logo();

  //calibration check histo
  TCanvas *cs=new TCanvas("cs","cs",1000,800);
  cs->Divide(2,2);
  for(Int_t j=0;j<4;j++){
    cs->cd(j+1);
    fHistCharge[j]->GetXaxis()->SetTitle("Charge [keV]");
    if(j==0) fHistCharge[j]->SetTitle("Drift inner layer");
    if(j==1) fHistCharge[j]->SetTitle("Drift outer layer");
    if(j==2) fHistCharge[j]->SetTitle("Strip inner layer");
    if(j==3) fHistCharge[j]->SetTitle("Strip inner layer");
    fHistCharge[j]->Draw();
    fHistCharge[j]->SetFillColor(7);
    fHistCharge[j]->Fit("gaus","QR","",70,100);
    fHistCharge[j]->Write();
  }

  //canvas MC spectra
  if(optMC){
    TCanvas *cspectraMC=new TCanvas("cspectraMC","cspectraMC",1000,800);
    cspectraMC->Divide(2,1);
    cspectraMC->cd(1);
    gPad->SetLogy();
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hSpectraPrimPosMC[i]->Draw("");
	hSpectraPrimPosMC[i]->SetTitle("ITS EvMC truth positive");
	hSpectraPrimPosMC[i]->SetMinimum(100);
	hSpectraPrimPosMC[i]->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy");
	hSpectraPrimPosMC[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hSpectraPrimPosMC[i]->Draw("same");
      hSpectraPrimPosMC[i]->SetLineColor(i+2);
      hSpectraPrimPosMC[i]->SetMarkerColor(i+2);
      hSpectraPrimPosMC[i]->SetMarkerStyle(21+i);
      hSpectraPrimPosMC[i]->Write();
    }
    cspectraMC->cd(2);
    gPad->SetLogy();
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hSpectraPrimNegMC[i]->Draw("");
	hSpectraPrimNegMC[i]->SetTitle("ITS MC truth negative");
	hSpectraPrimNegMC[i]->SetMinimum(100);
	hSpectraPrimNegMC[i]->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy");
	hSpectraPrimNegMC[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hSpectraPrimNegMC[i]->Draw("same");
      hSpectraPrimNegMC[i]->SetLineColor(i+2);
      hSpectraPrimNegMC[i]->SetMarkerColor(i+2);
      hSpectraPrimNegMC[i]->SetMarkerStyle(21+i);
      hSpectraPrimNegMC[i]->Write();
    }
    cspectraMC->Print(ps0.Data());
  }

  //dedx distribution to be fitted
  TCanvas *cgausPipos=new TCanvas("cgausPipos","PIONS pos",1000,800);
  cgausPipos->Divide(5,4,0.001,0.001);
  TCanvas *cgausKpos=new TCanvas("cgausKpos","KAONS pos",1000,800);
  cgausKpos->Divide(5,4,0.001,0.001);
  TCanvas *cgausPpos=new TCanvas("cgausPpos","PROTONS pos",1000,800);
  cgausPpos->Divide(5,4,0.001,0.001);
  TCanvas *cgausPineg=new TCanvas("cgausPineg","PIONS neg",1000,800);
  cgausPineg->Divide(5,4,0.001,0.001);
  TCanvas *cgausKneg=new TCanvas("cgausKneg","KAONS neg",1000,800);
  cgausKneg->Divide(5,4,0.001,0.001);
  TCanvas *cgausPneg=new TCanvas("cgausPneg","PROTONS neg",1000,800);
  cgausPneg->Divide(5,4,0.001,0.001);

  //binsize for one histo only because are all equal
  Float_t binsize = fHistPosPi[0]->GetBinWidth(1);
  Double_t fpar[5],efpar[5];
  for(Int_t i=0; i<nbins-2; i++){
    //------------------- positive particles
    //pions
    cgausPipos->cd(i+1);
    gPad->SetLogy();
    fHistPosPi[i]->Write();
    gres[0][i]=new TGraph();
    ITSsa->DoFit(fHistPosPi[i],i,211,gres[0][i]);
    ITSsa->FillHisto(hSpectraPos[0],i+1,binsize,211);
    if(optMC) ITSsa->FillHistoMC(hSpectraMCPIDPos[0],i+1,211,fHistMCPosPi[i]);
    if(ITSsa->IsGoodBin(i,211)){
      ITSsa->GetFitPar(fpar,efpar);
      hMeanPos[0]->SetBinContent(i+1,fpar[1]);
      hMeanPos[0]->SetBinError(i+1,efpar[1]);
      hSigmaPos[0]->SetBinContent(i+1,fpar[2]);
      hSigmaPos[0]->SetBinError(i+1,efpar[2]);
    }
    cgausPipos->Update();
    //kaons
    cgausKpos->cd(i+1);
    gPad->SetLogy();
    fHistPosK[i]->Write();
    gres[1][i]=new TGraph();
    ITSsa->DoFit(fHistPosK[i],i,321,gres[1][i]);
    ITSsa->FillHisto(hSpectraPos[1],i+1,binsize,321);
    if(optMC) ITSsa->FillHistoMC(hSpectraMCPIDPos[1],i+1,321,fHistMCPosK[i]);
    if(ITSsa->IsGoodBin(i,321)){
      ITSsa->GetFitPar(fpar,efpar);
      hMeanPos[1]->SetBinContent(i+1,fpar[1]);
      hMeanPos[1]->SetBinError(i+1,efpar[1]);
      hSigmaPos[1]->SetBinContent(i+1,fpar[2]);
      hSigmaPos[1]->SetBinError(i+1,efpar[2]);
    }
    cgausKpos->Update();
    //protons
    cgausPpos->cd(i+1);
    gPad->SetLogy();
    fHistPosP[i]->Write();
    fHistPosP[i]->SetFillColor(16);
    gres[2][i]=new TGraph();
    ITSsa->DoFitProton(fHistPosP[i],i,2212,gres[2][i]);
    ITSsa->FillHisto(hSpectraPos[2],i+1,binsize,2212);
    if(optMC) ITSsa->FillHistoMC(hSpectraMCPIDPos[2],i+1,2212,fHistMCPosP[i]);
    if(ITSsa->IsGoodBin(i,2212)){
      ITSsa->GetFitPar(fpar,efpar);
      hMeanPos[2]->SetBinContent(i+1,fpar[1]);
      hMeanPos[2]->SetBinError(i+1,efpar[1]);
      hSigmaPos[2]->SetBinContent(i+1,fpar[2]);
      hSigmaPos[2]->SetBinError(i+1,efpar[2]);
    }
    cgausPpos->Update();  

    //------------------- negative particles
    //pions
    cgausPineg->cd(i+1);
    gPad->SetLogy();
    fHistNegPi[i]->Write();
    gres[3][i]=new TGraph();
    ITSsa->DoFit(fHistNegPi[i],i,211,gres[3][i]);
    ITSsa->FillHisto(hSpectraNeg[0],i+1,binsize,211);
    if(optMC) ITSsa->FillHistoMC(hSpectraMCPIDNeg[0],i+1,211,fHistMCNegPi[i]);
    if(ITSsa->IsGoodBin(i,211)){
      ITSsa->GetFitPar(fpar,efpar);
      hMeanNeg[0]->SetBinContent(i+1,fpar[1]);
      hMeanNeg[0]->SetBinError(i+1,efpar[1]);
      hSigmaNeg[0]->SetBinContent(i+1,fpar[2]);
      hSigmaNeg[0]->SetBinError(i+1,efpar[2]);
    }
    cgausPineg->Update();
    //kaons
    cgausKneg->cd(i+1);
    gPad->SetLogy();
    fHistNegK[i]->Write();
    gres[4][i]=new TGraph();
    ITSsa->DoFit(fHistNegK[i],i,321,gres[4][i]);
    ITSsa->FillHisto(hSpectraNeg[1],i+1,binsize,321);
    if(optMC) ITSsa->FillHistoMC(hSpectraMCPIDNeg[1],i+1,321,fHistMCNegK[i]);
    if(ITSsa->IsGoodBin(i,321)){
      ITSsa->GetFitPar(fpar,efpar);
      hMeanNeg[1]->SetBinContent(i+1,fpar[1]);
      hMeanNeg[1]->SetBinError(i+1,efpar[1]);
      hSigmaNeg[1]->SetBinContent(i+1,fpar[2]);
      hSigmaNeg[1]->SetBinError(i+1,efpar[2]);
    }
    cgausKneg->Update();
    //protons
    cgausPneg->cd(i+1);
    gPad->SetLogy();
    fHistNegP[i]->Write();
    gres[5][i]=new TGraph();
    ITSsa->DoFitProton(fHistNegP[i],i,2212,gres[5][i]);
    ITSsa->FillHisto(hSpectraNeg[2],i+1,binsize,2212);
    if(optMC) ITSsa->FillHistoMC(hSpectraMCPIDNeg[2],i+1,2212,fHistMCNegP[i]);
    if(ITSsa->IsGoodBin(i,2212)){
      ITSsa->GetFitPar(fpar,efpar);
      hMeanNeg[2]->SetBinContent(i+1,fpar[1]);
      hMeanNeg[2]->SetBinError(i+1,efpar[1]);
      hSigmaNeg[2]->SetBinContent(i+1,fpar[2]);
      hSigmaNeg[2]->SetBinError(i+1,efpar[2]);
    }
    cgausPneg->Update(); 
  }

  //save histograms in the ps file
  cgausPipos->Print(ps1.Data());
  cgausKpos->Print(ps1.Data());
  cgausPpos->Print(ps1.Data());
  cgausPineg->Print(ps1.Data());
  cgausKneg->Print(ps1.Data());
  cgausPneg->Print(ps1.Data());

  //spectra REC
  TCanvas *cspecREC=new TCanvas("cspecREC","SPECTRA rec",1000,800);
  cspecREC->Divide(2,1);
  cspecREC->cd(1);
  gPad->SetLogy();
  gPad->SetGridy();
  for(Int_t i=0;i<3;i++){
    if(i==0){
      hSpectraPos[i]->Draw("text");
      hSpectraPos[i]->SetTitle("ITSsa RAW positive");
      hSpectraPos[i]->SetMinimum(10);
      hSpectraPos[i]->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy");
      hSpectraPos[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
    }
    hSpectraPos[i]->Draw("esame");
    hSpectraPos[i]->SetLineColor(i+2);
    hSpectraPos[i]->SetMarkerColor(i+2);
    hSpectraPos[i]->SetMarkerStyle(21+i);
    hSpectraPos[i]->Write();
  }
  TLegend *leg=new TLegend(0.51,0.11,0.84,0.35,NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  TLegendEntry *entry0=leg->AddEntry(hSpectraPos[0],"pions","p");
  entry0->SetTextColor(2);
  TLegendEntry *entry2=leg->AddEntry(hSpectraPos[1],"kaons","p");
  entry2->SetTextColor(3);
  TLegendEntry *entry4=leg->AddEntry(hSpectraPos[2],"protons","p");
  entry4->SetTextColor(4);
  leg->Draw("same");

  cspecREC->cd(2);
  gPad->SetLogy();
  gPad->SetGridy();
  for(Int_t i=0;i<3;i++){
    if(i==0){
      hSpectraNeg[i]->Draw("e");
      hSpectraNeg[i]->SetTitle("ITSsa RAW negative");
      hSpectraNeg[i]->SetMinimum(10);
      hSpectraNeg[i]->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy");
      hSpectraNeg[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
    }
    hSpectraNeg[i]->Draw("esame");
    hSpectraNeg[i]->SetLineColor(i+2);
    hSpectraNeg[i]->SetMarkerColor(i+2);
    hSpectraNeg[i]->SetMarkerStyle(21+i);
    hSpectraNeg[i]->Write();
  }


  //FitParameters
  TCanvas *cFitPar=new TCanvas("cFitPar","FitParameters",1000,800);
  cFitPar->Divide(2,3);
  TLatex** tplus=new TLatex*[3];
  TLatex** tminus=new TLatex*[3];
  tplus[0]=new TLatex(0.7,0.75,"#pi^{+}");
  tplus[1]=new TLatex(0.7,0.75,"K^{+}");
  tplus[2]=new TLatex(0.7,0.75,"p");
  tminus[0]=new TLatex(0.7,0.65,"#pi^{-}");
  tminus[1]=new TLatex(0.7,0.65,"K^{-}");
  tminus[2]=new TLatex(0.7,0.65,"#bar{p}");

  for(Int_t ipart=0; ipart<3; ipart++){
    tplus[ipart]->SetNDC();
    tplus[ipart]->SetTextColor(1);
    tplus[ipart]->SetTextSize(0.05);
    tminus[ipart]->SetNDC();
    tminus[ipart]->SetTextColor(4);
    tminus[ipart]->SetTextSize(0.05);
    cFitPar->cd(1+2*ipart);
    gPad->SetGridy();
    hMeanPos[ipart]->SetMarkerStyle(20);
    hMeanPos[ipart]->Draw("e");
    hMeanPos[ipart]->GetYaxis()->SetRangeUser(-1,1);
    hMeanNeg[ipart]->SetMarkerStyle(24);
    hMeanNeg[ipart]->SetMarkerColor(4);
    hMeanNeg[ipart]->SetLineColor(4);
    hMeanNeg[ipart]->Draw("esame");
    tplus[ipart]->Draw();
    tminus[ipart]->Draw();
    cFitPar->cd(2+2*ipart);
    gPad->SetGridy();
    hSigmaPos[ipart]->SetMarkerStyle(20);
    hSigmaPos[ipart]->Draw("e");
    hSigmaPos[ipart]->GetYaxis()->SetRangeUser(0.1,0.3);
    hSigmaNeg[ipart]->SetMarkerStyle(24);
    hSigmaNeg[ipart]->SetMarkerColor(4);
    hSigmaNeg[ipart]->SetLineColor(4);
    hSigmaNeg[ipart]->Draw("esame");
    hSigmaPos[ipart]->Write();
    hSigmaNeg[ipart]->Write();
    tplus[ipart]->Draw();
    tminus[ipart]->Draw();
  }

  //plus/minus ratio plots
  TH1F* hRatioPions=new TH1F(*hSpectraPos[0]);
  hRatioPions->Divide(hSpectraNeg[0]);
  TH1F* hRatioKaons=new TH1F(*hSpectraPos[1]);
  hRatioKaons->Divide(hSpectraNeg[1]);
  TH1F* hRatioProtons=new TH1F(*hSpectraPos[2]);
  hRatioProtons->Divide(hSpectraNeg[2]);

  TCanvas *cratios=new TCanvas("cratios","Ratios +/-",1000,800);
  cratios->Divide(1,3);
  cratios->cd(1);
  hRatioPions->SetMinimum(0.7);
  hRatioPions->SetMaximum(1.3);
  hRatioPions->GetXaxis()->SetTitle("p_{t} [GeV/c]");  
  hRatioPions->GetYaxis()->SetTitle("#pi^{+}/#pi^{-}");  
  hRatioPions->Draw("e");
  cratios->cd(2);
  hRatioKaons->SetMinimum(0.7);
  hRatioKaons->SetMaximum(1.3);
  hRatioKaons->GetXaxis()->SetTitle("p_{t} [GeV/c]");  
  hRatioKaons->GetYaxis()->SetTitle("K^{+}/K^{-}");  
  hRatioKaons->Draw("e");
  cratios->cd(3);
  hRatioProtons->SetMinimum(0.7);
  hRatioProtons->SetMaximum(1.3);
  hRatioProtons->GetXaxis()->SetTitle("p_{t} [GeV/c]");  
  hRatioProtons->GetYaxis()->SetTitle("p/#bar{p}");
  hRatioProtons->Draw("e");

  //efficiency and correction factor histograms
  if(optMC){
    for(Int_t i=0;i<3;i++){
      hEffPos[i] = (TH1F*)hSpectraPos[i]->Clone(Form("hEffPos%d",i));
      hEffPos[i]->SetTitle(Form("hEffPos%d",i));
      hEffPos[i]->Divide(hEffPos[i], hSpectraPrimPosMC[i], 1.0, 1.0, "B");//binomial errors
      hEffPos[i]->SetLineColor(i+2);
      hEffPos[i]->SetMarkerColor(i+2);
      hEffPos[i]->SetMarkerStyle(21+i);
      hEffPos[i]->Write();

      if(hSpectraPrimPosMCBefEvSel[i]){
	hCorrFacPos[i] = (TH1F*)hSpectraPos[i]->Clone(Form("hCorrFacPos%d",i));
	hCorrFacPos[i]->SetTitle(Form("hCorrFacPos%d",i));
	hCorrFacPos[i]->Divide(hCorrFacPos[i], hSpectraPrimPosMCBefEvSel[i], 1.0, 1.0, "B");//binomial errors
	hCorrFacPos[i]->SetLineColor(i+2);
	hCorrFacPos[i]->SetMarkerColor(i+2);
	hCorrFacPos[i]->SetMarkerStyle(25+i);     
	hCorrFacPos[i]->Write();
      }
      hEffMCPIDPos[i] = (TH1F*)hSpectraMCPIDPos[i]->Clone(Form("hEffMCPIDPos%d",i));
      hEffMCPIDPos[i]->SetTitle(Form("hEffMCPIDPos%d",i));
      hEffMCPIDPos[i]->Divide(hEffMCPIDPos[i], hSpectraPrimPosMC[i], 1.0, 1.0, "B");//binomial errors
      hEffMCPIDPos[i]->SetLineColor(i+2);
      hEffMCPIDPos[i]->SetMarkerColor(i+2);
      hEffMCPIDPos[i]->SetMarkerStyle(25+i);
      hEffMCPIDPos[i]->Write();

      if(hSpectraPrimPosMCBefEvSel[i]){
	hCorrFacMCPIDPos[i] = (TH1F*)hSpectraMCPIDPos[i]->Clone(Form("hCorrFacMCPIDPos%d",i));
	hCorrFacMCPIDPos[i]->SetTitle(Form("hCorrFacMCPIDPos%d",i));
	hCorrFacMCPIDPos[i]->Divide(hCorrFacMCPIDPos[i], hSpectraPrimPosMCBefEvSel[i], 1.0, 1.0, "B");//binomial errors
	hCorrFacMCPIDPos[i]->SetLineColor(i+2);
	hCorrFacMCPIDPos[i]->SetMarkerColor(i+2);
	hCorrFacMCPIDPos[i]->SetMarkerStyle(25+i);     
	hCorrFacMCPIDPos[i]->Write();
      }

      hEffNeg[i] = (TH1F*)hSpectraNeg[i]->Clone(Form("hEffNeg%d",i));
      hEffNeg[i]->SetTitle(Form("hEffNeg%d",i));
      hEffNeg[i]->Divide(hEffNeg[i], hSpectraPrimNegMC[i], 1.0, 1.0, "B");//binomial errors
      hEffNeg[i]->SetLineColor(i+2);
      hEffNeg[i]->SetMarkerColor(i+2);
      hEffNeg[i]->SetMarkerStyle(21+i);
      hEffNeg[i]->Write();

      if(hSpectraPrimNegMCBefEvSel[i]){
	hCorrFacNeg[i] = (TH1F*)hSpectraNeg[i]->Clone(Form("hCorrFacNeg%d",i));
	hCorrFacNeg[i]->SetTitle(Form("hCorrFacNeg%d",i));
	hCorrFacNeg[i]->Divide(hCorrFacNeg[i], hSpectraPrimNegMCBefEvSel[i], 1.0, 1.0, "B");//binomial errors
	hCorrFacNeg[i]->SetLineColor(i+2);
	hCorrFacNeg[i]->SetMarkerColor(i+2);
	hCorrFacNeg[i]->SetMarkerStyle(25+i);     
	hCorrFacNeg[i]->Write();
      }

      hEffMCPIDNeg[i] = (TH1F*)hSpectraMCPIDNeg[i]->Clone(Form("hEffMCPIDNeg%d",i));
      hEffMCPIDNeg[i]->SetTitle(Form("hEffMCPIDNeg%d",i));
      hEffMCPIDNeg[i]->Divide(hEffMCPIDNeg[i], hSpectraPrimNegMC[i], 1.0, 1.0, "B");//binomial errors
      hEffMCPIDNeg[i]->SetLineColor(i+2);
      hEffMCPIDNeg[i]->SetMarkerColor(i+2);
      hEffMCPIDNeg[i]->SetMarkerStyle(25+i);
      hEffMCPIDNeg[i]->Write();

      if(hSpectraPrimPosMCBefEvSel[i]){
	hCorrFacMCPIDNeg[i] = (TH1F*)hSpectraMCPIDPos[i]->Clone(Form("hCorrFacMCPIDNeg%d",i));
	hCorrFacMCPIDNeg[i]->SetTitle(Form("hCorrFacMCPIDNeg%d",i));
	hCorrFacMCPIDNeg[i]->Divide(hCorrFacMCPIDNeg[i], hSpectraPrimNegMCBefEvSel[i], 1.0, 1.0, "B");//binomial errors
	hCorrFacMCPIDNeg[i]->SetLineColor(i+2);
	hCorrFacMCPIDNeg[i]->SetMarkerColor(i+2);
	hCorrFacMCPIDNeg[i]->SetMarkerStyle(25+i);     
	hCorrFacMCPIDNeg[i]->Write();
      }

      hSpectraMCPIDPos[i]->SetLineColor(i+2);
      hSpectraMCPIDPos[i]->SetMarkerColor(i+2);
      hSpectraMCPIDPos[i]->SetMarkerStyle(25+i);
      hSpectraMCPIDPos[i]->Write();

      hSpectraMCPIDNeg[i]->SetLineColor(i+2);
      hSpectraMCPIDNeg[i]->SetMarkerColor(i+2);
      hSpectraMCPIDNeg[i]->SetMarkerStyle(25+i);
      hSpectraMCPIDNeg[i]->Write();
    }

    //spectra REC Ideal PID
    TCanvas *cspecREC2=new TCanvas("cspecIdPid","SPECTRA rec Ideal PID",1000,800);
    cspecREC2->Divide(2,1);
    cspecREC2->cd(1);
    gPad->SetLogy();
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hSpectraMCPIDPos[i]->Draw("text");
	hSpectraMCPIDPos[i]->SetTitle("ITSsa RAW positive - Ideal PID");
	hSpectraMCPIDPos[i]->SetMinimum(10);
	hSpectraMCPIDPos[i]->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy");
	hSpectraMCPIDPos[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hSpectraMCPIDPos[i]->Draw("esame");
    }
    leg->Draw("same");    
    cspecREC2->cd(2);
    gPad->SetLogy();
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hSpectraMCPIDNeg[i]->Draw("e");
	hSpectraMCPIDNeg[i]->SetTitle("ITSsa RAW negative -Ideal PID");
	hSpectraMCPIDNeg[i]->SetMinimum(10);
	hSpectraMCPIDNeg[i]->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy");
	hSpectraMCPIDNeg[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hSpectraMCPIDNeg[i]->Draw("esame");
    }

    // Efficiency
    TCanvas *ceff=new TCanvas("ceff","ceff",1000,800);
    ceff->Divide(2,1);
    ceff->cd(1);
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hEffPos[i]->SetTitle("ITSsa efficiency (positive)");
	hEffPos[i]->Draw("e");
	hEffPos[i]->GetYaxis()->SetRangeUser(0.1,0.9);
	hEffPos[i]->GetYaxis()->SetTitle("#epsilon");
	hEffPos[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hEffPos[i]->Draw("esame");
      hEffMCPIDPos[i]->Draw("esame");
    }
    TLegend *legeff=new TLegend(0.51,0.11,0.89,0.5,NULL,"brNDC");
    legeff->SetBorderSize(0);
    legeff->SetFillColor(0);
    TLegendEntry *entryeff0=legeff->AddEntry(hEffPos[0],"pions from fit","p");
    entryeff0->SetTextColor(2);
    TLegendEntry *entryeff2=legeff->AddEntry(hEffPos[1],"kaons from fit","p");
    entryeff2->SetTextColor(3);
    TLegendEntry *entryeff4=legeff->AddEntry(hEffPos[2],"protons from fit","p");
    entryeff4->SetTextColor(4);
    TLegendEntry *entryeff1=legeff->AddEntry(hEffMCPIDPos[0],"pions ideal PID","p");
    entryeff1->SetTextColor(2);
    TLegendEntry *entryeff3=legeff->AddEntry(hEffMCPIDPos[1],"kaons ideal PID","p");
    entryeff3->SetTextColor(3);
    TLegendEntry *entryeff5=legeff->AddEntry(hEffMCPIDPos[2],"protons ideal PID","p");
    entryeff5->SetTextColor(4);
    legeff->Draw("same");

    ceff->cd(2);
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hEffNeg[i]->SetTitle("ITSsa efficiency (negative)");
	hEffNeg[i]->Draw("e");
	hEffNeg[i]->GetYaxis()->SetRangeUser(0.1,0.9);
	hEffNeg[i]->GetYaxis()->SetTitle("#epsilon");
	hEffNeg[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hEffNeg[i]->Draw("esame");
      hEffMCPIDNeg[i]->Draw("esame");
    }
    ceff->Update();

    // Correction factor
    TCanvas *ccf=new TCanvas("ccf","Corr Factor",1000,800);
    ccf->Divide(2,1);
    ccf->cd(1);
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hEffPos[i]->SetTitle("ITSsa efficiency (positive)");
	hEffPos[i]->Draw("e");
	hEffPos[i]->GetYaxis()->SetRangeUser(0.1,0.9);
	hEffPos[i]->GetYaxis()->SetTitle("#epsilon");
	hEffPos[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hEffPos[i]->Draw("esame");
      if(hSpectraPrimPosMCBefEvSel[i]) hCorrFacPos[i]->Draw("esame");
    }
    TLegend *legcf=new TLegend(0.51,0.11,0.89,0.5,NULL,"brNDC");
    legcf->SetBorderSize(0);
    legcf->SetFillColor(0);
    TLegendEntry *entrycf0=legcf->AddEntry(hEffPos[0],"pions - eff.","p");
    entrycf0->SetTextColor(2);
    TLegendEntry *entrycf2=legcf->AddEntry(hEffPos[1],"kaons  - eff.","p");
    entrycf2->SetTextColor(3);
    TLegendEntry *entrycf4=legcf->AddEntry(hEffPos[2],"protons  - eff.","p");
    entrycf4->SetTextColor(4);
    TLegendEntry *entrycf1=legcf->AddEntry(hCorrFacPos[0],"pions - Corr. Fac.","p");
    entrycf1->SetTextColor(2);
    TLegendEntry *entrycf3=legcf->AddEntry(hCorrFacPos[1],"kaons - Corr. Fac.","p");
    entrycf3->SetTextColor(3);
    TLegendEntry *entrycf5=legcf->AddEntry(hCorrFacPos[2],"protons - Corr. Fac.","p");
    entrycf5->SetTextColor(4);
    legcf->Draw("same");

    ccf->cd(2);
    gPad->SetGridy();
    for(Int_t i=0;i<3;i++){
      if(i==0){
	hEffNeg[i]->SetTitle("ITSsa efficiency (negative)");
	hEffNeg[i]->Draw("e");
	hEffNeg[i]->GetYaxis()->SetRangeUser(0.1,0.9);
	hEffNeg[i]->GetYaxis()->SetTitle("#epsilon");
	hEffNeg[i]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
      }
      hEffNeg[i]->Draw("esame");
      if(hSpectraPrimNegMCBefEvSel[i]) hCorrFacNeg[i]->Draw("esame");
    }
    ccf->Update();

    cspecREC2->Print(ps1.Data());
    ceff->Print(ps1.Data());
    ccf->Print(ps1.Data());
  }

  //save canvas in the ps file
  cratios->Print(ps1.Data());
  cspecREC->Print(ps1.Data());
  cFitPar->Print(ps1.Data());
  cs->Print(ps1.Data());
  cDEDX->Print(ps1.Data());
  cDEDXdouble->Print(ps1.Data());
  cdummy->Print(ps2.Data());
  //	fout->Close();
  return;

}//end of the main 


//______________________________________________________________________
Double_t BetheBloch(Double_t *mom, Double_t *mass) {
  // BB PHobos parametrization fine tuned by G. Ortona
  // on data and MC in June 2010
  Double_t bg=mom[0]/mass[0];	
  Double_t par[5];
  Bool_t MC=kFALSE;
  if(MC){
    par[0]=1.39126e+02;//tuned on LHC10d (sim pass4)
    par[1]=2.33589e+01;
    par[2]=6.05219e-02;
    par[3]=2.04336e-01;
    par[4]=-4.99854e-04;
  }
  else{
    par[0]=5.33458e+04;//tuned on data pass6
    par[1]=1.65303e+01;
    par[2]=2.60065e-03;
    par[3]=3.59533e-04;
    par[4]=7.51168e-05;
  }
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t gamma=bg/beta;
  Double_t eff=1.0;
  if(bg<par[2]) eff=(bg-par[3])*(bg-par[3])+par[4];
  else eff=(par[2]-par[3])*(par[2]-par[3])+par[4];
  return (par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(beta*beta))*eff;
}

//______________________________________________________________________
void Logo(){
  //method to plot the ALICE logo in the propaganda plot
  TPaveText *tpave=new TPaveText(0.3,0.7,0.7,0.89,"brNDC");
  tpave->SetBorderSize(0);
  tpave->SetFillStyle(0);
  tpave->SetFillColor(0);
  TText *txt1=tpave->AddText("ALICE performance");
  TText *txt2=tpave->AddText("ITS stand-alone tracks");
  TText *txt3=tpave->AddText("pp @ #sqrt{s} = 7 TeV");
  txt1->SetTextFont(62);
  txt2->SetTextFont(62);
  txt3->SetTextFont(62);
  tpave->Draw();
  TPad *pad1=new TPad("pad1", "pad1", 0.75, 0.7, 0.89, 0.89);
  pad1->Draw();
  pad1->cd();
  TImage *img1 = TImage::Open("ALICElogo.gif");
  img1->FillRectangle(0);
  img1->Draw();
}

//EOF

