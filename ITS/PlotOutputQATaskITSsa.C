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
// dE/dx for ITS pure standalone tracks with parameterization
// dE/dx resolution for ITS pure standalone tracks with parameterization
// d0 resolution and bias with gaussian fit
// d0 resolution and bias with gaussian+tail fit
// Comparison of the two methods for d0 resolution
// Pt resolution (matching ITSèureSA with TPC+ITS tracks)
// General Plots: ratios between ITSsa, ITSpureSA and TPC+ITS tracks 
//                eta phi distributions of tracks
//                number of clusters per track
//
// Authors: Leonardo Milano, Francesco Prino
//
//-------------------------------------------------------

Int_t color[3]={1,2,4};
TString particle[3]={"Pion","Kaon","Proton"};


void PlotOutputQATaskITSsa(TString filename="LHC11h_QA95_MB_4cls.root",TString foutName="OutputQAITSstandaloneTracks_LHC11h_QA95_MB_4cls.root",Bool_t fMC=0){
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 

  if(filename.Contains("alien")) TGrid::Connect("alien:");
  TFile* f=TFile::Open(filename.Data());
  TDirectory *dirFile=(TDirectory*)f->Get("ITSsaTracks");
  TString lname="clistITSsaTracks";
  TList *li = (TList*)dirFile->Get(lname.Data());
  
  TFile *fout=new TFile(foutName.Data(),"recreate");
  fout->Close();
  delete fout;
  
  Bool_t useHybridITSsaParameterization=1; //BetheITSsaHybrid committed in Rev. 59830
  Bool_t PlotPIDBands=0; //Bands used for PID in the spectra analysis
  dedxPlot(li,fMC,foutName,useHybridITSsaParameterization,PlotPIDBands);
  dedxRes(li,foutName);
  ImpParRes(li,foutName);    
  ImpParResGausTail(li,foutName);    
  CompareImpParRes(foutName);
  Bool_t optFromMC=0;
  PlotPtResol(li,foutName,optFromMC);
  PlotITSsa(li,foutName);
}    



//----------------------------------------------------------------------------

void dedxPlot(TList *li,Bool_t fMC,TString foutName,Bool_t useHybridITSsaParameterization,Bool_t PlotPIDBands)    
{
  // Plot the dedx vs momentum from the ITSsaTracks QA output
  // Bethe Bloch parameterization is taken from AliITSPIDResponse
  // If useHybridITSsaParameterization the Hybrid ITS parameterization (Phobos + polinomial at low beta*gamma) is used. This requires AliRoot >  Rev. 59830
  
  TH2F *fHistDEDX=(TH2F*)li->FindObject("hdedxvsP2clsITSpureSA");
  fHistDEDX->Add((TH2F*)li->FindObject("hdedxvsP3clsITSpureSA"));
  fHistDEDX->Add((TH2F*)li->FindObject("hdedxvsP4clsITSpureSA"));
  
  TCanvas *c=new TCanvas("dEdx_vs_p","dEdx_vs_p",1);
  c->SetLogz();
  c->SetLogx();
  fHistDEDX->SetXTitle("#it{p} (GeV/#it{c})");
  fHistDEDX->SetYTitle("d#it{E}/dx (KeV/300#mum)");
  fHistDEDX->GetXaxis()->SetRangeUser(0.08,3);
  fHistDEDX->GetYaxis()->SetRangeUser(0.,700);
  fHistDEDX->GetYaxis()->SetTitleSize(0.05);
  fHistDEDX->GetXaxis()->SetTitleSize(0.05);
  //fHistDEDX->SetMinimum(2000);
  //fHistDEDX->SetMaximum(18000);
  fHistDEDX->DrawClone("col");
  
  Float_t pdgmass[5]={0.13957,0.493677,0.938272,1.875612762,0.00996}; //mass for pi, K, P, d (Gev/c^2)
  Double_t dedxvalue[5];
  Double_t funcvalue[7];
  TGraph *fdedx[7];
  AliITSPIDResponse *pidresp=new AliITSPIDResponse(fMC);
  
  for(Int_t i7=0;i7<7;i7++){
    fdedx[i7]=new TGraph();
    fdedx[i7]->SetName(Form("%d",i7));
    fdedx[i7]->SetTitle(Form("%d",i7));
    fdedx[i7]->SetLineWidth(2);
    fdedx[i7]->SetLineColor(1);
  }
  
  Float_t mom=0.07.;//low pt range
  Float_t step=(3-mom)/100;
  Double_t lr[7]={0.09,0.09,0.09,0.1,0.2,.25,.3};
  
  for(Int_t i=0;i<100;i++){
    if(useHybridITSsaParameterization)for(Int_t i2=0;i2<5;i2++)dedxvalue[i2]=pidresp->BetheITSsaHybrid(mom,pdgmass[i2]);
    else for(Int_t i2=0;i2<5;i2++)dedxvalue[i2]=pidresp->Bethe(mom,pdgmass[i2],1);
    funcvalue[0]=dedxvalue[0]-0.2*dedxvalue[0];
    funcvalue[1]=dedxvalue[0];
    funcvalue[2]=0.5*(dedxvalue[0]+dedxvalue[1]);
    funcvalue[3]=dedxvalue[1];
    funcvalue[4]=0.5*(dedxvalue[2]+dedxvalue[1]);
    funcvalue[5]=dedxvalue[2];
    funcvalue[6]=0.5*(dedxvalue[2]+dedxvalue[3]);
    
    for(Int_t i2=0;i2<7;i2++){
      fdedx[i2]->SetPoint(i,mom,funcvalue[i2]);
    }
    mom+=step;
  }
  fdedx[0]->SetLineStyle(2);
  fdedx[2]->SetLineStyle(2);
  fdedx[4]->SetLineStyle(2);
  fdedx[6]->SetLineStyle(2);
  for(Int_t i2=0;i2<7;i2++){
    fdedx[i2]->RemovePoint(0);
    fdedx[i2]->GetXaxis()->SetRangeUser(lr[i2],3);
    if(!PlotPIDBands)if(i2==0 || i2==2 || i2==4 || i2==6)continue;
    fdedx[i2]->DrawClone("Lsame");
  }
  TLatex *tex=new TLatex(0.55,630,"ITS standalone tracking");
  tex->SetTextSize(0.04);
  tex->DrawClone();
  
  TLatex *pi=new TLatex(0.1223724,176.0205,"#pi");
  TLatex *k=new TLatex(0.263512,291.1394,"K");
  TLatex *p=new TLatex(0.517675,296.7096,"p");
  TLatex *e=new TLatex(0.09856606,101.7502,"e");
  pi->DrawClone();
  k->DrawClone();
  p->DrawClone();
  e->DrawClone();
  
  //add electron line
  Float_t mom=0.09.;//low pt range
  Float_t step=(0.16-mom)/10;
  TGraph* fElectr=new TGraph();
  fElectr->SetName("fElectr");
  fElectr->SetTitle("fElectr");
  AliITSPIDResponse *pidresp=new AliITSPIDResponse(0);
  for(Int_t i=0;i<10;i++){
    if(useHybridITSsaParameterization)fElectr->SetPoint(i,mom,pidresp->BetheITSsaHybrid(mom,0.00051));
    else  fElectr->SetPoint(i,mom,pidresp->Bethe(mom,0.00051,1));
    mom+=step;
  }
  fElectr->SetLineWidth(2);
  fElectr->SetLineColor(1);
  fElectr->DrawClone("Lsame");
  
  TFile *fout=new TFile(foutName.Data(),"update");
  c->Write();
  fout->Close();
  delete fout;
  }

//----------------------------------------------------------------------------

void dedxRes(TList *li,TString foutName)    
{
  // plot the dedx resolution as a function of momuntum from the ITSsaTracks QA
  // Resolution is calculated from Gaussian fit to the pion peak for tracks with 2, 3 and 4 clusters in the ITS
  
  TCanvas *cdedxRes=new TCanvas("dedxRes_vs_p","dedxRes_vs_p",1000,800);
  cdedxRes->SetGridx();
  cdedxRes->SetLogx();
  
  for(Int_t ncls=2;ncls<5;ncls++){
    TH2F *hdedxvsPITSpureSA=(TH2F*)li->FindObject(Form("hdedxvsP%iclsITSpureSA",ncls));
    //binning
    const Int_t nbins = 29;
    Double_t xbins[nbins+1]={0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,
			     0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,
			     0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,1.90,2.00};
    TH1F *fHistdedxRes = new TH1F(Form("fHistdedxRes%icls",ncls),Form("fHistdedxRes%icls",ncls),nbins,xbins);
    fHistdedxRes->GetXaxis()->SetTitle("#font[52]{p} (GeV/#font[52]{c})");
    fHistdedxRes->GetYaxis()->SetTitle("de/dx relative resolution");
    fHistdedxRes->GetYaxis()->SetTitleOffset(1.5);
    TH1F *fHistDCA[nbins];
    TF1 *fPar = new TF1("fPar","gaus",0,1000);
    for(Int_t m=0;m<nbins;m++){
      fHistDCA[m]= (TH1F*)hdedxvsPITSpureSA->ProjectionY(Form("%i",m),hdedxvsPITSpureSA->GetXaxis()->FindBin(xbins[m]+0.0000001),hdedxvsPITSpureSA->GetXaxis()->FindBin(xbins[m+1]-0.0000001));
    }
    TCanvas *cgaus=new TCanvas(Form("cgausFitdEdxRes_%icls",ncls),Form("cgausFitdEdxRes_%icls",ncls),1000,800);
    cgaus->Divide(8,4,0.001,0.001);
    for(Int_t i=0; i<nbins; i++){
      cgaus->cd(i+1);
      fHistDCA[i]->SetFillColor(16);
      Float_t MaxPosition=fHistDCA[i]->GetBinCenter(fHistDCA[i]->GetMaximumBin());
      Float_t minFit=MaxPosition-0.1*MaxPosition;
      Float_t maxFit=MaxPosition+0.1*MaxPosition;
      fPar->SetParameter(1,MaxPosition);
      fHistDCA[i]->Fit(fPar,"NM","",minFit,maxFit);
      fHistDCA[i]->GetXaxis()->SetRangeUser(minFit-50,maxFit+50);
      fPar->SetLineColor(1);
      fHistDCA[i]->DrawClone();
      fPar->DrawClone("same");
      fHistdedxRes->Fill((xbins[i]+xbins[i+1])/2,fPar->GetParameter(2)/fPar->GetParameter(1));
      fHistdedxRes->SetBinError(fHistdedxRes->FindBin((xbins[i]+xbins[i+1])/2),fPar->GetParError(2)/fPar->GetParameter(1));
    }
    cdedxRes->cd();
    setDrawAtt(20,ncls,1,ncls,1,fHistdedxRes);
    fHistdedxRes->SetMinimum(0.0);
    fHistdedxRes->SetMaximum(0.2);
    fHistdedxRes->GetXaxis()->SetRangeUser(0.2,2);  
    if(ncls==2)fHistdedxRes->DrawCopy("p");
    else fHistdedxRes->DrawCopy("psame");
    
    TF1 *fFit = new TF1(Form("fFitdEdxRes_%d",ncls),"pol0",0.6,2);
    fFit->SetLineColor(ncls);
    fHistdedxRes->Fit(fFit,"NM","",0.6,2);
    fFit->DrawCopy("same");
    TPaveText *tpave1=new TPaveText(0.3,0.2-0.15*ncls,0.7,0.89,"brNDC");
    tpave1->SetBorderSize(0);
    tpave1->SetFillStyle(0);
    tpave1->SetFillColor(0);
    tpave1->SetTextColor(ncls);
    TText *txt2=tpave1->AddText(Form("pol 0: %.3f",fFit->GetParameter(0)));
    txt2->SetTextFont(62);
    tpave1->DrawClone();
    
    TFile *fout=new TFile(foutName.Data(),"update");
    fHistdedxRes->Write();
    fFit->Write();
    fout->Close();
    delete fout;
  }
  
  cdedxRes->BuildLegend()->SetFillColor(0);
  TPaveText *tpave=new TPaveText(0.5,0.9,0.8,0.99,"brNDC");
  tpave->SetBorderSize(0);
  tpave->SetFillStyle(0);
  tpave->SetFillColor(0);
  tpave->SetTextColor(2);
  TText *txt1=tpave->AddText("ITS standalone");
  txt1->SetTextFont(62);
  tpave->DrawClone();
  
  TFile *fout=new TFile(foutName.Data(),"update");
  cdedxRes->Write();
  fout->Close();
  delete fout;
}

//----------------------------------------------------------------------------

void ImpParRes(TList *li,TString foutName)    
{
  //Plot the resolution on the tranverse component of impact parmeter from the ITSsaTracks QA output as a function of pt
  // Gaussian fit is used

  TH2F *hd0rphiITSpureSA[3];
  //binning
  const Int_t nbins = 29;
  Double_t xbins[nbins+1]={0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,
			   0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,
			   0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,1.90,2.00};
  Int_t MinBin[3]={2,7,9};
  
  TH1F *fHistDCA[nbins];  
  TString particle[3]={"Pion","Kaon","Proton"};
  
  TCanvas *cImpPar=new TCanvas("ImpParRes_vs_pt_Gaus","ImpParRes_vs_pt_Gaus",1000,800);
  cImpPar->SetGridx();
  cImpPar->SetLogx();
  TCanvas *cImpParMean=new TCanvas("ImpParMean_vs_pt_Gaus","ImpParMean_vs_pt_Gaus",1000,800);
  cImpParMean->SetGridx();
  cImpParMean->SetLogx();
  
  for(Int_t iparticle=0;iparticle<3;iparticle++){
    hd0rphiITSpureSA[iparticle]=(TH2F*)li->FindObject(Form("hd0rphiITSpureSA%s",particle[iparticle].Data()));
    
    TH1F *fHistImpParRes = new TH1F(Form("fHistImpParResGaus%s",particle[iparticle].Data()),Form("%s",particle[iparticle].Data()),nbins,xbins);
    fHistImpParRes->GetXaxis()->SetTitle("#font[52]{p}_{T} (GeV/#font[52]{c})");
    fHistImpParRes->GetYaxis()->SetTitle("d0 r#phi resolution [#mum]");
    fHistImpParRes->GetYaxis()->SetTitleOffset(1.5);
    TH1F *fHistImpParMean = new TH1F(Form("fHistImpParMeanGaus%s",particle[iparticle].Data()),Form("%s",particle[iparticle].Data()),nbins,xbins);
    fHistImpParMean->GetXaxis()->SetTitle("#font[52]{p}_{T} (GeV/#font[52]{c})");
    fHistImpParMean->GetYaxis()->SetTitle("d0 r#phi mean [#mum]");
    fHistImpParMean->GetYaxis()->SetTitleOffset(1.5);
    
    TF1 *fPar = new TF1("fPar","gaus",-1,1);
    
    for(Int_t m=0;m<nbins;m++){
      fHistDCA[m]= (TH1F*)hd0rphiITSpureSA[iparticle]->ProjectionY(Form("%s%i",particle[iparticle].Data(),m),hd0rphiITSpureSA[iparticle]->GetXaxis()->FindBin(xbins[m]+0.000001),hd0rphiITSpureSA[iparticle]->GetXaxis()->FindBin(xbins[m+1]-0.000001));
    }
    
    TCanvas *cgaus=new TCanvas(Form("cgausFit_ImpParRes_vs_pt_%s",particle[iparticle].Data()),Form("cgausFit_ImpParRes_vs_pt_%s",particle[iparticle].Data()),1000,800);
    cgaus->Divide(8,4,0.001,0.001);
    Float_t sigma=0;
    
    for(Int_t i=0; i<nbins; i++){
      if(i<MinBin[iparticle])continue;
      
      cgaus->cd(i+1);
      gPad->SetLogy();
      fHistDCA[i]->SetLineColor(color[iparticle]);
      fHistDCA[i]->SetMarkerColor(color[iparticle]);
      fPar->SetLineColor(color[iparticle]);
      fHistDCA[i]->SetFillColor(16);
      Printf("first fit step [-1,1]");
      fHistDCA[i]->Fit(fPar,"NM","",-1,1);
      Printf("second fit step [-2sigma,2sigma]");
      Float_t nsigmas=2.;
      sigma=fPar->GetParameter(2);
      fHistDCA[i]->Fit(fPar,"NM","",fPar->GetParameter(1)-nsigmas*sigma,fPar->GetParameter(1)+nsigmas*sigma);
      Printf("Third fit step [-0.5 sigma,0.5 sigma]");
      nsigmas=.5; // narrow range, just want to fit the primaries
      fHistDCA[i]->Fit(fPar,"NM","",fPar->GetParameter(1)-nsigmas*sigma,fPar->GetParameter(1)+nsigmas*sigma);
      
      //set range to 3 sigmas (for the plot)
      sigma=fPar->GetParameter(2);
      fHistDCA[i]->GetXaxis()->SetRangeUser(fPar->GetParameter(1)-3*sigma,fPar->GetParameter(1)+3*sigma);
      fHistDCA[i]->DrawClone();
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
    
    setDrawAtt(iparticle+20,color[iparticle],1,color[iparticle],1,fHistImpParRes);
    fHistImpParRes->SetTitle(Form("%s - Gaus Fit",particle[iparticle].Data()));
    setDrawAtt(iparticle+20,color[iparticle],1,color[iparticle],1,fHistImpParMean);
    fHistImpParMean->SetTitle(Form("%s - Gaus Fit",particle[iparticle].Data()));
    
    cImpPar->cd();
    if(iparticle==0)fHistImpParRes->DrawCopy("p");
    else fHistImpParRes->DrawCopy("psame");
    
    cImpParMean->cd();
    if(iparticle==0)fHistImpParMean->DrawCopy("p");
    else fHistImpParMean->DrawCopy("psame");
    
    TFile *fout=new TFile(foutName.Data(),"update");
    cgaus->Write();
    fHistImpParRes->Write();
    fHistImpParMean->Write();
    fout->Close();
    delete fout;
  }
  
  cImpPar->BuildLegend()->SetFillColor(0);
  cImpParMean->BuildLegend()->SetFillColor(0);
  TFile *fout=new TFile(foutName.Data(),"update");
  cImpPar->Write();
  cImpParMean->Write();
  fout->Close();
  delete fout;
  
}

//----------------------------------------------------------------------------

void ImpParResGausTail(TList *li,TString foutName)    
{
  //Plot the resolution on the tranverse component of impact parmeter from the ITSsaTracks QA output as a function of pt
  // Gaussian + tail fit is used
 
  TH2F *hd0rphiITSpureSA[3];
  //binning
  const Int_t nbins = 29;
  Double_t xbins[nbins+1]={0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,
			   0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,
			   0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,1.90,2.00};
  TH1F *d0AllpointrphiSkip1_[nbins];  
  
  TCanvas *cImpPar=new TCanvas("ImpParResGausTail_vs_pt","ImpParResGausTail_vs_pt",1000,800);
  cImpPar->SetGridx();
  cImpPar->SetLogx();
  TCanvas *cImpParMean=new TCanvas("ImpParMeanGausTail_vs_pt","ImpParMeanGausTail_vs_pt",1000,800);
  cImpParMean->SetGridx();
  cImpParMean->SetLogx();
  
  
  //define the histgram
  TH1F **d0AllpointrphiSkipTail1_=new TH1F*[nbins]; 
  TH1F **d0AllpointrphiSkipGaus1_=new TH1F*[nbins];  
  TH1F **d0Pt_=new TH1F*[nbins];
  Float_t sigma=0;
  Double_t j =3.;  
  
  for(Int_t iparticle=0;iparticle<3;iparticle++){
    hd0rphiITSpureSA[iparticle]=(TH2F*)li->FindObject(Form("hd0rphiITSpureSA%s",particle[iparticle].Data()));
    
    TH1F *fHistImpParRes = new TH1F(Form("fHistImpParResGausTail%s",particle[iparticle].Data()),Form("%s",particle[iparticle].Data()),nbins,xbins);
    fHistImpParRes->GetXaxis()->SetTitle("#font[52]{p}_{T} (GeV/#font[52]{c})");
    fHistImpParRes->GetYaxis()->SetTitle("d0 r#phi resolution [#mum]");
    fHistImpParRes->GetYaxis()->SetTitleOffset(1.5);
    TH1F *fHistImpParMean = new TH1F(Form("fHistImpParMeanGausTail%s",particle[iparticle].Data()),Form("%s",particle[iparticle].Data()),nbins,xbins);
    fHistImpParMean->GetXaxis()->SetTitle("#font[52]{p}_{T} (GeV/#font[52]{c})");
    fHistImpParMean->GetYaxis()->SetTitle("d0 r#phi mean [#mum]");
    fHistImpParMean->GetYaxis()->SetTitleOffset(1.5);
    
    for(Int_t m=0;m<nbins;m++){
      d0AllpointrphiSkip1_[m]= (TH1F*)hd0rphiITSpureSA[iparticle]->ProjectionY(Form("%s%i",particle[iparticle].Data(),m),hd0rphiITSpureSA[iparticle]->GetXaxis()->FindBin(xbins[m]+0.000001),hd0rphiITSpureSA[iparticle]->GetXaxis()->FindBin(xbins[m+1]-0.000001));
    }
    
    TCanvas *cgaus=new TCanvas(Form("cgausTailFit_ImpParRes_vs_pt_%s",particle[iparticle].Data()),Form("cgausTailFit_ImpParRes_vs_pt_%s",particle[iparticle].Data()),1000,800);
    cgaus->Divide(8,4,0.001,0.001);
    
    for(Int_t i=0; i<nbins; i++){
      cgaus->cd(i+1);
      gPad->SetLogy();
      d0AllpointrphiSkip1_[i]->SetLineColor(1);
      d0AllpointrphiSkip1_[i]->SetMarkerColor(1);
      TF1 *h = new TF1("h","gaus",-1,1);
      d0AllpointrphiSkip1_[i]->Fit(h,"NM","",-1,1);
      Double_t d0rphirange_allpointSkip1 = h->GetParameter(2);   
      Double_t d0rphiMean_allpointSkip1 = h->GetParameter(1);   
      Double_t cutleft1= -j*(d0rphirange_allpointSkip1);
      Double_t cutright1 =j*d0rphirange_allpointSkip1;
      
      //fitting
      d0AllpointrphiSkipTail1_[i] = new TH1F(*d0AllpointrphiSkip1_[i]); 
      d0AllpointrphiSkipGaus1_[i] = new TH1F(*d0AllpointrphiSkip1_[i]); 
      d0AllpointrphiSkipTail1_[i]->Reset(0);
      d0AllpointrphiSkipGaus1_[i]->Reset(0);
      
      //Filling only with tail and Gaus
      for (Int_t bin=1;bin<d0AllpointrphiSkip1_[i]->GetNbinsX();bin++){
      	Float_t bincenter = d0AllpointrphiSkip1_[i]->GetBinCenter(bin);
	if(bincenter<cutleft1 || bincenter>cutright1) {
      	  d0AllpointrphiSkipTail1_[i]->SetBinContent(bin,d0AllpointrphiSkip1_[i]->GetBinContent(bin));
      	  d0AllpointrphiSkipTail1_[i]->SetBinError(bin,d0AllpointrphiSkip1_[i]->GetBinError(bin));
      	  d0AllpointrphiSkipGaus1_[i]->SetBinContent(bin,0.);  
      	  //This sentence is very important,otherwise we will get the information the data is empty when we fit it .
      	}
      	else if(bincenter>=cutleft1 && bincenter<=cutright1){
      	  d0AllpointrphiSkipTail1_[i]->SetBinContent(bin,0);
      	  d0AllpointrphiSkipGaus1_[i]->SetBinContent(bin,d0AllpointrphiSkip1_[i]->GetBinContent(bin));
      	  d0AllpointrphiSkipGaus1_[i]->SetBinError(bin,d0AllpointrphiSkip1_[i]->GetBinError(bin));  
      	}
      }
      d0AllpointrphiSkipGaus1_[i]->SetLineColor(2);
      d0AllpointrphiSkipGaus1_[i]->SetMarkerColor(2);
      d0AllpointrphiSkipTail1_[i]->SetLineColor(4);
      d0AllpointrphiSkipTail1_[i]->SetMarkerColor(4);
      
      TF1 *hh;
      hh =CreateFuncTail(d0AllpointrphiSkipTail1_[i],"hh");
      hh->SetLineColor(d0AllpointrphiSkipTail1_[i]->GetLineColor());
      d0AllpointrphiSkipTail1_[i]->Fit(hh,"NM","",-1,1);  
      Double_t Sigmatail_allpointSkip1 = hh->GetParameter(2);
      d0AllpointrphiSkipGaus1_[i]->Fit(h,"NM","",d0rphiMean_allpointSkip1-d0rphirange_allpointSkip1,d0rphiMean_allpointSkip1+d0rphirange_allpointSkip1);//narrow around the peak
      //d0AllpointrphiSkipGaus1_[i]->Fit(h,"NM","",-1,1);
      h->SetLineColor(d0AllpointrphiSkipGaus1_[i]->GetLineColor());
      Double_t Sigmagaus_allpointSkip1 = h->GetParameter(2);
      Double_t Meangaus_allpointSkip1 = h->GetParameter(1);
      Double_t Constgaus_allpointSkip1 = h->GetParameter(0);
      d0AllpointrphiSkipGaus1_[i]->DrawClone("");
      d0AllpointrphiSkipTail1_[i]->DrawClone("same");
      hh->DrawClone("same");
      h->DrawClone("same");
      
      TF1 * fPar=CreateFuncGaussTail(d0AllpointrphiSkip1_[i],"allpointSkip1",Constgaus_allpointSkip1,Meangaus_allpointSkip1,Sigmagaus_allpointSkip1,Sigmatail_allpointSkip1);
      d0AllpointrphiSkip1_[i]->Fit(fPar,"NM","",-1,1);
      fPar->SetLineColor(d0AllpointrphiSkip1_[i]->GetLineColor());
      fPar->DrawClone("same");
      
      sigma=fPar->GetParameter(2);
      fHistImpParRes->Fill((xbins[i]+xbins[i+1])/2,10000*fPar->GetParameter(2));
      fHistImpParRes->SetBinError(fHistImpParRes->FindBin((xbins[i]+xbins[i+1])/2),10000*fPar->GetParError(2));
      fHistImpParMean->Fill((xbins[i]+xbins[i+1])/2,10000*fPar->GetParameter(1));
      fHistImpParMean->SetBinError(fHistImpParMean->FindBin((xbins[i]+xbins[i+1])/2),10000*fPar->GetParError(1));
    }
    
    fHistImpParRes->SetMaximum(1000);
    fHistImpParRes->SetMinimum(0);
    fHistImpParMean->SetMaximum(80);
    fHistImpParMean->SetMinimum(-80);
    
    setDrawAtt(iparticle+20,color[iparticle],1,color[iparticle],1,fHistImpParRes);
    fHistImpParRes->SetTitle(Form("%s - GausTail Fit",particle[iparticle].Data()));
    setDrawAtt(iparticle+20,color[iparticle],1,color[iparticle],1,fHistImpParMean);
    fHistImpParMean->SetTitle(Form("%s - GausTail Fit",particle[iparticle].Data()));
   
    cImpPar->cd();
    if(iparticle==0)fHistImpParRes->DrawCopy("p");
    else fHistImpParRes->DrawCopy("psame");
    
    cImpParMean->cd();
    if(iparticle==0)fHistImpParMean->DrawCopy("p");
    else fHistImpParMean->DrawCopy("psame");
    
    
    TFile *fout=new TFile(foutName.Data(),"update");
    cgaus->Write();
    fHistImpParRes->Write();
    fHistImpParMean->Write();
    fout->Close();
    delete fout;
  }
  
  cImpPar->BuildLegend()->SetFillColor(0);
  cImpParMean->BuildLegend()->SetFillColor(0);
  TFile *fout=new TFile(foutName.Data(),"update");
  cImpPar->Write();
  cImpParMean->Write();
  fout->Close();
  delete fout;
  
}
//----------------------------------------------------------------------------

TF1 *CreateFuncTail(TH1F *hh,TString funcname,Double_t wholeRangeInt=-1.)
{
  TF1 *tail=new TF1(funcname.Data(),"[0]*(1./(2.*[2])*TMath::Exp(-TMath::Abs(x-[1])/[2]))",-1,1);
  Double_t binwidth=hh->GetBinWidth(10);
  Double_t integral=hh->Integral();
  if(wholeRangeInt!=-1.)tail->SetParLimits(0,(0.2)*wholeRangeInt*binwidth,(0.5)*wholeRangeInt*binwidth);
  Double_t RMS1=TMath::Abs(hh->GetRMS());
  Double_t firstvalue1=binwidth*integral;
  tail->SetParameters(1.,0,100.);//Set the initial value of parameter
  return tail;

} 
//----------------------------------------------------------------------------
TF1 *CreateFuncGaussTail(TH1F *h,TString funcname,Double_t Norm,Double_t parzero,Double_t parone,Double_t partwo)
{
  TF1 *gaustail=new TF1(funcname.Data(),"[0]*([4]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))+(1.-[4])/(2.*[3])*TMath::Exp(-TMath::Abs(x-[1])/[3]))",-1,1);
  gaustail->SetParameters(Norm,parzero,parone,partwo,0.82);//Set the initial value of parameter
  return gaustail;

} 

//----------------------------------------------------------------------------

void CompareImpParRes(TString foutName){
  //Copare the resolution obtained from Gaus and GausTail fit!
  // BE CAREFUL!! look at the fits "by eye", GausTail sometimes fail, especially pions at low pt
  
  TFile *fout=new TFile(foutName.Data(),"read");
  
  TCanvas *c=new TCanvas("CompareImpParRes","CompareImpParRes",1);
  c->SetLogx();
  
  for(Int_t iparticle=0;iparticle<3;iparticle++){
    TH1F *fHistImpParResGaus = (TH1F*)fout->Get(Form("fHistImpParResGaus%s",particle[iparticle].Data()));
    TH1F *fHistImpParResGausTail = (TH1F*)fout->Get(Form("fHistImpParResGausTail%s",particle[iparticle].Data()));
    setDrawAtt(iparticle+24,color[iparticle],1,color[iparticle],1,fHistImpParResGaus);
    if(iparticle==0)fHistImpParResGaus->DrawClone();
    else fHistImpParResGaus->DrawClone("same");
    fHistImpParResGausTail->DrawClone("same");
  }
  gPad->BuildLegend()->SetFillColor(0);
  
  TFile *fout=new TFile(foutName.Data(),"update");
  c->Write();
  fout->Close();
  delete fout;
  
}


//-----------------------------------------------------

void PlotPtResol(TList* l,TString foutName, Bool_t optFromMC){
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
	h1ptrelres[iSpec][iptbin]->Fit("gaus","","",-0.1,0.1);//,"L");
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
 
  TFile *fout=new TFile(foutName.Data(),"update");
  for(Int_t iSpec=0; iSpec<3; iSpec++){
    c2d[iSpec]->Write();
    c1dA[iSpec]->Write();
    c1dR[iSpec]->Write();
  }
  cb->Write();
  cr->Write();
  fout->Close();
  delete fout;
}
//----------------------------------------------------------------------------

void PlotITSsa(TList* li,TString foutName){

  TH1F* hPtTPCITS=(TH1F*)li->FindObject("hPtTPCITS");
  TH1F* hPtITSsa=(TH1F*)li->FindObject("hPtITSsa");
  TH1F* hPtITSpureSA=(TH1F*)li->FindObject("hPtITSpureSA");

  TH2F* hEtaPhiTPCITS=(TH2F*)li->FindObject("hEtaPhiTPCITS");
  TH2F* hEtaPhiITSsa=(TH2F*)li->FindObject("hEtaPhiITSsa");
  TH2F* hEtaPhiITSpureSA=(TH2F*)li->FindObject("hEtaPhiITSpureSA");

  TH1F* hChi2TPCITS=(TH1F*)li->FindObject("hChi2TPCITS");
  TH1F* hChi2ITSsa=(TH1F*)li->FindObject("hChi2ITSsa");
  TH1F* hChi2ITSpureSA=(TH1F*)li->FindObject("hChi2ITSpureSA");

  TH1F* hNcluTPCITS=(TH1F*)li->FindObject("hNcluTPCITS");
  TH1F* hNcluITSsa=(TH1F*)li->FindObject("hNcluITSsa");
  TH1F* hNcluITSpureSA=(TH1F*)li->FindObject("hNcluITSpureSA");


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
  
  TFile *fout=new TFile(foutName.Data(),"update");
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  fout->Close();
  delete fout;
 
}


//----------------------------------------------------------------------------

void setDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1)
{ 
  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
}
