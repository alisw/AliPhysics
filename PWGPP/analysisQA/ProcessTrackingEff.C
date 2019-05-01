#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TFile.h>
#include <TStyle.h>
#include <TMath.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "AliCFContainer.h"
#endif

void SetHistoStyle(TH1* h,  Int_t col, Int_t mar);
TH1D* GetEff(TDirectoryFile* d, TString typeEff, Int_t var, Double_t ptmin, Double_t ptmax);

Bool_t firstCall=kTRUE;

void ProcessTrackingEff(TString filname="AnalysisResults.root",
			TString pionDirName="PionFiltBit4",
			TString kaonDirName="KaonFiltBit4",
			TString protonDirName="ProtonFiltBit4",
			TString electronDirName="ElectronFiltBit4",
			TString deuteronDirName="DeuteronFiltBit4"
			){

  const Int_t totSteps=9;
  const Int_t nPtBins=3;
  Double_t ptLow[nPtBins]={-1.,0.3,3.};
  Double_t ptHigh[nPtBins]={999999.,0.5,10.};
  TString varname[6]={"Pt","Eta","Phi","Theta","Zvert","Mult"};
  Double_t theVar[totSteps]={0,2,2,1,1,4,4,5,5};
  Int_t thePtBin[totSteps]={0,1,2,1,2,1,2,1,2};
  Int_t colorPion[totSteps]={1,kGray+1,1,kGray+1,1,kGray+1,1,kGray+1,1};
  Int_t markerPion[totSteps]={20,20,24,20,24,20,24,20,24};
  Int_t colorKaon[totSteps]={2,kRed-7,2,kRed-7,2,kRed-7,2,kRed-7,2};
  Int_t markerKaon[totSteps]={21,21,25,21,25,21,25,21,25};
  Int_t colorProton[totSteps]={4,kBlue-7,4,kBlue-7,4,kBlue-7,4,kBlue-7,4};
  Int_t markerProton[totSteps]={33,33,27,33,27,33,27,33,27};
  Int_t colorElectron[totSteps]={kGreen+2,kSpring-5,kGreen+2,kSpring-5,kGreen+2,kSpring-5,kGreen+2,kSpring-5,kGreen+2};
  Int_t markerElectron[totSteps]={22,22,26,22,26,22,26,22,26};
  Int_t colorDeuteron[totSteps]={6,kMagenta+1,6,kMagenta+1,6,kMagenta+1,6,kMagenta+1,6};
  Int_t markerDeuteron[totSteps]={23,23,32,23,32,23,32,23,32};

  
  TFile* f = new TFile(filname.Data());
  TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TH1D* hEffPion[totSteps];
  TH1D* hEffKaon[totSteps];
  TH1D* hEffProton[totSteps];
  TH1D* hEffElectron[totSteps];
  TH1D* hEffDeuteron[totSteps];

  for(Int_t iStep=0; iStep<totSteps; iStep++){
    Int_t iVar=theVar[iStep];
    Int_t iPt=thePtBin[iStep];
    Double_t ptmin=ptLow[iPt];
    Double_t ptmax=ptHigh[iPt];
    TString ptrange="";
    if(iPt==1) ptrange="_LowPt";
    if(iPt==2) ptrange="_HighPt";
    hEffPion[iStep]=GetEff(d,pionDirName.Data(),iVar,ptmin,ptmax);
    if(hEffPion[iStep]){
      hEffPion[iStep]->SetName(Form("hEffPion%s%s",varname[iVar].Data(),ptrange.Data()));
      SetHistoStyle(hEffPion[iStep],colorPion[iStep],markerPion[iStep]);
    }
    hEffKaon[iStep]=GetEff(d,kaonDirName.Data(),iVar,ptmin,ptmax);
    if(hEffKaon[iStep]){
      hEffKaon[iStep]->SetName(Form("hEffKaon%s%s",varname[iVar].Data(),ptrange.Data()));
      SetHistoStyle(hEffKaon[iStep],colorKaon[iStep],markerKaon[iStep]);
    }
    hEffProton[iStep]=GetEff(d,protonDirName.Data(),iVar,ptmin,ptmax);
    if(hEffProton[iStep]){
      hEffProton[iStep]->SetName(Form("hEffProton%s%s",varname[iVar].Data(),ptrange.Data()));
      SetHistoStyle(hEffProton[iStep],colorProton[iStep],markerProton[iStep]);
    }
    hEffElectron[iStep]=GetEff(d,electronDirName.Data(),iVar,ptmin,ptmax);
    if(hEffElectron[iStep]){
      hEffElectron[iStep]->SetName(Form("hEffElectron%s%s",varname[iVar].Data(),ptrange.Data()));
      SetHistoStyle(hEffElectron[iStep],colorElectron[iStep],markerElectron[iStep]);
    }
    hEffDeuteron[iStep]=GetEff(d,deuteronDirName.Data(),iVar,ptmin,ptmax);
    if(hEffDeuteron[iStep]){
      hEffDeuteron[iStep]->SetName(Form("hEffDeuteron%s%s",varname[iVar].Data(),ptrange.Data()));
      SetHistoStyle(hEffDeuteron[iStep],colorDeuteron[iStep],markerDeuteron[iStep]);
    }
  }

  
  TH2F* hFramePt=new TH2F("hFramePt","",100,0.,30.,100.,0.,1.2);
  hFramePt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hFramePt->GetYaxis()->SetTitle("Efficiency");
  hFramePt->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFramePhi=new TH2F("hFramePhi","",100,0.,2.*TMath::Pi(),100.,0.,1.2);
  hFramePhi->GetXaxis()->SetTitle("#varphi (rad)");
  hFramePhi->GetYaxis()->SetTitle("Efficiency");
  hFramePhi->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameEta=new TH2F("hFrameEta","",100,-1.,1.,100.,0.,1.2);
  hFrameEta->GetXaxis()->SetTitle("#eta");
  hFrameEta->GetYaxis()->SetTitle("Efficiency");
  hFrameEta->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameZvert=new TH2F("hFrameZvert","",100,-12.,12.,100.,0.,1.2);
  hFrameZvert->GetXaxis()->SetTitle("z_{vertex} (cm)");
  hFrameZvert->GetYaxis()->SetTitle("Efficiency");
  hFrameZvert->GetYaxis()->SetTitleOffset(1.2);

  TH2F* hFrameMult=new TH2F("hFrameMult","",10000,0.,10000.,100.,0.,1.2);
  hFrameMult->GetXaxis()->SetTitle("Multiplicity");
  hFrameMult->GetYaxis()->SetTitle("Efficiency");
  hFrameMult->GetYaxis()->SetTitleOffset(1.2);

  TLatex* tcuts=new TLatex(0.7,0.2,"#splitline{|#eta|<0.8}{|z_{vert}|<10 cm}");
  tcuts->SetNDC();
  tcuts->SetTextFont(43);
  tcuts->SetTextSize(26);
  
  TCanvas* cept=new TCanvas("cept","EffVsPt",1600,800);
  cept->Divide(2,1);
  cept->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hFramePt->Draw();
  TLegend* legpt=new TLegend(0.7,0.3,0.89,0.55);
  cept->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetLogx();
  hFramePt->Draw();
  for(Int_t iStep=0; iStep<totSteps; iStep++){
    Int_t iPt=thePtBin[iStep];
    Double_t ptmin=ptLow[iPt];
    Double_t ptmax=ptHigh[iPt];
    if(theVar[iStep]==0){
      if(hEffPion[iStep]){
	cept->cd(1);
	hEffPion[iStep]->Draw("same");
	cept->cd(2);
	hEffPion[iStep]->Draw("same");
	legpt->AddEntry(hEffPion[iStep],"#pi","P");
      }
      if(hEffKaon[iStep]){
	cept->cd(1);
	hEffKaon[iStep]->Draw("same");
	cept->cd(2);
	hEffKaon[iStep]->Draw("same");
 	legpt->AddEntry(hEffKaon[iStep],"K","P");
      }
      if(hEffProton[iStep]){
	cept->cd(1);
	hEffProton[iStep]->Draw("same");
	cept->cd(2);
	hEffProton[iStep]->Draw("same");
 	legpt->AddEntry(hEffProton[iStep],"p","P");
      }
      if(hEffElectron[iStep]){
	cept->cd(1);
	hEffElectron[iStep]->Draw("same");
 	cept->cd(2);
	hEffElectron[iStep]->Draw("same");
 	legpt->AddEntry(hEffElectron[iStep],"e","P");
      }
      if(hEffDeuteron[iStep]){
	cept->cd(1);
	hEffDeuteron[iStep]->Draw("same");
	cept->cd(2);
	hEffDeuteron[iStep]->Draw("same");
 	legpt->AddEntry(hEffDeuteron[iStep],"d","P");
      }
    }
  }
  cept->cd(1);
  legpt->Draw();
  tcuts->Draw();
  cept->cd(2);
  legpt->Draw();
  tcuts->Draw();

  cept->SaveAs("EfficVsPt.png");

  TCanvas** cevar=new TCanvas*[4];
  cevar[0]=new TCanvas("cephi","EffVsPhi",1600,800);
  cevar[1]=new TCanvas("ceeta","EffVsEta",1600,800);
  cevar[2]=new TCanvas("cezv","EffVsZvert",1600,800);
  cevar[3]=new TCanvas("cemult","EffVsMult",1600,800);

  TLegend** legvar=new TLegend*[8];
  for(Int_t iv=0; iv<4; iv++){
    cevar[iv]->Divide(2,1);
    for(Int_t ip=1; ip<=2; ip++){
      cevar[iv]->cd(ip);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.08);
      gPad->SetTickx();
      gPad->SetTicky();
      if(iv==0) hFramePhi->Draw();
      else if(iv==1) hFrameEta->Draw();
      else if(iv==2) hFrameZvert->Draw();
      else if(iv==3) hFrameMult->Draw();
      legvar[iv*2+(ip-1)]=new TLegend(0.16,0.73,0.87,0.89);
      legvar[iv*2+(ip-1)]->SetNColumns(2);
      legvar[iv*2+(ip-1)]->SetColumnSeparation(0.15);
      legvar[iv*2+(ip-1)]->SetMargin(0.15);
    }
  }

  for(Int_t iStep=0; iStep<totSteps; iStep++){
    Int_t iPt=thePtBin[iStep];
    Double_t ptmin=ptLow[iPt];
    Double_t ptmax=ptHigh[iPt];
    Int_t iCanv=-1;
    Int_t iPad=-1;
    if(theVar[iStep]==2) iCanv=0;
    else if(theVar[iStep]==1) iCanv=1;
    else if(theVar[iStep]==4) iCanv=2;
    else if(theVar[iStep]==5) iCanv=3;
    if(thePtBin[iStep]==1) iPad=1;
    else if(thePtBin[iStep]==2) iPad=2;
    if(iCanv<0 || iPad<0) continue;
    cevar[iCanv]->cd(iPad);
    if(hEffPion[iStep]){
      hEffPion[iStep]->Draw("same");
      legvar[2*iCanv+iPad-1]->AddEntry(hEffPion[iStep],Form("#pi, %.1f<p_{T}<%.1f GeV/c",ptmin,ptmax),"P");	
      if(theVar[iStep]==5) hFrameMult->GetXaxis()->SetRangeUser(0.,hEffPion[iStep]->GetXaxis()->GetXmax());
    }
    if(hEffKaon[iStep]){
      hEffKaon[iStep]->Draw("same");
      legvar[2*iCanv+iPad-1]->AddEntry(hEffKaon[iStep],Form("K, %.1f<p_{T}<%.1f GeV/c",ptmin,ptmax),"P");
    }
    if(hEffProton[iStep]){
      hEffProton[iStep]->Draw("same");
      legvar[2*iCanv+iPad-1]->AddEntry(hEffProton[iStep],Form("p, %.1f<p_{T}<%.1f GeV/c",ptmin,ptmax),"P");
    }
    if(hEffDeuteron[iStep]){
      hEffDeuteron[iStep]->Draw("same");
      legvar[2*iCanv+iPad-1]->AddEntry(hEffDeuteron[iStep],Form("d, %.1f<p_{T}<%.1f GeV/c",ptmin,ptmax),"P");
    }
    if(hEffElectron[iStep]){
      hEffElectron[iStep]->Draw("same");
      legvar[2*iCanv+iPad-1]->AddEntry(hEffElectron[iStep],Form("e, %.1f<p_{T}<%.1f GeV/c",ptmin,ptmax),"P");
    }
  }
  for(Int_t iv=0; iv<4; iv++){
    cevar[iv]->cd(1);
    legvar[iv*2]->Draw();
    tcuts->Draw();
    cevar[iv]->cd(2);
    legvar[iv*2+1]->Draw();
    tcuts->Draw();
  }
  cevar[0]->SaveAs("EfficVsPhi.png");
  cevar[1]->SaveAs("EfficVsEta.png");
  cevar[2]->SaveAs("EfficVsZvert.png");
  cevar[3]->SaveAs("EfficVsMult.png");

  TFile* filout=new TFile("TrackingEffs.root","recreate");
  for(Int_t iStep=0; iStep<totSteps; iStep++){
    if(hEffPion[iStep]) hEffPion[iStep]->Write();
    if(hEffKaon[iStep]) hEffKaon[iStep]->Write();
    if(hEffProton[iStep]) hEffProton[iStep]->Write();
    if(hEffElectron[iStep]) hEffElectron[iStep]->Write();
    if(hEffDeuteron[iStep]) hEffDeuteron[iStep]->Write();
  }
  filout->Close();
}

TH1D* GetEff(TDirectoryFile* d, TString typeEff, Int_t var, Double_t ptmin, Double_t ptmax){

  if(typeEff=="") return 0x0;
  AliCFContainer *data = (AliCFContainer*) (d->Get(Form("container%s",typeEff.Data())));
  if(!data) return 0x0;
  if(firstCall) data->Print("");
  firstCall=kFALSE;

  AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(6); // Reco
  THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
  AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
  THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();
  Int_t ndim=numData->GetNdimensions();


  printf("---- Project CF %s vs. %s ----\n",typeEff.Data(),numData->GetAxis(var)->GetTitle());

  for(Int_t j=0; j<ndim; j++){
    TAxis* ax=numData->GetAxis(j);
    Int_t nbins=ax->GetNbins();
    if(j!=var) ax->SetRange(0,nbins+1);
    else ax->SetRange(1,nbins);
    ax=denData->GetAxis(j);
    if(j!=var) ax->SetRange(0,nbins+1);
    else ax->SetRange(1,nbins);
  }

  if(var!=0 && ptmin>0 && ptmax<100){
    TAxis* ax=numData->GetAxis(0);
    Int_t bl=ax->FindBin(ptmin*1.0001);
    Int_t bu=ax->FindBin(ptmax*0.9999);
    printf("Pt bins considered: %d-%d\n",bl,bu);
    ax->SetRange(bl,bu);
    TAxis* ax2=denData->GetAxis(0);
    ax2->SetRange(bl,bu);
  }

  //cut at |eta|<0.8 if projecting vs. other variables
  if(var!=1 && var!=3){
    TAxis* ax1=numData->GetAxis(1);
    Int_t bl=ax1->FindBin(-0.7999);
    Int_t bu=ax1->FindBin(0.7999);
    printf("Eta bins considered: %d-%d\n",bl,bu);
    ax1->SetRange(bl,bu);
    TAxis* ax1d=denData->GetAxis(1);
    ax1d->SetRange(bl,bu);
  }

  // cut zvertex at 10 cm
  TAxis* ax4=numData->GetAxis(4);
  Int_t bl=ax4->FindBin(-9.9999);
  Int_t bu=ax4->FindBin(9.9999);
  ax4->SetRange(bl,bu);
  TAxis* ax4d=denData->GetAxis(4);
  ax4d->SetRange(bl,bu);
  printf("zvertex bins considered: %d-%d\n",bl,bu);



  TH1D* hNumer=numData->Projection(var);
  TH1D* hDenom=denData->Projection(var);
  TH1D* hEff=(TH1D*)hNumer->Clone("hEff");
  hEff->Divide(hNumer,hDenom,1,1,"B");
  delete hNumer;
  delete hDenom;
  hEff->GetYaxis()->SetTitle("Efficiency");
  if(var==0) hEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  else if(var==1) hEff->GetXaxis()->SetTitle("#eta");
  else if(var==2) hEff->GetXaxis()->SetTitle("#varphi (rad)");
  else if(var==3) hEff->GetXaxis()->SetTitle("#theta (rad)");
  else if(var==4) hEff->GetXaxis()->SetTitle("z_{Vert} (cm)");
  else if(var==5) hEff->GetXaxis()->SetTitle("Multiplicity");
  else if(var==6) hEff->GetXaxis()->SetTitle("Centrality");
  hEff->SetMinimum(0.);
  hEff->SetMaximum(1.2);
  hEff->GetYaxis()->SetTitleOffset(1.1);

  return hEff;
}


void SetHistoStyle(TH1* h,  Int_t col, Int_t mar){
  h->SetMarkerStyle(mar);
  h->SetMarkerColor(col);
  h->SetLineColor(col);
}
