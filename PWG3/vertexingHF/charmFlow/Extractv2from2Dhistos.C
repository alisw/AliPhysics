#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TDatabasePDG.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>


#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#include "AliEventPlaneResolution.h"
#endif


// Common variables: to be configured by the user
TString filname="AnalysisResultsgoodruns.root";
TString listname="coutputv2D0Std";//"coutputv2D0RAACuts"; //"coutputv2D0Std";
TString filcutsname="AnalysisResultsgoodruns.root";
Int_t minCent=30;
Int_t maxCent=50;
Int_t mesonPDG=421;

const Int_t nFinalPtBins=3;
Double_t ptlims[nFinalPtBins+1]={2.,5.,8.,12.};
Double_t sigmaRangeForSig=2.5;
Double_t sigmaRangeForBkg=4.5;
Int_t rebinHistoSideBands[nFinalPtBins]={2,2,2};
Bool_t useConstantvsBkgVsMass=kFALSE;
Int_t rebinHistov2Mass[nFinalPtBins]={2,2,5};
Int_t factor4refl=0;
Int_t minPtBin[nFinalPtBins]={-1,-1,-1};
Int_t maxPtBin[nFinalPtBins]={-1,-1,-1};

// systematic errors for 2-5, 5-8 and 8-12 (no run-by-run weights)
/*
Double_t systErrMeth1[nFinalPtBins]={
  (0.308-0.169)/2.,
  (0.14-0.1)/2.,
  (0.04+0.02)/2.
};
Double_t systErrMeth2[nFinalPtBins]={
  (0.305-0.252)/2.,
  (0.129-0.020)/2.,
  (0.101+0.06)/2.
};
*/

// systematic errors for 2-5, 5-8 and 8-12 (93 runs with run-by-run weights)
Double_t systErrMeth1[nFinalPtBins]={
  (0.23-0.10)/2.,
  (0.152-0.078)/2.,
  (0.161-0.097)/2.
};
Double_t systErrMeth2[nFinalPtBins]={
  (0.164-0.097)/2.,
  (0.110-0.012)/2.,
  (0.131-0.036)/2.
};


// systematic errors for 2-5, 5-8 and 8-12 (93 runs with run-by-run weights, RAA cuts)
/*
Double_t systErrMeth1[nFinalPtBins]={
  (0.265-0.122)/2.,
  (0.165-0.117)/2.,
  (0.238-0.169)/2.
};
Double_t systErrMeth2[nFinalPtBins]={
  (0.174-0.135)/2.,
  (0.18-0.11)/2.,
  (0.311-0.28)/2.
};
*/

// output of mass fitter
Int_t hMinBin;
Int_t hMaxBin;
Double_t signalFromFit,esignalFromFit;
Double_t massFromFit,sigmaFromFit;
TF1* fB1=0x0;
TF1* fB2=0x0;
TF1* fSB=0x0;

void LoadMassHistos(TList* lst, TH2F** hMassDphi);
Bool_t DefinePtBins(TDirectoryFile* df);
Double_t v2vsMass(Double_t *x, Double_t *par);
Double_t GetEPResolution(TList* lst, Double_t &rcflow, Double_t &rcfhigh);
Bool_t FitMassSpectrum(TH1F* hMass, TPad* pad);
TH1F* DoSideBands(Int_t iFinalPtBin, TH2F* hMassDphi, TH1F* hMass, Int_t rebin, TCanvas* c1, Int_t optBkg=0);
TF1* DoFitv2vsMass(Int_t iFinalPtBin, TH2F* hMassDphi, TH1F* hMass, Int_t rebin, TCanvas* c2, Int_t optErrors=0, Bool_t useConst=kTRUE);


void Extractv2from2Dhistos(){

  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");

  TFile* filcuts=new TFile(filcutsname.Data());
  TDirectoryFile* dfcuts=(TDirectoryFile*)filcuts->Get("PWG3_D2H_HFv2");
  Bool_t binOK=DefinePtBins(dfcuts);
  if(!binOK) return;

  TFile* fil=new TFile(filname.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("PWG3_D2H_HFv2");
  TList* lst=(TList*)df->Get(listname.Data());
  if(!lst){
    printf("ERROR: list %s not found in file\n",listname.Data());
    return;
  }
  Double_t rcfmin,rcfmax;
  Double_t resolFull=GetEPResolution(lst,rcfmin,rcfmax);
  Double_t resolSyst=(rcfmax-rcfmin)/2./resolFull;
  printf("Relative Systematic Error on RCF=%f\n",resolSyst);

  TH2F** hMassDphi=new TH2F*[nFinalPtBins];
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    hMassDphi[iFinalPtBin]=0x0;
  }
  LoadMassHistos(lst, hMassDphi);

  TCanvas** c1=new TCanvas*[nFinalPtBins];
  TCanvas** c2=new TCanvas*[nFinalPtBins];

  Double_t v2M1[nFinalPtBins],errv2M1[nFinalPtBins];
  Double_t v2M2[nFinalPtBins],errv2M2[nFinalPtBins];

  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);

  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("**************** BIN %d ******************\n",iFinalPtBin);
    printf("\n--------- Method 1: Side Bands ----------\n");
    TH1F* hMass=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionY();

    c1[iFinalPtBin]=new TCanvas(Form("cMeth1PtBin%d",iFinalPtBin),Form("Meth1PtBin%d",iFinalPtBin));
    c1[iFinalPtBin]->Divide(2,2);
    c1[iFinalPtBin]->cd(1);
    hMassDphi[iFinalPtBin]->Draw("colz");
    c1[iFinalPtBin]->cd(2);
    Bool_t out=FitMassSpectrum(hMass,(TPad*)gPad);
    if(!out) continue;

    TH1F* hCos2PhiSig=DoSideBands(iFinalPtBin,hMassDphi[iFinalPtBin],hMass,rebinHistoSideBands[iFinalPtBin],c1[iFinalPtBin]);
    Double_t v2obsM1=hCos2PhiSig->GetMean();
    Double_t errv2obsM1=hCos2PhiSig->GetMeanError();  
    printf("v2obs = %f +- %f\n",v2obsM1,errv2obsM1);
    v2M1[iFinalPtBin]=v2obsM1/resolFull;
    errv2M1[iFinalPtBin]=errv2obsM1/resolFull;
    printf("v2 = %f +- %f\n",v2M1[iFinalPtBin],errv2M1[iFinalPtBin]);
    
    printf("\n--------- Method 2: S/S+B ----------\n");
    c2[iFinalPtBin]=new TCanvas(Form("cMeth2Bin%d",iFinalPtBin),Form("Meth2Bin%d",iFinalPtBin));
    TF1* fv2=DoFitv2vsMass(iFinalPtBin,hMassDphi[iFinalPtBin],hMass,rebinHistov2Mass[iFinalPtBin],c2[iFinalPtBin],0,useConstantvsBkgVsMass);

    Double_t v2obsM2=fv2->GetParameter(3);
    Double_t errv2obsM2=fv2->GetParError(3);
    printf("v2obs = %f +- %f\n",v2obsM2,errv2obsM2);
    v2M2[iFinalPtBin]=v2obsM2/resolFull;
    errv2M2[iFinalPtBin]=errv2obsM2/resolFull;
    printf("v2 = %f +- %f\n",v2M2[iFinalPtBin],errv2M2[iFinalPtBin]);
    c1[iFinalPtBin]->SaveAs(Form("cMeth1Bin%d.root",iFinalPtBin));
    c2[iFinalPtBin]->SaveAs(Form("cMeth2Bin%d.root",iFinalPtBin));
  }

  printf("\n--------- Summary ------------\n");

  TH1F* hv2m1=new TH1F("hv2m1","Side Band subtraction",nFinalPtBins,ptlims);
  TH1F* hv2m2=new TH1F("hv2m2","Fit to v2 vs. mass",nFinalPtBins,ptlims);
   for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("PtBin %d   v2method1 = %f +- %f   v2method2 = %f +-%f\n",iFinalPtBin,
	   v2M1[iFinalPtBin],errv2M1[iFinalPtBin],
	   v2M2[iFinalPtBin],errv2M2[iFinalPtBin]
	   );
    hv2m1->SetBinContent(iFinalPtBin+1,v2M1[iFinalPtBin]);
    hv2m1->SetBinError(iFinalPtBin+1,errv2M1[iFinalPtBin]);
    hv2m2->SetBinContent(iFinalPtBin+1,v2M2[iFinalPtBin]);
    hv2m2->SetBinError(iFinalPtBin+1,errv2M2[iFinalPtBin]);
  }
    
   TH1F* hSystErr1=(TH1F*)hv2m1->Clone("hSystErr1");
   TH1F* hSystErr2=(TH1F*)hv2m2->Clone("hSystErr2");
   for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
     Double_t systRes1=resolSyst*v2M1[iFinalPtBin];
     Double_t systRes2=resolSyst*v2M2[iFinalPtBin];
     printf("%f\n",systRes1);
     Double_t syste1=TMath::Sqrt(systErrMeth1[iFinalPtBin]*systErrMeth1[iFinalPtBin]+systRes1*systRes1);
     Double_t syste2=TMath::Sqrt(systErrMeth2[iFinalPtBin]*systErrMeth2[iFinalPtBin]+systRes2*systRes2);
     hSystErr1->SetBinError(iFinalPtBin+1,syste1);
     hSystErr2->SetBinError(iFinalPtBin+1,syste2);
   }


   Double_t maxy=TMath::Max(hv2m2->GetMaximum(),hv2m1->GetMaximum())+0.1;
   Double_t miny=TMath::Min(hv2m2->GetMinimum(),hv2m1->GetMinimum())-0.1;
   TH2F* hempty=new TH2F("hempty","",10,0.,hv2m1->GetXaxis()->GetXmax()+2.,10,miny,maxy);
   hempty->GetXaxis()->SetTitle("p_{t} (GeV/c)");
   hempty->GetYaxis()->SetTitle("v_{2}");

   TCanvas* cv2=new TCanvas("cv2","v2");
   hempty->Draw();
   hv2m1->SetMarkerStyle(26);
   hSystErr1->SetFillColor(kGray);
   hSystErr1->SetFillStyle(3002);
   hSystErr1->Draw("e2same");

   hv2m2->SetLineColor(kRed+1);
   hv2m2->SetMarkerColor(kRed+1);
   hv2m2->SetMarkerStyle(20);
   hSystErr2->SetFillColor(kRed-9);
   hSystErr2->SetFillStyle(3005);
   hSystErr2->Draw("e2same");

   hv2m1->Draw("same");
   hv2m2->Draw("same");

   TLegend* leg2=new TLegend(0.5,0.7,0.89,0.89);
   leg2->SetFillStyle(0);
   TLegendEntry* ent=leg2->AddEntry(hv2m1,"Side Band subtraction","P");
   ent->SetTextColor(hv2m1->GetMarkerColor());
   ent=leg2->AddEntry(hv2m2,"Fit to v2 vs. mass","P");
   ent->SetTextColor(hv2m2->GetMarkerColor());
   leg2->Draw();
   cv2->Update();
   cv2->SaveAs("Dzero-v2-2Dmethods.gif");

   TFile* outfil=new TFile("Dzero-v2-2Dmethods.root","recreate");
   outfil->cd();
   hv2m1->Write();
   hv2m2->Write();
   hSystErr1->Write();
   hSystErr2->Write();
   outfil->Close();
}


Double_t v2vsMass(Double_t *x, Double_t *par){
  // Fit function for signal+background
  // par[0] = S/B at the mass peak
  // par[1] = mass
  // par[2] = sigma
  // par[3] = v2sig
  // par[4] = v2back

  Double_t fracsig=par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
  Double_t fracbkg=1-fracsig;
  return par[3]*fracsig+par[4]*fracbkg;
}

Double_t v2vsMassLin(Double_t *x, Double_t *par){
  // Fit function for signal+background
  // par[0] = S/B at the mass peak
  // par[1] = mass
  // par[2] = sigma
  // par[3] = v2sig
  // par[4] = v2back constant
  // par[5] = v2back slope

  Double_t fracsig=par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
  Double_t fracbkg=1-fracsig;
  return par[3]*fracsig+(par[4]+par[5]*x[0])*fracbkg;
}

Bool_t FitMassSpectrum(TH1F* hMass, TPad* pad){

  Int_t nMassBins=hMass->GetNbinsX();
  hMinBin=3;
  hMaxBin=nMassBins-2;
  Double_t hmin=hMass->GetBinLowEdge(hMinBin);
  Double_t hmax=hMass->GetBinLowEdge(hMaxBin)+hMass->GetBinWidth(hMaxBin);
  Float_t massD=TDatabasePDG::Instance()->GetParticle(mesonPDG)->Mass();
    
  AliHFMassFitter* fitter=new AliHFMassFitter(hMass,hmin,hmax,2,0,0);
  fitter->SetReflectionSigmaFactor(factor4refl);
  fitter->SetInitialGaussianMean(massD);
  Bool_t out=fitter->MassFitter(0);
  if(!out) return kFALSE;
  fitter->Signal(sigmaRangeForSig, signalFromFit,esignalFromFit);
  massFromFit=fitter->GetMean();
  sigmaFromFit=fitter->GetSigma();
  fB1=fitter->GetBackgroundFullRangeFunc();
  fB2=fitter->GetBackgroundRecalcFunc();
  fSB=fitter->GetMassFunc();
  if(!fB1) return kFALSE;
  if(!fB2) return kFALSE;
  if(!fSB) return kFALSE;
  if(pad){
    fitter->DrawHere(gPad,3.,0.);
  }
  return kTRUE;
}

TH1F* DoSideBands(Int_t iFinalPtBin,
		  TH2F* hMassDphi, 
		  TH1F* hMass, 
		  Int_t rebin,
		  TCanvas* c1,
		  Int_t optBkg){

  // Build histo with cos(2*deltaphi) distribution for signal

  Double_t hmin=hMass->GetBinLowEdge(hMinBin);
  Double_t hmax=hMass->GetBinLowEdge(hMaxBin)+hMass->GetBinWidth(hMaxBin);

  Double_t minMassSig=massFromFit-sigmaRangeForSig*sigmaFromFit;
  Double_t maxMassSig=massFromFit+sigmaRangeForSig*sigmaFromFit;
  Int_t minBinSig=hMass->FindBin(minMassSig);
  Int_t maxBinSig=hMass->FindBin(maxMassSig);
  Double_t minMassSigBin=hMass->GetBinLowEdge(minBinSig);
  Double_t maxMassSigBin=hMass->GetBinLowEdge(maxBinSig)+hMass->GetBinWidth(maxBinSig);
  printf("Signal Fit Limits = %f %f\n",minMassSigBin,maxMassSigBin);
  Double_t maxMassBkgLow=massFromFit-sigmaRangeForBkg*sigmaFromFit;
  Int_t minBinBkgLow=hMinBin;
  Int_t maxBinBkgLow=hMass->FindBin(maxMassBkgLow);
  Double_t minMassBkgLowBin=hmin;
  Double_t maxMassBkgLowBin=hMass->GetBinLowEdge(maxBinBkgLow)+hMass->GetBinWidth(maxBinBkgLow);
  Double_t minMassBkgHi=massFromFit+sigmaRangeForBkg*sigmaFromFit;
  Int_t minBinBkgHi=hMass->FindBin(minMassBkgHi);
  Int_t maxBinBkgHi=hMaxBin;
  Double_t minMassBkgHiBin=hMass->GetBinLowEdge(minBinBkgHi);
  Double_t maxMassBkgHiBin=hmax;
  printf("BKG Fit Limits = %f %f  && %f %f\n",minMassBkgLowBin,maxMassBkgLowBin,minMassBkgHiBin,maxMassBkgHiBin);
  Double_t bkgSig=fB2->Integral(minMassSigBin,maxMassSigBin);
  Double_t bkgLow=fB2->Integral(minMassBkgLowBin,maxMassBkgLowBin);
  Double_t bkgHi=fB2->Integral(minMassBkgHiBin,maxMassBkgHiBin);
  printf("Background integrals = %f %f %f\n",bkgLow,bkgSig,bkgHi);
  if(c1){
    TBox* bleft=new TBox(minMassBkgLowBin,0.,maxMassBkgLowBin,hMass->GetMaximum());
    bleft->SetFillColor(kRed+1);
    bleft->SetFillStyle(3002);
    bleft->Draw();
    TBox* bright=new TBox(minMassBkgHiBin,0.,maxMassBkgHiBin,hMass->GetMaximum());
    bright->SetFillColor(kBlue+1);
    bright->SetFillStyle(3002);
    bright->Draw();
    TBox* bsig=new TBox(minMassSigBin,0.,maxMassSigBin,hMass->GetMaximum()*2);
    bsig->SetFillColor(1);
    bsig->SetFillStyle(3002);
    bsig->Draw();
  }
  TH1F* hCos2PhiBkgLo=(TH1F*)hMassDphi->ProjectionX(Form("hCos2PhiBkgLoBin%d",iFinalPtBin),minBinBkgLow,maxBinBkgLow);
  TH1F* hCos2PhiBkgHi=(TH1F*)hMassDphi->ProjectionX(Form("hCos2PhiBkgHiBin%d",iFinalPtBin),minBinBkgHi,maxBinBkgHi);
  TH1F* hCos2PhiSigReg=(TH1F*)hMassDphi->ProjectionX(Form("hCos2PhiBkgSigBin%d",iFinalPtBin),minBinSig,maxBinSig);
  
  hCos2PhiBkgLo->Rebin(rebin);
  hCos2PhiBkgHi->Rebin(rebin);
  hCos2PhiSigReg->Rebin(rebin);
  hCos2PhiSigReg->SetLineWidth(2);
  hCos2PhiBkgLo->SetLineWidth(2);
  hCos2PhiBkgHi->SetLineWidth(2);
  hCos2PhiBkgLo->SetLineColor(kRed+1);
  hCos2PhiBkgHi->SetLineColor(kBlue+1);
  
  TH1F* hCos2PhiBkgLoScal=(TH1F*)hCos2PhiBkgLo->Clone(Form("hCos2PhiBkgLoScalBin%d",iFinalPtBin));
  hCos2PhiBkgLoScal->Scale(bkgSig/bkgLow);
  TH1F* hCos2PhiBkgHiScal=(TH1F*)hCos2PhiBkgHi->Clone(Form("hCos2PhiBkgHiScalBin%d",iFinalPtBin));
  hCos2PhiBkgHiScal->Scale(bkgSig/bkgHi);
  hCos2PhiBkgLoScal->SetLineWidth(2);
  hCos2PhiBkgHiScal->SetLineWidth(2);
  hCos2PhiBkgLoScal->SetLineColor(kRed+1);
  hCos2PhiBkgHiScal->SetLineColor(kBlue+1);
  TH1F* hCos2PhiBkgAver=0x0;
  if(optBkg==0){
    hCos2PhiBkgAver=(TH1F*)hCos2PhiBkgLoScal->Clone(Form("hCos2PhiBkgAverBin%d",iFinalPtBin));
    hCos2PhiBkgAver->Add(hCos2PhiBkgHiScal);
    hCos2PhiBkgAver->Scale(0.5);
  }else if(optBkg==-1){
    hCos2PhiBkgAver=(TH1F*)hCos2PhiBkgLoScal->Clone(Form("hCos2PhiBkgAverBin%d",iFinalPtBin));
  }else{
    hCos2PhiBkgAver=(TH1F*)hCos2PhiBkgHiScal->Clone(Form("hCos2PhiBkgAverBin%d",iFinalPtBin));
  }
  hCos2PhiBkgAver->SetLineWidth(2);
  hCos2PhiBkgAver->SetLineColor(kGreen+1);
  TH1F* hCos2PhiSig=(TH1F*)hCos2PhiSigReg->Clone(Form("hCos2PhiSigBin%d",iFinalPtBin));
  hCos2PhiSig->Add(hCos2PhiBkgAver,-1.);   
  
  TLegend* leg0=new TLegend(0.3,0.6,0.75,0.89);
  TPaveText* t0= new TPaveText(0.15,0.70,0.45,0.89,"NDC");
  if(c1){
    c1->cd(3);
    hCos2PhiSigReg->Draw();
    hCos2PhiBkgLoScal->Draw("same");
    hCos2PhiBkgHiScal->Draw("same");
    hCos2PhiBkgAver->Draw("same");
    leg0->SetFillColor(0);
    TLegendEntry* ent=leg0->AddEntry(hCos2PhiSigReg,"Signal region","L");
    ent->SetTextColor(hCos2PhiSigReg->GetLineColor());
    ent=leg0->AddEntry(hCos2PhiBkgLoScal,"Left side band","L");
    ent->SetTextColor(hCos2PhiBkgLoScal->GetLineColor());
    ent=leg0->AddEntry(hCos2PhiBkgHiScal,"Right side band","L");
    ent->SetTextColor(hCos2PhiBkgHiScal->GetLineColor());
    ent=leg0->AddEntry(hCos2PhiBkgAver,"Average of side bands","L");
    ent->SetTextColor(hCos2PhiBkgAver->GetLineColor());
    leg0->Draw();
    c1->cd(4);
    hCos2PhiSig->Draw("EP");
    t0->SetFillColor(0);
    t0->AddText(Form("v2=%.3f+-%.3f\n",hCos2PhiSig->GetMean(),hCos2PhiSig->GetMeanError()));
    t0->Draw();
  
    c1->Update();
  }
  if(!c1){
    delete leg0;
    delete t0;
    delete hCos2PhiBkgLo;
    delete hCos2PhiBkgHi;
    delete hCos2PhiSigReg;
    delete hCos2PhiBkgLoScal;
    delete hCos2PhiBkgHiScal;
    delete hCos2PhiBkgAver;
  }
  printf("Signal from mass fitter = %f  Signal from subracted histo= %f\n",
	 signalFromFit,hCos2PhiSig->Integral());
  return hCos2PhiSig;
}


TF1* DoFitv2vsMass(Int_t iFinalPtBin, TH2F* hMassDphi, TH1F* hMass, Int_t rebin, TCanvas* c2, Int_t optErrors, Bool_t useConst){

  Int_t npars=fSB->GetNpar();
  Double_t sigmaSB=fSB->GetParameter(npars-1);
  Double_t massSB=fSB->GetParameter(npars-2);
  Double_t integr=fSB->GetParameter(npars-3);
  Double_t sOverAll=(fSB->Eval(massSB)-fB2->Eval(massSB))/fSB->Eval(massSB);
  printf("mass=%f  S+B=%f   bkg=%f S/(S+B)=%f\n",massSB,fSB->Eval(massSB),fB2->Eval(massSB),sOverAll);
  printf("Number of parameters: %d. Signal params: %f %f %f\n",npars,massSB,sigmaSB,integr);
  if(optErrors==1){
    sOverAll+=fSB->GetParError(npars-3)/fSB->GetParameter(npars-3);
  }
  else if(optErrors==-1){
    sOverAll-=fSB->GetParError(npars-3)/fSB->GetParameter(npars-3);
  }

  Int_t nbinsmass=hMassDphi->GetNbinsY();
  Double_t minmass=hMassDphi->GetYaxis()->GetXmin();
  Double_t maxmass=hMassDphi->GetYaxis()->GetXmax();

  TF1* fSig=new TF1(Form("fSig%d",iFinalPtBin),"[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./[2]/[2])",minmass,maxmass);
  fSig->SetParameters(integr,massSB,sigmaSB);
 
  TH1F* hAverCos2Phi=new TH1F(Form("hAverCos2PhiBin%d",iFinalPtBin),"",nbinsmass,minmass,maxmass);
  TH1F* hFractionSig=new TH1F(Form("hFractionSigBin%d",iFinalPtBin),"",nbinsmass,minmass,maxmass);
  TH1F* hFractionBkg=new TH1F(Form("hFractionBkgBin%d",iFinalPtBin),"",nbinsmass,minmass,maxmass);

  for(Int_t iBin=1; iBin<=hMassDphi->GetNbinsY(); iBin++){
    TH1F* htemp=(TH1F*)hMassDphi->ProjectionX("htemp",iBin,iBin);
    hAverCos2Phi->SetBinContent(iBin,htemp->GetMean());
    hAverCos2Phi->SetBinError(iBin,htemp->GetMeanError());
    Double_t sig=fSig->Eval(hFractionSig->GetBinCenter(iBin));
    Double_t bkg=fB2->Eval(hFractionSig->GetBinCenter(iBin));
    if(bkg<1 && sig<1){
      hFractionSig->SetBinContent(iBin,0.);
      hFractionSig->SetBinError(iBin,0.);
      hFractionBkg->SetBinContent(iBin,1.);
      hFractionBkg->SetBinError(iBin,0.);
    }else{
      Double_t fracs=sig/(sig+bkg);
      Double_t fracb=bkg/(sig+bkg);
      Double_t efracs=0.;//TMath::Sqrt(fracs*(1.-fracs)/(sig+bkg));
      Double_t efracb=0.;//TMath::Sqrt(fracb*(1.-fracb)/(sig+bkg));
      
      hFractionSig->SetBinContent(iBin,fracs);
      hFractionSig->SetBinError(iBin,efracs);
      hFractionBkg->SetBinContent(iBin,fracb);      
      hFractionBkg->SetBinError(iBin,efracb);
    }
    delete htemp;
  }
  
  TF1* fv2=0x0;
  if(useConst){
    fv2=new TF1(Form("fv2Bin%d",iFinalPtBin),v2vsMass,minmass,maxmass,5);
  }else{
    fv2=new TF1(Form("fv2Bin%d",iFinalPtBin),v2vsMassLin,minmass,maxmass,6);
    fv2->SetParameter(5,0.);
  }
  fv2->SetParameter(0,sOverAll);
  fv2->SetParameter(1,massSB);
  fv2->SetParameter(2,sigmaSB);
  fv2->SetParameter(3,0.2);
  fv2->SetParameter(4,0.2);
  fv2->FixParameter(0,sOverAll);
  fv2->FixParameter(1,massSB);
  fv2->FixParameter(2,sigmaSB);

  if((hAverCos2Phi->GetNbinsX()%rebin)==0){
    hAverCos2Phi->Rebin(rebin);
    hAverCos2Phi->Scale(1./(Double_t)rebin);
  }

  if(c2){
    c2->Divide(2,2);
    c2->cd(1);
    hMassDphi->Draw("colz");
    c2->cd(2);
    hMass->Rebin(2);
    hMass->SetMinimum(0.);
    hMass->SetMarkerStyle(20);
    hMass->Draw("E");
    fSB->Draw("same");
    fSig->Draw("same");
    fB2->Draw("same");
    c2->cd(3);
    hFractionSig->SetMaximum(1.2);
    hFractionSig->Draw();
    hFractionSig->GetXaxis()->SetTitle("Mass (GeV/c^2)");
    hFractionSig->GetYaxis()->SetTitle("Fraction");
    hFractionBkg->SetLineColor(2);
    hFractionBkg->Draw("same");
    TLegend* leg1=new TLegend(0.15,0.15,0.35,0.35);
    leg1->SetFillColor(0);
    TLegendEntry* ent=leg1->AddEntry(hFractionSig,"S/(S+B)","L");
    ent->SetTextColor(hFractionSig->GetLineColor());
    ent=leg1->AddEntry(hFractionBkg,"B/(S+B)","L");
    ent->SetTextColor(hFractionBkg->GetLineColor());
    leg1->Draw();
    c2->cd(4);
    hAverCos2Phi->Fit(fv2);
    fv2->DrawCopy("same");
    hAverCos2Phi->GetXaxis()->SetTitle("Mass (GeV/c^2)");
    hAverCos2Phi->GetYaxis()->SetTitle("v_2^{obs}");
    TPaveText* t1= new TPaveText(0.55,0.70,0.89,0.89,"NDC");
    t1->SetFillColor(0);
    t1->AddText(Form("v2sig=%.3f+-%.3f\n",fv2->GetParameter(3),fv2->GetParError(3)));
    if(useConst){
      t1->AddText(Form("v2bkg=%.3f+-%.3f\n",fv2->GetParameter(4),fv2->GetParError(4)));
    }else{
      t1->AddText(Form("v2bkg=(%.3f+-%.3f) + (%.3g+-%.3g)*mass\n",fv2->GetParameter(4),fv2->GetParError(4),fv2->GetParameter(5),fv2->GetParError(5)));
    }
    t1->Draw();
    c2->Update();
  }else{
    hAverCos2Phi->Fit(fv2);
  }
  return fv2;
}


Double_t GetEPResolution(TList* lst, Double_t &rcflow, Double_t &rcfhigh){
  // Event plane resolution syst err (from wide centrality bin
  TH1F* hResolSubAB=0x0;
  Double_t xmin=1.;
  Double_t xmax=-1.;
  TGraphAsymmErrors* grSingle=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* grInteg=new TGraphAsymmErrors(0);
  Int_t iPt=0;
  for(Int_t iHisC=minCent; iHisC<=maxCent-5; iHisC+=5){    
    TString hisnameEP=Form("hEvPlaneResocentr%d_%d",iHisC,iHisC+5);
    TH1F* hResolSubABsing=(TH1F*)lst->FindObject(hisnameEP.Data());
    Double_t resolFull=AliEventPlaneResolution::GetFullEvResol(hResolSubABsing,1);
    Double_t resolFullmin=AliEventPlaneResolution::GetFullEvResolLowLim(hResolSubABsing,1);
    Double_t resolFullmax=AliEventPlaneResolution::GetFullEvResolHighLim(hResolSubABsing,1);
    grSingle->SetPoint(iPt,iHisC+2.5,resolFull);
    grSingle->SetPointEXlow(iPt,2.5);
    grSingle->SetPointEXhigh(iPt,2.5);
    grSingle->SetPointEYlow(iPt,resolFullmax-resolFull);
    grSingle->SetPointEYhigh(iPt,resolFull-resolFullmin);
    ++iPt;
    if(resolFullmin<xmin) xmin=resolFullmin;
    if(resolFullmax>xmax) xmax=resolFullmax;
    if(iHisC==minCent){
      hResolSubAB=(TH1F*)hResolSubABsing->Clone("hResolSubAB");
    }else{
      hResolSubAB->Add(hResolSubABsing);
    }
    printf("Adding histogram %s\n",hisnameEP.Data());
  }
  rcflow=xmin;
  rcfhigh=xmax;

  Double_t resolSub=AliEventPlaneResolution::GetSubEvResol(hResolSubAB);
  printf("Resolution on sub event  = %.4f\n",resolSub);
  Double_t chisub=AliEventPlaneResolution::FindChi(resolSub,1);
  printf("Chi Subevent             = %.4f\n",chisub);
  Double_t chifull=chisub*TMath::Sqrt(2);
  printf("Chi Full Event           = %.4f\n",chifull);
  Double_t resolFull=AliEventPlaneResolution::GetFullEvResol(hResolSubAB,1);
  Double_t resolFullmin=AliEventPlaneResolution::GetFullEvResolLowLim(hResolSubAB,1);
  Double_t resolFullmax=AliEventPlaneResolution::GetFullEvResolHighLim(hResolSubAB,1);

  AliEventPlaneResolution *resol=new AliEventPlaneResolution(1);
  resol->SetSubEventHisto(hResolSubAB);  
  Double_t resolFull2=resol->GetFullEvResol();
  printf("Resolution on full event = %.4f %.4f\n",resolFull,resolFull2);

  grInteg->SetPoint(0,40.,resolFull);
  grInteg->SetPointEXlow(0,10);
  grInteg->SetPointEXhigh(0,10.);
  grInteg->SetPointEYlow(0,resolFullmax-resolFull);
  grInteg->SetPointEYhigh(0,resolFull-resolFullmin);
  
 
  TCanvas* cEP=new TCanvas("cEP","EvPlaneResol");
  cEP->Divide(1,2);
  cEP->cd(1);
  hResolSubAB->Draw();
  TLatex* tres=new TLatex(0.15,0.7,Form("Resolution on full event = %.4f\n",resolFull));
  tres->SetNDC();
  tres->Draw();
  cEP->cd(2);
  grSingle->SetMarkerStyle(20);
  grInteg->SetMarkerColor(kRed+1);
  grInteg->SetLineColor(kRed+1);
  grInteg->SetMarkerStyle(29);
  grSingle->Draw("AP");
  grSingle->GetXaxis()->SetTitle("Centrality Percentile");
  grSingle->GetYaxis()->SetTitle("Resolution Correction Factor");
  grInteg->Draw("Psame");
 
  return resolFull;
  
}
  

void LoadMassHistos(TList* lst, TH2F** hMassDphi){

  for(Int_t iHisC=minCent; iHisC<=maxCent-5; iHisC+=5){    
    for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
      for(Int_t iPtBin=minPtBin[iFinalPtBin]; iPtBin<=maxPtBin[iFinalPtBin]; iPtBin++){
    	TString hisname=Form("hMc2phi_pt%dcentr%d_%d",iPtBin,iHisC,iHisC+5);
    	TH2F* htmp=(TH2F*)lst->FindObject(hisname.Data());
	if(hMassDphi[iFinalPtBin]==0x0){
	  hMassDphi[iFinalPtBin]=(TH2F*)htmp->Clone(Form("hMassCos2DphiBin%d",iFinalPtBin));
	}else{ 
	  hMassDphi[iFinalPtBin]->Add(htmp);
	}
    	printf("Adding histogram %s to PtBin %d\n",hisname.Data(),iFinalPtBin);
      }
    }
  }
}

Bool_t DefinePtBins(TDirectoryFile* df){
  AliRDHFCutsD0toKpi *d0cuts=(AliRDHFCutsD0toKpi*)df->Get("D0toKpiCuts");
  Int_t nPtBinsCuts=d0cuts->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBinsCuts);
  Float_t *ptlimsCuts=d0cuts->GetPtBinLimits();
  for(Int_t iPt=0; iPt<nPtBinsCuts; iPt++) printf(" %d %f \n",iPt,ptlimsCuts[iPt]);
  for(Int_t iPtCuts=0; iPtCuts<nPtBinsCuts; iPtCuts++){
    for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){  
      if(TMath::Abs(ptlimsCuts[iPtCuts]-ptlims[iFinalPtBin])<0.0001){ 
	minPtBin[iFinalPtBin]=iPtCuts;
	if(iFinalPtBin>0) maxPtBin[iFinalPtBin-1]=iPtCuts;
      }
    }
    if(TMath::Abs(ptlimsCuts[iPtCuts]-ptlims[nFinalPtBins])<0.0001) maxPtBin[nFinalPtBins-1]=iPtCuts;
 }

  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("Pt bins to be merged: %d %d\n",minPtBin[iFinalPtBin],maxPtBin[iFinalPtBin]);
    if(minPtBin[iFinalPtBin]<0 || maxPtBin[iFinalPtBin]<0) return kFALSE;
  }
  return kTRUE;
}


void SystForSideBands(){
  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");

  TFile* filcuts=new TFile(filcutsname.Data());
  TDirectoryFile* dfcuts=(TDirectoryFile*)filcuts->Get("PWG3_D2H_HFv2");
  Bool_t binOK=DefinePtBins(dfcuts);
  if(!binOK) return;

  TFile* fil=new TFile(filname.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("PWG3_D2H_HFv2");
  TList* lst=(TList*)df->Get(listname.Data());
  if(!lst){
    printf("ERROR: list %s not found in file\n",listname.Data());
    return;
  }

  Double_t rcfmin,rcfmax;
  Double_t resolFull=GetEPResolution(lst,rcfmin,rcfmax);
  
  TH2F** hMassDphi=new TH2F*[nFinalPtBins];
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    hMassDphi[iFinalPtBin]=0x0;
  }
  LoadMassHistos(lst, hMassDphi);

  Int_t nSteps=21;

  TGraphErrors** gSystSigRange=new TGraphErrors*[nFinalPtBins];
  TGraphErrors** gSystBkgRange=new TGraphErrors*[nFinalPtBins];
  TGraphErrors** gSystRebin=new TGraphErrors*[nFinalPtBins];
  TGraphErrors** gSystWhichSide=new TGraphErrors*[nFinalPtBins];

  Double_t min1[nFinalPtBins],max1[nFinalPtBins];
  Double_t min123[nFinalPtBins],max123[nFinalPtBins];
  Double_t min2[nFinalPtBins],max2[nFinalPtBins];
  Double_t min3[nFinalPtBins],max3[nFinalPtBins];
  Double_t min4[nFinalPtBins],max4[nFinalPtBins];

  Double_t sigmaRangeForSigOrig=sigmaRangeForSig;
  Double_t sigmaRangeForBkgOrig=sigmaRangeForBkg;


  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("**************** BIN %d ******************\n",iFinalPtBin);

    Int_t rebinHistoSideBandsOrig=rebinHistoSideBands[iFinalPtBin];

    gSystSigRange[iFinalPtBin]=new TGraphErrors(0);
    gSystSigRange[iFinalPtBin]->SetTitle(Form("v2 vs. nSigma Signal Region Ptbin %d",iFinalPtBin));

    TH1F* hMass=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionY();

    Bool_t out=FitMassSpectrum(hMass,0x0);
    if(!out) continue;

    min1[iFinalPtBin]=99999.;
    max1[iFinalPtBin]=-99999.;
    min123[iFinalPtBin]=99999.;
    max123[iFinalPtBin]=-99999.;
    for(Int_t iStep=0; iStep<nSteps; iStep++){
      Int_t index=iFinalPtBin*nSteps+iStep;
      sigmaRangeForSig=1.5+0.1*iStep;
      sigmaRangeForBkg=sigmaRangeForBkgOrig;
      rebinHistoSideBands[iFinalPtBin]=rebinHistoSideBandsOrig;
      TH1F* hCos2PhiSig=DoSideBands(index,hMassDphi[iFinalPtBin],hMass,rebinHistoSideBands[iFinalPtBin],0x0);
      Double_t v2obsM1=hCos2PhiSig->GetMean();
      Double_t errv2obsM1=hCos2PhiSig->GetMeanError();  
      delete hCos2PhiSig;
      Double_t v2M1=v2obsM1/resolFull;
      Double_t errv2M1=errv2obsM1/resolFull;
      gSystSigRange[iFinalPtBin]->SetPoint(iStep,sigmaRangeForSig,v2M1);
      gSystSigRange[iFinalPtBin]->SetPointError(iStep,0.,errv2M1);
      if(v2M1>max1[iFinalPtBin]) max1[iFinalPtBin]=v2M1;
      if(v2M1<min1[iFinalPtBin]) min1[iFinalPtBin]=v2M1;
      if(sigmaRangeForSig>=2. && sigmaRangeForSig<=3){
	if(v2M1>max123[iFinalPtBin]) max123[iFinalPtBin]=v2M1;
	if(v2M1<min123[iFinalPtBin]) min123[iFinalPtBin]=v2M1;
      }
    }
 
    min2[iFinalPtBin]=99999.;
    max2[iFinalPtBin]=-99999.;
    gSystBkgRange[iFinalPtBin]=new TGraphErrors(0);
    gSystBkgRange[iFinalPtBin]->SetTitle(Form("v2 vs. nSigma Bkg Region Ptbin %d",iFinalPtBin));
    for(Int_t iStep=0; iStep<nSteps; iStep++){
      Int_t index=nSteps*nFinalPtBins+iFinalPtBin*nSteps+iStep;
      sigmaRangeForSig=sigmaRangeForSigOrig;
      sigmaRangeForBkg=4.+0.1*iStep;
      rebinHistoSideBands[iFinalPtBin]=rebinHistoSideBandsOrig;
      TH1F* hCos2PhiSig=DoSideBands(index,hMassDphi[iFinalPtBin],hMass,rebinHistoSideBands[iFinalPtBin],0x0);
      Double_t v2obsM1=hCos2PhiSig->GetMean();
      Double_t errv2obsM1=hCos2PhiSig->GetMeanError();  
      delete hCos2PhiSig;
      Double_t v2M1=v2obsM1/resolFull;
      Double_t errv2M1=errv2obsM1/resolFull;
      gSystBkgRange[iFinalPtBin]->SetPoint(iStep,sigmaRangeForBkg,v2M1);
      gSystBkgRange[iFinalPtBin]->SetPointError(iStep,0.,errv2M1);
      if(v2M1>max2[iFinalPtBin]) max2[iFinalPtBin]=v2M1;
      if(v2M1<min2[iFinalPtBin]) min2[iFinalPtBin]=v2M1;
    }

    min3[iFinalPtBin]=99999.;
    max3[iFinalPtBin]=-99999.;
    gSystRebin[iFinalPtBin]=new TGraphErrors(0);
    gSystRebin[iFinalPtBin]->SetTitle(Form("v2 vs. Rebin Ptbin %d",iFinalPtBin));
    Int_t nPts=0;
    for(Int_t iStep=0; iStep<nSteps; iStep++){
      Int_t index=2*nSteps*nFinalPtBins+iFinalPtBin*nSteps+iStep;
      sigmaRangeForSig=sigmaRangeForSigOrig;
      sigmaRangeForBkg=sigmaRangeForBkgOrig;
      rebinHistoSideBands[iFinalPtBin]=iStep+1;
      if((hMassDphi[iFinalPtBin]->GetNbinsY())%rebinHistoSideBands[iFinalPtBin]!=0) continue;
      TH1F* hCos2PhiSig=DoSideBands(index,hMassDphi[iFinalPtBin],hMass,rebinHistoSideBands[iFinalPtBin],0x0);
      Double_t v2obsM1=hCos2PhiSig->GetMean();
      Double_t errv2obsM1=hCos2PhiSig->GetMeanError();  
      delete hCos2PhiSig;
      Double_t v2M1=v2obsM1/resolFull;
      Double_t errv2M1=errv2obsM1/resolFull;
      gSystRebin[iFinalPtBin]->SetPoint(nPts,(Double_t)rebinHistoSideBands[iFinalPtBin],v2M1);
      gSystRebin[iFinalPtBin]->SetPointError(nPts,0.,errv2M1);
      if(v2M1>max3[iFinalPtBin]) max3[iFinalPtBin]=v2M1;
      if(v2M1<min3[iFinalPtBin]) min3[iFinalPtBin]=v2M1;
      ++nPts;
    }

    min4[iFinalPtBin]=99999.;
    max4[iFinalPtBin]=-99999.;
    gSystWhichSide[iFinalPtBin]=new TGraphErrors(0);
    gSystWhichSide[iFinalPtBin]->SetTitle(Form("v2 vs. WhichSide Ptbin %d",iFinalPtBin));
    nPts=0;
    for(Int_t iStep=-1; iStep<=1; iStep++){
      Int_t index=3*nSteps*nFinalPtBins+2+iStep;
      sigmaRangeForSig=sigmaRangeForSigOrig;
      sigmaRangeForBkg=sigmaRangeForBkgOrig;
      rebinHistoSideBands[iFinalPtBin]=rebinHistoSideBandsOrig;
      TH1F* hCos2PhiSig=DoSideBands(index,hMassDphi[iFinalPtBin],hMass,rebinHistoSideBands[iFinalPtBin],0x0,iStep);
      Double_t v2obsM1=hCos2PhiSig->GetMean();
      Double_t errv2obsM1=hCos2PhiSig->GetMeanError();  
      delete hCos2PhiSig;
      Double_t v2M1=v2obsM1/resolFull;
      Double_t errv2M1=errv2obsM1/resolFull;
      gSystWhichSide[iFinalPtBin]->SetPoint(nPts,(Double_t)iStep,v2M1);
      gSystWhichSide[iFinalPtBin]->SetPointError(nPts,0.,errv2M1);
      if(v2M1>max4[iFinalPtBin]) max4[iFinalPtBin]=v2M1;
      if(v2M1<min4[iFinalPtBin]) min4[iFinalPtBin]=v2M1;
      ++nPts;
    }
    rebinHistoSideBands[iFinalPtBin]=rebinHistoSideBandsOrig;
  }

  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("------ Pt Bin %d ------\n",iFinalPtBin);
    printf("Range of values for sig variation = %f %f\n",min1[iFinalPtBin],max1[iFinalPtBin]);
    printf("           (limited to 2-3 sigma) = %f %f\n",min123[iFinalPtBin],max123[iFinalPtBin]);
    printf("Range of values for bkg variation = %f %f\n",min2[iFinalPtBin],max2[iFinalPtBin]);
    printf("Range of values for rebin = %f %f\n",min3[iFinalPtBin],max3[iFinalPtBin]);
    printf("Range of values for whichside = %f %f\n",min4[iFinalPtBin],max4[iFinalPtBin]);
    Float_t minenv=TMath::Min(min123[iFinalPtBin],TMath::Min(min2[iFinalPtBin],min4[iFinalPtBin]));
    Float_t maxenv=TMath::Max(max123[iFinalPtBin],TMath::Max(max2[iFinalPtBin],max4[iFinalPtBin]));
    printf(" ---> Envelope                = %f %f\n",minenv,maxenv);
  }

  TCanvas* cs1=new TCanvas("cs1");
  cs1->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs1->cd(iFinalPtBin+1);
    gSystSigRange[iFinalPtBin]->SetMarkerStyle(20);
    gSystSigRange[iFinalPtBin]->Draw("AP");
    gSystSigRange[iFinalPtBin]->GetXaxis()->SetTitle("nSigmaSignal");
    gSystSigRange[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }

  TCanvas* cs2=new TCanvas("cs2");
  cs2->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs2->cd(iFinalPtBin+1);
    gSystBkgRange[iFinalPtBin]->SetMarkerStyle(20);
    gSystBkgRange[iFinalPtBin]->Draw("AP");
    gSystBkgRange[iFinalPtBin]->GetXaxis()->SetTitle("nSigmaBackground");
    gSystBkgRange[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }

  TCanvas* cs3=new TCanvas("cs3");
  cs3->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs3->cd(iFinalPtBin+1);
    gSystRebin[iFinalPtBin]->SetMarkerStyle(20);
    gSystRebin[iFinalPtBin]->Draw("AP");
    gSystRebin[iFinalPtBin]->GetXaxis()->SetTitle("Rebin factor");
    gSystRebin[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }

  TCanvas* cs4=new TCanvas("cs4");
  cs4->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs4->cd(iFinalPtBin+1);
    gSystWhichSide[iFinalPtBin]->SetMarkerStyle(20);
    gSystWhichSide[iFinalPtBin]->Draw("AP");
    gSystWhichSide[iFinalPtBin]->GetXaxis()->SetTitle("Side band used");
    gSystWhichSide[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }
}

void SystForFitv2Mass(){
  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");

  TFile* filcuts=new TFile(filcutsname.Data());
  TDirectoryFile* dfcuts=(TDirectoryFile*)filcuts->Get("PWG3_D2H_HFv2");
  Bool_t binOK=DefinePtBins(dfcuts);
  if(!binOK) return;

  TFile* fil=new TFile(filname.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("PWG3_D2H_HFv2");
  TList* lst=(TList*)df->Get(listname.Data());
  if(!lst){
    printf("ERROR: list %s not found in file\n",listname.Data());
    return;
  }

  Double_t rcfmin,rcfmax;
  Double_t resolFull=GetEPResolution(lst,rcfmin,rcfmax);
  
  TH2F** hMassDphi=new TH2F*[nFinalPtBins];
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    hMassDphi[iFinalPtBin]=0x0;
  }
  LoadMassHistos(lst, hMassDphi);

  Int_t nSteps=11;

  TGraphErrors** gSystParamErr=new TGraphErrors*[nFinalPtBins];
  TGraphErrors** gSystRebin=new TGraphErrors*[nFinalPtBins];
  TGraphErrors** gSystLinConst=new TGraphErrors*[nFinalPtBins];

  Double_t min1[nFinalPtBins],max1[nFinalPtBins];
  Double_t min2[nFinalPtBins],max2[nFinalPtBins];
  Double_t min3[nFinalPtBins],max3[nFinalPtBins];


  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("**************** BIN %d ******************\n",iFinalPtBin);

    Int_t rebinHistov2MassOrig=rebinHistov2Mass[iFinalPtBin];
    TH1F* hMass=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionY();

    Bool_t out=FitMassSpectrum(hMass,0x0);
    if(!out) continue;

    min1[iFinalPtBin]=99999.;
    max1[iFinalPtBin]=-99999.;
    gSystRebin[iFinalPtBin]=new TGraphErrors(0);
    gSystRebin[iFinalPtBin]->SetTitle(Form("v2 vs. Rebin PtBin %d",iFinalPtBin));
    Int_t nPts=0;
    for(Int_t iStep=0; iStep<nSteps; iStep++){
      Int_t index=iStep;
      rebinHistov2Mass[iFinalPtBin]=iStep+1;
      if((hMassDphi[iFinalPtBin]->GetNbinsY())%rebinHistov2Mass[iFinalPtBin]!=0) continue;
      TF1* fv2=DoFitv2vsMass(index,hMassDphi[iFinalPtBin],hMass,rebinHistov2Mass[iFinalPtBin],0x0,0,useConstantvsBkgVsMass);
      Double_t v2obsM2=fv2->GetParameter(3);
      Double_t errv2obsM2=fv2->GetParError(3);
      delete fv2;
      Double_t v2M2=v2obsM2/resolFull;
      Double_t errv2M2=errv2obsM2/resolFull;
      gSystRebin[iFinalPtBin]->SetPoint(nPts,(Double_t)rebinHistov2Mass[iFinalPtBin],v2M2);
      gSystRebin[iFinalPtBin]->SetPointError(nPts,0.,errv2M2);
      if(v2M2>max1[iFinalPtBin]) max1[iFinalPtBin]=v2M2;
      if(v2M2<min1[iFinalPtBin]) min1[iFinalPtBin]=v2M2;
      ++nPts;
    }
    rebinHistov2Mass[iFinalPtBin]=rebinHistov2MassOrig;

    min2[iFinalPtBin]=99999.;
    max2[iFinalPtBin]=-99999.;
    gSystParamErr[iFinalPtBin]=new TGraphErrors(0);
    gSystParamErr[iFinalPtBin]->SetTitle(Form("v2 vs. ParamErr PtBin %d",iFinalPtBin));
    nPts=0;
    for(Int_t iStep=-1; iStep<=1; iStep++){
      Int_t index=nSteps*2+iStep;
      TF1* fv2=DoFitv2vsMass(index,hMassDphi[iFinalPtBin],hMass,rebinHistov2Mass[iFinalPtBin],0x0,iStep,useConstantvsBkgVsMass);
      Double_t v2obsM2=fv2->GetParameter(3);
      Double_t errv2obsM2=fv2->GetParError(3);
      delete fv2;
      Double_t v2M2=v2obsM2/resolFull;
      Double_t errv2M2=errv2obsM2/resolFull;
      gSystParamErr[iFinalPtBin]->SetPoint(nPts,(Double_t)iStep,v2M2);
      gSystParamErr[iFinalPtBin]->SetPointError(nPts,0.,errv2M2);
      if(v2M2>max2[iFinalPtBin]) max2[iFinalPtBin]=v2M2;
      if(v2M2<min2[iFinalPtBin]) min2[iFinalPtBin]=v2M2;
      ++nPts;
    }

    min3[iFinalPtBin]=99999.;
    max3[iFinalPtBin]=-99999.;
    gSystLinConst[iFinalPtBin]=new TGraphErrors(0);
    gSystLinConst[iFinalPtBin]->SetTitle(Form("v2 LinVsConst Bin%d",iFinalPtBin));
    nPts=0;
    for(Int_t iStep=0; iStep<=1; iStep++){
      Int_t index=nSteps*3+iStep;
      TF1* fv2=DoFitv2vsMass(index,hMassDphi[iFinalPtBin],hMass,rebinHistov2Mass[iFinalPtBin],0x0,0,iStep);
      Double_t v2obsM2=fv2->GetParameter(3);
      Double_t errv2obsM2=fv2->GetParError(3);
      delete fv2;
      Double_t v2M2=v2obsM2/resolFull;
      Double_t errv2M2=errv2obsM2/resolFull;
      gSystLinConst[iFinalPtBin]->SetPoint(nPts,1.-(Double_t)iStep,v2M2);
      gSystLinConst[iFinalPtBin]->SetPointError(nPts,0.,errv2M2);
      if(v2M2>max3[iFinalPtBin]) max3[iFinalPtBin]=v2M2;
      if(v2M2<min3[iFinalPtBin]) min3[iFinalPtBin]=v2M2;
      ++nPts;
    }


  }

  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("------ Pt Bin %d ------\n",iFinalPtBin);
    printf("Range of values for rebin = %f %f\n",min1[iFinalPtBin],max1[iFinalPtBin]);
    printf("Range of values for par err = %f %f\n",min2[iFinalPtBin],max2[iFinalPtBin]);
    printf("Range of values for lin const = %f %f\n",min3[iFinalPtBin],max3[iFinalPtBin]);
    Float_t minenv=TMath::Min(min2[iFinalPtBin],min3[iFinalPtBin]);
    Float_t maxenv=TMath::Max(max2[iFinalPtBin],max3[iFinalPtBin]);
    printf(" ---> Envelope                = %f %f\n",minenv,maxenv);
  }

  TCanvas* cs1=new TCanvas("cs1");
  cs1->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs1->cd(iFinalPtBin+1);
    gSystRebin[iFinalPtBin]->SetMarkerStyle(20);
    gSystRebin[iFinalPtBin]->Draw("AP");
    gSystRebin[iFinalPtBin]->GetXaxis()->SetTitle("Rebin factor");
    gSystRebin[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }

  TCanvas* cs2=new TCanvas("cs2");
  cs2->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs2->cd(iFinalPtBin+1);
    gSystParamErr[iFinalPtBin]->SetMarkerStyle(20);
    gSystParamErr[iFinalPtBin]->Draw("AP");
    gSystParamErr[iFinalPtBin]->GetXaxis()->SetTitle("Error on Signal yield");
    gSystParamErr[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }

  TCanvas* cs3=new TCanvas("cs3");
  cs3->Divide(2,2);
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    cs3->cd(iFinalPtBin+1);
    gSystLinConst[iFinalPtBin]->SetMarkerStyle(20);
    gSystLinConst[iFinalPtBin]->Draw("AP");
    gSystLinConst[iFinalPtBin]->GetXaxis()->SetLimits(-0.1,1.1);
    gSystLinConst[iFinalPtBin]->GetXaxis()->SetTitle("Const/Linear v2 of background");
    gSystLinConst[iFinalPtBin]->GetYaxis()->SetTitle("v2");
  }
}
