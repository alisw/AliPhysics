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

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#include "AliEventPlaneResolution.h"
#endif

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

void Extractv2from2Dhistos(){

  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");

  TFile* fil=new TFile("AnalysisResults_ptbins.root");
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("PWG3_D2H_HFv2");
  TList* lst=(TList*)df->Get("coutputv2D0Std");

  AliRDHFCutsD0toKpi *d0cuts=(AliRDHFCutsD0toKpi*)df->Get("D0toKpiCuts");
  Int_t nPtBinsCuts=d0cuts->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBinsCuts);
  Float_t *ptlimsCuts=d0cuts->GetPtBinLimits();
  for(Int_t iPt=0; iPt<nPtBinsCuts; iPt++) printf(" %d %f \n",iPt,ptlimsCuts[iPt]);

  Int_t minCent=30;
  Int_t maxCent=50;


  const Int_t nFinalPtBins=3;
  Double_t ptlims[nFinalPtBins+1]={2.,4.,5.,12.};

  Int_t minPtBin[nFinalPtBins]={-1,-1,-1};
  Int_t maxPtBin[nFinalPtBins]={-1,-1,-1};
  for(Int_t iPtCuts=0; iPtCuts<nPtBinsCuts; iPtCuts++){
    for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){  
      if(TMath::Abs(ptlimsCuts[iPtCuts]-ptlims[iFinalPtBin])<0.0001){ 
	minPtBin[iFinalPtBin]=iPtCuts;
	if(iFinalPtBin>0) maxPtBin[iFinalPtBin-1]=iPtCuts;
      }
    }
    if(TMath::Abs(ptlimsCuts[iPtCuts]-ptlims[nFinalPtBins])<0.0001) maxPtBin[nFinalPtBins-1]=iPtCuts;
 }
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++) printf("Pt bins to be merged: %d %d\n",minPtBin[iFinalPtBin],maxPtBin[iFinalPtBin]);

  TH1F* hResolSubAB=0x0;
  TH2F** hMassDphi=new TH2F*[nFinalPtBins];
  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    hMassDphi[iFinalPtBin]=0x0;
  }


  for(Int_t iHisC=minCent; iHisC<=maxCent-5; iHisC+=5){    
    TString hisnameEP=Form("hEvPlaneResocentr%d_%d",iHisC,iHisC+5);
    TH2F* htmpEP=(TH2F*)lst->FindObject(hisnameEP.Data());

    if(iHisC==minCent){ 
      hResolSubAB=(TH1F*)htmpEP->Clone("hResolSubAB");
    }else{
      hResolSubAB->Add(htmpEP);
    }
    printf("Adding histogram %s\n",hisnameEP.Data());

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

  AliEventPlaneResolution *resol=new AliEventPlaneResolution(1);
  resol->SetSubEventHisto(hResolSubAB);  

  Double_t resolSub=AliEventPlaneResolution::GetSubEvResol(hResolSubAB);
  printf("Resolution on sub event  = %.4f\n",resolSub);
  Double_t chisub=AliEventPlaneResolution::FindChi(resolSub,1);
  printf("Chi Subevent             = %.4f\n",chisub);
  Double_t chifull=chisub*TMath::Sqrt(2);
  printf("Chi Full Event           = %.4f\n",chifull);
  Double_t resolFull=resol->GetFullEvResol();
  Double_t resolFull2=AliEventPlaneResolution::GetFullEvResol(hResolSubAB,1);
  printf("Resolution on full event = %.4f %.4f\n",resolFull,resolFull2);
 
  TCanvas* cEP=new TCanvas("cEP","EvPlaneResol");
  hResolSubAB->Draw();
  TLatex* tres=new TLatex(0.15,0.7,Form("Resolution on full event = %.4f\n",resolFull));
  tres->SetNDC();
  tres->Draw();
  cEP->Update();

  TCanvas** c1=new TCanvas*[nFinalPtBins];
  TCanvas** c2=new TCanvas*[nFinalPtBins];
  TF1** fSig=new TF1*[nFinalPtBins];
  TF1** fv2=new TF1*[nFinalPtBins];
  TH1F** hAverCos2Phi=new TH1F*[nFinalPtBins];
  TH1F** hFractionSig=new TH1F*[nFinalPtBins];
  TH1F** hFractionBkg=new TH1F*[nFinalPtBins];

  TLegend* leg0=0x0;
  TLegend* leg1=0x0;
  TLegendEntry* ent;

  TH1F* hMass=0x0;
  TH1F* hCos2PhiBkgLo=0x0;
  TH1F* hCos2PhiBkgHi=0x0;
  TH1F* hCos2PhiBkgLoScal=0x0;
  TH1F* hCos2PhiBkgHiScal=0x0;
  TH1F* hCos2PhiBkgAver=0x0;
  TH1F* hCos2PhiSigReg=0x0;
  TH1F* hCos2PhiSig=0x0;
 
  Double_t v2M1[nFinalPtBins],errv2M1[nFinalPtBins];
  Double_t v2M2[nFinalPtBins],errv2M2[nFinalPtBins];

  for(Int_t iFinalPtBin=0; iFinalPtBin<nFinalPtBins; iFinalPtBin++){
    printf("**************** BIN %d ******************\n",iFinalPtBin);
    printf("\n--------- Method 1: Side Bands ----------\n");
    hMass=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionY();

    Double_t sigmaRangeForSig=3.;
    Double_t sigmaRangeForBkg=4.5;

    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    c1[iFinalPtBin]=new TCanvas(Form("cMeth1PtBin%d",iFinalPtBin),Form("Meth1PtBin%d",iFinalPtBin));
    c1[iFinalPtBin]->Divide(2,2);
    c1[iFinalPtBin]->cd(1);
    hMassDphi[iFinalPtBin]->Draw("colz");
    c1[iFinalPtBin]->cd(2);
    
    Int_t nMassBins=hMass->GetNbinsX();
    Int_t hMinBin=3;
    Int_t hMaxBin=nMassBins-2;
    Double_t hmin=hMass->GetBinLowEdge(hMinBin);
    Double_t hmax=hMass->GetBinLowEdge(hMaxBin)+hMass->GetBinWidth(hMaxBin);
    Int_t factor4refl=0;
    Float_t massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    
    AliHFMassFitter* fitter=new AliHFMassFitter(hMass,hmin,hmax,2,0,0);
    fitter->SetReflectionSigmaFactor(factor4refl);
    fitter->SetInitialGaussianMean(massD);
    Bool_t out=fitter->MassFitter(0);
    if(!out) continue;
    fitter->DrawHere(gPad);
    Double_t sigfitter,esigfitter;
    fitter->Signal(sigmaRangeForSig, sigfitter,esigfitter);
    Double_t mass=fitter->GetMean();
    Double_t sigma=fitter->GetSigma();
    TF1* fB1=fitter->GetBackgroundFullRangeFunc();
    TF1* fB2=fitter->GetBackgroundRecalcFunc();
    TF1* fSB=fitter->GetMassFunc();
    Double_t minMassSig=mass-sigmaRangeForSig*sigma;
    Double_t maxMassSig=mass+sigmaRangeForSig*sigma;
    Int_t minBinSig=hMass->FindBin(minMassSig);
    Int_t maxBinSig=hMass->FindBin(maxMassSig);
    Double_t minMassSigBin=hMass->GetBinLowEdge(minBinSig);
    Double_t maxMassSigBin=hMass->GetBinLowEdge(maxBinSig)+hMass->GetBinWidth(maxBinSig);
    printf("Signal Fit Limits = %f %f\n",minMassSigBin,maxMassSigBin);
    Double_t maxMassBkgLow=mass-sigmaRangeForBkg*sigma;
    Int_t minBinBkgLow=hMinBin;
    Int_t maxBinBkgLow=hMass->FindBin(maxMassBkgLow);
    Double_t minMassBkgLowBin=hmin;
    Double_t maxMassBkgLowBin=hMass->GetBinLowEdge(maxBinBkgLow)+hMass->GetBinWidth(maxBinBkgLow);
    Double_t minMassBkgHi=mass+sigmaRangeForBkg*sigma;
    Int_t minBinBkgHi=hMass->FindBin(minMassBkgHi);
    Int_t maxBinBkgHi=hMaxBin;
    Double_t minMassBkgHiBin=hMass->GetBinLowEdge(minBinBkgHi);
    Double_t maxMassBkgHiBin=hmax;
    printf("BKG Fit Limits = %f %f  && %f %f\n",minMassBkgLowBin,maxMassBkgLowBin,minMassBkgHiBin,maxMassBkgHiBin);
    Double_t bkgSig=fB2->Integral(minMassSigBin,maxMassSigBin);
    Double_t bkgLow=fB2->Integral(minMassBkgLowBin,maxMassBkgLowBin);
    Double_t bkgHi=fB2->Integral(minMassBkgHiBin,maxMassBkgHiBin);
    printf("Background integrals = %f %f %f\n",bkgLow,bkgSig,bkgHi);
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
    
    hCos2PhiBkgLo=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionX(Form("hCos2PhiBkgLoBin%d",iFinalPtBin),minBinBkgLow,maxBinBkgLow);
    hCos2PhiBkgHi=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionX(Form("hCos2PhiBkgHiBin%d",iFinalPtBin),minBinBkgHi,maxBinBkgHi);
    hCos2PhiSigReg=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionX(Form("hCos2PhiBkgSigBin%d",iFinalPtBin),minBinSig,maxBinSig);

    hCos2PhiBkgLo->Rebin(4);
    hCos2PhiBkgHi->Rebin(4);
    hCos2PhiSigReg->Rebin(4);
    hCos2PhiSigReg->SetLineWidth(2);
    hCos2PhiBkgLo->SetLineWidth(2);
    hCos2PhiBkgHi->SetLineWidth(2);
    hCos2PhiBkgLo->SetLineColor(kRed+1);
    hCos2PhiBkgHi->SetLineColor(kBlue+1);

    hCos2PhiBkgLoScal=(TH1F*)hCos2PhiBkgLo->Clone(Form("hCos2PhiBkgLoScalBin%d",iFinalPtBin));
    hCos2PhiBkgLoScal->Scale(bkgSig/bkgLow);
    hCos2PhiBkgHiScal=(TH1F*)hCos2PhiBkgHi->Clone(Form("hCos2PhiBkgHiScalBin%d",iFinalPtBin));
    hCos2PhiBkgHiScal->Scale(bkgSig/bkgHi);
    hCos2PhiBkgLoScal->SetLineWidth(2);
    hCos2PhiBkgHiScal->SetLineWidth(2);
    hCos2PhiBkgLoScal->SetLineColor(kRed+1);
    hCos2PhiBkgHiScal->SetLineColor(kBlue+1);
    hCos2PhiBkgAver=(TH1F*)hCos2PhiBkgLoScal->Clone(Form("hCos2PhiBkgAverBin%d",iFinalPtBin));
    hCos2PhiBkgAver->Add(hCos2PhiBkgHiScal);
    hCos2PhiBkgAver->Scale(0.5);
    hCos2PhiBkgAver->SetLineWidth(2);
    hCos2PhiBkgAver->SetLineColor(kGreen+1);
    hCos2PhiSig=(TH1F*)hCos2PhiSigReg->Clone(Form("hCos2PhiSigBin%d",iFinalPtBin));
    hCos2PhiSig->Add(hCos2PhiBkgAver,-1.);   
  
    c1[iFinalPtBin]->cd(3);
    hCos2PhiSigReg->Draw();
    hCos2PhiBkgLoScal->Draw("same");
    hCos2PhiBkgHiScal->Draw("same");
    hCos2PhiBkgAver->Draw("same");
    if(iFinalPtBin==0){
      leg0=new TLegend(0.3,0.6,0.75,0.89);
      leg0->SetFillColor(0);
      ent=leg0->AddEntry(hCos2PhiSigReg,"Signal region","L");
      ent->SetTextColor(hCos2PhiSigReg->GetLineColor());
      ent=leg0->AddEntry(hCos2PhiBkgLoScal,"Left side band","L");
      ent->SetTextColor(hCos2PhiBkgLoScal->GetLineColor());
      ent=leg0->AddEntry(hCos2PhiBkgHiScal,"Right side band","L");
      ent->SetTextColor(hCos2PhiBkgHiScal->GetLineColor());
      ent=leg0->AddEntry(hCos2PhiBkgAver,"Average of side bands","L");
      ent->SetTextColor(hCos2PhiBkgAver->GetLineColor());
    }
    leg0->Draw();
    c1[iFinalPtBin]->cd(4);
    hCos2PhiSig->Draw("EP");
    TPaveText* t0= new TPaveText(0.15,0.70,0.45,0.89,"NDC");
    t0->SetFillColor(0);
    t0->AddText(Form("v2=%.3f+-%.3f\n",hCos2PhiSig->GetMean(),hCos2PhiSig->GetMeanError()));
    t0->Draw();

    printf("Signal from mass fitter = %f  Signal from subracted histo= %f\n",
	   sigfitter,hCos2PhiSig->Integral());
    Double_t v2obsM1=hCos2PhiSig->GetMean();
    Double_t errv2obsM1=hCos2PhiSig->GetMeanError();  
    printf("v2obs = %f +- %f\n",v2obsM1,errv2obsM1);
    v2M1[iFinalPtBin]=v2obsM1/resolFull;
    errv2M1[iFinalPtBin]=errv2obsM1/resolFull;
    printf("v2 = %f +- %f\n",v2M1[iFinalPtBin],errv2M1[iFinalPtBin]);
    
    printf("\n--------- Method 2: S/S+B ----------\n");
    Int_t npars=fSB->GetNpar();
    Double_t sigmaSB=fSB->GetParameter(npars-1);
    Double_t massSB=fSB->GetParameter(npars-2);
    Double_t integr=fSB->GetParameter(npars-3);
    Double_t sOverAll=(fSB->Eval(massSB)-fB2->Eval(massSB))/fSB->Eval(massSB);
    printf("mass=%f  S+B=%f   bkg=%f S/(S+B)=%f\n",massSB,fSB->Eval(massSB),fB2->Eval(massSB),sOverAll);
    printf("Number of parameters: %d. Signal params: %f %f %f\n",npars,massSB,sigmaSB,integr);

    Int_t nbinsmass=hMassDphi[iFinalPtBin]->GetNbinsY();
    Double_t minmass=hMassDphi[iFinalPtBin]->GetYaxis()->GetXmin();
    Double_t maxmass=hMassDphi[iFinalPtBin]->GetYaxis()->GetXmax();

    fSig[iFinalPtBin]=new TF1(Form("fSig%d",iFinalPtBin),"[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./[2]/[2])",minmass,maxmass);
    fSig[iFinalPtBin]->SetParameters(integr,massSB,sigmaSB);
 
    hAverCos2Phi[iFinalPtBin]=new TH1F(Form("hAverCos2PhiBin%d",iFinalPtBin),"",nbinsmass,minmass,maxmass);
    hFractionSig[iFinalPtBin]=new TH1F(Form("hFractionSigBin%d",iFinalPtBin),"",nbinsmass,minmass,maxmass);
    hFractionBkg[iFinalPtBin]=new TH1F(Form("hFractionBkgBin%d",iFinalPtBin),"",nbinsmass,minmass,maxmass);

    for(Int_t iBin=1; iBin<=hMassDphi[iFinalPtBin]->GetNbinsY(); iBin++){
      TH1F* htemp=(TH1F*)hMassDphi[iFinalPtBin]->ProjectionX("htemp",iBin,iBin);
      hAverCos2Phi[iFinalPtBin]->SetBinContent(iBin,htemp->GetMean());
      hAverCos2Phi[iFinalPtBin]->SetBinError(iBin,htemp->GetMeanError());
      Double_t sig=fSig[iFinalPtBin]->Eval(hFractionSig[iFinalPtBin]->GetBinCenter(iBin));
      Double_t bkg=fB2->Eval(hFractionSig[iFinalPtBin]->GetBinCenter(iBin));
      if(bkg<1 && sig<1){
	hFractionSig[iFinalPtBin]->SetBinContent(iBin,0.);
	hFractionSig[iFinalPtBin]->SetBinError(iBin,0.);
	hFractionBkg[iFinalPtBin]->SetBinContent(iBin,1.);
	hFractionBkg[iFinalPtBin]->SetBinError(iBin,0.);
      }else{
	Double_t fracs=sig/(sig+bkg);
	Double_t fracb=bkg/(sig+bkg);
	Double_t efracs=0.;//TMath::Sqrt(fracs*(1.-fracs)/(sig+bkg));
	Double_t efracb=0.;//TMath::Sqrt(fracb*(1.-fracb)/(sig+bkg));
       
	hFractionSig[iFinalPtBin]->SetBinContent(iBin,fracs);
	hFractionSig[iFinalPtBin]->SetBinError(iBin,efracs);
	hFractionBkg[iFinalPtBin]->SetBinContent(iBin,fracb);      
	hFractionBkg[iFinalPtBin]->SetBinError(iBin,efracb);
      }
      delete htemp;
    }
  
    fv2[iFinalPtBin]=new TF1(Form("fv2Bin%d",iFinalPtBin),v2vsMass,minmass,maxmass,5);
    fv2[iFinalPtBin]->SetParameter(0,sOverAll);
    fv2[iFinalPtBin]->SetParameter(1,massSB);
    fv2[iFinalPtBin]->SetParameter(2,sigmaSB);
    fv2[iFinalPtBin]->SetParameter(3,0.2);
    fv2[iFinalPtBin]->SetParameter(4,0.2);
    fv2[iFinalPtBin]->FixParameter(0,sOverAll);
    fv2[iFinalPtBin]->FixParameter(1,massSB);
    fv2[iFinalPtBin]->FixParameter(2,sigmaSB);

    hAverCos2Phi[iFinalPtBin]->Rebin(2);
    hAverCos2Phi[iFinalPtBin]->Scale(0.5);

    c2[iFinalPtBin]=new TCanvas(Form("cMeth2Bin%d",iFinalPtBin),Form("cMeth2Bin%d",iFinalPtBin));
    c2[iFinalPtBin]->Divide(2,2);
    c2[iFinalPtBin]->cd(1);
    hMassDphi[iFinalPtBin]->Draw("colz");
    c2[iFinalPtBin]->cd(2);
    hMass->Rebin(2);
    hMass->SetMinimum(0.);
    hMass->SetMarkerStyle(20);
    hMass->Draw("E");
    fSB->Draw("same");
    fSig[iFinalPtBin]->Draw("same");
    fB2->Draw("same");
    c2[iFinalPtBin]->cd(3);
    hFractionSig[iFinalPtBin]->SetMaximum(1.2);
    hFractionSig[iFinalPtBin]->Draw();
    hFractionSig[iFinalPtBin]->GetXaxis()->SetTitle("Mass (GeV/c^2)");
    hFractionSig[iFinalPtBin]->GetYaxis()->SetTitle("Fraction");
    hFractionBkg[iFinalPtBin]->SetLineColor(2);
    hFractionBkg[iFinalPtBin]->Draw("same");
    if(iFinalPtBin==0){
      leg1=new TLegend(0.15,0.15,0.35,0.35);
      leg1->SetFillColor(0);
      ent=leg1->AddEntry(hFractionSig[iFinalPtBin],"S/(S+B)","L");
      ent->SetTextColor(hFractionSig[iFinalPtBin]->GetLineColor());
      ent=leg1->AddEntry(hFractionBkg[iFinalPtBin],"B/(S+B)","L");
      ent->SetTextColor(hFractionBkg[iFinalPtBin]->GetLineColor());
    }
    leg1->Draw();
    c2[iFinalPtBin]->cd(4);
    hAverCos2Phi[iFinalPtBin]->Fit(fv2[iFinalPtBin]);
    hAverCos2Phi[iFinalPtBin]->GetXaxis()->SetTitle("Mass (GeV/c^2)");
    hAverCos2Phi[iFinalPtBin]->GetYaxis()->SetTitle("v_2^{obs}");
    TPaveText* t1= new TPaveText(0.55,0.70,0.89,0.89,"NDC");
    t1->SetFillColor(0);
    t1->AddText(Form("v2sig=%.3f+-%.3f\n",fv2[iFinalPtBin]->GetParameter(3),fv2[iFinalPtBin]->GetParError(3)));
    t1->AddText(Form("v2bkg=%.3f+-%.3f\n",fv2[iFinalPtBin]->GetParameter(4),fv2[iFinalPtBin]->GetParError(4)));
    t1->Draw();

    Double_t v2obsM2=fv2[iFinalPtBin]->GetParameter(3);
    Double_t errv2obsM2=fv2[iFinalPtBin]->GetParError(3);
    printf("v2obs = %f +- %f\n",v2obsM2,errv2obsM2);
    v2M2[iFinalPtBin]=v2obsM2/resolFull;
    errv2M2[iFinalPtBin]=errv2obsM2/resolFull;
    printf("v2 = %f +- %f\n",v2M2[iFinalPtBin],errv2M2[iFinalPtBin]);
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
    
   Double_t maxy=TMath::Max(hv2m2->GetMaximum(),hv2m1->GetMaximum())+0.1;
   Double_t miny=TMath::Min(hv2m2->GetMinimum(),hv2m1->GetMinimum())-0.1;
   TH2F* hempty=new TH2F("hempty","",10,0.,hv2m1->GetXaxis()->GetXmax()+2.,10,miny,maxy);
   hempty->GetXaxis()->SetTitle("p_{t} (GeV/c)");
   hempty->GetYaxis()->SetTitle("v_{2}");

   TCanvas* cv2=new TCanvas("cv2","v2");
   hempty->Draw();
   hv2m1->SetMarkerStyle(26);
   hv2m1->Draw("same");
   hv2m2->SetLineColor(2);
   hv2m2->SetMarkerColor(2);
   hv2m2->SetMarkerStyle(20);
   hv2m2->Draw("same");
   TLegend* leg2=new TLegend(0.5,0.7,0.89,0.89);
   leg2->SetFillStyle(0);
   ent=leg2->AddEntry(hv2m1,"Side Band subtraction","P");
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
   outfil->Close();
}
