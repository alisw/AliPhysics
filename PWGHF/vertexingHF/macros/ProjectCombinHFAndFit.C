#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include "AliHFMassFitter.h"
#include "AliHFMassFitterVAR.h"
#include "AliNormalizationCounter.h"
#endif

enum Method{kME,kRot,kLS,kSB};

// input files and pt binning
TString fileName="AnalysisResults_train543544.root";
TString fileNameMC="AnalysisResultsMC_train471472.root";
TString suffix="c3SigPID_Pt400_EM1";
TString meson="Dzero";
const Int_t nPtBins=8;
Double_t binLims[nPtBins+1]={0.,1.,2.,3.,4.,5.,6.,8.,12.};
Double_t sigmas[nPtBins]={0.0079,0.0084,0.0092,0.0101,
			  0.0109,0.0118,0.0129,0.0145}; // MC values up to 12 (p-Pb pt(track)>400 MeV/c)

// outputfiles
Bool_t saveCanvasAsRoot=kTRUE;
Int_t saveCanvasAsEps=2;   //0=none, 1=main ones, 2=all

// fit configuration
Int_t rebin[nPtBins]={5,5,5,5,5,5,5,5};
Bool_t fixSigma=kTRUE;
Bool_t fixMean=kFALSE;
Double_t minMass=1.74;
Double_t maxMass=2.04;
Int_t optForNorm=1;
Double_t rangeForNorm=0.05;
TString fitoption="E";
Bool_t useEMwithLS=kTRUE;
Int_t typeb=2;
Double_t nsigmaBinCounting=4.;      // defines the mass interval over which the signal is bin counted
Double_t massD;

Int_t smoothLS=0;

// reflection option
TString reflopt="2gaus";
Bool_t correctForRefl=kFALSE;
Double_t rOverSmodif=1;

// objects and options related to side-band fit method
Bool_t tryDirectFit=kTRUE;

Int_t nparback=0;
Double_t fitrangelow[nPtBins]={1.74,1.74,1.74,1.72,1.72,1.72,1.72,1.72};
Double_t fitrangeup[nPtBins]={2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};
Int_t nDegreeBackPol[nPtBins]={4,4,4,2,2,2,2,2};             // degree of polynomial function describing the background

TH1D* hMCReflPtBin;
TH1D* hMCSigPtBin;

void WriteFitInfo(AliHFMassFitter *fitter, TH1D* histo);
TH1* AdaptRangesAndBinning(const TH1 *hRef,TH1 *h1);

AliHFMassFitter* ConfigureFitter(TH1D* histo, Int_t iPtBin){
  TH1F* histof=(TH1F*)histo->Clone(Form("%s_Fl",histo->GetName()));
  AliHFMassFitter* fitter=new AliHFMassFitter(histof,minMass,maxMass,1,typeb,0);
  fitter->SetReflectionSigmaFactor(0);
  fitter->SetUseChi2Fit();
  fitter->SetInitialGaussianMean(massD);
  fitter->SetInitialGaussianSigma(sigmas[iPtBin]);
  if(fixSigma) fitter->SetFixGaussianSigma(sigmas[iPtBin]);
  if(fixMean) fitter->SetFixGaussianMean(massD);
  return fitter;
}

AliHFMassFitterVAR* ConfigureFitterVAR(TH1D* histo, Int_t iPtBin,Int_t typeback=6,Bool_t saveReflModifPlot=kTRUE){
  
  TString hname=histo->GetName();
  Double_t xfmin=fitrangelow[iPtBin];
  Double_t xfmax=fitrangeup[iPtBin];
  if(hname.Contains("Rot") || hname.Contains("ME") || hname.Contains("LS")){
    xfmin=minMass;
    xfmax=maxMass;
  }
  
  TH1F* histof=(TH1F*)histo->Clone(Form("%s_Fl",histo->GetName()));
  AliHFMassFitterVAR* fitter;
  if(!correctForRefl){
    fitter=new AliHFMassFitterVAR(histof,xfmin,xfmax,1,typeback,0); 
    if(typeback==6)fitter->SetBackHighPolDegree(nDegreeBackPol[iPtBin]);
  }
  else {
    fitter=new AliHFMassFitterVAR(histof,xfmin,xfmax,1,typeback,2); 
    if(typeback==6)fitter->SetBackHighPolDegree(nDegreeBackPol[iPtBin]);


    TCanvas *cTest=new TCanvas("cTest","cTest",800,800);    
    TH1F *hmasstemp=fitter->GetHistoClone();
    TH1F *hReflModif=(TH1F*)AdaptRangesAndBinning(hmasstemp,hMCReflPtBin);
    TH1F *hSigModif=(TH1F*)AdaptRangesAndBinning(hmasstemp,hMCSigPtBin);
    hReflModif->SetLineColor(kRed);
    hSigModif->SetLineColor(kBlue);
    hSigModif->Draw();
    hReflModif->Draw("same");
    if(saveReflModifPlot){  
      cTest->SaveAs(Form("cTest%d.root",iPtBin));
    }
 
    delete hmasstemp;
    Double_t fixSoverRefAt=rOverSmodif*(hReflModif->Integral(hReflModif->FindBin(fitrangelow[iPtBin]*1.0001),hReflModif->FindBin(fitrangeup[iPtBin]*0.999))/hSigModif->Integral(hSigModif->FindBin(fitrangelow[iPtBin]*1.0001),hSigModif->FindBin(fitrangeup[iPtBin]*0.999)));
    if(!fitter->SetTemplateReflections(hReflModif,reflopt,fitrangelow[iPtBin],fitrangeup[iPtBin])){
      Printf("SOMETHING WENT WRONG WHILE SETTINGS REFLECTIONS TEMPLATE");
      delete hReflModif;
      delete hSigModif;
      delete cTest; 
      return 0x0;
    }
    fitter->SetFixReflOverS(fixSoverRefAt,kTRUE);//fixSoverRefAt[iBin]);
    delete cTest; 
    delete hReflModif;
    delete hSigModif;
  }
  fitter->SetReflectionSigmaFactor(0);
  fitter->SetInitialGaussianMean(massD);
  fitter->SetInitialGaussianSigma(sigmas[iPtBin]);
  if(fixSigma) fitter->SetFixGaussianSigma(sigmas[iPtBin]);
  if(fixMean) fitter->SetFixGaussianMean(massD);
  fitter->SetFitOption(fitoption.Data());
  return fitter;
}

Double_t GetSignalBinCounting(TH1 *h,TF1 *fbackground,Double_t &err,Double_t nsigmaBC=4.,Double_t sigmafit=0.010,Double_t minx=-999,Double_t maxx=-999){// default range for counting is pdg mass +- nsigmaBC x sigmafit, unless minx and maxx are specified
  
  Int_t binL,binR;
  if(minx>=0.){
    binL=h->FindBin(minx);
    if(binL<1){
      Printf("Left range for bin counting smaller than allowed by histogram axis, setting it to the minimum ");
      binL=1;
    }
  }
  else {
    Double_t lr=massD-sigmafit*nsigmaBC;
    if(minx<0&&minx>-999.){
      if(lr<-1.*minx){
	Printf("Left range for counting set to min of value defined");
	lr=-1.*minx;
      }
    }
    binL=h->FindBin(lr);

    if(binL<1){
      Printf("Left range for bin counting smaller than allowed by histogram axis, setting it to the minimum ");
      binL=1;
    }
  }

  if(maxx>=0.){
    binR=h->FindBin(maxx);
    if(binR>h->GetNbinsX()){
      Printf("Right range for bin counting larger than allowed by histogram axis, setting it to the maximum allowed ");
      binR=h->GetNbinsX();
    }
  }
  else {
    Double_t rr=massD+sigmafit*nsigmaBC;
    if(maxx<0. && maxx>-999){
      if(rr>-1.*maxx){
	Printf("Right range for counting set to max of value defined");
	rr=-1.*maxx;
      }
    }
    binR=h->FindBin(rr);
    if(binR>h->GetNbinsX()){
      Printf("Right range for bin counting larger than allowed by histogram axis, setting it to the maximum allowed ");
      binR=h->GetNbinsX();
    }
  }

  Double_t sign=0;
  err=0;
  for(Int_t j=binL;j<=binR;j++){
    sign+=h->GetBinContent(j)-fbackground->Integral(h->GetBinLowEdge(j),h->GetBinLowEdge(j)+h->GetBinWidth(j))/h->GetBinWidth(j);
    err+=h->GetBinError(j)*h->GetBinError(j);
  }
  err=TMath::Sqrt(err);
  return sign; 
}

void PrintGausParams(TH1F* hPulls){
  TF1* fgfit=(TF1*)hPulls->GetListOfFunctions()->FindObject("gaus");
  TLatex* tg1=new TLatex(0.2,0.8,Form("Gauss mean = %.2f#pm%.2f",fgfit->GetParameter(1),fgfit->GetParError(1)));
  tg1->SetNDC();
  tg1->Draw();
  TLatex* tg2=new TLatex(0.2,0.7,Form("Gauss #sigma = %.2f#pm%.2f",fgfit->GetParameter(2),fgfit->GetParError(2)));
  tg2->SetNDC();
  tg2->Draw();

}

TH1* AdaptRangesAndBinning(const TH1 *hRef,TH1 *h1){
  
  Double_t min=hRef->GetBinLowEdge(1);
  Double_t max=hRef->GetBinLowEdge(hRef->GetNbinsX())+hRef->GetBinWidth(hRef->GetNbinsX());

  Int_t binmin=h1->FindBin(min*1.00001);
  Int_t binmax=h1->FindBin(max*0.99999);

  if(binmin<1||binmax>h1->GetNbinsX()){
    Printf("Cannot adapt range: different bins: Ref histo min, max = %f, %f, hist 1 min %f -> bin min %d , max %f -> bin max %d, nbins= %d",min,max,h1->GetBinLowEdge(1),binmin,h1->GetBinLowEdge(h1->GetNbinsX())+h1->GetBinWidth(1),binmax,h1->GetNbinsX());
    return 0x0;
  }
  if(TMath::Abs(h1->GetBinLowEdge(binmin)-min)>0.0001*min||TMath::Abs(h1->GetBinLowEdge(binmax)+h1->GetBinWidth(binmax)-max)>0.0001*max){
    Printf("Cannot adapt range and rebin histo");
    return 0x0;
  }
  Double_t *binsNew=new Double_t[binmax-binmin+2];
  for(Int_t j=binmin;j<=binmax;j++){
    binsNew[j-binmin]=h1->GetBinLowEdge(j);
  }
  binsNew[binmax-binmin+1]=binsNew[binmax-binmin]+h1->GetBinWidth(binmax-binmin);
  TH1 *hOut;
  TString stype=h1->ClassName();
  if(!stype.Contains("TH1F")){
    hOut=(TH1F*)h1->Rebin(binmax-binmin+1,Form("%s_reb",h1->GetName()),binsNew);
  }
  else if(!stype.Contains("TH1D")){
    hOut=(TH1D*)h1->Rebin(binmax-binmin+1,Form("%s_reb",h1->GetName()),binsNew);
  }
  

  Double_t width=h1->GetBinWidth(binmin);
  Double_t widthRef=hRef->GetBinWidth(1);
  Int_t rebinh=1;
  if((TMath::Abs(width-widthRef)<0.00001*widthRef)){
    rebinh=1;
  }
  else{
    Double_t wr;
    if(width>widthRef){
      Printf("Should rebin the reference histo by %d",(Int_t)(width/widthRef));
      return 0x0;
    }
    else{
      wr=widthRef/width;
      if(TMath::Abs(wr-(Int_t)wr)>0.00001){
	Printf("Cannot adapt histo: rebin issue, width ref = %f, width h1=%f",widthRef,width);
	return 0x0;
      }
      else rebinh=(Int_t)wr;
    }
  }
  hOut->Rebin(rebinh);
   
  return hOut;
}

Bool_t QuadraticSmooth(TH1 *h,Int_t ntimes=1){// quadratic fit of 5 points
  ntimes--;
  TF1 *fp2=new TF1("fp2","pol2",h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  TH1D *htmp=(TH1D*)h->Clone("htmp");
  for(Int_t i=3;i<=h->GetNbinsX()-3;i++){ // first intermediate bins
    htmp->Fit(fp2,"RLEMN0","",h->GetXaxis()->GetBinLowEdge(i-2),h->GetXaxis()->GetBinUpEdge(i+2));
    h->SetBinContent(i,(fp2->Integral(h->GetXaxis()->GetBinLowEdge(i),h->GetXaxis()->GetBinUpEdge(i)))/h->GetBinWidth(i));// change only content, errors unchanged
  }
  // now initial and final bins
  htmp->Fit(fp2,"RLEMN0","",h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(5));// 5 pts, asymmetric fit region
  h->SetBinContent(2,(fp2->Integral(h->GetXaxis()->GetBinLowEdge(2),h->GetXaxis()->GetBinUpEdge(2)))/h->GetBinWidth(2));// change only content, errors unchanged

  htmp->Fit(fp2,"RLEMN0","",h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(4));//  use only 4 pts
  h->SetBinContent(1,(fp2->Integral(h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(1)))/h->GetBinWidth(1));// change only content, errors unchanged


  htmp->Fit(fp2,"RLEMN0","",h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()-4),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));// 5 pts, asymmetric fit region
  h->SetBinContent(h->GetNbinsX()-1,(fp2->Integral(h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()-1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()-1)))/h->GetBinWidth(h->GetNbinsX()-1));// change only content, errors unchanged

  htmp->Fit(fp2,"RLEMN0","",h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()-3),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  h->SetBinContent(h->GetNbinsX(),(fp2->Integral(h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX())))/h->GetBinWidth(h->GetNbinsX()));// change only content, errors unchanged
  if(ntimes>0)QuadraticSmooth(h,ntimes);
  delete htmp;
  return kTRUE; 
}


void SetStyleHisto(TH1 *h,Int_t method,Int_t isXpt=-1){
  if(isXpt==1)h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");// just for convenience, most of the time it's pt :)

  if(method==kRot){
    h->SetMarkerStyle(21);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(kBlack);
    h->SetMarkerSize(1.0);
    h->SetLineColor(kBlack);
    return;
  }

  if(method==kME){
    h->SetMarkerStyle(25);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(4);
    h->SetMarkerSize(1.0);
    h->SetLineColor(4);
    return;
  }

  if(method==kLS){
    h->SetMarkerStyle(22);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(kGreen+2);
    h->SetMarkerSize(1.0);
    h->SetLineColor(kGreen+2);
    return;
  }
  
  if(method==kSB){
    h->SetMarkerStyle(27);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(6);
    h->SetMarkerSize(1.0);
    h->SetLineColor(6);
    return;
  }

  return;
}

void DivideCanvas(TCanvas *c,Int_t ndivisions){
  if(ndivisions>1){
    if(ndivisions<3){
      c->Divide(2,1);
    }
    else if(ndivisions<5){
      c->Divide(2,2);
    }
    else if(ndivisions<7){
      c->Divide(3,2);

    }
    else if(ndivisions<=8){
      c->Divide(4,2);

    }
    else if(ndivisions<10){
      c->Divide(3,3);

    }
    else if(ndivisions<13){
      c->Divide(4,3);

    }
    else if(ndivisions<17){
      c->Divide(4,4);

    }
    else {
      c->Divide(5,5);

    }
  }
  return;
}



TF1 *GausPlusLine(Double_t minRange=1.72,Double_t maxRange=2.05){
  TF1 *f=new TF1("GausPlusLine","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]+[4]*x",minRange,maxRange);

  f->SetParameter(1,massD);
  f->SetParLimits(0,0.,1.e+7);
  f->SetParLimits(1,1.8,1.92);
  f->SetParameter(2,0.010);
  f->SetParLimits(2,0,0.1);
  f->SetParameter(3,0.);
  //  f->SetParLimits(3,-10000,10000);
  return f;

}


Double_t GetBackgroundNormalizationFactor(TH1D* hRatio){
  //
  Double_t norm=hRatio->GetMaximum();

  if(optForNorm==1){
      norm=0.0001;
      for(Int_t iMassBin=1; iMassBin<hRatio->GetNbinsX(); iMassBin++){
	Double_t bce=hRatio->GetBinCenter(iMassBin);
	if(bce>minMass && bce<maxMass){
	  Double_t bco=hRatio->GetBinContent(iMassBin);
	  if(bco>norm) norm=bco;
	}
      }
  }else if(optForNorm==2){ 
    hRatio->Fit("pol0","","",minMass,minMass+rangeForNorm);
    TF1* func0=(TF1*)hRatio->GetListOfFunctions()->FindObject("pol0");
    Double_t norml=func0->GetParameter(0);
    hRatio->Fit("pol0","","",maxMass-rangeForNorm,maxMass);
    func0=(TF1*)hRatio->GetListOfFunctions()->FindObject("pol0");
    Double_t normh=func0->GetParameter(0);
    norm=TMath::Max(norml,normh);
  }
  return norm;
}

void ProjectCombinHFAndFit(){

  
  
  TString dirName=Form("PWG3_D2H_InvMass%sLowPt%s",meson.Data(),suffix.Data());
  TString lstName=Form("coutput%s%s",meson.Data(),suffix.Data());

  if(meson=="Dplus") massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  else if(meson=="Dzero") massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  TFile* fil=new TFile(fileName.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dirName.Data());
  
  TFile* filMC;
  TDirectoryFile* dfMC;
  TList* lMC;
  TH3F* h3drefl;
  TH3F* h3dmcsig;
  
  AliNormalizationCounter *nC=(AliNormalizationCounter*)df->Get("NormalizationCounter");
  TH1F* hnEv=new TH1F("hEvForNorm","events for normalization",1,0,1);
  printf("Number of Ev. for norm=%d\n",(Int_t)nC->GetNEventsForNorm());
  hnEv->SetBinContent(1,nC->GetNEventsForNorm());

  TList* l=(TList*)df->Get(lstName.Data());
  l->ls();

  TH3F* h3d=(TH3F*)l->FindObject("hMassVsPtVsY");
  TH3F* h3dr=(TH3F*)l->FindObject("hMassVsPtVsYRot");
  TH3F* h3dme=(TH3F*)l->FindObject("hMassVsPtVsYME");
  TH3F* h3dmepp=(TH3F*)l->FindObject("hMassVsPtVsYMELSpp");
  TH3F* h3dmemm=(TH3F*)l->FindObject("hMassVsPtVsYMELSmm");
  TH2F* hpoolMix=(TH2F*)l->FindObject("hMixingsPerPool");
  TH2F* hpoolEv=(TH2F*)l->FindObject("hEventsPerPool");
  TH3F* h3dlsp=(TH3F*)l->FindObject("hMassVsPtVsYLSpp");
  TH3F* h3dlsm=(TH3F*)l->FindObject("hMassVsPtVsYLSmm");

  if(correctForRefl){
    filMC=new TFile(fileNameMC.Data());
    dfMC=(TDirectoryFile*)filMC->Get(dirName.Data());

    lMC=(TList*)dfMC->Get(lstName.Data());
    h3drefl=(TH3F*)lMC->FindObject("hMassVsPtVsYRefl");
    h3dmcsig=(TH3F*)lMC->FindObject("hMassVsPtVsYSig");
  }


  TCanvas* cem=new TCanvas("cem","Pools",1200,600);
  cem->Divide(2,1);
  cem->cd(1);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  hpoolMix->GetXaxis()->SetTitle("zVertex");
  hpoolMix->GetYaxis()->SetTitle("Multiplicity");
  hpoolMix->GetYaxis()->SetTitleOffset(1.4);
  hpoolMix->Draw("colztext");
  cem->cd(2);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  hpoolEv->Draw("colztext");
  hpoolEv->GetXaxis()->SetTitle("zVertex");
  hpoolEv->GetYaxis()->SetTitle("Multiplicity");
  hpoolEv->GetYaxis()->SetTitleOffset(1.4);


  TCanvas* c1=new TCanvas("c1","Mass",1200,800);
  DivideCanvas(c1,nPtBins);

  TCanvas* c2=new TCanvas("c2","Mass-Bkg Rot",1200,800);
  DivideCanvas(c2,nPtBins);
  TCanvas* c2pulls=new TCanvas("c2pulls","Mass-Bkg Rot pulls",1200,800);
  DivideCanvas(c2pulls,nPtBins);
  TCanvas* c2residuals=new TCanvas("c2residuals","Mass-Bkg Rot residuals",1200,800);
  DivideCanvas(c2residuals,nPtBins);
  TCanvas* c2residualTrend=new TCanvas("c2residualTrend","Mass-Bkg Rot residual trend vs. mass",1200,800);
  DivideCanvas(c2residualTrend,nPtBins);
  TCanvas* c2pullTrend=new TCanvas("c2pullTrend","Mass-Bkg Rot pull trend vs. mass",1200,800);
  DivideCanvas(c2pullTrend,nPtBins);

  TCanvas* c3=new TCanvas("c3","Mass-Bkg LS",1200,800);
  DivideCanvas(c3,nPtBins);
  TCanvas* c3pulls=new TCanvas("c3pulls","Mass-Bkg LS pulls",1200,800);
  DivideCanvas(c3pulls,nPtBins);
  TCanvas* c3residuals=new TCanvas("c3residuals","Mass-Bkg LS residuals",1200,800);
  DivideCanvas(c3residuals,nPtBins);
  TCanvas* c3residualTrend=new TCanvas("c3residualTrend","Mass-Bkg LS residual trend vs. mass",1200,800);
  DivideCanvas(c3residualTrend,nPtBins);
  TCanvas* c3pullTrend=new TCanvas("c3pullTrend","Mass-Bkg LS pull trend vs. mass",1200,800);
  DivideCanvas(c3pullTrend,nPtBins);

  TCanvas* c4=new TCanvas("c4","Mass-Bkg EM",1200,800);
  DivideCanvas(c4,nPtBins);
  TCanvas* c4pulls=new TCanvas("c4pulls","Mass-Bkg EM pulls",1200,800);
  DivideCanvas(c4pulls,nPtBins);
  TCanvas* c4residuals=new TCanvas("c4residuals","Mass-Bkg EM residuals",1200,800);
  DivideCanvas(c4residuals,nPtBins);
  TCanvas* c4residualTrend=new TCanvas("c4residualTrend","Mass-Bkg EM residual trend vs. mass",1200,800);
  DivideCanvas(c4residualTrend,nPtBins);
  TCanvas* c4pullTrend=new TCanvas("c4pullTrend","Mass-Bkg EM pull trend vs. mass",1200,800);
  DivideCanvas(c4pullTrend,nPtBins);

  TCanvas *c5=0x0;
  TCanvas *c5sub=0x0;
  TCanvas *c5pulls=0x0;
  TCanvas *c5residuals=0x0;
  TCanvas *c5residualTrend=0x0;
  TCanvas *c5pullTrend=0x0;

  if(tryDirectFit){
    c5=new TCanvas("c5","Mass SB Fit",1200,800);
    DivideCanvas(c5,nPtBins);
    c5sub=new TCanvas("c5sub","Mass-Bkg SB",1200,800);
    DivideCanvas(c5sub,nPtBins);
    c5pulls=new TCanvas("c5pulls","Mass-Bkg SB pulls",1200,800);
    DivideCanvas(c5pulls,nPtBins);
    c5residuals=new TCanvas("c5residuals","Mass-Bkg SB residuals",1200,800);
    DivideCanvas(c5residuals,nPtBins);
    c5residualTrend=new TCanvas("c5residualTrend","Mass-Bkg SB residual trend vs. mass",1200,800);
    DivideCanvas(c5residualTrend,nPtBins);
    c5pullTrend=new TCanvas("c5pullTrend","Mass-Bkg SB pull trend vs. mass",1200,800);
    DivideCanvas(c5pullTrend,nPtBins);    
  }

  TCanvas *cCompareResidualTrends=new TCanvas("cCompareResidualTrends","cCompareResidualTrends",1200,800);
  DivideCanvas(cCompareResidualTrends,nPtBins);  


  AliHFMassFitter* fitterRot[nPtBins];
  AliHFMassFitter* fitterLS[nPtBins];
  AliHFMassFitter* fitterME[nPtBins];
  AliHFMassFitterVAR* fitterSB[nPtBins];

  TH1F* hRawYieldRot=new TH1F("hRawYieldRot","",nPtBins,binLims);
  TH1F* hRawYieldLS=new TH1F("hRawYieldLS","",nPtBins,binLims);
  TH1F* hRawYieldME=new TH1F("hRawYieldME","",nPtBins,binLims);
  TH1F* hRawYieldSBfit=new TH1F("hRawYieldSBfit","",nPtBins,binLims);

  TH1F* hChiSqRot=new TH1F("hChiSqRot","",nPtBins,binLims);
  TH1F* hChiSqLS=new TH1F("hChiSqLS","",nPtBins,binLims);
  TH1F* hChiSqME=new TH1F("hChiSqME","",nPtBins,binLims);
  TH1F* hChiSqSBfit=new TH1F("hChiSqSBfit","",nPtBins,binLims);

  TH1F* hNdfRot=new TH1F("hNdfRot","",nPtBins,binLims);
  TH1F* hNdfLS=new TH1F("hNdfLS","",nPtBins,binLims);
  TH1F* hNdfME=new TH1F("hNdfME","",nPtBins,binLims);
  TH1F* hNdfSBfit=new TH1F("hNdfSBfit","",nPtBins,binLims);

  TH1F* hRawYieldRotBC=new TH1F("hRawYieldRotBC","BC yield (rotational background)",nPtBins,binLims);
  TH1F* hRawYieldLSBC=new TH1F("hRawYieldLSBC","BC yield (like-sign background)",nPtBins,binLims);
  TH1F* hRawYieldMEBC=new TH1F("hRawYieldMEBC","BC yield (mixed-event background)",nPtBins,binLims);
  TH1F* hRawYieldSBfitBC=new TH1F("hRawYieldSBfitBC","BC yield (direct fit background)",nPtBins,binLims);

  TLatex* tME=new TLatex(0.65,0.82,"MixEv +- pairs");
  tME->SetNDC();
  TLatex* tMEpp=new TLatex(0.65,0.75,"MixEv ++ pairs");
  tMEpp->SetNDC();
  TLatex* tMEmm=new TLatex(0.65,0.68,"MixEv -- pairs");
  tMEmm->SetNDC();

  TF1 *fpeak=new TF1("fpeak","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",minMass,maxMass);

  TDirectory *current = gDirectory;
  TFile* fout=new TFile(Form("outputMassFits_%s.root",suffix.Data()),"recreate");
  current->cd();


  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){

    Int_t bin1=h3d->GetYaxis()->FindBin(binLims[iPtBin]);
    Int_t bin2=h3d->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
    printf("Bin %d   Pt range=%f %f\n",iPtBin,h3d->GetYaxis()->GetBinLowEdge(bin1),h3d->GetYaxis()->GetBinUpEdge(bin2));
    TH1D* hMassPtBin=h3d->ProjectionX(Form("hMassPtBin%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinr=h3dr->ProjectionX(Form("hMassPtBinr%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinme=h3dme->ProjectionX(Form("hMassPtBinme%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinmeLSpp=h3dmepp->ProjectionX(Form("hMassPtBinmeLSpp%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinmeLSmm=h3dmemm->ProjectionX(Form("hMassPtBinmeLSmm%d",iPtBin),bin1,bin2);

    if(correctForRefl){
      Int_t bin1MC=h3drefl->GetYaxis()->FindBin(binLims[iPtBin]);
      Int_t bin2MC=h3drefl->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
      hMCReflPtBin=h3drefl->ProjectionX(Form("hMCReflPtBin%d",iPtBin),bin1MC,bin2MC);
      hMCSigPtBin=h3dmcsig->ProjectionX(Form("hMCMassPtBin%d",iPtBin),bin1MC,bin2MC);
    }

    TH1D* hMassPtBinlsp=0x0;
    TH1D* hMassPtBinlsm=0x0;
    TH1D* hMassPtBinls=0x0;
    if(h3dlsp){ 
      hMassPtBinlsp=h3dlsp->ProjectionX(Form("hMassPtBinlsp%d",iPtBin),bin1,bin2);
      hMassPtBinlsm=h3dlsm->ProjectionX(Form("hMassPtBinlsm%d",iPtBin),bin1,bin2);
      hMassPtBinls=(TH1D*)hMassPtBinlsp->Clone(Form("hMassPtBinls%d",iPtBin));
      hMassPtBinls->Reset();
      for(Int_t iBin=1; iBin<=hMassPtBinlsp->GetNbinsX(); iBin++){
	Double_t np=hMassPtBinlsp->GetBinContent(iBin);
	Double_t nm=hMassPtBinlsm->GetBinContent(iBin);
	Double_t tt=2*TMath::Sqrt(np*nm);
	Double_t enp=hMassPtBinlsp->GetBinError(iBin);
	Double_t enm=hMassPtBinlsm->GetBinError(iBin);
	Double_t ett=0;
	if(tt>0) ett=2./tt*TMath::Sqrt(np*np*enm*enm+nm*nm*enp*enp);
	hMassPtBinls->SetBinContent(iBin,tt);
	hMassPtBinls->SetBinError(iBin,ett);
      }
      if(smoothLS<-0.5)hMassPtBinls->Smooth(-1*smoothLS);
      else if(smoothLS>0.5)QuadraticSmooth(hMassPtBinls,smoothLS);
      hMassPtBinls->SetLineColor(kGreen+1);
    }

    // hMassPtBin->Sumw2();
    // hMassPtBinr->Sumw2();
    // hMassPtBinme->Sumw2();
    // hMassPtBinlsp->Sumw2();
    // hMassPtBinls->Sumw2();
    hMassPtBin->SetLineColor(1);
    hMassPtBinr->SetLineColor(4);
    hMassPtBinme->SetLineColor(kOrange+1);
    hMassPtBinmeLSpp->SetLineColor(kGreen+2);
    hMassPtBinmeLSmm->SetLineColor(kCyan);
    hMassPtBin->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    hMassPtBinr->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    hMassPtBinme->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    TH1D* hRatioMEpp=(TH1D*)hMassPtBinmeLSpp->Clone(Form("hRatioPtBinmeLSpp%d",iPtBin));
    hRatioMEpp->Divide(hMassPtBinme);
    TH1D* hRatioMEmm=(TH1D*)hMassPtBinmeLSmm->Clone(Form("hRatioPtBinmeLSmm%d",iPtBin));
    hRatioMEmm->Divide(hMassPtBinme);
    TH1D* hMassPtBinmeAll=(TH1D*)hMassPtBinme->Clone(Form("hRatioPtBinmeAll%d",iPtBin));
    hMassPtBinmeAll->Add(hMassPtBinmeLSpp);
    hMassPtBinmeAll->Add(hMassPtBinmeLSmm);
    hMassPtBinmeAll->SetLineColor(kRed);
    TH1D* hRatioME=(TH1D*)hMassPtBinme->Clone(Form("hRatioME%d",iPtBin));
    hRatioME->Divide(hMassPtBin);
    TH1D* hRatioMEAll=(TH1D*)hMassPtBinmeAll->Clone(Form("hRatioMEAll%d",iPtBin));
    hRatioMEAll->Divide(hMassPtBin);


    TCanvas* c0=new TCanvas(Form("CBin%d",iPtBin),Form("Bin%d norm",iPtBin),1000,700);
    c0->Divide(2,2);
    c0->cd(1);
    TH1D* hRatio=(TH1D*)hMassPtBinr->Clone("hRatio");
    hRatio->Divide(hMassPtBin);
    hRatio->Draw();
    hRatio->GetYaxis()->SetTitle("Rotational/All");
    hRatio->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    Double_t normRot=GetBackgroundNormalizationFactor(hRatio);
    hMassPtBinr->Scale(1./normRot);
    c0->cd(2);
    hMassPtBinme->GetYaxis()->SetTitle("Entries (EvMix)");
    hMassPtBinme->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hMassPtBinme->DrawCopy();
    hMassPtBinmeLSpp->Draw("same");
    hMassPtBinmeLSmm->Draw("same");
    tME->SetTextColor(hMassPtBinme->GetLineColor());
    tMEpp->SetTextColor(hMassPtBinmeLSpp->GetLineColor());
    tMEmm->SetTextColor(hMassPtBinmeLSmm->GetLineColor());
    tME->Draw();
    tMEpp->Draw();
    tMEmm->Draw();
    c0->cd(4);
    hRatioMEpp->Draw();
    hRatioMEpp->SetMinimum(0.4);
    hRatioMEpp->SetMaximum(0.6);
    hRatioMEpp->GetYaxis()->SetTitle("ME with LS / ME with OS");
    hRatioMEpp->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hRatioMEmm->Draw("same");
    c0->cd(3);
    hRatioME->Draw();
    hRatioME->SetMaximum(hRatioMEAll->GetMaximum()*1.05);
    hRatioMEAll->Draw("same");
    hRatioME->GetYaxis()->SetTitle("EvMix/All");
    hRatioME->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
 
    Double_t normME=GetBackgroundNormalizationFactor(hRatioME);
    Double_t normMEAll=GetBackgroundNormalizationFactor(hRatioMEAll);
    hMassPtBinme->Scale(1./normME);
    hMassPtBinmeAll->Scale(1./normMEAll);

    c1->cd(iPtBin+1);
    hMassPtBin->Draw();
    hMassPtBin->GetYaxis()->SetTitle("Counts");
    hMassPtBin->GetYaxis()->SetTitleOffset(2.);
    hMassPtBin->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hMassPtBinr->Draw("same");
    if(!useEMwithLS) hMassPtBinme->Draw("same");
    else  hMassPtBinmeAll->Draw("same");
    if(hMassPtBinls) hMassPtBinls->Draw("same");
    if(iPtBin==0){
      TLegend* leg=new TLegend(0.5,0.6,0.89,0.89);
      leg->SetFillStyle(0);
      leg->AddEntry(hMassPtBin,"All candidates","L")->SetTextColor(hMassPtBin->GetLineColor());
      leg->AddEntry(hMassPtBinr,"Background (rotations)","L")->SetTextColor(hMassPtBinr->GetLineColor());
      if(!useEMwithLS) leg->AddEntry(hMassPtBinme,"Background (ME)","L")->SetTextColor(hMassPtBinme->GetLineColor());
      else leg->AddEntry(hMassPtBinmeAll,"Background (ME)","L")->SetTextColor(hMassPtBinmeAll->GetLineColor());
      if(hMassPtBinls) leg->AddEntry(hMassPtBinls,"Like-sign","L")->SetTextColor(hMassPtBinls->GetLineColor());
      leg->Draw();
    }
  
    TH1D* hMassSubRot=(TH1D*)hMassPtBin->Clone(Form("hMassSubRot_bin%d",iPtBin));
    hMassSubRot->Add(hMassPtBinr,-1);
    hMassSubRot->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c -- Rotational",binLims[iPtBin],binLims[iPtBin+1]));
    TH1D* hMassSubME=(TH1D*)hMassPtBin->Clone(Form("hMassSubME_bin%d",iPtBin));
    if(!useEMwithLS) hMassSubME->Add(hMassPtBinme,-1);
    else hMassSubME->Add(hMassPtBinmeAll,-1);
    hMassSubME->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c -- Mixed Ev",binLims[iPtBin],binLims[iPtBin+1]));
    TH1D* hMassSubLS=0x0;
    if(hMassPtBinls){
      hMassSubLS=(TH1D*)hMassPtBin->Clone(Form("hMassSubLS_bin%d",iPtBin));
      hMassSubLS->Add(hMassPtBinls,-1);
      hMassSubLS->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c  -- Like Sign",binLims[iPtBin],binLims[iPtBin+1]));
    }
 
    fout->cd();
    hMassPtBin->Write();
    hMassSubRot->Write();
    hMassSubME->Write();
    if(hMassPtBinls) hMassSubLS->Write();
    current->cd();

    hMassSubRot->Rebin(rebin[iPtBin]);
    hMassSubME->Rebin(rebin[iPtBin]);
    if(hMassPtBinls) hMassSubLS->Rebin(rebin[iPtBin]);
 
    if(!correctForRefl){
      fitterRot[iPtBin]=ConfigureFitter(hMassSubRot,iPtBin);
      if(hMassPtBinls) fitterLS[iPtBin]=ConfigureFitter(hMassSubLS,iPtBin);
      fitterME[iPtBin]=ConfigureFitter(hMassSubME,iPtBin);
    }
    else {
      fitterRot[iPtBin]=ConfigureFitterVAR(hMassSubRot,iPtBin,typeb);
      if(hMassPtBinls) fitterLS[iPtBin]=ConfigureFitterVAR(hMassSubLS,iPtBin,typeb);
      fitterME[iPtBin]=ConfigureFitterVAR(hMassSubME,iPtBin,typeb);  
    }

    Bool_t out1=fitterRot[iPtBin]->MassFitter(0);
    Bool_t out2=kFALSE;
    if(hMassPtBinls) out2=fitterLS[iPtBin]->MassFitter(0);
    Bool_t out3=fitterME[iPtBin]->MassFitter(0);
    
    Bool_t out4=kFALSE;
    if(tryDirectFit){
      TH1D *hMassDirectFit=(TH1D*)hMassPtBin->Clone(Form("hMassDirectFit_bin%d",iPtBin));
      hMassDirectFit->Rebin(rebin[iPtBin]);
      fitterSB[iPtBin]=ConfigureFitterVAR(hMassDirectFit,iPtBin);
      out4=fitterSB[iPtBin]->MassFitter(0);//DirectFit(hMassDirectFit,iPtBin,hRawYieldSBfit);

      if(out4){
	c5->cd(iPtBin+1);
	if(!correctForRefl)fitterSB[iPtBin]->DrawHere(gPad,3,0);
	else fitterSB[iPtBin]->DrawHere(gPad,3,1);
	hRawYieldSBfit->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield());
	hRawYieldSBfit->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetRawYieldError());
	hChiSqSBfit->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetReducedChiSquare());
	hChiSqSBfit->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	hNdfSBfit->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetMassFunc()->GetNDF());
	hNdfSBfit->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	//	if(!correctForRefl)
	WriteFitInfo(fitterSB[iPtBin],hMassPtBin);
	

	c5sub->cd(iPtBin+1);
	TH1F* hsubTemp=(TH1F*)hMassDirectFit->Clone(Form("%sSubBack%d",hMassDirectFit->GetName(),iPtBin));
	TH1F* hsubTempAllRange=(TH1F*)hMassDirectFit->Clone(Form("%sSubBackAllRange%d",hMassDirectFit->GetName(),iPtBin));
	TF1* funcAll=fitterSB[iPtBin]->GetMassFunc();
	TF1* funcBkg=fitterSB[iPtBin]->GetBackgroundRecalcFunc();
	//	TF1* funcBkg2=fitterSB[iPtBin]->GetBackgroundFullRangeFunc();	
	for(Int_t jst=1;jst<=hsubTemp->GetNbinsX();jst++){
	  Double_t backg=funcBkg->Integral(hsubTemp->GetBinLowEdge(jst),hsubTemp->GetBinLowEdge(jst)+hsubTemp->GetBinWidth(jst))/hsubTemp->GetBinWidth(jst);
	  Double_t tot=funcAll->Integral(hsubTempAllRange->GetBinLowEdge(jst),hsubTempAllRange->GetBinLowEdge(jst)+hsubTempAllRange->GetBinWidth(jst))/hsubTempAllRange->GetBinWidth(jst);
	  hsubTemp->SetBinContent(jst,hsubTemp->GetBinContent(jst)-backg);
	  hsubTempAllRange->SetBinContent(jst,hsubTempAllRange->GetBinContent(jst)-tot);
	}
	hsubTemp->SetLineColor(kBlue);
	hsubTempAllRange->SetLineColor(kGray+2);
	Double_t ymin=0;
	Double_t ymax=1;
	for(Int_t ibs=1; ibs<hsubTemp->GetNbinsX(); ibs++){
	  Double_t binc=hsubTemp->GetBinCenter(ibs);
	  if(binc>fitrangelow[iPtBin] && binc<fitrangeup[iPtBin]){
	    Double_t yl=hsubTemp->GetBinContent(ibs)-hsubTemp->GetBinError(ibs);
	    Double_t yu=hsubTemp->GetBinContent(ibs)+hsubTemp->GetBinError(ibs);
	    if(yl<ymin) ymin=yl;
	    if(yu>ymax) ymax=yu;
	  }
	}
	if(ymax>0) ymax*=1.2;
	else ymax*=0.8;
	if(ymin<0) ymin*=1.2;
	else ymin*=0.8;
	hsubTemp->GetXaxis()->SetRangeUser(fitrangelow[iPtBin],fitrangeup[iPtBin]);
	hsubTemp->SetMinimum(ymin);
	hsubTemp->SetMaximum(ymax);
	hsubTemp->SetMarkerStyle(20);
	hsubTemp->SetMarkerColor(hsubTemp->GetLineColor());
	hsubTemp->DrawCopy();
	hsubTempAllRange->DrawCopy("same");
	hsubTemp->DrawCopy("same");
	fpeak->SetParameter(0,funcAll->GetParameter(nDegreeBackPol[iPtBin]+1));
	fpeak->SetParameter(1,funcAll->GetParameter(nDegreeBackPol[iPtBin]+2));
	fpeak->SetParameter(2,funcAll->GetParameter(nDegreeBackPol[iPtBin]+3));
	fpeak->DrawCopy("same");

	Double_t errbc;
	Double_t bc=GetSignalBinCounting(hMassDirectFit,fitterSB[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterSB[iPtBin]->GetSigma(),-1.*fitrangelow[iPtBin],-1.*fitrangeup[iPtBin]);
	hRawYieldSBfitBC->SetBinContent(iPtBin+1,bc);
	hRawYieldSBfitBC->SetBinError(iPtBin+1,errbc);
      
	c5pulls->cd(iPtBin+1);
	TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
	TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
	TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
	TH1F *hResiduals=fitterSB[iPtBin]->GetResidualsAndPulls(hMassDirectFit,fitterSB[iPtBin]->GetMassFunc(),fitrangelow[iPtBin],fitrangeup[iPtBin],hPulls,hResidualTrend,hPullsTrend);

	hPulls->Draw();
	PrintGausParams(hPulls);
	
	c5residuals->cd(iPtBin+1);
	hResiduals->Draw();
	
	SetStyleHisto(hResidualTrend,kSB);
	cCompareResidualTrends->cd(iPtBin+1);
	hResidualTrend->Draw();

	c5residualTrend->cd(iPtBin+1);
	hResidualTrend->Draw();
	
	c5pullTrend->cd(iPtBin+1);
	hPullsTrend->Draw();
	//TPaveStats *tp=(TPaveStats*)hPulls->FindObject("stats");
	//      tp->SetOptStat("remnpcev");
	delete hMassDirectFit;
      }
    }

    c2->cd(iPtBin+1);
    if(out1){
      if(!correctForRefl)fitterRot[iPtBin]->DrawHere(gPad,3,0);
      else static_cast<AliHFMassFitterVAR*>(fitterRot[iPtBin])->AliHFMassFitterVAR::DrawHere(gPad,3,1);
      hRawYieldRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetRawYield());
      hRawYieldRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetRawYieldError());
      hChiSqRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetReducedChiSquare());
      hChiSqRot->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      hNdfRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetMassFunc()->GetNDF());
      hNdfRot->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      //      if(!correctForRefl)
      WriteFitInfo(fitterRot[iPtBin],hMassPtBin);
      
      Double_t errbc;
      Double_t bc=GetSignalBinCounting(hMassSubRot,fitterRot[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterRot[iPtBin]->GetSigma(),-1.*minMass,-1.*maxMass);
      hRawYieldRotBC->SetBinContent(iPtBin+1,bc);
      hRawYieldRotBC->SetBinError(iPtBin+1,errbc);

      c2pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterRot[iPtBin]->GetResidualsAndPulls(hMassSubRot,fitterRot[iPtBin]->GetMassFunc(),minMass,maxMass,hPulls,hResidualTrend,hPullsTrend);

      hPulls->Draw();
      PrintGausParams(hPulls);

      c2residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kRot);
      c2residualTrend->cd(iPtBin+1);
      hResidualTrend->Draw();
      cCompareResidualTrends->cd(iPtBin+1);
      if(out4)hResidualTrend->Draw("same");
      else hResidualTrend->Draw();

      c2pullTrend->cd(iPtBin+1);
      hPullsTrend->Draw();

    }
    else hMassSubRot->Draw("");
    

    c3->cd(iPtBin+1);
    if(out2){
      if(!correctForRefl)fitterLS[iPtBin]->DrawHere(gPad,3,0); 
      else static_cast<AliHFMassFitterVAR*>(fitterLS[iPtBin])->AliHFMassFitterVAR::DrawHere(gPad,3,1); 
      hRawYieldLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetRawYield());
      hRawYieldLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetRawYieldError());
      hChiSqLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetReducedChiSquare());
      hChiSqLS->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      hNdfLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetMassFunc()->GetNDF());
      hNdfLS->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      //      if(!correctForRefl)
WriteFitInfo(fitterLS[iPtBin],hMassPtBin);

      Double_t errbc;
      Double_t bc=GetSignalBinCounting(hMassSubLS,fitterLS[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterLS[iPtBin]->GetSigma(),-1.*minMass,-1.*maxMass);
      hRawYieldLSBC->SetBinContent(iPtBin+1,bc);
      hRawYieldLSBC->SetBinError(iPtBin+1,errbc);

      c3pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterLS[iPtBin]->GetResidualsAndPulls(hMassSubLS,fitterLS[iPtBin]->GetMassFunc(),minMass,maxMass,hPulls,hResidualTrend,hPullsTrend);
  
      hPulls->Draw();
      PrintGausParams(hPulls);

      c3residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kLS);
      c3residualTrend->cd(iPtBin+1);
      hResidualTrend->Draw();
      
      cCompareResidualTrends->cd(iPtBin+1);
      if(out4||out1)hResidualTrend->Draw("same");
      else hResidualTrend->Draw();

      c3pullTrend->cd(iPtBin+1);
      hPullsTrend->Draw();
    }
    else if(hMassPtBinls) hMassSubLS->Draw("");

    c4->cd(iPtBin+1);
    if(out3){ 
      if(!correctForRefl)fitterME[iPtBin]->DrawHere(gPad,3,0); 
      else static_cast<AliHFMassFitterVAR*>(fitterME[iPtBin])->AliHFMassFitterVAR::DrawHere(gPad,3,1); 
      hRawYieldME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetRawYield());
      hRawYieldME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetRawYieldError());
      hChiSqME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetReducedChiSquare());
      hChiSqME->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      hNdfME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetMassFunc()->GetNDF());
      hNdfME->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      //      if(!correctForRefl)
      WriteFitInfo(fitterME[iPtBin],hMassPtBin);

      Double_t errbc;
      Double_t bc=GetSignalBinCounting(hMassSubME,fitterME[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterME[iPtBin]->GetSigma(),-1.*minMass,-1.*maxMass);
      hRawYieldMEBC->SetBinContent(iPtBin+1,bc);
      hRawYieldMEBC->SetBinError(iPtBin+1,errbc);


      c4pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterME[iPtBin]->GetResidualsAndPulls(hMassSubME,fitterME[iPtBin]->GetMassFunc(),minMass,maxMass,hPulls,hResidualTrend,hPullsTrend);
      hPulls->Draw();
      PrintGausParams(hPulls);

      c4residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kME);
      c4residualTrend->cd(iPtBin+1);
      hResidualTrend->Draw();

      cCompareResidualTrends->cd(iPtBin+1);
      if(out4||out1||out2)hResidualTrend->Draw("same");
      else hResidualTrend->Draw();

      c4pullTrend->cd(iPtBin+1);
      hPullsTrend->Draw();
    }
    else hMassSubME->Draw("");
    if(correctForRefl){
      delete hMCReflPtBin;
      delete hMCSigPtBin;
    }

  }
  TString path(gSystem->pwd());
  path.Append("/figures");
  if(gSystem->AccessPathName(path.Data())){ 
    gROOT->ProcessLine(Form(".!mkdir -p %s",path.Data()));  
  }

  if(saveCanvasAsEps>0){
    c2->SaveAs(Form("figures/InvMassSpectra_%s_Rot.eps",suffix.Data()));
    c3->SaveAs(Form("figures/InvMassSpectra_%s_LS.eps",suffix.Data()));
    c4->SaveAs(Form("figures/InvMassSpectra_%s_EM.eps",suffix.Data()));
    if(tryDirectFit){
      c5->SaveAs(Form("figures/InvMassSpectra_%s_SB.eps",suffix.Data()));
      c5sub->SaveAs(Form("figures/InvMassSpectra_%s_SBsub.eps",suffix.Data()));
      c5sub->SaveAs(Form("figures/InvMassSpectra_%s_SBsub.root",suffix.Data()));
    }

    if(saveCanvasAsEps>1){
      c2residuals->SaveAs(Form("figures/ResidualDistribution_%s_Rot.eps",suffix.Data()));
      c3residuals->SaveAs(Form("figures/ResidualDistribution_%s_LS.eps",suffix.Data()));
      c4residuals->SaveAs(Form("figures/ResidualDistribution_%s_EM.eps",suffix.Data()));

      c2residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_Rot.eps",suffix.Data()));
      c3residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_LS.eps",suffix.Data()));
      c4residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_EM.eps",suffix.Data()));

      c2pulls->SaveAs(Form("figures/PullDistribution_%s_Rot.eps",suffix.Data()));
      c3pulls->SaveAs(Form("figures/PullDistribution_%s_LS.eps",suffix.Data()));
      c4pulls->SaveAs(Form("figures/PullDistribution_%s_EM.eps",suffix.Data()));

      c2pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_Rot.eps",suffix.Data()));
      c3pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_LS.eps",suffix.Data()));
      c4pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_EM.eps",suffix.Data()));

      if(tryDirectFit){
	c5residuals->SaveAs(Form("figures/ResidualDistribution_%s_SB.eps",suffix.Data()));
	c5residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_SB.eps",suffix.Data()));
	c5pulls->SaveAs(Form("figures/PullDistribution_%s_SB.eps",suffix.Data()));
	c5pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_SB.eps",suffix.Data()));
      }
    }
  }

  // save also .root
  if(saveCanvasAsRoot){
    c2->SaveAs(Form("figures/InvMassSpectra_%s_Rot.root",suffix.Data()));
    c3->SaveAs(Form("figures/InvMassSpectra_%s_LS.root",suffix.Data()));
    c4->SaveAs(Form("figures/InvMassSpectra_%s_EM.root",suffix.Data()));

    c2residuals->SaveAs(Form("figures/ResidualDistribution_%s_Rot.root",suffix.Data()));
    c3residuals->SaveAs(Form("figures/ResidualDistribution_%s_LS.root",suffix.Data()));
    c4residuals->SaveAs(Form("figures/ResidualDistribution_%s_EM.root",suffix.Data()));

    c2residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_Rot.root",suffix.Data()));
    c3residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_LS.root",suffix.Data()));
    c4residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_EM.root",suffix.Data()));

    c2pulls->SaveAs(Form("figures/PullDistribution_%s_Rot.root",suffix.Data()));
    c3pulls->SaveAs(Form("figures/PullDistribution_%s_LS.root",suffix.Data()));
    c4pulls->SaveAs(Form("figures/PullDistribution_%s_EM.root",suffix.Data()));

    c2pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_Rot.root",suffix.Data()));
    c3pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_LS.root",suffix.Data()));
    c4pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_EM.root",suffix.Data()));

    if(tryDirectFit){
      c5->SaveAs(Form("figures/InvMassSpectra_%s_SB.root",suffix.Data()));
      c5residuals->SaveAs(Form("figures/ResidualDistribution_%s_SB.root",suffix.Data()));
      c5residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_SB.root",suffix.Data()));
      c5pulls->SaveAs(Form("figures/PullDistribution_%s_SB.root",suffix.Data()));
      c5pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_SB.root",suffix.Data()));
    }
  }

  TCanvas* cry=new TCanvas("cry","RawYield",800,700);
  cry->SetLeftMargin(0.15);
  hRawYieldRot->SetMarkerStyle(21);
  hRawYieldRot->Draw("P");
  hRawYieldRot->SetMinimum(0);
  hRawYieldRot->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRawYieldRot->GetYaxis()->SetTitle("Raw Yield");
  hRawYieldRot->GetYaxis()->SetTitleOffset(1.8);
  Double_t max=hRawYieldRot->GetMaximum();
  if(hRawYieldLS->GetMaximum()>max)max=hRawYieldLS->GetMaximum();
  if(hRawYieldME->GetMaximum()>max)max=hRawYieldME->GetMaximum();
  if(tryDirectFit){
    if(hRawYieldSBfit->GetMaximum()>max)max=hRawYieldSBfit->GetMaximum();
  }
  hRawYieldRot->SetMaximum(max*1.2);

  hRawYieldLS->SetMarkerStyle(22);
  hRawYieldLS->SetMarkerColor(kGreen+2);
  hRawYieldLS->SetLineColor(kGreen+2);
  hRawYieldLS->Draw("PZSAME");
  hRawYieldME->SetMarkerStyle(25);
  hRawYieldME->SetMarkerColor(4);
  hRawYieldME->SetLineColor(4);
  hRawYieldME->Draw("PSAME");
  if(tryDirectFit){
    hRawYieldSBfit->SetMarkerStyle(27);
    hRawYieldSBfit->SetMarkerColor(6);
    hRawYieldSBfit->SetLineColor(6);
    hRawYieldSBfit->Draw("PSAME");
  }
  TLegend* legry=new TLegend(0.7,0.7,0.89,0.89);
  legry->SetFillStyle(0);
  legry->SetBorderSize(0);
  legry->AddEntry(hRawYieldRot,"Rotational","PL")->SetTextColor(1);
  legry->AddEntry(hRawYieldLS,"Like Sign","PL")->SetTextColor(kGreen+2);
  legry->AddEntry(hRawYieldME,"Ev Mix","PL")->SetTextColor(4);
  if(tryDirectFit){
    legry->AddEntry(hRawYieldSBfit,"Direct Fit","PL")->SetTextColor(6);
  }
  legry->Draw();
  if(saveCanvasAsEps>0) cry->SaveAs(Form("figures/RawYield_%s.eps",suffix.Data()));
  if(saveCanvasAsRoot) cry->SaveAs(Form("figures/RawYield_%s.root",suffix.Data()));

  TCanvas* cch2=new TCanvas("cch2","Chi2",800,700);
  cch2->SetLeftMargin(0.15);
  hChiSqRot->SetMarkerStyle(21);
  hChiSqRot->Draw("P");
  hChiSqRot->SetMinimum(0);
  hChiSqRot->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hChiSqRot->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hChiSqRot->GetYaxis()->SetTitleOffset(1.8);
  Double_t maxc=hChiSqRot->GetMaximum();
  if(hChiSqLS->GetMaximum()>maxc)maxc=hChiSqLS->GetMaximum();
  if(hChiSqME->GetMaximum()>maxc)maxc=hChiSqME->GetMaximum();
  if(tryDirectFit){
    if(hChiSqSBfit->GetMaximum()>maxc)maxc=hChiSqSBfit->GetMaximum();
  }
  hChiSqRot->SetMaximum(maxc*1.2);

  hChiSqLS->SetMarkerStyle(22);
  hChiSqLS->SetMarkerColor(kGreen+2);
  hChiSqLS->SetLineColor(kGreen+2);
  hChiSqLS->Draw("PZSAME");
  hChiSqME->SetMarkerStyle(25);
  hChiSqME->SetMarkerColor(4);
  hChiSqME->SetLineColor(4);
  hChiSqME->Draw("PSAME");
  if(tryDirectFit){
    hChiSqSBfit->SetMarkerStyle(27);
    hChiSqSBfit->SetMarkerColor(6);
    hChiSqSBfit->SetLineColor(6);
    hChiSqSBfit->Draw("PSAME");
  }
  legry->Draw();
  if(saveCanvasAsEps>0) cch2->SaveAs(Form("figures/ChiSq_%s.eps",suffix.Data()));
  if(saveCanvasAsRoot) cch2->SaveAs(Form("figures/ChiSq_%s.root",suffix.Data()));




  TH1F* hRatioLSToME=(TH1F*)hRawYieldLS->Clone("hRatioLStoME");
  TH1F* hRatioRotToME=(TH1F*)hRawYieldRot->Clone("hRatioRottoME");
  TH1F* hRatioMEToME=(TH1F*)hRawYieldME->Clone("hRatioMEtoME");
  TH1F* hRatioSBToME=(TH1F*)hRawYieldSBfit->Clone("hRatioSBtoME");
  for(Int_t ib=1; ib<=hRawYieldME->GetNbinsX(); ib++){
    Double_t yme=hRawYieldME->GetBinContent(ib);
    if(yme>0.){
      hRatioLSToME->SetBinContent(ib,hRawYieldLS->GetBinContent(ib)/yme);
      hRatioLSToME->SetBinError(ib,hRawYieldLS->GetBinError(ib)/yme);
      hRatioMEToME->SetBinContent(ib,hRawYieldME->GetBinContent(ib)/yme);
      hRatioMEToME->SetBinError(ib,hRawYieldME->GetBinError(ib)/yme);
      hRatioRotToME->SetBinContent(ib,hRawYieldRot->GetBinContent(ib)/yme);
      hRatioRotToME->SetBinError(ib,hRawYieldRot->GetBinError(ib)/yme);
      hRatioSBToME->SetBinContent(ib,hRawYieldSBfit->GetBinContent(ib)/yme);
      hRatioSBToME->SetBinError(ib,hRawYieldSBfit->GetBinError(ib)/yme);
    }
  }

  TCanvas* cry2=new TCanvas("cry2","RawYield+Ratios",1400,700);
  cry2->Divide(2,1);
  cry2->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);    
  hRawYieldRot->Draw("P");
  hRawYieldLS->Draw("PZSAME");
  hRawYieldME->Draw("PSAME");
  if(tryDirectFit){
    hRawYieldSBfit->Draw("PSAME");
  }
  legry->Draw();
  cry2->cd(2);
  hRatioLSToME->SetStats(0);
  hRatioLSToME->SetMinimum(0.3);
  hRatioLSToME->SetMaximum(1.7);
  hRatioLSToME->GetYaxis()->SetTitle("Ratio To EvMix");
  hRatioLSToME->Draw("same");
  hRatioRotToME->Draw("same");
  hRatioMEToME->Draw("same");
  hRatioSBToME->Draw("same");
  if(saveCanvasAsEps>0) cry2->SaveAs(Form("figures/RawYieldAndRatios_%s.eps",suffix.Data()));
  if(saveCanvasAsRoot) cry2->SaveAs(Form("figures/RawYieldAndRatios_%s.root",suffix.Data()));


  cCompareResidualTrends->cd(1);
  TLegend *legRT=new TLegend(*legry);
  legRT->SetX1NDC(0.12);
  legRT->SetY1NDC(0.7);
  legRT->SetX2NDC(0.3);
  legRT->SetY2NDC(0.9);
  legRT->Draw();

  // NOW COMPARE YIELDS OBTAINED WITH BC FROM DIFFERENT APPROACHES
  TCanvas* cryBC=new TCanvas("cryBC","RawYield with BC",800,700);
  cryBC->SetLeftMargin(0.15);
  hRawYieldRot->Draw("P");
  hRawYieldME->Draw("PSAME");
  hRawYieldLS->Draw("PZSAME");
  if(tryDirectFit){
    hRawYieldSBfit->Draw("PSAME");
  }

  if(hRawYieldRotBC->GetMaximum()>max)max=hRawYieldRotBC->GetMaximum();
  if(hRawYieldLSBC->GetMaximum()>max)max=hRawYieldLSBC->GetMaximum();
  if(hRawYieldMEBC->GetMaximum()>max)max=hRawYieldMEBC->GetMaximum();
  if(tryDirectFit){
    if(hRawYieldSBfitBC->GetMaximum()>max)max=hRawYieldSBfitBC->GetMaximum();
  }
  hRawYieldRot->SetMaximum(max*1.2);
  
  hRawYieldRotBC->SetMarkerStyle(21);
  hRawYieldRotBC->SetLineStyle(2);
  hRawYieldRotBC->Draw("PSAME");


  hRawYieldLSBC->SetMarkerStyle(22);
  hRawYieldLSBC->SetMarkerColor(kGreen+2);
  hRawYieldLSBC->SetLineColor(kGreen+2);
  hRawYieldLSBC->SetLineStyle(2);
  hRawYieldLSBC->Draw("PZSAME");

  hRawYieldMEBC->SetMarkerStyle(25);
  hRawYieldMEBC->SetMarkerColor(4);
  hRawYieldMEBC->SetLineColor(4);
  hRawYieldMEBC->SetLineStyle(2);
  hRawYieldMEBC->Draw("PSAME");
  if(tryDirectFit){
    hRawYieldSBfitBC->SetMarkerStyle(27);
    hRawYieldSBfitBC->SetMarkerColor(6);
    hRawYieldSBfitBC->SetLineColor(6);
    hRawYieldSBfitBC->SetLineStyle(2);
    hRawYieldSBfitBC->Draw("PSAME");
  }

  TLegend* legryBC=new TLegend(0.7,0.7,0.89,0.89);
  legryBC->SetFillStyle(0);
  legryBC->SetBorderSize(0);
  legryBC->AddEntry(hRawYieldRot,"Rotational","PL")->SetTextColor(1);
  legryBC->AddEntry(hRawYieldRotBC,"Rotational BC","PL")->SetTextColor(1);
  legryBC->AddEntry(hRawYieldLS,"Like Sign ","PL")->SetTextColor(kGreen+2);
  legryBC->AddEntry(hRawYieldLSBC,"Like Sign BC","PL")->SetTextColor(kGreen+2);
  legryBC->AddEntry(hRawYieldME,"Ev Mix","PL")->SetTextColor(4);
  legryBC->AddEntry(hRawYieldMEBC,"Ev Mix BC","PL")->SetTextColor(4);
  if(tryDirectFit){
    legryBC->AddEntry(hRawYieldSBfit,"Direct Fit","PL")->SetTextColor(6);
    legryBC->AddEntry(hRawYieldSBfitBC,"Direct Fit BC","PL")->SetTextColor(6);
  }
  legryBC->Draw();
  if(saveCanvasAsEps>0) cryBC->SaveAs(Form("figures/RawYieldBC_%s.eps",suffix.Data()));
  if(saveCanvasAsRoot) cryBC->SaveAs(Form("figures/RawYieldBC_%s.root",suffix.Data()));




  fout->cd();
  hRawYieldRot->Write();
  hRawYieldLS->Write();
  hRawYieldME->Write();
  hRawYieldSBfit->Write();
  hChiSqRot->Write();
  hChiSqLS->Write();
  hChiSqME->Write();
  hChiSqSBfit->Write();
  hNdfRot->Write();
  hNdfLS->Write();
  hNdfME->Write();
  hNdfSBfit->Write();
  hnEv->Write();
  fout->Close();

  return;
}


void WriteFitInfo(AliHFMassFitter *fitter, TH1D* histo){
  Double_t sig,esig;
  fitter->Signal(3,sig,esig);
  Double_t mean=fitter->GetMean();
  Double_t emean=fitter->GetMeanUncertainty();
  Double_t sigma=fitter->GetSigma();
  Double_t esigma=fitter->GetSigmaUncertainty();
  Double_t minBin=histo->FindBin(mean-3.*sigma);
  Double_t maxBin=histo->FindBin(mean+3.*sigma);
  Double_t back=histo->Integral(minBin,maxBin);
  TPaveText* tpar=new TPaveText(0.5,0.7,0.89,.87,"NDC");
  tpar->SetBorderSize(0);
  tpar->SetFillStyle(0);
  tpar->AddText(Form("Mean = %.3f #pm %.3f",mean,emean));
  tpar->AddText(Form("Sigma = %.3f #pm %.3f",sigma,esigma));
  tpar->SetTextColor(4);
  tpar->Draw();

  TPaveText* tss=new TPaveText(0.15,0.15,0.5,0.4,"NDC");
  tss->SetBorderSize(0);
  tss->SetFillStyle(0);
  tss->AddText(Form("S(3#sigma) = %.0f #pm %.0f",sig,esig));
  tss->AddText(Form("B(3#sigma) = %.3g",back));
  tss->AddText(Form("S/B (3#sigma) = %.4f",sig/back));
  tss->AddText(Form("Significance(3#sigma) = %.2f",sig/TMath::Sqrt(back+sig)));
  tss->SetTextColor(1);
  tss->Draw();

}



