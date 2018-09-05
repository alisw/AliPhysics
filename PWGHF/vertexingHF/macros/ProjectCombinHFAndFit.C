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
#include <TStyle.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include "AliHFInvMassFitter.h"
#include "AliVertexingHFUtils.h"
#include "AliNormalizationCounter.h"
#endif

enum Method{kME,kRot,kLS,kSB};

// input files and pt binning
TString fileName="DataTrains/AnalysisResults_17pq_FAST_wSDD_train2201.root";
TString suffix="3SigPID_Pt300_FidY_PilMV_EM1";
TString fileNameMC="MCTrains/AnalysisResults_LHC17pq_FAST_CENTwSDD_G3_train890-889.root";
//TString fileNameMC="MCTrains/AnalysisResults_LHC17pq_FAST_CENTwSDD_G4_train892-891.root";
TString suffixMC="_Prompt_3SigPID_Pt300_FidY_PilSPD5_EM1";

// TString fileName="DataTrains/AnalysisResults_17pq_FAST_wSDD_train2104.root";
// TString suffix="2SigPID_Pt300_FidY_PilSPD5_EM1";
// //TString fileNameMC="MCTrains/AnalysisResults_LHC17pq_FAST_CENTwSDD_G4_train870-869.root";
// TString fileNameMC="MCTrains/AnalysisResults_LHC17pq_FAST_CENTwSDD_G3_train868-867.root";
// TString suffixMC="_Prompt_2SigPID_Pt300_FidY_PilSPD5_EM1";

TString meson="Dzero";
const Int_t nPtBins=8;
Double_t binLims[nPtBins+1]={0.,1.,2.,3.,4.,5.,6.,8.,12.};
//Double_t binLims[nPtBins+1]={0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.};
Double_t sigmas[nPtBins]={0.006,0.008,0.009,0.010,0.011,0.012,0.013,0.013};

// outputfiles
Bool_t saveCanvasAsRoot=kTRUE;
Int_t saveCanvasAsEps=1;   //0=none, 1=main ones, 2=all

// fit configuration
Int_t rebin[nPtBins]={4,6,7,8,9,10,10,12};//8,9,10,12,12,12,12};
Int_t fixSigmaConf=0; // 0= all free, 1=all fixed, 2=use values per pt bin
Bool_t fixSigma[nPtBins]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
Double_t tuneSigmaOnData=-1.00; // scaling factor for the Gaussian sigma, if -1. use MC out of the box
Int_t fixMeanConf=0; // 0= all free, 1=all fixed, 2=use values per pt bin
Bool_t fixMean[nPtBins]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
Double_t minMass=1.72;
Double_t maxMass=2.04;
Int_t optForNorm=1;
Double_t rangeForNorm=0.05;
TString fitoption="E";
Bool_t useEMwithLS=kTRUE;
Int_t typeb=2;
Double_t nsigmaBinCounting=4.;      // defines the mass interval over which the signal is bin counted
Int_t optBkgBinCount=1;
Double_t massD;

Int_t smoothLS=0;

// reflection option
TString reflopt="2gaus";
Bool_t correctForRefl=kTRUE;
Double_t rOverSmodif=1;

// objects and options related to side-band fit method
Bool_t tryDirectFit=kTRUE;

Int_t nparback=0;
Double_t fitrangelow[nPtBins]={1.74,1.74,1.74,1.72,1.72,1.72,1.72,1.72};
Double_t fitrangeup[nPtBins]={2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.04};
Int_t nDegreeBackPol[nPtBins]={4,4,4,2,2,2,2,2};             // degree of polynomial function describing the background

TH1D* hMCReflPtBin;
TH1D* hMCSigPtBin;

void WriteFitInfo(AliHFInvMassFitter *fitter, TH1D* histo);
void WriteFitFunctionsToFile(AliHFInvMassFitter *fitter, TString meth, Int_t iPtBin);
TH1F* FitMCInvMassSpectra(TList* lMC);

AliHFInvMassFitter* ConfigureFitter(TH1D* histo, Int_t iPtBin, Int_t backcase, Double_t minFit, Double_t maxFit){
  TH1F* histof=(TH1F*)histo->Clone(Form("%s_Fl",histo->GetName()));


  AliHFInvMassFitter* fitter=new AliHFInvMassFitter(histof,minFit,maxFit,backcase,0);
  if(backcase==6)fitter->SetPolDegreeForBackgroundFit(nDegreeBackPol[iPtBin]);
  fitter->SetUseChi2Fit();
  fitter->SetFitOption(fitoption.Data());
  fitter->SetInitialGaussianMean(massD);
  fitter->SetInitialGaussianSigma(sigmas[iPtBin]);
  if(fixSigmaConf==1 || (fixSigmaConf==2 && fixSigma[iPtBin])) fitter->SetFixGaussianSigma(sigmas[iPtBin]);
  if(fixMeanConf==1 || (fixMeanConf==2 && fixMean[iPtBin])) fitter->SetFixGaussianMean(massD);
  if(correctForRefl){
    TCanvas *cTest=new TCanvas("cTest","cTest",800,800);    
    TH1F *hmasstemp=fitter->GetHistoClone();
    TH1F *hReflModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(hMCReflPtBin,hmasstemp,minFit,maxFit);
    TH1F *hSigModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(hMCSigPtBin,hmasstemp,minFit,maxFit);
    hReflModif->SetLineColor(kRed);
    hSigModif->SetLineColor(kBlue);
    hSigModif->Draw();
    hReflModif->Draw("same");
    cTest->SaveAs(Form("figures/cTest%d.eps",iPtBin));
    delete hmasstemp;
    Double_t fixSoverRefAt=rOverSmodif*(hReflModif->Integral(hReflModif->FindBin(minFit*1.0001),hReflModif->FindBin(maxFit*0.999))/hSigModif->Integral(hSigModif->FindBin(minFit*1.0001),hSigModif->FindBin(maxFit*0.999)));
    TH1F* hrfl=fitter->SetTemplateReflections(hReflModif,reflopt,minFit,maxFit);
    if(!hrfl){
      Printf("SOMETHING WENT WRONG WHILE SETTINGS REFLECTIONS TEMPLATE");
      delete hReflModif;
      delete hSigModif;
      delete cTest; 
      return 0x0;
    }
    fitter->SetFixReflOverS(fixSoverRefAt);
    delete cTest; 
    delete hReflModif;
    delete hSigModif;
  }
  return fitter;
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
  }else{
    c->Divide(1,1);
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
  TString dirNameMC=Form("PWG3_D2H_InvMass%sLowPt%s",meson.Data(),suffixMC.Data());
  TString lstNameMC=Form("coutput%s%s",meson.Data(),suffixMC.Data());

  if(correctForRefl) suffix.Prepend("Refl_");
  if(fileName.Contains("FAST") && !fileName.Contains("wSDD")){
    suffix.Prepend("FAST_");
  }else if(!fileName.Contains("FAST") && fileName.Contains("wSDD")){
    suffix.Prepend("wSDD_");
  }
  if(fileNameMC.Contains("_G3")) suffix.Append("_Geant3MC");
  else if(fileNameMC.Contains("_G4")) suffix.Append("_Geant4MC");
  if(smoothLS!=0) suffix.Append(Form("_smoothLS%d",smoothLS));

  if(meson=="Dplus") massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  else if(meson=="Dzero") massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  Int_t nBinsWithFixSig=0;
  Int_t nBinsWithFixMean=0;
  Int_t patSig=0;
  Int_t patMean=0;
  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    if(fixSigmaConf==1 || (fixSigmaConf==2 && fixSigma[iPtBin])){ 
      nBinsWithFixSig++;
      patSig+=1<<iPtBin;
    }
    if(fixMeanConf==1 || (fixMeanConf==2 && fixMean[iPtBin])){
      nBinsWithFixMean++;
      patMean+=1<<iPtBin;
    }
  }

  TFile* fil=new TFile(fileName.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dirName.Data());
  
  
  AliNormalizationCounter *nC=(AliNormalizationCounter*)df->Get("NormalizationCounter");
  TH1F* hnEv=new TH1F("hEvForNorm","events for normalization",1,0,1);
  printf("Number of Ev. for norm=%d\n",(Int_t)nC->GetNEventsForNorm());
  hnEv->SetBinContent(1,nC->GetNEventsForNorm());

  TH1F* hRebin=new TH1F("hRebin","",nPtBins,binLims);
  TH1F* hBkgFitFunc=new TH1F("hBkgFitFunc","",nPtBins,binLims);
  TH1F* hBkgFitFuncSB=new TH1F("hBkgFitFuncSB","",nPtBins,binLims);

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

  TH3F* h3drefl=0x0;
  TH3F* h3dmcsig=0x0;
  TH1F* hSigmaMC=0x0;
  TFile* filMC=new TFile(fileNameMC.Data());
  if(filMC && filMC->IsOpen()){
    TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get(dirNameMC.Data());
    TList* lMC=(TList*)dfMC->Get(lstNameMC.Data());
    hSigmaMC=FitMCInvMassSpectra(lMC);
    if(nBinsWithFixSig>0 && !hSigmaMC){
      printf("Fit to MC inv. mass spectra failed\n");
      return;
    }
    if(correctForRefl){
      h3drefl=(TH3F*)lMC->FindObject("hMassVsPtVsYRefl");
      h3dmcsig=(TH3F*)lMC->FindObject("hMassVsPtVsYSig");
    }
  }

  TString sigConf="FixedSigma";
  if(nBinsWithFixSig==0) sigConf="FreeSigma";
  else if(nBinsWithFixSig==nPtBins) sigConf="FixedSigmaAll";
  else sigConf=Form("FixedSigma%d",patSig);
  if(nBinsWithFixSig>0 && tuneSigmaOnData>0.) sigConf+=Form("%d",TMath::Nint(tuneSigmaOnData*100.));
  if(nBinsWithFixMean==nPtBins) sigConf+="FixedMeanAll";
  else if(nBinsWithFixMean>0) sigConf+=Form("FixedMean%d",patMean);

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


  AliHFInvMassFitter* fitterRot[nPtBins];
  AliHFInvMassFitter* fitterLS[nPtBins];
  AliHFInvMassFitter* fitterME[nPtBins];
  AliHFInvMassFitter* fitterSB[nPtBins];

  TH1F* hRawYieldRot=new TH1F("hRawYieldRot"," ; p_{T} (GeV/c) ; Raw yield",nPtBins,binLims);
  TH1F* hRawYieldLS=new TH1F("hRawYieldLS"," ; p_{T} (GeV/c) ; Raw yield",nPtBins,binLims);
  TH1F* hRawYieldME=new TH1F("hRawYieldME"," ; p_{T} (GeV/c) ; Raw yield",nPtBins,binLims);
  TH1F* hRawYieldSB=new TH1F("hRawYieldSB"," ; p_{T} (GeV/c) ; Raw yield",nPtBins,binLims);

  TH1F* hRelStatRot=new TH1F("hRelStatRot"," ; p_{T} (GeV/c) ; Relative stat. unc.",nPtBins,binLims);
  TH1F* hRelStatLS=new TH1F("hRelStatLS"," ; p_{T} (GeV/c) ; Relative stat. unc.",nPtBins,binLims);
  TH1F* hRelStatME=new TH1F("hRelStatME"," ; p_{T} (GeV/c) ; Relative stat. unc.",nPtBins,binLims);
  TH1F* hRelStatSB=new TH1F("hRelStatSB"," ; p_{T} (GeV/c) ; Relative stat. unc.",nPtBins,binLims);

  TH1F* hSignifRot=new TH1F("hSignifRot","",nPtBins,binLims);
  TH1F* hSignifLS=new TH1F("hSignifLS","",nPtBins,binLims);
  TH1F* hSignifME=new TH1F("hSignifME","",nPtBins,binLims);
  TH1F* hSignifSB=new TH1F("hSignifSB","",nPtBins,binLims);

  TH1F* hSoverBRot=new TH1F("hSoverBRot","",nPtBins,binLims);
  TH1F* hSoverBLS=new TH1F("hSoverBLS","",nPtBins,binLims);
  TH1F* hSoverBME=new TH1F("hSoverBME","",nPtBins,binLims);
  TH1F* hSoverBSB=new TH1F("hSoverBSB","",nPtBins,binLims);

  TH1F* hGausMeanRot=new TH1F("hGausMeanRot","",nPtBins,binLims);
  TH1F* hGausMeanLS=new TH1F("hGausMeanLS","",nPtBins,binLims);
  TH1F* hGausMeanME=new TH1F("hGausMeanME","",nPtBins,binLims);
  TH1F* hGausMeanSB=new TH1F("hGausMeanSB","",nPtBins,binLims);

  TH1F* hGausSigmaRot=new TH1F("hGausSigmaRot","",nPtBins,binLims);
  TH1F* hGausSigmaLS=new TH1F("hGausSigmaLS","",nPtBins,binLims);
  TH1F* hGausSigmaME=new TH1F("hGausSigmaME","",nPtBins,binLims);
  TH1F* hGausSigmaSB=new TH1F("hGausSigmaSB","",nPtBins,binLims);

  TH1F* hChiSqRot=new TH1F("hChiSqRot","",nPtBins,binLims);
  TH1F* hChiSqLS=new TH1F("hChiSqLS","",nPtBins,binLims);
  TH1F* hChiSqME=new TH1F("hChiSqME","",nPtBins,binLims);
  TH1F* hChiSqSB=new TH1F("hChiSqSB","",nPtBins,binLims);

  TH1F* hNdfRot=new TH1F("hNdfRot","",nPtBins,binLims);
  TH1F* hNdfLS=new TH1F("hNdfLS","",nPtBins,binLims);
  TH1F* hNdfME=new TH1F("hNdfME","",nPtBins,binLims);
  TH1F* hNdfSB=new TH1F("hNdfSB","",nPtBins,binLims);

  TH1F* hRawYieldRotBC=new TH1F("hRawYieldRotBC","BC yield (rotational background)",nPtBins,binLims);
  TH1F* hRawYieldLSBC=new TH1F("hRawYieldLSBC","BC yield (like-sign background)",nPtBins,binLims);
  TH1F* hRawYieldMEBC=new TH1F("hRawYieldMEBC","BC yield (mixed-event background)",nPtBins,binLims);
  TH1F* hRawYieldSBBC=new TH1F("hRawYieldSBBC","BC yield (side-band fit background)",nPtBins,binLims);

  TH1F* hInvMassHistoBinWidthRot=new TH1F("hInvMassHistoBinWidthRot","",nPtBins,binLims);
  TH1F* hInvMassHistoBinWidthLS=new TH1F("hInvMassHistoBinWidthLS","",nPtBins,binLims);
  TH1F* hInvMassHistoBinWidthME=new TH1F("hInvMassHistoBinWidthME","",nPtBins,binLims);
  TH1F* hInvMassHistoBinWidthSB=new TH1F("hInvMassHistoBinWidthSB","",nPtBins,binLims);

  TLatex* tME=new TLatex(0.65,0.82,"MixEv +- pairs");
  tME->SetNDC();
  TLatex* tMEpp=new TLatex(0.65,0.75,"MixEv ++ pairs");
  tMEpp->SetNDC();
  TLatex* tMEmm=new TLatex(0.65,0.68,"MixEv -- pairs");
  tMEmm->SetNDC();

  TF1 *fpeak=new TF1("fpeak","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",minMass,maxMass);

 
  TDirectory *current = gDirectory;
  TFile* fout=new TFile(Form("outputMassFits_%s_%s.root",sigConf.Data(),suffix.Data()),"recreate");
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
      hMCSigPtBin=h3dmcsig->ProjectionX(Form("hMCSigPtBin%d",iPtBin),bin1MC,bin2MC);
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
    printf("Background norm bin %d DONE\n",iPtBin);

    c1->cd(iPtBin+1);
    hMassPtBin->GetXaxis()->SetRangeUser(minMass,maxMass);
    hMassPtBin->SetMinimum(0);
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
    gPad->Update();
    
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
    if(hMCReflPtBin) hMCReflPtBin->Write();
    if(hMCSigPtBin) hMCSigPtBin->Write();
    current->cd();

    hMassSubRot=AliVertexingHFUtils::RebinHisto(hMassSubRot,rebin[iPtBin]);
    hMassSubME=AliVertexingHFUtils::RebinHisto(hMassSubME,rebin[iPtBin]);
    if(hMassPtBinls) hMassSubLS=AliVertexingHFUtils::RebinHisto(hMassSubLS,rebin[iPtBin]);
 
    hRebin->SetBinContent(iPtBin+1,rebin[iPtBin]);
    Int_t bkgToFill=typeb;
    if(typeb==6) bkgToFill=typeb+nDegreeBackPol[iPtBin];
    Int_t bkgToFillSB=6+nDegreeBackPol[iPtBin];

    hBkgFitFunc->SetBinContent(iPtBin+1,bkgToFill);
    hBkgFitFuncSB->SetBinContent(iPtBin+1,bkgToFillSB);

    fitterRot[iPtBin]=ConfigureFitter(hMassSubRot,iPtBin,typeb,minMass,maxMass);
    if(hMassPtBinls) fitterLS[iPtBin]=ConfigureFitter(hMassSubLS,iPtBin,typeb,minMass,maxMass);
    fitterME[iPtBin]=ConfigureFitter(hMassSubME,iPtBin,typeb,minMass,maxMass);

    Bool_t out1=fitterRot[iPtBin]->MassFitter(0);
    Bool_t out2=kFALSE;
    if(hMassPtBinls) out2=fitterLS[iPtBin]->MassFitter(0);
    Bool_t out3=fitterME[iPtBin]->MassFitter(0);
    
    Double_t background=999999999.;
    Bool_t out4=kFALSE;
    if(tryDirectFit){
      TH1D *hMassDirectFit=(TH1D*)hMassPtBin->Clone(Form("hMassDirectFit_bin%d",iPtBin));
      hMassDirectFit=AliVertexingHFUtils::RebinHisto(hMassDirectFit,rebin[iPtBin]);
      fitterSB[iPtBin]=ConfigureFitter(hMassDirectFit,iPtBin,6,fitrangelow[iPtBin],fitrangeup[iPtBin]);
      out4=fitterSB[iPtBin]->MassFitter(0);//DirectFit(hMassDirectFit,iPtBin,hRawYieldSB);

      Double_t background,ebkg;

      if(out4 && fitterSB[iPtBin]->GetMassFunc()){
	c5->cd(iPtBin+1);
	fitterSB[iPtBin]->DrawHere(gPad,3,0);
	hRawYieldSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield());
	hRawYieldSB->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetRawYieldError());
	if(fitterSB[iPtBin]->GetRawYield()>0){
	  hRelStatSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYieldError()/fitterSB[iPtBin]->GetRawYield());
	  hRelStatSB->SetBinError(iPtBin+1,0.00000001);
	}
	fitterSB[iPtBin]->Background(3.,background,ebkg);
	hSignifSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield()/TMath::Sqrt(background+fitterSB[iPtBin]->GetRawYield()));
	hSignifSB->SetBinError(iPtBin+1,0.00000001);
	hSoverBSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield()/background);
	hSoverBSB->SetBinError(iPtBin+1,0.00000001);
	hInvMassHistoBinWidthSB->SetBinContent(iPtBin+1,hMassDirectFit->GetBinWidth(1));
	hGausMeanSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetMean());
	hGausMeanSB->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetMeanUncertainty());
	hGausSigmaSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetSigma());
	hGausSigmaSB->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetSigmaUncertainty());
	hChiSqSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetReducedChiSquare());
	hChiSqSB->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	hNdfSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetMassFunc()->GetNDF());
	hNdfSB->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	//	if(!correctForRefl)
	WriteFitInfo(fitterSB[iPtBin],hMassPtBin);
	fout->cd();
	WriteFitFunctionsToFile(fitterSB[iPtBin],"SB",iPtBin);
	current->cd();

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
	fpeak->SetRange(fitrangelow[iPtBin],fitrangeup[iPtBin]);
	fpeak->SetParameter(0,funcAll->GetParameter(nDegreeBackPol[iPtBin]+1));
	fpeak->SetParameter(1,funcAll->GetParameter(nDegreeBackPol[iPtBin]+2));
	fpeak->SetParameter(2,funcAll->GetParameter(nDegreeBackPol[iPtBin]+3));
	fpeak->DrawCopy("same");

	Double_t errbc;
	Double_t bc=fitterSB[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
	hRawYieldSBBC->SetBinContent(iPtBin+1,bc);
	hRawYieldSBBC->SetBinError(iPtBin+1,errbc);
      
	c5pulls->cd(iPtBin+1);
	TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
	TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
	TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
	TH1F *hResiduals=fitterSB[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);

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
    if(out1 && fitterRot[iPtBin]->GetMassFunc()){
      fitterRot[iPtBin]->DrawHere(gPad,3,0);
      gPad->Update();
      hRawYieldRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetRawYield());
      hRawYieldRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetRawYieldError());
      Double_t minBinBkg=hMassPtBin->FindBin(fitterRot[iPtBin]->GetMean()-3.*fitterRot[iPtBin]->GetSigma());
      Double_t maxBinBkg=hMassPtBin->FindBin(fitterRot[iPtBin]->GetMean()+3.*fitterRot[iPtBin]->GetSigma());
      background=hMassPtBin->Integral(minBinBkg,maxBinBkg);
      hSoverBRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetRawYield()/background);
      hSoverBRot->SetBinError(iPtBin+1,0.000001);
      hRelStatRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetRawYieldError()/fitterRot[iPtBin]->GetRawYield());
      hRelStatRot->SetBinError(iPtBin+1,0.000001);
      hSignifRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetRawYield()/TMath::Sqrt(background+fitterRot[iPtBin]->GetRawYield()));
      hSignifRot->SetBinError(iPtBin+1,0.00000001);
      hInvMassHistoBinWidthRot->SetBinContent(iPtBin+1,hMassSubRot->GetBinWidth(1));
      hGausMeanRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetMean());
      hGausMeanRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetMeanUncertainty());
      hGausSigmaRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetSigma());
      hGausSigmaRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetSigmaUncertainty());
      hChiSqRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetReducedChiSquare());
      hChiSqRot->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      hNdfRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetMassFunc()->GetNDF());
      hNdfRot->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      //      if(!correctForRefl)
      WriteFitInfo(fitterRot[iPtBin],hMassPtBin);
      fout->cd();
      WriteFitFunctionsToFile(fitterRot[iPtBin],"Rot",iPtBin);
      current->cd();
      
      Double_t errbc;
      Double_t bc=fitterRot[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
      hRawYieldRotBC->SetBinContent(iPtBin+1,bc);
      hRawYieldRotBC->SetBinError(iPtBin+1,errbc);

      c2pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterRot[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);
 
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
    if(out2 && fitterLS[iPtBin]->GetMassFunc()){
      fitterLS[iPtBin]->DrawHere(gPad,3,0); 
      hRawYieldLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetRawYield());
      hRawYieldLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetRawYieldError());
      Double_t minBinBkg=hMassPtBin->FindBin(fitterLS[iPtBin]->GetMean()-3.*fitterLS[iPtBin]->GetSigma());
      Double_t maxBinBkg=hMassPtBin->FindBin(fitterLS[iPtBin]->GetMean()+3.*fitterLS[iPtBin]->GetSigma());
      background=hMassPtBin->Integral(minBinBkg,maxBinBkg);
      hSoverBLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetRawYield()/background);
      hSoverBLS->SetBinError(iPtBin+1,0.000001);
      hRelStatLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetRawYieldError()/fitterLS[iPtBin]->GetRawYield());
      hRelStatLS->SetBinError(iPtBin+1,0.0000001);
      hSignifLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetRawYield()/TMath::Sqrt(background+fitterLS[iPtBin]->GetRawYield()));
      hSignifLS->SetBinError(iPtBin+1,0.00000001);
      hInvMassHistoBinWidthLS->SetBinContent(iPtBin+1,hMassSubLS->GetBinWidth(1));
      hGausMeanLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetMean());
      hGausMeanLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetMeanUncertainty());
      hGausSigmaLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetSigma());
      hGausSigmaLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetSigmaUncertainty());
      hChiSqLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetReducedChiSquare());
      hChiSqLS->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      hNdfLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetMassFunc()->GetNDF());
      hNdfLS->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      //      if(!correctForRefl)
      WriteFitInfo(fitterLS[iPtBin],hMassPtBin);
      fout->cd();
      WriteFitFunctionsToFile(fitterLS[iPtBin],"LS",iPtBin);
      current->cd();


      Double_t errbc;
      Double_t bc=fitterLS[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
      hRawYieldLSBC->SetBinContent(iPtBin+1,bc);
      hRawYieldLSBC->SetBinError(iPtBin+1,errbc);

      c3pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterLS[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);
  
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
    if(out3 && fitterME[iPtBin]->GetMassFunc()){ 
      fitterME[iPtBin]->DrawHere(gPad,3,0); 
      hRawYieldME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetRawYield());
      hRawYieldME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetRawYieldError());
      Double_t minBinBkg=hMassPtBin->FindBin(fitterME[iPtBin]->GetMean()-3.*fitterME[iPtBin]->GetSigma());
      Double_t maxBinBkg=hMassPtBin->FindBin(fitterME[iPtBin]->GetMean()+3.*fitterME[iPtBin]->GetSigma());
      background=hMassPtBin->Integral(minBinBkg,maxBinBkg);
      hSoverBME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetRawYield()/background);
      hSoverBME->SetBinError(iPtBin+1,0.000001);
      hRelStatME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetRawYieldError()/fitterME[iPtBin]->GetRawYield());
      hRelStatME->SetBinError(iPtBin+1,0.000001);
      hSignifME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetRawYield()/TMath::Sqrt(background+fitterME[iPtBin]->GetRawYield()));
      hSignifME->SetBinError(iPtBin+1,0.00000001);
      hInvMassHistoBinWidthME->SetBinContent(iPtBin+1,hMassSubME->GetBinWidth(1));
      hGausMeanME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetMean());
      hGausMeanME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetMeanUncertainty());
      hGausSigmaME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetSigma());
      hGausSigmaME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetSigmaUncertainty());
      hChiSqME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetReducedChiSquare());
      hChiSqME->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      hNdfME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetMassFunc()->GetNDF());
      hNdfME->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      //      if(!correctForRefl)
      WriteFitInfo(fitterME[iPtBin],hMassPtBin);
      fout->cd();
      WriteFitFunctionsToFile(fitterME[iPtBin],"ME",iPtBin);
      current->cd();

      Double_t errbc;
      Double_t bc=fitterME[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
      hRawYieldMEBC->SetBinContent(iPtBin+1,bc);
      hRawYieldMEBC->SetBinError(iPtBin+1,errbc);


      c4pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterME[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);
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
    c1->SaveAs(Form("figures/InvMassSpectra_%s_%s_NoBkgSub.eps",sigConf.Data(),suffix.Data()));
    c2->SaveAs(Form("figures/InvMassSpectra_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
    c3->SaveAs(Form("figures/InvMassSpectra_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
    c4->SaveAs(Form("figures/InvMassSpectra_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));
    if(tryDirectFit){
      c5->SaveAs(Form("figures/InvMassSpectra_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
      c5sub->SaveAs(Form("figures/InvMassSpectra_%s_%s_SBsub.eps",sigConf.Data(),suffix.Data()));
    }

    if(saveCanvasAsEps>1){
      c2residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
      c3residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
      c4residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));

      c2residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
      c3residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
      c4residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));

      c2pulls->SaveAs(Form("figures/PullDistribution_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
      c3pulls->SaveAs(Form("figures/PullDistribution_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
      c4pulls->SaveAs(Form("figures/PullDistribution_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));

      c2pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
      c3pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
      c4pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));

      if(tryDirectFit){
	c5residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
	c5residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
	c5pulls->SaveAs(Form("figures/PullDistribution_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
	c5pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
      }
    }
  }

  // save also .root
  if(saveCanvasAsRoot){
    c2->SaveAs(Form("figures/InvMassSpectra_%s_%s_Rot.root",sigConf.Data(),suffix.Data()));
    c3->SaveAs(Form("figures/InvMassSpectra_%s_%s_LS.root",sigConf.Data(),suffix.Data()));
    c4->SaveAs(Form("figures/InvMassSpectra_%s_%s_EM.root",sigConf.Data(),suffix.Data()));

    c2residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_Rot.root",sigConf.Data(),suffix.Data()));
    c3residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_LS.root",sigConf.Data(),suffix.Data()));
    c4residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_EM.root",sigConf.Data(),suffix.Data()));

    c2residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_Rot.root",sigConf.Data(),suffix.Data()));
    c3residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_LS.root",sigConf.Data(),suffix.Data()));
    c4residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_EM.root",sigConf.Data(),suffix.Data()));

    c2pulls->SaveAs(Form("figures/PullDistribution_%s_%s_Rot.root",sigConf.Data(),suffix.Data()));
    c3pulls->SaveAs(Form("figures/PullDistribution_%s_%s_LS.root",sigConf.Data(),suffix.Data()));
    c4pulls->SaveAs(Form("figures/PullDistribution_%s_%s_EM.root",sigConf.Data(),suffix.Data()));

    c2pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_Rot.root",sigConf.Data(),suffix.Data()));
    c3pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_LS.root",sigConf.Data(),suffix.Data()));
    c4pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_EM.root",sigConf.Data(),suffix.Data()));

    if(tryDirectFit){
      c5->SaveAs(Form("figures/InvMassSpectra_%s_%s_SB.root",sigConf.Data(),suffix.Data()));
      c5sub->SaveAs(Form("figures/InvMassSpectra_%s_%s_SBsub.root",sigConf.Data(),suffix.Data()));
      c5residuals->SaveAs(Form("figures/ResidualDistribution_%s_%s_SB.root",sigConf.Data(),suffix.Data()));
      c5residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_%s_SB.root",sigConf.Data(),suffix.Data()));
      c5pulls->SaveAs(Form("figures/PullDistribution_%s_%s_SB.root",sigConf.Data(),suffix.Data()));
      c5pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_%s_SB.root",sigConf.Data(),suffix.Data()));
    }
  }


  hRawYieldRot->SetMarkerStyle(21);
  hRawYieldLS->SetMarkerStyle(22);
  hRawYieldLS->SetMarkerColor(kGreen+2);
  hRawYieldLS->SetLineColor(kGreen+2);
  hRawYieldME->SetMarkerStyle(25);
  hRawYieldME->SetMarkerColor(4);
  hRawYieldME->SetLineColor(4);
  hRawYieldSB->SetMarkerStyle(27);
  hRawYieldSB->SetMarkerColor(6);
  hRawYieldSB->SetLineColor(6);


  TCanvas* cry=new TCanvas("cry","RawYield",800,700);
  cry->SetLeftMargin(0.15);
  hRawYieldRot->Draw("P");
  hRawYieldRot->SetMinimum(0);
  hRawYieldRot->GetYaxis()->SetTitleOffset(1.8);
  Double_t max=hRawYieldRot->GetMaximum();
  if(hRawYieldLS->GetMaximum()>max)max=hRawYieldLS->GetMaximum();
  if(hRawYieldME->GetMaximum()>max)max=hRawYieldME->GetMaximum();
  if(tryDirectFit){
    if(hRawYieldSB->GetMaximum()>max)max=hRawYieldSB->GetMaximum();
  }
  hRawYieldRot->SetMaximum(max*1.2);
  hRawYieldLS->Draw("PZSAME");
  hRawYieldME->Draw("PSAME");
  if(tryDirectFit) hRawYieldSB->Draw("PSAME");
  TLegend* legry=new TLegend(0.7,0.7,0.89,0.89);
  legry->SetFillStyle(0);
  legry->SetBorderSize(0);
  legry->AddEntry(hRawYieldRot,"Rotational","PL")->SetTextColor(1);
  legry->AddEntry(hRawYieldLS,"Like Sign","PL")->SetTextColor(kGreen+2);
  legry->AddEntry(hRawYieldME,"Ev Mix","PL")->SetTextColor(4);
  if(tryDirectFit){
    legry->AddEntry(hRawYieldSB,"Side-Band Fit","PL")->SetTextColor(6);
  }
  legry->Draw();
  if(saveCanvasAsEps>0) cry->SaveAs(Form("figures/RawYield_%s_%s.eps",sigConf.Data(),suffix.Data()));
  if(saveCanvasAsRoot) cry->SaveAs(Form("figures/RawYield_%s_%s.root",sigConf.Data(),suffix.Data()));



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
    if(hChiSqSB->GetMaximum()>maxc)maxc=hChiSqSB->GetMaximum();
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
    hChiSqSB->SetMarkerStyle(27);
    hChiSqSB->SetMarkerColor(6);
    hChiSqSB->SetLineColor(6);
    hChiSqSB->Draw("PSAME");
  }
  legry->Draw();
  if(saveCanvasAsEps>0) cch2->SaveAs(Form("figures/ChiSq_%s_%s.eps",sigConf.Data(),suffix.Data()));
  if(saveCanvasAsRoot) cch2->SaveAs(Form("figures/ChiSq_%s_%s.root",sigConf.Data(),suffix.Data()));

  TH1F* hRatioLSToME=(TH1F*)hRawYieldLS->Clone("hRatioLStoME");
  TH1F* hRatioRotToME=(TH1F*)hRawYieldRot->Clone("hRatioRottoME");
  TH1F* hRatioMEToME=(TH1F*)hRawYieldME->Clone("hRatioMEtoME");
  TH1F* hRatioSBToME=(TH1F*)hRawYieldSB->Clone("hRatioSBtoME");
  for(Int_t ib=1; ib<=hRawYieldME->GetNbinsX(); ib++){
    Double_t yme=hRawYieldME->GetBinContent(ib);
    if(yme>0.){
      hRatioLSToME->SetBinContent(ib,hRawYieldLS->GetBinContent(ib)/yme);
      hRatioLSToME->SetBinError(ib,hRawYieldLS->GetBinError(ib)/yme);
      hRatioMEToME->SetBinContent(ib,hRawYieldME->GetBinContent(ib)/yme);
      hRatioMEToME->SetBinError(ib,hRawYieldME->GetBinError(ib)/yme);
      hRatioRotToME->SetBinContent(ib,hRawYieldRot->GetBinContent(ib)/yme);
      hRatioRotToME->SetBinError(ib,hRawYieldRot->GetBinError(ib)/yme);
      hRatioSBToME->SetBinContent(ib,hRawYieldSB->GetBinContent(ib)/yme);
      hRatioSBToME->SetBinError(ib,hRawYieldSB->GetBinError(ib)/yme);
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
    hRawYieldSB->Draw("PSAME");
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
  if(saveCanvasAsEps>0) cry2->SaveAs(Form("figures/RawYieldAndRatios_%s_%s.eps",sigConf.Data(),suffix.Data()));
  if(saveCanvasAsRoot) cry2->SaveAs(Form("figures/RawYieldAndRatios_%s_%s.root",sigConf.Data(),suffix.Data()));

  hRelStatRot->GetYaxis()->SetTitleOffset(1.8);
  hRelStatRot->SetMarkerStyle(21);
  hRelStatLS->SetMarkerStyle(22);
  hRelStatLS->SetMarkerColor(kGreen+2);
  hRelStatLS->SetLineColor(kGreen+2);
  hRelStatME->SetMarkerStyle(25);
  hRelStatME->SetMarkerColor(4);
  hRelStatME->SetLineColor(4);
  hRelStatSB->SetMarkerStyle(27);
  hRelStatSB->SetMarkerColor(6);
  hRelStatSB->SetLineColor(6);

  TCanvas* cry3=new TCanvas("cry3","RawYield+Ratios+Unc",1800,600);
  cry3->Divide(3,1);
  cry3->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);    
  hRawYieldRot->Draw("P");
  hRawYieldLS->Draw("PZSAME");
  hRawYieldME->Draw("PSAME");
  if(tryDirectFit){
    hRawYieldSB->Draw("PSAME");
  }
  legry->Draw();
  cry3->cd(2);
  hRatioLSToME->SetStats(0);
  hRatioLSToME->SetMinimum(0.3);
  hRatioLSToME->SetMaximum(1.7);
  hRatioLSToME->GetYaxis()->SetTitle("Ratio To EvMix");
  hRatioLSToME->Draw("same");
  hRatioRotToME->Draw("same");
  hRatioMEToME->Draw("same");
  hRatioSBToME->Draw("same");
  cry3->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);    
  hRelStatRot->SetStats(0);
  hRelStatRot->SetMinimum(0.04);
  hRelStatRot->SetMaximum(0.2);
  hRelStatRot->Draw();
  hRelStatLS->Draw("same");
  hRelStatME->Draw("same");
  if(tryDirectFit) hRelStatSB->Draw("same");

  if(saveCanvasAsEps>0) cry3->SaveAs(Form("figures/RawYieldRatiosAndUnc_%s_%s.eps",sigConf.Data(),suffix.Data()));
  if(saveCanvasAsRoot) cry3->SaveAs(Form("figures/RawYieldRatiosAndUnc_%s_%s.root",sigConf.Data(),suffix.Data()));




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
    hRawYieldSB->Draw("PSAME");
  }

  if(hRawYieldRotBC->GetMaximum()>max)max=hRawYieldRotBC->GetMaximum();
  if(hRawYieldLSBC->GetMaximum()>max)max=hRawYieldLSBC->GetMaximum();
  if(hRawYieldMEBC->GetMaximum()>max)max=hRawYieldMEBC->GetMaximum();
  if(tryDirectFit){
    if(hRawYieldSBBC->GetMaximum()>max)max=hRawYieldSBBC->GetMaximum();
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
    hRawYieldSBBC->SetMarkerStyle(27);
    hRawYieldSBBC->SetMarkerColor(6);
    hRawYieldSBBC->SetLineColor(6);
    hRawYieldSBBC->SetLineStyle(2);
    hRawYieldSBBC->Draw("PSAME");
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
    legryBC->AddEntry(hRawYieldSB,"Side-Band Fit","PL")->SetTextColor(6);
    legryBC->AddEntry(hRawYieldSBBC,"Side-Band Fit BC","PL")->SetTextColor(6);
  }
  legryBC->Draw();
  if(saveCanvasAsEps>0) cryBC->SaveAs(Form("figures/RawYieldBC_%s_%s.eps",sigConf.Data(),suffix.Data()));
  if(saveCanvasAsRoot) cryBC->SaveAs(Form("figures/RawYieldBC_%s_%s.root",sigConf.Data(),suffix.Data()));




  fout->cd();
  hRawYieldRot->Write();
  hRawYieldLS->Write();
  hRawYieldME->Write();
  hRawYieldSB->Write();
  hRelStatRot->Write();
  hRelStatLS->Write();
  hRelStatME->Write();
  hRelStatSB->Write();
  hSignifRot->Write();
  hSignifLS->Write();
  hSignifME->Write();
  hSignifSB->Write();
  hSoverBRot->Write();
  hSoverBLS->Write();
  hSoverBME->Write();
  hSoverBSB->Write();
  hGausMeanRot->Write();
  hGausMeanLS->Write();
  hGausMeanME->Write();
  hGausMeanSB->Write();
  hGausSigmaRot->Write();
  hGausSigmaLS->Write();
  hGausSigmaME->Write();
  hGausSigmaSB->Write();
  hChiSqRot->Write();
  hChiSqLS->Write();
  hChiSqME->Write();
  hChiSqSB->Write();
  hNdfRot->Write();
  hNdfLS->Write();
  hNdfME->Write();
  hNdfSB->Write();
  hInvMassHistoBinWidthRot->Write();
  hInvMassHistoBinWidthLS->Write();
  hInvMassHistoBinWidthME->Write();
  hInvMassHistoBinWidthSB->Write();
  hRebin->Write();
  hBkgFitFunc->Write();
  hBkgFitFuncSB->Write();
  hnEv->Write();
  if(hSigmaMC) hSigmaMC->Write();
  fout->Close();

  return;
}

void WriteFitFunctionsToFile(AliHFInvMassFitter *fitter, TString meth, Int_t iPtBin){
  TString nameh;
  TF1* fTot=fitter->GetMassFunc();
  nameh.Form("FitFuncTot_%s_PtBin%d",meth.Data(),iPtBin);
  fTot->SetRange(1.6,2.2);
  fTot->SetNpx(500);
  fTot->Write(nameh.Data());
  TF1* fSig=fitter->GetSignalFunc();
  nameh.Form("FitFuncSig_%s_PtBin%d",meth.Data(),iPtBin);
  fSig->SetRange(1.6,2.2);
  fSig->SetNpx(500);
  fSig->Write(nameh.Data());
  TF1* fBkg=fitter->GetBackgroundRecalcFunc();
  nameh.Form("FitFuncBkg_%s_PtBin%d",meth.Data(),iPtBin);
  fBkg->SetRange(1.6,2.2);
  fBkg->SetNpx(500);
  fBkg->Write(nameh.Data());  
  if(meson=="Dzero"){
    TF1* fBkgR=fitter->GetBkgPlusReflFunc();
    nameh.Form("FitFuncBkgRefl_%s_PtBin%d",meth.Data(),iPtBin);
    fBkgR->SetRange(1.6,2.2);
    fBkgR->SetNpx(500);
    fBkgR->Write(nameh.Data());  
  }
  return;
}

void WriteFitInfo(AliHFInvMassFitter *fitter, TH1D* histo){
  Double_t sig=fitter->GetRawYield();
  Double_t esig=fitter->GetRawYieldError();
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
  tss->AddText(Form("S = %.0f #pm %.0f",sig,esig));
  tss->AddText(Form("B(3#sigma) = %.3g",back));
  tss->AddText(Form("S/B (3#sigma) = %.4f",sig/back));
  if(correctForRefl) tss->AddText(Form("Refl/Sig =  %.3f #pm %.3f ",fitter->GetReflOverSig(),fitter->GetReflOverSigUncertainty()));
  tss->AddText(Form("Significance(3#sigma) = %.2f",sig/TMath::Sqrt(back+sig)));
  tss->SetTextColor(1);
  tss->Draw();

}

TH1F* FitMCInvMassSpectra(TList* lMC){
  TH3F* h3d=(TH3F*)lMC->FindObject("hMassVsPtVsYSig");
  if(!h3d){
    printf("hMassVsPtVsYSig not found\n");
    return 0x0;
  }
  TH3F* h3dr=(TH3F*)lMC->FindObject("hMassVsPtVsYRefl");
  if(!h3dr){
    printf("hMassVsPtVsYRefl not found\n");
  }
  TH1F* hSigmaMC=new TH1F("hSigmaMC","",nPtBins,binLims);

  TCanvas* cmc1=new TCanvas("InvMassMC","InvMassMC",1200,800);
  cmc1->Divide(4,2);

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    Int_t bin1=h3d->GetYaxis()->FindBin(binLims[iPtBin]);
    Int_t bin2=h3d->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
    Double_t ptmed=(binLims[iPtBin]+binLims[iPtBin+1])*0.5;
    Double_t pthalfwid=(binLims[iPtBin+1]-binLims[iPtBin])*0.5;
    printf("Bin %d   Pt range=%f %f\n",iPtBin,h3d->GetYaxis()->GetBinLowEdge(bin1),h3d->GetYaxis()->GetBinUpEdge(bin2));
    TH1D* hMassMCPtBin=h3d->ProjectionX(Form("hMassMCPtBin%d",iPtBin),bin1,bin2);
    hMassMCPtBin->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    hMassMCPtBin->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hMassMCPtBin->GetXaxis()->SetRangeUser(1.65,2.06);
    cmc1->cd(iPtBin+1);
    gPad->SetLogy();
    hMassMCPtBin->Fit("gaus","");
    TF1* fg=(TF1*)hMassMCPtBin->GetListOfFunctions()->FindObject("gaus");
    TLatex* tsig=new TLatex(0.65,0.82,"Signal");
    tsig->SetNDC();
    tsig->SetTextColor(kBlue+1);
    tsig->Draw();

    TPaveText* t1=new TPaveText(0.15,0.55,0.4,0.88,"ndc");
    t1->SetFillStyle(0);
    t1->SetBorderSize(0);
    t1->AddText(Form("#mu = %.4f GeV/c^{2}",fg->GetParameter(1)));
    t1->AddText(Form("#sigma = %.4f GeV/c^{2}",fg->GetParameter(2)));
    t1->AddText(Form("Integral  = %.0f",fg->Integral(1.7,2.1)/hMassMCPtBin->GetBinWidth(1)));
    t1->AddText(Form("Entries  = %.0f",hMassMCPtBin->GetEntries()));
 
    if(h3dr){
      TH1D* hReflPtBin=h3dr->ProjectionX(Form("hReflPtBin%d",iPtBin),bin1,bin2);
      hReflPtBin->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
      hReflPtBin->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
      hReflPtBin->SetLineColor(2);      
      Int_t bin1=hReflPtBin->FindBin(fg->GetParameter(1)-3.*fg->GetParameter(2));
      Int_t bin2=hReflPtBin->FindBin(fg->GetParameter(1)+3.*fg->GetParameter(2));
      TLatex* tref=new TLatex(0.65,0.75,"Reflections");
      tref->SetNDC();
      tref->SetTextColor(2);
      tref->Draw();
      t1->AddText(Form("Refl(3#sigma) = %.0f",hReflPtBin->Integral(bin1,bin2)));
      hReflPtBin->Draw("same");
      t1->Draw();
    }
    sigmas[iPtBin]=fg->GetParameter(2);
    Double_t errMCsigma=fg->GetParError(2);
    if(tuneSigmaOnData>0.){
      sigmas[iPtBin]*=tuneSigmaOnData;
      errMCsigma*=tuneSigmaOnData;
    }
    hSigmaMC->SetBinContent(iPtBin+1,sigmas[iPtBin]);
    hSigmaMC->SetBinError(iPtBin+1,errMCsigma);
  }
  return hSigmaMC;
}
