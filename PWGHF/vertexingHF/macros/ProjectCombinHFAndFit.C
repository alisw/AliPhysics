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

enum Method{kME,kRot,kLikeS,kSB};

TString configFileName="configfile4lowptanalysis.txt";

// input files and pt binning
TString fileName="";
TString suffix="";
TString fileNameMC="";
TString suffixMC="";

TString meson="Dzero";
const Int_t maxPtBins=30;
Int_t nPtBins=8;
Double_t binLims[maxPtBins+1]={0.,1.,2.,3.,4.,5.,6.,8.,12.};

// fit configuration
Int_t rebin[maxPtBins]={4,6,7,8,9,10,10,12};
Double_t minMass4Fit[maxPtBins]={1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72};
Double_t maxMass4Fit[maxPtBins]={2.04,2.04,2.04,2.04,2.04,2.04,2.04,2.04};
Int_t fixSigmaConf=1; // 0= all free, 1=all fixed, 2=use values per pt bin
Bool_t fixSigma[maxPtBins]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
Double_t tuneSigmaOnData=-1.; // scaling factor for the Gaussian sigma, if -1. use MC out of the box
Bool_t useSigmaManual=kFALSE;
Double_t sigmas[maxPtBins]={0.006,0.008,0.009,0.010,0.011,0.012,0.013,0.013};
Int_t fixMeanConf=0; // 0= all free, 1=all fixed, 2=use values per pt bin
Bool_t fixMean[maxPtBins]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
Double_t tuneMeanOnData=1.868; // if -1. use PDG value
// objects and options related to side-band fit method
Bool_t tryDirectFit=kTRUE;
Double_t fitSBrangelow[maxPtBins]={1.74,1.74,1.74,1.72,1.72,1.72,1.72,1.72};
Double_t fitSBrangeup[maxPtBins]={2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.04};
Int_t nDegreeBackPolSB[maxPtBins]={4,4,4,2,2,2,2,2};             // degree of polynomial function describing the background

Int_t optForNorm=1;
Double_t rangeForNorm=0.05;
Bool_t useEMwithLS=kTRUE;
Bool_t useGeomMeanLS=kTRUE;
Bool_t renormLS=kTRUE;
Int_t smoothLS=0;

TString fitoption="E";
Int_t typeb=2;
Int_t nDegreeBackPol[maxPtBins]={2,2,2,2,2,2,2,2};             // degree of polynomial function describing the background
// reflection option
Bool_t correctForRefl=kTRUE;
TString reflopt="2gaus";
Double_t rOverSmodif=1;

Int_t optBkgBinCount=1;
Double_t nsigmaBinCounting=4.;      // defines the mass interval over which the signal is bin counted
Double_t costhstcut=1.1;

// outputfiles
Bool_t saveCanvasAsRoot=kTRUE;
Int_t saveCanvasAsEps=1;   //0=none, 1=main ones, 2=all


TH1D* hMCReflPtBin;
TH1D* hMCSigPtBin;
Double_t minMass=1.72;
Double_t maxMass=2.04;
Double_t massD;


void WriteFitInfo(AliHFInvMassFitter *fitter, TH1D* histo);
void WriteFitFunctionsToFile(AliHFInvMassFitter *fitter, TString meth, Int_t iPtBin);
TH1F* FitMCInvMassSpectra(TList* lMC, TString var);
Bool_t ReadConfig(TString configName);
void PrintConfig();
void CheckMCLineShapes();

AliHFInvMassFitter* ConfigureFitter(TH1D* histo, Int_t iPtBin, Int_t backcase, Double_t minFit, Double_t maxFit, TCanvas* crf, Bool_t isDirect=kFALSE){
  TH1F* histof=(TH1F*)histo->Clone(Form("%s_Fl",histo->GetName()));


  AliHFInvMassFitter* fitter=new AliHFInvMassFitter(histof,minFit,maxFit,backcase,0);
  if(backcase==6){
    fitter->SetPolDegreeForBackgroundFit(nDegreeBackPol[iPtBin]);
    if(isDirect) fitter->SetPolDegreeForBackgroundFit(nDegreeBackPolSB[iPtBin]);
  }
  fitter->SetUseChi2Fit();
  fitter->SetFitOption(fitoption.Data());
  fitter->SetInitialGaussianMean(massD);
  fitter->SetInitialGaussianSigma(sigmas[iPtBin]);
  if(fixSigmaConf==1 || (fixSigmaConf==2 && fixSigma[iPtBin])) fitter->SetFixGaussianSigma(sigmas[iPtBin]);
  if(fixMeanConf==1 || (fixMeanConf==2 && fixMean[iPtBin])) fitter->SetFixGaussianMean(massD);
  if(correctForRefl){
    TCanvas *cTest=new TCanvas("cTest","cTest",1600,800);
    cTest->Divide(2,1);
    cTest->cd(1);
    TH1F *hmasstemp=fitter->GetHistoClone();
    TH1F *hReflModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(hMCReflPtBin,hmasstemp,minFit,maxFit);
    TH1F *hSigModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(hMCSigPtBin,hmasstemp,minFit,maxFit);
    hReflModif->SetLineColor(kRed);
    hSigModif->SetLineColor(kBlue);
    hSigModif->Draw();
    hMCSigPtBin->SetLineColor(kRed-9);
    hMCSigPtBin->Draw("same");
    hMCReflPtBin->SetLineColor(kGray+1);
    hMCReflPtBin->Draw("same");
    hReflModif->Draw("same");
    delete hmasstemp;
    Double_t fixSoverRefAt=rOverSmodif*(hReflModif->Integral(hReflModif->FindBin(minFit*1.0001),hReflModif->FindBin(maxFit*0.999))/hSigModif->Integral(hSigModif->FindBin(minFit*1.0001),hSigModif->FindBin(maxFit*0.999)));
    TH1F* hrfl=fitter->SetTemplateReflections(hReflModif,reflopt,-1.,-1.);
    cTest->cd(2);
    hReflModif->Draw();
    hrfl->SetLineColor(kBlue+1);
    hrfl->Draw("same");
    cTest->SaveAs(Form("figures/ReflectionConfig_PtBin%d.eps",iPtBin));
    if(crf){
      crf->cd(iPtBin+1);
      hReflModif->DrawCopy();
      hrfl->DrawCopy("same");
    }
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

  if(method==kLikeS){
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



Double_t GetBackgroundNormalizationFactor(TH1D* hRatio, Int_t reb=1){
  //
  Double_t norm=hRatio->GetMaximum();
  Double_t massForNorm=0;

  if(optForNorm==0){
    Int_t binl=hRatio->GetXaxis()->FindBin(minMass);
    Int_t binh=hRatio->GetXaxis()->FindBin(maxMass);
    Double_t norml=hRatio->GetBinContent(1);
    if(binl>0 && binl<=hRatio->GetNbinsX()) norml=hRatio->GetBinContent(binl);
    Double_t normh=hRatio->GetBinContent(hRatio->GetNbinsX());
    if(binh>0 && binh<=hRatio->GetNbinsX()) normh=hRatio->GetBinContent(binh);
    if(norml>normh){
      norm=norml;
      massForNorm=hRatio->GetBinCenter(binl);
    }else{
      norm=normh;
      massForNorm=hRatio->GetBinCenter(binh);
    }
  }else if(optForNorm==1){
    norm=0.0001;
    for(Int_t iMassBin=1; iMassBin<hRatio->GetNbinsX(); iMassBin++){
      Double_t bce=hRatio->GetBinCenter(iMassBin);
      if(bce>minMass && bce<maxMass){
	Double_t bco=hRatio->GetBinContent(iMassBin);
	if(bco>norm){ 
	  norm=bco;
	  massForNorm=bce;
	}
      }
    }
  }else if(optForNorm==2){ 
    norm=0.0001;
    for(Int_t iMassBin=1; iMassBin<hRatio->GetNbinsX(); iMassBin++){
      Double_t bce=hRatio->GetBinCenter(iMassBin);
      Double_t bco=hRatio->GetBinContent(iMassBin);
      if(bco>norm){
	norm=bco;
	massForNorm=bce;
      }
    }
  }else if(optForNorm==3){ 
    Int_t binl=hRatio->GetXaxis()->FindBin(minMass);
    Int_t binh=hRatio->GetXaxis()->FindBin(maxMass);
    Double_t norml=0;
    Double_t normReb=0;
    Int_t iFirstBin=TMath::Max(binl-reb/2,1);
    Int_t iLastBin=iFirstBin+reb;
    if(iLastBin>hRatio->GetNbinsX()){
      iLastBin=hRatio->GetNbinsX();
      iFirstBin=iLastBin-reb;
    }
    for(Int_t jjj=iFirstBin; jjj<=iLastBin; jjj++){
      norml+=hRatio->GetBinContent(jjj);
      normReb+=1;
    }
    if(normReb>=1) norml/=normReb;
    Double_t normh=0;
    normReb=0;
    iFirstBin=TMath::Max(binh-reb/2,1);
    iLastBin=iFirstBin+reb;
    if(iLastBin>hRatio->GetNbinsX()){
      iLastBin=hRatio->GetNbinsX();
      iFirstBin=iLastBin-reb;
    }
    for(Int_t jjj=iFirstBin; jjj<=iLastBin; jjj++){
      normh+=hRatio->GetBinContent(jjj);
      normReb+=1;
    }
    if(normReb>=1) normh/=normReb;
    if(norml>normh){
      norm=norml;
      massForNorm=hRatio->GetBinCenter(binl);
    }else{
      norm=normh;
      massForNorm=hRatio->GetBinCenter(binh);
    }
  }else if(optForNorm==4){ 
    norm=0.0001;
    for(Int_t iMassBin=1; iMassBin<hRatio->GetNbinsX(); iMassBin++){
      Double_t bce=hRatio->GetBinCenter(iMassBin);
      if(bce>minMass && bce<maxMass){
	Double_t bco=0;
	Double_t normReb=0;
	Int_t iFirstBin=TMath::Max(iMassBin-reb/2,1);
	Int_t iLastBin=iFirstBin+reb;
	if(iLastBin>hRatio->GetNbinsX()){
	  iLastBin=hRatio->GetNbinsX();
	  iFirstBin=iLastBin-reb;
	}
	for(Int_t jjj=iFirstBin; jjj<=iLastBin; jjj++){
	  bco+=hRatio->GetBinContent(jjj);
	  normReb+=1;
	}
	if(normReb>=1){
	  bco/=normReb;
	  if(bco>norm){
	    norm=bco;
	    massForNorm=bce;
	  }
	}
      }
    }
  }else if(optForNorm==5){ 
    hRatio->Fit("pol0","","",minMass,minMass+rangeForNorm);
    TF1* func0=(TF1*)hRatio->GetListOfFunctions()->FindObject("pol0");
    Double_t norml=func0->GetParameter(0);
    hRatio->Fit("pol0","","",maxMass-rangeForNorm,maxMass);
    func0=(TF1*)hRatio->GetListOfFunctions()->FindObject("pol0");
    Double_t normh=func0->GetParameter(0);
    if(norml>normh){
      norm=norml;
      massForNorm=maxMass-0.5*rangeForNorm;
    }else{
      norm=normh;
      massForNorm=minMass+0.5*rangeForNorm;
    }
  }

  printf("Normalization factor = %f --> taken at mass=%f\n",norm,massForNorm);

  return norm;
}

void CheckMCLineShapes(){
  
  if(configFileName.Length()>0){
    if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",configFileName.Data()))==0){
      printf("Read configuration from file %s\n",configFileName.Data());
      Bool_t readOK=ReadConfig(configFileName);
      if(!readOK){
	printf("Error in reading configuration file\n");
	return;
      }
    }
  }
  printf("***************************************************\n");
  printf("*** This is the configuration that will be used ***\n");
  PrintConfig();
  printf("***                                             ***\n");
  printf("***************************************************\n");
  
  TString dirNameMC=Form("PWG3_D2H_InvMass%sLowPt%s",meson.Data(),suffixMC.Data());
  TString lstNameMC=Form("coutput%s%s",meson.Data(),suffixMC.Data());
  TFile* filMC=new TFile(fileNameMC.Data());
  if(filMC && filMC->IsOpen()){
    TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get(dirNameMC.Data());
    if(!dfMC){
      printf("TDirectoryFile %s not found in TFile for MC\n",dirNameMC.Data());
      filMC->ls();
      return;
    }
    TList* lMC=(TList*)dfMC->Get(lstNameMC.Data());
    TH1F* hSigmaMC=FitMCInvMassSpectra(lMC,"Y");
    TCanvas* csigma=new TCanvas("csigma","",800,700);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);
    hSigmaMC->SetMarkerStyle(24);
    hSigmaMC->SetLineWidth(2);
    hSigmaMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hSigmaMC->GetYaxis()->SetTitle("Gaussian #sigma (GeV/c^{2})");
    hSigmaMC->GetYaxis()->SetTitleOffset(1.85);
    hSigmaMC->GetXaxis()->SetTitleOffset(1.2);
    hSigmaMC->Draw();
  }
}


void ProjectCombinHFAndFit(TString configInput=""){

  if(configInput!="") configFileName=configInput.Data();

  if(configFileName.Length()>0){
    if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",configFileName.Data()))==0){
      printf("Read configuration from file %s\n",configFileName.Data());
      Bool_t readOK=ReadConfig(configFileName);
      if(!readOK){
	printf("Error in reading configuration file\n");
	return;
      }
    }
  }
  printf("***************************************************\n");
  printf("*** This is the configuration that will be used ***\n");
  PrintConfig();
  printf("***                                             ***\n");
  printf("***************************************************\n");
  
  TString dirName=Form("PWG3_D2H_InvMass%sLowPt%s",meson.Data(),suffix.Data());
  TString lstName=Form("coutput%s%s",meson.Data(),suffix.Data());
  TString dirNameMC=Form("PWG3_D2H_InvMass%sLowPt%s",meson.Data(),suffixMC.Data());
  TString lstNameMC=Form("coutput%s%s",meson.Data(),suffixMC.Data());

  if(correctForRefl){
    if(reflopt=="template") suffix.Prepend("TemplRefl_");
    else suffix.Prepend("Refl_");
    if(TMath::Abs(rOverSmodif-1)>0.01) suffix.ReplaceAll("Refl_",Form("Refl%02d_",(Int_t)(rOverSmodif*10)));
  }
  if(fileName.Contains("FAST") && !fileName.Contains("wSDD")){
    suffix.Prepend("FAST_");
  }else if(!fileName.Contains("FAST") && fileName.Contains("wSDD")){
    suffix.Prepend("wSDD_");
  }
  if(fileNameMC.Contains("_G3")) suffix.Append("_Geant3MC");
  else if(fileNameMC.Contains("_G4")) suffix.Append("_Geant4MC");
  if(smoothLS!=0) suffix.Append(Form("_smoothLS%d",smoothLS));
  if(costhstcut<1.) suffix.Append(Form("_CosthSt%d",TMath::Nint(costhstcut*100)));
  if(configFileName.Contains("coarse")) suffix.Append("_CoarsePt");
  if(useEMwithLS) suffix.Append("_EMwithLS");
  
  if(tuneMeanOnData<0){
    if(meson=="Dplus") massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
    else if(meson=="Dzero") massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }else{
    massD=tuneMeanOnData;
  }

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

  if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",fileName.Data())) !=0){
    printf("File %s with raw data results does not exist -> exiting\n",fileName.Data());
    return;
  }
  TFile* fil=new TFile(fileName.Data());
  if(!fil){
    printf("File %s with raw data not opened -> exiting\n",fileName.Data());
    return;
  }
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dirName.Data());
  if(!df){
    printf("TDirectoryFile %s not found in TFile\n",dirName.Data());
    fil->ls();
    return;
  }
  
  AliNormalizationCounter *nC=(AliNormalizationCounter*)df->Get("NormalizationCounter");
  if(!nC){
    printf("AliNormalizationCounter object missing -> exiting\n");
    return;
  }
  TH1F* hnEv=new TH1F("hEvForNorm","events for normalization",1,0,1);
  printf("Number of Ev. for norm=%d\n",(Int_t)nC->GetNEventsForNorm());
  hnEv->SetBinContent(1,nC->GetNEventsForNorm());

  TH1F* hRebin=new TH1F("hRebin","",nPtBins,binLims);
  TH1F* hBkgFitFunc=new TH1F("hBkgFitFunc","",nPtBins,binLims);
  TH1F* hBkgFitFuncSB=new TH1F("hBkgFitFuncSB","",nPtBins,binLims);

  TList* l=(TList*)df->Get(lstName.Data());

  TString var="Y";
  if(costhstcut<1.) var="CosthSt";

  TH3F* h3d=(TH3F*)l->FindObject(Form("hMassVsPtVs%s",var.Data()));
  if(!h3d){
    printf("Histogram hMassVsPtVsCosthSt does not exist->  use Y\n");
    var="Y";
    h3d=(TH3F*)l->FindObject(Form("hMassVsPtVs%s",var.Data()));
    if(!h3d){
      printf("Histogram hMassVsPtVsY does not exist->  exit\n");
    }
  }
  Int_t zbin1=0;
  Int_t zbin2=h3d->GetZaxis()->GetNbins()+1;
  if(var=="CosthSt") zbin2=h3d->GetZaxis()->FindBin(costhstcut-0.000001);
  
  printf("Binning in %s: %d %d\n",var.Data(),zbin1,zbin2);

  TH3F* h3dr=(TH3F*)l->FindObject(Form("hMassVsPtVs%sRot",var.Data()));
  TH3F* h3dme=(TH3F*)l->FindObject(Form("hMassVsPtVs%sME",var.Data()));
  TH3F* h3dmepp=(TH3F*)l->FindObject(Form("hMassVsPtVs%sMELSpp",var.Data()));
  TH3F* h3dmemm=(TH3F*)l->FindObject(Form("hMassVsPtVs%sMELSmm",var.Data()));
  TH2F* hpoolMix=(TH2F*)l->FindObject("hMixingsPerPool");
  TH2F* hpoolEv=(TH2F*)l->FindObject("hEventsPerPool");
  TH3F* h3dlsp=(TH3F*)l->FindObject(Form("hMassVsPtVs%sLSpp",var.Data()));
  TH3F* h3dlsm=(TH3F*)l->FindObject(Form("hMassVsPtVs%sLSmm",var.Data()));

  TH3F* h3drefl=0x0;
  TH3F* h3dmcsig=0x0;
  TH1F* hSigmaMC=0x0;
  if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",fileNameMC.Data())) !=0){
    printf("File %s with MC results does not exist -> exiting\n",fileNameMC.Data());
    return;
  }  
  TFile* filMC=new TFile(fileNameMC.Data());
  if(filMC && filMC->IsOpen()){
    TDirectoryFile* dfMC=(TDirectoryFile*)filMC->Get(dirNameMC.Data());
    if(!dfMC){
      printf("TDirectoryFile %s not found in TFile for MC\n",dirNameMC.Data());
      filMC->ls();
      return;
    }
    TList* lMC=(TList*)dfMC->Get(lstNameMC.Data());
    hSigmaMC=FitMCInvMassSpectra(lMC,var);
    if(nBinsWithFixSig>0 && !hSigmaMC){
      printf("Fit to MC inv. mass spectra failed\n");
      return;
    }
    if(hSigmaMC && !useSigmaManual){
      for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++) sigmas[iPtBin]=hSigmaMC->GetBinContent(iPtBin+1);
    }
    if(correctForRefl){
      h3drefl=(TH3F*)lMC->FindObject(Form("hMassVsPtVs%sRefl",var.Data()));
      h3dmcsig=(TH3F*)lMC->FindObject(Form("hMassVsPtVs%sSig",var.Data()));
    }
  }
  
  TString sigConf="FixedSigma";
  if(nBinsWithFixSig==0) sigConf="FreeSigma";
  else if(nBinsWithFixSig==nPtBins) sigConf="FixedSigmaAll";
  else sigConf=Form("FixedSigma%d",patSig);
  if(useSigmaManual) sigConf.ReplaceAll("Sigma","SigmaManual");
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


  TCanvas* crf=new TCanvas("crf","Reflections",1200,800);
  DivideCanvas(crf,nPtBins);

  TCanvas* c1=new TCanvas("c1","Mass",1200,800);
  DivideCanvas(c1,nPtBins);

  TCanvas* c2=new TCanvas("c2","Mass-Bkg Rot",1200,800);
  DivideCanvas(c2,nPtBins);
  TCanvas* c2sub=new TCanvas("c2sub","Mass-Bkg-Fit Rot",1200,800);
  DivideCanvas(c2sub,nPtBins);
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
  TCanvas* c3sub=new TCanvas("c3sub","Mass-Bkg-Fit LS",1200,800);
  DivideCanvas(c3sub,nPtBins);
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
  TCanvas* c4sub=new TCanvas("c4sub","Mass-Bkg-Fit EM",1200,800);
  DivideCanvas(c4sub,nPtBins);
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

  TH1F* hIsSigmaFixed=new TH1F("hIsSigmaFixed","",nPtBins,binLims);
  TH1F* hSigmaFixedVal=new TH1F("hSigmaFixedVal","",nPtBins,binLims);

  TLatex* tME=new TLatex(0.65,0.82,"MixEv +- pairs");
  tME->SetNDC();
  TLatex* tMEpp=new TLatex(0.65,0.75,"MixEv ++ pairs");
  tMEpp->SetNDC();
  TLatex* tMEmm=new TLatex(0.65,0.68,"MixEv -- pairs");
  tMEmm->SetNDC();

  //  TF1 *fpeak=new TF1("fpeak","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",minMass,maxMass);

 
  TDirectory *current = gDirectory;
  TFile* fout=new TFile(Form("outputMassFits_%s_%s.root",sigConf.Data(),suffix.Data()),"recreate");
  current->cd();


  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    printf("\n---------- pt interval %d (%.1f-%.1f)\n",iPtBin,binLims[iPtBin],binLims[iPtBin+1]);
    minMass=minMass4Fit[iPtBin];
    maxMass=maxMass4Fit[iPtBin];
    Int_t bin1=h3d->GetYaxis()->FindBin(binLims[iPtBin]);
    Int_t bin2=h3d->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
    printf("Bin %d   Pt range=%f %f\n",iPtBin,h3d->GetYaxis()->GetBinLowEdge(bin1),h3d->GetYaxis()->GetBinUpEdge(bin2));
    TH1D* hMassPtBin=h3d->ProjectionX(Form("hMassPtBin%d",iPtBin),bin1,bin2,zbin1,zbin2);
    TH1D* hMassPtBinr=h3dr->ProjectionX(Form("hMassPtBinr%d",iPtBin),bin1,bin2,zbin1,zbin2);
    TH1D* hMassPtBinme=h3dme->ProjectionX(Form("hMassPtBinme%d",iPtBin),bin1,bin2,zbin1,zbin2);
    TH1D* hMassPtBinmeLSpp=h3dmepp->ProjectionX(Form("hMassPtBinmeLSpp%d",iPtBin),bin1,bin2,zbin1,zbin2);
    TH1D* hMassPtBinmeLSmm=h3dmemm->ProjectionX(Form("hMassPtBinmeLSmm%d",iPtBin),bin1,bin2,zbin1,zbin2);

    if(correctForRefl){
      Int_t bin1MC=h3drefl->GetYaxis()->FindBin(binLims[iPtBin]);
      Int_t bin2MC=h3drefl->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
      hMCReflPtBin=h3drefl->ProjectionX(Form("hMCReflPtBin%d",iPtBin),bin1MC,bin2MC,zbin1,zbin2);
      hMCReflPtBin->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
      hMCSigPtBin=h3dmcsig->ProjectionX(Form("hMCSigPtBin%d",iPtBin),bin1MC,bin2MC,zbin1,zbin2);
      hMCSigPtBin->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    }

    TH1D* hMassPtBinlsp=0x0;
    TH1D* hMassPtBinlsm=0x0;
    TH1D* hMassPtBinls=0x0;
    if(h3dlsp){
      hMassPtBinlsp=h3dlsp->ProjectionX(Form("hMassPtBinlsp%d",iPtBin),bin1,bin2,zbin1,zbin2);
      hMassPtBinlsm=h3dlsm->ProjectionX(Form("hMassPtBinlsm%d",iPtBin),bin1,bin2,zbin1,zbin2);
      hMassPtBinls=(TH1D*)hMassPtBinlsp->Clone(Form("hMassPtBinls%d",iPtBin));
      hMassPtBinls->Reset();
      for(Int_t iBin=1; iBin<=hMassPtBinlsp->GetNbinsX(); iBin++){
	Double_t np=hMassPtBinlsp->GetBinContent(iBin);
	Double_t nm=hMassPtBinlsm->GetBinContent(iBin);
	Double_t enp=hMassPtBinlsp->GetBinError(iBin);
	Double_t enm=hMassPtBinlsm->GetBinError(iBin);
	Double_t tt=0;
	Double_t ett=0;
	if(useGeomMeanLS){
	  tt=2*TMath::Sqrt(np*nm);
	  if(tt>0) ett=2./tt*TMath::Sqrt(np*np*enm*enm+nm*nm*enp*enp);
	}else{
	  tt=0.5*(np+nm);
	  ett=0.5*TMath::Sqrt(enm*enm+enp*enp);
	}
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


    TCanvas* c0=new TCanvas(Form("CBin%d",iPtBin),Form("Bin%d norm",iPtBin),1300,700);
    c0->Divide(3,2);
    c0->cd(1);
    TH1D* hRatio=(TH1D*)hMassPtBinr->Clone(Form("hRatioFormNorm%d",iPtBin));
    hRatio->Divide(hMassPtBin);
    hRatio->Draw();
    hRatio->GetYaxis()->SetTitle("Rotational/All");
    hRatio->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    Double_t normRot=GetBackgroundNormalizationFactor(hRatio,rebin[iPtBin]);
    hMassPtBinr->Scale(1./normRot);
    TLatex* tnr=new TLatex(0.2,0.2,Form("Normaliz. factor = %f",normRot));
    tnr->SetNDC();
    tnr->Draw();
    c0->cd(2);
    TH1D* hRatioLS=(TH1D*)hMassPtBinls->Clone(Form("hRatioFormNormLS%d",iPtBin));
    hRatioLS->Divide(hMassPtBin);
    hRatioLS->Draw();
    hRatioLS->GetYaxis()->SetTitle("LS/All");
    hRatioLS->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    if(renormLS){
      Double_t normLS=GetBackgroundNormalizationFactor(hRatioLS,rebin[iPtBin]);
      hMassPtBinls->Scale(1./normLS);
      TLatex* tnl=new TLatex(0.2,0.2,Form("Normaliz. factor = %f",normLS));
      tnl->SetNDC();
      tnl->Draw();
    }
    c0->cd(3);
    hMassPtBinme->GetYaxis()->SetTitle("Entries (EvMix)");
    hMassPtBinme->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hMassPtBinme->SetMinimum(0);
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
    hRatioME->Draw();
    hRatioME->GetYaxis()->SetTitle("EvMix (+-)/All");
    hRatioME->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    c0->cd(5);
    hRatioMEAll->Draw();
    hRatioMEAll->GetYaxis()->SetTitle("EvMix (+-,++,--)/All");
    hRatioMEAll->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    c0->cd(6);
    hRatioMEpp->Draw();
    hRatioMEpp->SetMinimum(0.49);
    hRatioMEpp->SetMaximum(0.51);
    hRatioMEpp->GetYaxis()->SetTitle("ME with LS / ME with OS");
    hRatioMEpp->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hRatioMEmm->Draw("same");
 
    Double_t normME=GetBackgroundNormalizationFactor(hRatioME,rebin[iPtBin]);
    Double_t normMEAll=GetBackgroundNormalizationFactor(hRatioMEAll,rebin[iPtBin]);
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
    if(fixSigmaConf==1 || (fixSigmaConf==2 && fixSigma[iPtBin])) hIsSigmaFixed->SetBinContent(iPtBin+1,1);
    else hIsSigmaFixed->SetBinContent(iPtBin+1,0);
    hSigmaFixedVal->SetBinContent(iPtBin+1,sigmas[iPtBin]);
    hRebin->SetBinContent(iPtBin+1,rebin[iPtBin]);
    Int_t bkgToFill=typeb;
    if(typeb==6) bkgToFill=typeb+nDegreeBackPol[iPtBin];
    Int_t bkgToFillSB=6+nDegreeBackPol[iPtBin];

    hBkgFitFunc->SetBinContent(iPtBin+1,bkgToFill);
    hBkgFitFuncSB->SetBinContent(iPtBin+1,bkgToFillSB);

    fitterRot[iPtBin]=ConfigureFitter(hMassSubRot,iPtBin,typeb,minMass,maxMass,0x0);
    if(hMassPtBinls) fitterLS[iPtBin]=ConfigureFitter(hMassSubLS,iPtBin,typeb,minMass,maxMass,0x0);
    fitterME[iPtBin]=ConfigureFitter(hMassSubME,iPtBin,typeb,minMass,maxMass,crf);

    Bool_t out1=fitterRot[iPtBin]->MassFitter(0);
    Bool_t out2=kFALSE;
    if(hMassPtBinls) out2=fitterLS[iPtBin]->MassFitter(0);
    Bool_t out3=fitterME[iPtBin]->MassFitter(0);
    
    Double_t background=999999999.;
    Bool_t out4=kFALSE;
    if(tryDirectFit){
      TH1D *hMassDirectFit=(TH1D*)hMassPtBin->Clone(Form("hMassDirectFit_bin%d",iPtBin));
      hMassDirectFit=AliVertexingHFUtils::RebinHisto(hMassDirectFit,rebin[iPtBin]);
      fitterSB[iPtBin]=ConfigureFitter(hMassDirectFit,iPtBin,6,fitSBrangelow[iPtBin],fitSBrangeup[iPtBin],0x0,kTRUE);
      out4=fitterSB[iPtBin]->MassFitter(0);//DirectFit(hMassDirectFit,iPtBin,hRawYieldSB);

      Double_t background,ebkg;

      if(out4 && fitterSB[iPtBin]->GetMassFunc()){
	c5->cd(iPtBin+1);
	fitterSB[iPtBin]->DrawHere(gPad,3,0);
	if(fitterSB[iPtBin]->GetRawYield()>0 && fitterSB[iPtBin]->GetReducedChiSquare()>0 && fitterSB[iPtBin]->GetReducedChiSquare()<5){
	  hRawYieldSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield());
	  hRawYieldSB->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetRawYieldError());
	  hRelStatSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYieldError()/fitterSB[iPtBin]->GetRawYield());
	  hRelStatSB->SetBinError(iPtBin+1,0.00000001);
	  fitterSB[iPtBin]->Background(3.,background,ebkg);
	  hSignifSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield()/TMath::Sqrt(background+fitterSB[iPtBin]->GetRawYield()));
	  hSignifSB->SetBinError(iPtBin+1,0.00000001);
	  hSoverBSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetRawYield()/background);
	  hSoverBSB->SetBinError(iPtBin+1,0.00000001);
	  hGausMeanSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetMean());
	  hGausMeanSB->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetMeanUncertainty());
	  hGausSigmaSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetSigma());
	  hGausSigmaSB->SetBinError(iPtBin+1,fitterSB[iPtBin]->GetSigmaUncertainty());
	  hChiSqSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetReducedChiSquare());
	  hChiSqSB->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	  hNdfSB->SetBinContent(iPtBin+1,fitterSB[iPtBin]->GetMassFunc()->GetNDF());
	  hNdfSB->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	}
	hInvMassHistoBinWidthSB->SetBinContent(iPtBin+1,hMassDirectFit->GetBinWidth(1));
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
	  if(binc>fitSBrangelow[iPtBin] && binc<fitSBrangeup[iPtBin]){
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
	hsubTemp->GetXaxis()->SetRangeUser(fitSBrangelow[iPtBin],fitSBrangeup[iPtBin]);
	hsubTemp->SetMinimum(ymin);
	hsubTemp->SetMaximum(ymax);
	hsubTemp->SetMarkerStyle(20);
	hsubTemp->SetMarkerColor(hsubTemp->GetLineColor());
	hsubTemp->DrawCopy();
	hsubTempAllRange->DrawCopy("same");
	hsubTemp->DrawCopy("same");
	TF1 *fpeak=fitterSB[iPtBin]->GetSignalFunc();
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
	hResiduals->SetName(Form("hResidualsSB_PtBin%d",iPtBin));

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
      if(fitterRot[iPtBin]->GetRawYield()>0 && fitterRot[iPtBin]->GetReducedChiSquare()>0 && fitterRot[iPtBin]->GetReducedChiSquare()<5){
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
	hGausMeanRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetMean());
	hGausMeanRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetMeanUncertainty());
	hGausSigmaRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetSigma());
	hGausSigmaRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetSigmaUncertainty());
	hChiSqRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetReducedChiSquare());
	hChiSqRot->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	hNdfRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetMassFunc()->GetNDF());
	hNdfRot->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      }
      hInvMassHistoBinWidthRot->SetBinContent(iPtBin+1,hMassSubRot->GetBinWidth(1));
      //      if(!correctForRefl)
      WriteFitInfo(fitterRot[iPtBin],hMassPtBin);
      fout->cd();
      WriteFitFunctionsToFile(fitterRot[iPtBin],"Rot",iPtBin);
      current->cd();
      
      Double_t errbc;
      Double_t bc=fitterRot[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
      hRawYieldRotBC->SetBinContent(iPtBin+1,bc);
      hRawYieldRotBC->SetBinError(iPtBin+1,errbc);

      c2sub->cd(iPtBin+1);
      fitterRot[iPtBin]->DrawHistoMinusFit(gPad);

      c2pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterRot[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);
      hResiduals->SetName(Form("hResidualsRot_PtBin%d",iPtBin));

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
      if(fitterLS[iPtBin]->GetRawYield()>0 && fitterLS[iPtBin]->GetReducedChiSquare()>0 && fitterLS[iPtBin]->GetReducedChiSquare()<5){
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
	hGausMeanLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetMean());
	hGausMeanLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetMeanUncertainty());
	hGausSigmaLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetSigma());
	hGausSigmaLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetSigmaUncertainty());
	hChiSqLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetReducedChiSquare());
	hChiSqLS->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	hNdfLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetMassFunc()->GetNDF());
	hNdfLS->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      }
      hInvMassHistoBinWidthLS->SetBinContent(iPtBin+1,hMassSubLS->GetBinWidth(1));
      //      if(!correctForRefl)
      WriteFitInfo(fitterLS[iPtBin],hMassPtBin);
      fout->cd();
      WriteFitFunctionsToFile(fitterLS[iPtBin],"LS",iPtBin);
      current->cd();


      Double_t errbc;
      Double_t bc=fitterLS[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
      hRawYieldLSBC->SetBinContent(iPtBin+1,bc);
      hRawYieldLSBC->SetBinError(iPtBin+1,errbc);

      c3sub->cd(iPtBin+1);
      fitterLS[iPtBin]->DrawHistoMinusFit(gPad);
      
      c3pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterLS[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);
      hResiduals->SetName(Form("hResidualsLS_PtBin%d",iPtBin));
 
      hPulls->Draw();
      PrintGausParams(hPulls);

      c3residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kLikeS);
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
      if(fitterME[iPtBin]->GetRawYield()>0 && fitterME[iPtBin]->GetReducedChiSquare()>0 && fitterME[iPtBin]->GetReducedChiSquare()<5){
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
	hGausMeanME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetMean());
	hGausMeanME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetMeanUncertainty());
	hGausSigmaME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetSigma());
	hGausSigmaME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetSigmaUncertainty());
	hChiSqME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetReducedChiSquare());
	hChiSqME->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
	hNdfME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetMassFunc()->GetNDF());
	hNdfME->SetBinError(iPtBin+1,0.00001); // very small number, for graphics
      }
      hInvMassHistoBinWidthME->SetBinContent(iPtBin+1,hMassSubME->GetBinWidth(1));
      //      if(!correctForRefl)
      WriteFitInfo(fitterME[iPtBin],hMassPtBin);
      fout->cd();
      WriteFitFunctionsToFile(fitterME[iPtBin],"ME",iPtBin);
      current->cd();

      Double_t errbc;
      Double_t bc=fitterME[iPtBin]->GetRawYieldBinCounting(errbc,nsigmaBinCounting,optBkgBinCount);
      hRawYieldMEBC->SetBinContent(iPtBin+1,bc);
      hRawYieldMEBC->SetBinError(iPtBin+1,errbc);

      c4sub->cd(iPtBin+1);
      fitterME[iPtBin]->DrawHistoMinusFit(gPad);

      c4pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=fitterME[iPtBin]->GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend);
      hResiduals->SetName(Form("hResidualsME_PtBin%d",iPtBin));
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
    crf->SaveAs(Form("figures/ReflTemplates_%s_%s.eps",sigConf.Data(),suffix.Data()));
    c1->SaveAs(Form("figures/InvMassSpectra_%s_%s_NoBkgSub.eps",sigConf.Data(),suffix.Data()));
    c2->SaveAs(Form("figures/InvMassSpectra_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
    c3->SaveAs(Form("figures/InvMassSpectra_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
    c4->SaveAs(Form("figures/InvMassSpectra_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));
    c2sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_Rot.eps",sigConf.Data(),suffix.Data()));
    c3sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_LS.eps",sigConf.Data(),suffix.Data()));
    c4sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_EM.eps",sigConf.Data(),suffix.Data()));
    if(tryDirectFit){
      c5->SaveAs(Form("figures/InvMassSpectra_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
      c5sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_SB.eps",sigConf.Data(),suffix.Data()));
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
    c2sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_Rot.root",sigConf.Data(),suffix.Data()));
    c3sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_LS.root",sigConf.Data(),suffix.Data()));
    c4sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_EM.root",sigConf.Data(),suffix.Data()));

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
      c5sub->SaveAs(Form("figures/InvMassSpectraFitSubtr_%s_%s_SB.root",sigConf.Data(),suffix.Data()));
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
  hRelStatRot->SetMaximum(0.48);
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
  hIsSigmaFixed->Write();
  hSigmaFixedVal->Write();
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
  if(fTot){
    nameh.Form("FitFuncTot_%s_PtBin%d",meth.Data(),iPtBin);
    fTot->SetRange(1.6,2.2);
    fTot->SetNpx(500);
    fTot->Write(nameh.Data());
  }
  TF1* fSig=fitter->GetSignalFunc();
  nameh.Form("FitFuncSig_%s_PtBin%d",meth.Data(),iPtBin);
  if(fSig){
    fSig->SetRange(1.6,2.2);
    fSig->SetNpx(500);
    fSig->Write(nameh.Data());
  }
  TF1* fBkg=fitter->GetBackgroundRecalcFunc();
  nameh.Form("FitFuncBkg_%s_PtBin%d",meth.Data(),iPtBin);
  if(fBkg){
    fBkg->SetRange(1.6,2.2);
    fBkg->SetNpx(500);
    fBkg->Write(nameh.Data());
  }
  if(meson=="Dzero"){
    TF1* fBkgR=fitter->GetBkgPlusReflFunc();
    nameh.Form("FitFuncBkgRefl_%s_PtBin%d",meth.Data(),iPtBin);
    if(fBkgR){
      fBkgR->SetRange(1.6,2.2);
      fBkgR->SetNpx(500);
      fBkgR->Write(nameh.Data());  
    }
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
  tss->AddText(Form("S/B (3#sigma) = %.4g",sig/back));
  if(correctForRefl) tss->AddText(Form("Refl/Sig =  %.3f #pm %.3f ",fitter->GetReflOverSig(),fitter->GetReflOverSigUncertainty()));
  tss->AddText(Form("Significance(3#sigma) = %.2f",sig/TMath::Sqrt(back+sig)));
  tss->SetTextColor(1);
  tss->Draw();

}

TH1F* FitMCInvMassSpectra(TList* lMC, TString var){
  TH3F* h3d=(TH3F*)lMC->FindObject(Form("hMassVsPtVs%sSig",var.Data()));
  if(!h3d){
    printf("hMassVsPtVsYSig not found\n");
    return 0x0;
  }
  TH3F* h3dr=(TH3F*)lMC->FindObject(Form("hMassVsPtVs%sRefl",var.Data()));
  if(!h3dr){
    printf("hMassVsPtVsYRefl not found\n");
  }
  TH1F* hSigmaMC=new TH1F("hSigmaMC","",nPtBins,binLims);
  Int_t zbin1=0;
  Int_t zbin2=h3d->GetZaxis()->GetNbins()+1;
  if(var=="CosthSt" && costhstcut<1.){
    zbin2=h3d->GetZaxis()->FindBin(costhstcut-0.000001);
  }

  TCanvas* cmc1=new TCanvas("InvMassMC","InvMassMC",1200,800);
  DivideCanvas(cmc1,nPtBins);

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    Int_t bin1=h3d->GetYaxis()->FindBin(binLims[iPtBin]);
    Int_t bin2=h3d->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
    Double_t ptmed=(binLims[iPtBin]+binLims[iPtBin+1])*0.5;
    Double_t pthalfwid=(binLims[iPtBin+1]-binLims[iPtBin])*0.5;
    printf("Bin %d   Pt range=%f %f\n",iPtBin,h3d->GetYaxis()->GetBinLowEdge(bin1),h3d->GetYaxis()->GetBinUpEdge(bin2));
    TH1D* hMassMCPtBin=h3d->ProjectionX(Form("hMassMCPtBin%d",iPtBin),bin1,bin2,zbin1,zbin2);
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
      TH1D* hReflPtBin=h3dr->ProjectionX(Form("hReflPtBin%d",iPtBin),bin1,bin2,zbin1,zbin2);
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
    Double_t mcSigma=fg->GetParameter(2);
    Double_t errMCsigma=fg->GetParError(2);
    if(tuneSigmaOnData>0.){
      mcSigma*=tuneSigmaOnData;
      errMCsigma*=tuneSigmaOnData;
    }
    hSigmaMC->SetBinContent(iPtBin+1,mcSigma);
    hSigmaMC->SetBinError(iPtBin+1,errMCsigma);
  }
  return hSigmaMC;
}

void PrintConfig(){
  printf("Meson: %s\n",meson.Data());
  printf("Data file: %s\n", fileName.Data());
  printf("Suffix Data: %s\n", suffix.Data());
  printf("MC file: %s\n", fileNameMC.Data());
  printf("Suffix MC: %s\n", suffixMC.Data());
  printf("Number of pt bins = %d\n",nPtBins);
  for(Int_t j=0; j<nPtBins+1; j++) printf(" %.1f",binLims[j]);
  printf("\n");
  for(Int_t j=0; j<nPtBins; j++) printf(" %d",rebin[j]);
  printf("\n");
  for(Int_t j=0; j<nPtBins; j++) printf(" %.2f",minMass4Fit[j]);
  printf("\n");
  for(Int_t j=0; j<nPtBins; j++) printf(" %.2f",maxMass4Fit[j]);
  printf("\n");
  printf("Gaussian sigma option: %d (tune on data = %f) (use manual sigma = %d)\n",fixSigmaConf,tuneSigmaOnData,useSigmaManual);
  for(Int_t j=0; j<nPtBins; j++) printf(" %d",fixSigma[j]);
  printf("\n");
  for(Int_t j=0; j<nPtBins; j++) printf(" %.3f",sigmas[j]);
  printf("\n");
  printf("Gaussian mean option: %d (tune on data = %f)\n",fixMeanConf,tuneMeanOnData);
  for(Int_t j=0; j<nPtBins; j++) printf(" %d",fixMean[j]);
  printf("\n");
  printf("Background options: Norm=%d  NormRange=%.2f  UseLSinEvMex=%d  SmoothLS=%d\n",optForNorm,rangeForNorm,useEMwithLS,smoothLS);
  printf("Direct fit: %d\n",tryDirectFit);
  for(Int_t j=0; j<nPtBins; j++) printf(" %d",nDegreeBackPol[j]);
  printf("\n");
  for(Int_t j=0; j<nPtBins; j++) printf(" %.2f",fitSBrangelow[j]);
  printf("\n");
  for(Int_t j=0; j<nPtBins; j++) printf(" %.2f",fitSBrangeup[j]);
  printf("\n");
  printf("Fit Configuration: Option=%s  Background Func=%d\n",fitoption.Data(),typeb);
  printf("Reflections: Use=%d  Func=%s  rOverSmodif=%f\n",correctForRefl,reflopt.Data(),rOverSmodif);
  printf("BinCounting: Option=%d  nsigma=%.2f\n",optBkgBinCount,nsigmaBinCounting);
}

Bool_t ReadConfig(TString configName){
  FILE* confFil=fopen(configName.Data(),"r");
  char line[50];
  char name[200];
  int n;
  float x;
  bool readok;
  while(!feof(confFil)){
    readok=fscanf(confFil,"%s:",line);
    if(strstr(line,"DataFile")){
      readok=fscanf(confFil,"%s",name);
      fileName=name;
    }
    else if(strstr(line,"SuffixData")){
      readok=fscanf(confFil,"%s",name);
      suffix=name;
    }
    else if(strstr(line,"MCFile")){
      readok=fscanf(confFil,"%s",name);
      fileNameMC=name;
    }
    else if(strstr(line,"SuffixMC")){
      readok=fscanf(confFil,"%s",name);
      suffixMC=name;
    }
    else if(strstr(line,"Meson")){
      readok=fscanf(confFil,"%s",name);
      meson=name;
    }
    else if(strstr(line,"NumOfPtBins")){
      readok=fscanf(confFil,"%d",&n);
      nPtBins=n;
    }
    else if(strstr(line,"BinLimits")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	readok=fscanf(confFil,"%f,",&x);
	binLims[j]=x;
	if(j>0 && binLims[j]<=binLims[j-1]){
	  printf("ERROR in array of pt bin limits\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil,"%f",&x);
      binLims[nPtBins]=x;
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"Rebin")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%d,",&n);
	else readok=fscanf(confFil,"%d",&n);
	rebin[j]=n;
	if(rebin[j]<=0 || rebin[j]>20){
	  printf("ERROR in array of rebin values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"MinFit")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%f,",&x);
	else readok=fscanf(confFil,"%f",&x);
	minMass4Fit[j]=x;
	if(minMass4Fit[j]<1.5 || minMass4Fit[j]>1.86){
	  printf("ERROR in array of min mass values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"MaxFit")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%f,",&x);
	else readok=fscanf(confFil,"%f",&x);
	maxMass4Fit[j]=x;
	if(maxMass4Fit[j]<1.86 || maxMass4Fit[j]>2.2){
	  printf("ERROR in array of max mass values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"FixSigmaConf")){
      readok=fscanf(confFil,"%d",&n);
      fixSigmaConf=n;
    }
    else if(strstr(line,"FixSigmaPerBin")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%d,",&n);
	else readok=fscanf(confFil,"%d",&n);
	fixSigma[j]=n;
	if(n<0 || n>1){
	  printf("ERROR in array of fix-sigma settings\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"TuneSigmaOnData")){
      readok=fscanf(confFil,"%f",&x);
      tuneSigmaOnData=x;
    }
    else if(strstr(line,"UseSigmaManual")){
      readok=fscanf(confFil,"%d",&n);
      useSigmaManual=n;
    }
    else if(strstr(line,"SigmaManual")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%f,",&x);
	else readok=fscanf(confFil,"%f",&x);
	sigmas[j]=x;
	if(sigmas[j]<0.001 || sigmas[j]>0.04){
	  printf("ERROR in array of sigma values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"FixMeanConf")){
      readok=fscanf(confFil,"%d",&n);
      fixMeanConf=n;
    }
    else if(strstr(line,"FixMeanPerBin")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%d,",&n);
	else readok=fscanf(confFil,"%d",&n);
	fixMean[j]=n;
	if(n<0 || n>1){
	  printf("ERROR in array of fix-sigma settings\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"TuneMeanOnData")){
      readok=fscanf(confFil,"%f",&x);
      tuneMeanOnData=x;
    }
    else if(strstr(line,"NormalizationOption")){
      readok=fscanf(confFil,"%d",&n);
      optForNorm=n;
    }
    else if(strstr(line,"RangeForNormalization")){
      readok=fscanf(confFil,"%f",&x);
      rangeForNorm=x;
    }
    else if(strstr(line,"UseEvMixWithLS")){
      readok=fscanf(confFil,"%d",&n);
      if(n<0 || n>1){
	printf("ERROR in UseEvMixWithLS setting\n");
	return kFALSE;
      }
      useEMwithLS=n;
    }
    else if(strstr(line,"UseGeomMeanLS")){
      readok=fscanf(confFil,"%d",&n);
      if(n<0 || n>1){
	printf("ERROR in UseGeomMeanLS setting\n");
	return kFALSE;
      }
      useGeomMeanLS=n;
    }
    else if(strstr(line,"RenormalizeLS")){
      readok=fscanf(confFil,"%d",&n);
      if(n<0 || n>1){
	printf("ERROR in RenormalizeLS setting\n");
	return kFALSE;
      }
      renormLS=n;
    }
    else if(strstr(line,"SmoothLS")){
      readok=fscanf(confFil,"%d",&n);
      smoothLS=n;
    }
    else if(strstr(line,"TryDirectFit")){
      readok=fscanf(confFil,"%d",&n);
      if(n<0 || n>1){
	printf("ERROR in TryDirectFit setting\n");
	return kFALSE;
      }
      tryDirectFit=n;
    }
    else if(strstr(line,"MinSBFit")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%f,",&x);
	else readok=fscanf(confFil,"%f",&x);
	fitSBrangelow[j]=x;
	if(fitSBrangelow[j]<1.5 || fitSBrangelow[j]>1.86){
	  printf("ERROR in array of min mass SB values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"MaxSBFit")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%f,",&x);
	else readok=fscanf(confFil,"%f",&x);
	fitSBrangeup[j]=x;
	if(fitSBrangeup[j]<1.86 || fitSBrangeup[j]>2.2){
	  printf("ERROR in array of max mass SB values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"PolDegreeSB")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%d,",&n);
	else readok=fscanf(confFil,"%d",&n);
	nDegreeBackPolSB[j]=n;
	if(nDegreeBackPolSB[j]<=0 || nDegreeBackPolSB[j]>10){
	  printf("ERROR in array of polynomial degree values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"FitOption")){
      readok=fscanf(confFil,"%s",name);
      fitoption=name;
    }
    else if(strstr(line,"BackgroundFunction")){
      readok=fscanf(confFil,"%d",&n);
      typeb=n;
    }
    else if(strstr(line,"PolDegreeFit")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	if(j<nPtBins-1) readok=fscanf(confFil,"%d,",&n);
	else readok=fscanf(confFil,"%d",&n);
	nDegreeBackPol[j]=n;
	if(nDegreeBackPol[j]<=0 || nDegreeBackPol[j]>10){
	  printf("ERROR in array of polynomial degree values\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil," ] ");
    }
    else if(strstr(line,"CorrectForReflections")){
      readok=fscanf(confFil,"%d",&n);
      correctForRefl=n;
    }
    else if(strstr(line,"ReflectionOption")){
      readok=fscanf(confFil,"%s",name);
      reflopt=name;
    }
    else if(strstr(line,"ReflOverSignalModif")){
      readok=fscanf(confFil,"%f",&x);
      rOverSmodif=x;
    }
    else if(strstr(line,"BackgroundFuncOptionForBinCount")){
      readok=fscanf(confFil,"%d",&n);
      optBkgBinCount=n;
    }
    else if(strstr(line,"NumOfSigmaForBinCount")){
      readok=fscanf(confFil,"%f",&x);
      nsigmaBinCounting=x;
    }
    else if(strstr(line,"CosThetaStar")){
      readok=fscanf(confFil,"%f",&x);
      costhstcut=x;
    }
  }
  return kTRUE;
}
