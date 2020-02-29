#ifdef __CLING__
// Tell  ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#endif
#if !defined (__CINT__) || defined (__CLING__)
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TF1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TInterpreter.h>
#include <TString.h>
#include <TLatex.h>
#include "AliNormalizationCounter.h"
#endif

#include <vector>

/* ROUGH MACRO TO ANALYSE SigmaC TASK OUTPUT
   A. Rossi, M. Faggin

   Goal: bash script that reproduces all plots, studies and results
   - yield extraction:
           - study of Lc SB for background shape, rotational background
	          - chi2 shape selection for functional forms?
		  - real Lc with background pion? 
      	               - check that Lc SB background reproduce real Sc background
		       --> would require EM but still not expected to give ideal answer (better if EM with events with Lc candidate with m in Lc peak)
           - multitrial 
   
   - efficiency, including period dependence, resonant decay-channel, steps
         - report validation with Lc efficiency in Analysis Note
   - cross section and final plots (ratios as well)

Other open points:
- in MC in the task: include resonant channels + forsee acceptance correction from outside assumptions
*/

Int_t version=1;
Double_t minMassScFit=0.148,maxMassScFit=0.190;// so-far 0.2
Double_t fixSigma=0.0012;
Double_t fixSigmaBW=0.0007;
Double_t widthSigmaBW=0.00189;
enum ECandStatus {kGenLimAcc=1,kGenAccMother,kGenAcc,kReco=6,kRecoCuts,kRecoPID,kRecoLc=13,kRecoLcCuts,kRecoLcPID};

//const Int_t nbinsGlobal=6;
//Double_t ptbinsGlobal[nbinsGlobal+1]={1.,2.,3.,5.,8,16,24.};
//Int_t skipBinGlobal[nbinsGlobal]={0,0,0,0,0,0};

const Int_t nbinsGlobal=6;
Double_t ptbinsGlobal[nbinsGlobal+1]={1.,2.,4.,6.,8,12,24.};
Int_t skipBinGlobal[nbinsGlobal]={0,0,0,0,0,0};

// pT to consider for the Sc analysis
//  0: pT of Lc(<-Sc)
//  10: pT of SigmaC
//
//  TODO
//

TH1D *hEffRecoLc;
TH1D *hEffRecoLcSc;
//
//  create other histos in order to consider the 
//  efficiency separately for the different decay channels (direct, resonant)
//
TH1D *hEffRecoLc_direct, *hEffRecoLc_res2, *hEffRecoLc_res3, *hEffRecoLc_res4;
TH1D *hEffRecoLcSc_direct, *hEffRecoLcSc_res2, *hEffRecoLcSc_res3, *hEffRecoLcSc_res4;
double br_direct = 3.50/100;  double err_br_direct = 0.40/100;  // decay in pK-pi+ (direct)
double br_res_2  = 1.96/100;  double err_br_res_2  = 0.27/100;  // decay in pK*(892) case
double br_res_3  = 1.08/100;  double err_br_res_3  = 0.25/100;  // decay in Δ(1232)++ K-
double br_res_4  = 2.20/100;  double err_br_res_4  = 0.50/100;  // decay in Lambda(1520) pi+

TH1D *hRawYieldLc;
TH1D *hRawYieldLcSc;
AliNormalizationCounter *normCount;

void DivideCanvas(TCanvas *c,Int_t nbins){
  if(nbins<=4)c->Divide(2,2);
  else if(nbins<=6)c->Divide(2,3);
  else if(nbins<=9)c->Divide(3,3);
  else if(nbins<=12)c->Divide(4,3);
  else if(nbins<=16)c->Divide(4,4);
  else c->Divide(5,5);    
}

void WriteFitResultsAndAxisOnPad(TH1 *h,TF1 *f,TPad *p,Double_t signal,Double_t errsignal,Double_t signifRange=-1/*used to define the range around the function mean that will be used to calculate the significance and S/B*/,Bool_t isPerformance=kFALSE,TString strPtPerf=""){
  p->cd();
  TLatex *tlS=new TLatex(0.16,0.25,Form("S = %.1f #pm %.1f",signal,TMath::Abs(errsignal)));
  tlS->SetNDC();
  tlS->SetTextSizePixels(21.);
  TLatex *tlMean=new TLatex(0.16,0.21,Form("#mu = (%.4f #pm %.4f) GeV/#it{c}^{2}",f->GetParameter(1),f->GetParError(1)));
  tlMean->SetNDC();
  tlMean->SetTextSizePixels(21.);
  TLatex *tlSigma=new TLatex(0.16,0.17,Form("#sigma = (%.4f #pm %.4f) GeV/#it{c}^{2}",f->GetParameter(2),f->GetParError(2)));
  tlSigma->SetNDC();
  tlSigma->SetTextSizePixels(21.);
  TLatex *tlSignifErrSign=new TLatex(0.16,0.12,Form("S/#sigma(S) = %.1f",signal/TMath::Abs(errsignal)));
  tlSignifErrSign->SetNDC();
  tlSignifErrSign->SetTextSizePixels(21.);
    
  TLatex *tlSignif,*tlSoverB;
  Double_t sOverB;
  Double_t background;
  Double_t signalReduceRange=signal;  
  if(errsignal>0){

    if(signifRange<0){
      signalReduceRange=signal*0.9973;// assumed Gaussian is used for signal
      background=f->Integral(f->GetParameter(1)-3.*f->GetParameter(2),f->GetParameter(1)+3.*f->GetParameter(2))/h->GetBinWidth(1)-signalReduceRange;
      tlSignif=new TLatex(0.16,0.07,Form("S/#sqrt{S+B}(3#sigma) = %.1f",signalReduceRange/TMath::Sqrt(signalReduceRange+background)));// assumed Gaussian is used for signal  // ([0]/binwidth) / sqrt(([0]+[back])/binwidth) -> ([0]) / sqrt(([0]+[back])*binwidth)       
      tlSoverB=new TLatex(0.25,0.07,Form("S/B(3#sigma) = %.3f",signalReduceRange/background));
    }
    else{
      TF1 *ftemp=f->DrawCopy();
      ftemp->SetParameter(0,0);
      background=ftemp->Integral(ftemp->GetParameter(1)-signifRange,ftemp->GetParameter(1)+signifRange)/h->GetBinWidth(1);
      ftemp->SetParameter(0,f->GetParameter(0));
      ftemp->SetParameter(3,0);
      signalReduceRange=ftemp->Integral(ftemp->GetParameter(1)-signifRange,ftemp->GetParameter(1)+signifRange)/h->GetBinWidth(1);//ftemp->Integral(f->GetParameter(1)-signifRange,f->GetParameter(1)+signifRange)/h->GetBinWidth(1)-signalReduceRange;
      Printf("s(FWHM) %f, b (FWHM) %f",signalReduceRange,background);
      delete ftemp;
      if(isPerformance){
	tlSignif=new TLatex(0.16,0.07,Form("S/#sqrt{S+B}(2 FWHM) = %.1f",signalReduceRange/TMath::Sqrt(signalReduceRange+background)));
	tlSoverB=new TLatex(0.25,0.07,Form("S/B(2 FWHM) = %.3f",signalReduceRange/background));
      }
      else {
	tlSignif=new TLatex(0.16,0.07,Form("S/#sqrt{S+B}(|#it{M}-#mu_{#Sigma_{c}^{++}}| < %.1f MeV/#it{c}^{2}) = %.1f",signifRange,signalReduceRange/TMath::Sqrt(signalReduceRange+background)));
	tlSoverB=new TLatex(0.25,0.07,Form("S/B(|#it{M}-#mu_{#Sigma_{c}^{++}}| < %.1f MeV/#it{c}^{2}) = %.3f",signifRange,signalReduceRange/background));
      }
    }
  }
  else {
    Int_t binm=0,binM=-1;
    if(signifRange<0){
      // assumes h is background subtracted histo and f a fit to it --> must recalculate the background
      binm=h->FindBin(f->GetParameter(1)-3.*f->GetParameter(2)*1.001);
      binM=h->FindBin(f->GetParameter(1)+3.*f->GetParameter(2)*0.999);
      signalReduceRange=signal*0.9973; // assumed Gaussian is used for signal
    }
    else{
      binm=h->FindBin(f->GetParameter(1)-signifRange);
      binM=h->FindBin(f->GetParameter(1)+signifRange);
      TF1 *ftemp=new TF1();
      ftemp->Copy(*f);
      ftemp->SetParameter(0,signal*h->GetBinWidth(1));
      ftemp->SetParameter(3,0);
      signalReduceRange=ftemp->Integral(ftemp->GetParameter(1)-signifRange,ftemp->GetParameter(1)+signifRange)/h->GetBinWidth(1);
      delete ftemp;
    }
    Double_t sqrtSplusB=0;
    for(Int_t ib=binm;ib<=binM;ib++){
      sqrtSplusB+=h->GetBinError(ib)*h->GetBinError(ib);
    }
    sqrtSplusB=TMath::Sqrt(sqrtSplusB);
    
    tlSignif=new TLatex(0.16,0.07,Form("S/#sqrt{S+B}(3#sigma) = %.1f",signalReduceRange/TMath::Sqrt(signalReduceRange+background)));// assumed Gaussian is used for signal  // ([0]/binwidth) / sqrt(([0]+[back])/binwidth) -> ([0]) / sqrt(([0]+[back])*binwidth)       
    tlSoverB=new TLatex(0.25,0.07,Form("S/B(3#sigma) = %.3f",signalReduceRange/background));
    if(signifRange<0){
      if(isPerformance){
	tlSignif->SetTitle(Form("S/#sqrt{S+B}(2 FWHM) = %.1f",signalReduceRange/TMath::Sqrt(signalReduceRange+background)));
	tlSoverB->SetTitle(Form("S/B(2FWHM) = %.3f",signalReduceRange/background));
      }
      else {
	tlSignif->SetTitle(Form("S/#sqrt{S+B}(|#it{M}-#mu_{#Sigma_{c}^{++}}| < %.1f MeV/#it{c}^{2}) = %.1f",signifRange,signalReduceRange/TMath::Sqrt(signalReduceRange+background)));
	tlSoverB->SetTitle(Form("S/B(|#it{M}-#mu_{#Sigma_{c}^{++}}| < %.1f MeV/#it{c}^{2}) = %.3f",signifRange,signalReduceRange/background));
      }
    }

  }
  
  sOverB=signalReduceRange/background;
  tlSignif->SetNDC();
  tlSignif->SetTextSizePixels(21.);
  tlSoverB->SetNDC();
  tlSoverB->SetTextSizePixels(21.);
  //      h->GetXaxis()->SetRangeUser(0.135,0.23);
  h->SetXTitle("#it{M}_{pK#pi#pi}#font[122]\{-}#it{M}_{pK#pi} (GeV/#it{c}^{2})");
  //      h->SetYTitle(Form("Counts/%.1f MeV/#it{c}^{2}",h->GetBinWidth(1)*1000));
  h->SetYTitle(TMath::Abs((Int_t)(h->GetBinWidth(1)*1000)-h->GetBinWidth(1)*1000)<1.e-5 ? Form("Counts/%.0f MeV/#it{c}^{2}",h->GetBinWidth(1)*1000) : Form("Counts/%.1f MeV/#it{c}^{2}",h->GetBinWidth(1)*1000));
  
  if(errsignal>0 && signifRange<0){
    signalReduceRange=signal*0.9973;// assumed Gaussian is used for signal
    background=f->Integral(f->GetParameter(1)-3.*f->GetParameter(2),f->GetParameter(1)+3.*f->GetParameter(2))/h->GetBinWidth(1)-signalReduceRange;
    sOverB=signalReduceRange/background;
    tlSignif=new TLatex(0.16,0.07,Form("S/#sqrt{S+B}(3#sigma) = %.1f",signal/TMath::Sqrt(signalReduceRange+background)));// assumed Gaussian is used for signal  // ([0]/binwidth) / sqrt(([0]+[back])/binwidth) -> ([0]) / sqrt(([0]+[back])*binwidth) 
  }
  else if(signifRange<0){// assumes h is background subtracted histo and f a fit to it --> must recalculate the background
    Int_t binm=h->FindBin(f->GetParameter(1)-3.*f->GetParameter(2)*1.001);
    Int_t binM=h->FindBin(f->GetParameter(1)+3.*f->GetParameter(2)*0.999);
    Double_t sqrtSplusB=0;
    for(Int_t ib=binm;ib<=binM;ib++){
      sqrtSplusB+=h->GetBinError(ib)*h->GetBinError(ib);
    }
    sqrtSplusB=TMath::Sqrt(sqrtSplusB);
    
    tlSignif=new TLatex(0.16,0.07,Form("S/#sqrt{S+B} = %.1f",signal/sqrtSplusB));// ([0]/binwidth) / sqrt(([0]+[back])/binwidth) -> ([0]) / sqrt(([0]+[back])*binwidth)      
    }
    if(!isPerformance){
      tlS->Draw();
      tlMean->Draw();
      tlSigma->Draw();
      tlSignifErrSign->Draw();
      tlSignif->Draw();  
    }
    else{
      tlMean->SetTitle(Form("#mu_{#Sigma_{c}^{++}} = (%.4f #pm %.4f) GeV/#it{c}^{2}",f->GetParameter(1),f->GetParError(1)));
      tlMean->SetX(0.29);
      tlMean->SetY(0.37);
      tlMean->Draw();    
      
      TLatex *tlDeltaMean=new TLatex(0.29,0.32,"#mu_{#Sigma_{c}^{++}} #font[122]\{-} #mu_{#Sigma_{c}^{0}} = 0.22 MeV/#it{c}^{2} (fixed)");
      tlDeltaMean->SetNDC();
      tlDeltaMean->SetTextSizePixels(21.);
      tlDeltaMean->Draw();

      if(fixSigmaBW>0.)tlSigma->SetTitle(Form("#Gamma_{#Sigma_{c}^{0,++}} = %.2f MeV/#it{c}^{2}, #sigma = %.2f MeV/#it{c}^{2} (fixed)",widthSigmaBW*1000.,fixSigmaBW*1000.));
      else tlSigma->SetTitle(Form("#Gamma = %.2f (fixed) MeV/#it{c^}{2}, #sigma = (%.2f #pm %.2f) MeV/#it{c}^{2}",widthSigmaBW*1000,f->GetParameter(2),f->GetParError(2)*1000.));
      tlSigma->SetX(0.29);
      tlSigma->SetY(0.26);
      tlSigma->Draw();
      
      tlS->SetX(0.29);
      tlS->SetY(0.21);
      tlS->Draw();
      
      //      tlSignifErrSign->SetX(0.29);
      //      tlSignifErrSign->SetY(0.16);
      //      tlSignifErrSign->Draw();
      //      tlSignif->SetX(0.29);
      //      tlSignif->SetY(0.16);
      //      tlSignif->Draw();

      tlSoverB->SetX(0.29);
      tlSoverB->SetY(0.16);
      tlSoverB->Draw();


      TLatex *tlALICEperf=new TLatex(0.60,0.88,"ALICE Performance");
      tlALICEperf->SetNDC();
      tlALICEperf->SetTextSizePixels(24.);
      TLatex *tlChannel = new TLatex(0.18,0.83,"#Sigma_{c}^{0,++}#rightarrow #Lambda_{c}^{+}#pi^{#font[122]\{-},+},#Lambda_{c}^{+}#rightarrow pK^{#font[122]\{-}}#pi^{+}");
      tlChannel->SetNDC();
      tlChannel->SetTextSizePixels(24.);
      TLatex *tlChConj = new TLatex(0.22,0.78," + charge conj.");
      tlChConj->SetNDC();
      tlChConj->SetTextSizePixels(24.);
      TLatex *tlCollEnergy=new TLatex(0.65,0.83,"pp, #sqrt{#it{s}} = 13 TeV");
      tlCollEnergy->SetNDC();
      tlCollEnergy->SetTextSizePixels(24.);
      TLatex *tptrange= new TLatex(0.18,0.73,strPtPerf.Data());
      tptrange->SetNDC();
      tptrange->SetTextSizePixels(24.);
      tlALICEperf->Draw();
      tlChannel->Draw();
      tlChConj->Draw();
      tlCollEnergy->Draw();
      tptrange->Draw();
    }
}

Double_t fitSigmaBWGauss(Double_t *x,Double_t *par){
  return par[0]*0.5*(TMath::Voigt(x[0]-par[1],par[2],widthSigmaBW)+TMath::Voigt(x[0]-par[1]+0.00022,par[2],widthSigmaBW))+par[3]*TMath::Sqrt(TMath::Power(x[0]-0.139,par[4]));
}

Double_t fitSigmaBWGaussLinBack(Double_t *x,Double_t *par){
  return par[0]*0.5*(TMath::Voigt(x[0]-par[1],par[2],widthSigmaBW)+TMath::Voigt(x[0]-par[1]+0.00022,par[2],widthSigmaBW))+par[3]*(1.+par[4]*x[0]);
}

Double_t fitSigmaBWGaussPol2Back(Double_t *x,Double_t *par){
  return par[0]*0.5*(TMath::Voigt(x[0]-par[1],par[2],widthSigmaBW)+TMath::Voigt(x[0]-par[1]+0.00022,par[2],widthSigmaBW))+par[3]*(1.+par[4]*x[0]+par[5]*x[0]*x[0]);
}

Double_t fitSigmaBWGaussPol3Back(Double_t *x,Double_t *par){
  return par[0]*0.5*(TMath::Voigt(x[0]-par[1],par[2],widthSigmaBW)+TMath::Voigt(x[0]-par[1]+0.00022,par[2],widthSigmaBW))+par[3]*(1.+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]);
}

Double_t fitSigmaBWGaussPol4Back(Double_t *x,Double_t *par){
  return par[0]*0.5*(TMath::Voigt(x[0]-par[1],par[2],widthSigmaBW)+TMath::Voigt(x[0]-par[1]+0.00022,par[2],widthSigmaBW))+par[3]*(1.+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0]*x[0]);
}



TF1* FitHisto(TH1D* h, Int_t particle/*Lc=0, Sigma_c=1, Xic=2*/,Int_t backOpt,Int_t signalOpt,Double_t &signal,Double_t &errsignal,TString strFitOption=""){
  if(strFitOption.IsNull()){
    strFitOption="RLEMI";
  }
  TF1 *f;
  if(particle==1){ 
    if(signalOpt==1&&backOpt==-1){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*TMath::Sqrt(TMath::Power(x-0.139,[4]))",0.139,0.22);
    }
    else if(signalOpt==2&&backOpt==-2){
      f=new TF1("fFit",&fitSigmaBWGauss,0.139,0.22,4);
    }
    else if(signalOpt==2&&(backOpt==-3 || backOpt==-4)){
      f=new TF1("fFit",&fitSigmaBWGaussLinBack,0.139,0.22,5);
    }
    else if(signalOpt==2&&backOpt==2){
      f=new TF1("fFit",&fitSigmaBWGaussPol2Back,0.139,0.22,6);
    }
    else if(signalOpt==2&&backOpt==3){
      f=new TF1("fFit",&fitSigmaBWGaussPol3Back,0.139,0.22,7);
    }
    else if(signalOpt==2&&backOpt==4){
      f=new TF1("fFit",&fitSigmaBWGaussPol4Back,0.139,0.22,8);
    }
    else if(signalOpt==1&&backOpt==0){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]",0.139,0.22);
    } 
    else if(signalOpt==1&&backOpt==1){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*(1.+[4]*x)",0.139,0.22);
    } 
    else if(signalOpt==1&&backOpt==2){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*(1.+[4]*x+[5]*x*x)",0.139,0.22);
    }
    else if(signalOpt==1&&backOpt==3){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*(1.+[4]*x+[5]*x*x+[6]*x*x*x)",0.139,0.22);
    }
    else if(signalOpt==1&&backOpt==4){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*(1.+[4]*x+[5]*x*x+[6]*x*x*x+[7]*x*x*x*x)",0.139,0.22);
    }
    else if(signalOpt==1&&backOpt==5){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*TMath::Sqrt(x-0.139)*TMath::Exp([4]*(x-0.139))",0.139,0.22);
    }
    else if(signalOpt==1&&backOpt==6){
      f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*TMath::Log([4]*(x-0.139))",0.139,0.22);
    }
    else f=new TF1("fFit","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]*(1.+[4]*x+[5]*x*x)",0.,3.);
  
    // SETTING PARAMETERS

    f->SetParLimits(0,0.001,h->Integral(h->FindBin(0.1390),h->FindBin(0.220))/(h->FindBin(0.220)-(h->FindBin(0.1390))));
    f->SetParLimits(1,0.1390,0.190);
    f->SetParLimits(2,0.0001,0.004);
    if(backOpt<0){
      f->SetParLimits(3,0.00001,h->GetEntries());
    }
    if(backOpt==-1)f->SetParameters(h->Integral(h->FindBin(0.160),h->FindBin(0.180))/(h->FindBin(0.180)-(h->FindBin(0.160))),0.168,TMath::Abs(fixSigma));
    if(backOpt==-2 || backOpt==-3 || signalOpt==2)f->SetParameters(h->Integral(h->FindBin(0.160),h->FindBin(0.180))/(h->FindBin(0.180)-(h->FindBin(0.160))),0.168,TMath::Abs(fixSigmaBW));
    if(backOpt==5)f->SetParLimits(4,-1.e+5,0.);
    if(backOpt==6)f->SetParLimits(4,1.e-6,999.);
    // now fit
    if(backOpt<0||(backOpt>=2&&backOpt<=5))f->SetParameter(3,h->Integral(h->FindBin(0.1390),h->FindBin(0.220))/(h->FindBin(0.220)-(h->FindBin(0.1390))));
    else if(backOpt==0)f->SetParameter(3,0);
    //    Printf("%f",h->Integral(h->FindBin(0.1390),h->FindBin(0.220))/(h->FindBin(0.220)-(h->FindBin(0.1390))));
    if(fixSigma>0. && backOpt>=-1 && signalOpt!=2){//overwritten later for BW signal fit
      f->FixParameter(2,fixSigma);
    }
    if(backOpt==-4){
      f->FixParameter(3,0);
      //      f->FixParameter(4,0);
    }
    if(TMath::Abs(signal)<1.e-9){
      f->FixParameter(0,0);
      f->FixParameter(1,0.167);// this will not be effective since the yield is fixed to 0, but need to avoid numerical issues as 0/0 divisions
      f->FixParameter(2,0.0012);// this will not be effective since the yield is fixed to 0, but need to avoid numerical issues as 0/0 divisions
    }
    if(fixSigmaBW>0. && (backOpt==-2 || backOpt==-3 || signalOpt==2 ))f->FixParameter(2,fixSigmaBW);
    if(backOpt==3 || backOpt==4){
      f->SetParameter(3,h->Integral(h->FindBin(0.1390),h->FindBin(0.220))/(h->FindBin(0.220)-(h->FindBin(0.1390))));
      f->FixParameter(6,0);
      if(backOpt==4)f->FixParameter(7,0);
      h->Fit(f,"REM","",minMassScFit*1.2,maxMassScFit*0.8);   
      f->ReleaseParameter(6);
      if(backOpt==4){
	h->Fit(f,"REM","",minMassScFit*1.1,maxMassScFit*0.9);   
	f->ReleaseParameter(7);
      }
    }
    
    h->Fit(f,strFitOption.Data(),"",minMassScFit,maxMassScFit);   

  }


  if(particle==0){
    //    f->FixParameter(0,0);
    f=new TF1("f","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]+[4]*x+[5]*x*x",0.,3.);
    Double_t b=(h->GetBinContent(h->FindBin(2.416))-h->GetBinContent(h->FindBin(2.315)))/(2.416-2.315);
    Double_t a=h->GetBinContent(h->FindBin(2.315))-b*2.315;
    Double_t c=(h->GetBinContent(h->FindBin(2.180))-a-b*2.18)/(2.18*2.18);
    c=((h->GetBinContent(h->FindBin(2.380))-a-b*2.380)/(2.38*2.38));
    f->SetParameter(3,a);
    f->SetParameter(4,b);
    f->SetParameter(5,c);
    f->SetParameters(100,2.286,0.005,a,b,c);
    //f->SetParLimits(5,f->GetParameter(5)*0.3,f->GetParameter(5)*2);
    //    signal=10;errsignal=1;
    //    f->SetParameters(100,2.286,0.005,a,b,c);
    //    return f;
    f->SetParLimits(3,f->GetParameter(3)*0.3,f->GetParameter(3)*3);
    //    f->SetParLimits(4,f->GetParameter(4)*0.3,f->GetParameter(4)*3.);
    //    f->SetParLimits(5,f->GetParameter(5)*0.3,f->GetParameter(5)*3.);

    h->Fit(f,"REM","0",2.225,2.346);
    f->ReleaseParameter(0);
    f->SetParLimits(0,0.001,1.e+9);
    //    f->SetParLimits(1,2.390,2.430);//    
    f->SetParLimits(1,2.25,2.31);
    f->SetParLimits(2,0.0001,0.030);
    f->SetParLimits(3,f->GetParameter(3)*0.3,f->GetParameter(3)*3);
    //    f->SetParLimits(4,f->GetParameter(4)*0.3,f->GetParameter(4)*3.);
    //f->SetParLimits(5,f->GetParameter(5)*0.3,f->GetParameter(5)*3.);

    f->SetParameters(100,2.286,0.005,f->GetParameter(3),f->GetParameter(4),f->GetParameter(5));
    //f->SetParameters(100,2.415,0.005,f->GetParameter(3),f->GetParameter(4),f->GetParameter(5));
  }
  else if(particle==1){
    // nothing to do here
  }
  else if(particle==2){
    f->SetParLimits(0,0.001,1.e+9);
    f->SetParLimits(1,2.418,2.508);
    f->SetParLimits(2,0.0001,0.080);
    f->SetParLimits(3,0.00001,h->GetEntries());
  }
  else {
    Printf("Wrong particle option for fitting");
    return 0x0;
  }

  if(particle==0){
    h->Fit(f,strFitOption.Data(),"",2.225,2.346);   
    //h->Fit(f,"RLEMI","",2.330,2.470);   
  }
  if(particle==2){
    h->Fit(f,strFitOption.Data(),"",2.350,2.590);   
  }
  
  if(!(TMath::Abs(signal)<1.e-9))signal=f->GetParameter(0)/h->GetBinWidth(1);
  if(!(TMath::Abs(signal)<1.e-9))errsignal=f->GetParError(0)/h->GetBinWidth(1);
  return f;
  
}

TF1* FitHisto(TH1D* h, Int_t particle/*Lc=0, Sigma_c=1, Xic=2*/,Int_t linPol){

  Double_t signal,errsignal;
  return FitHisto(h,particle,linPol,1,signal,errsignal);
}

void FitHistoFromFile(TString strfile,Int_t noCut=0,Int_t particle=1,Int_t nptbins=6,Int_t ptbinUnique=-1,Int_t var1fix=-1,Int_t var2fix=-2,Bool_t pickupPad=kFALSE,Bool_t refit=kFALSE){
  TFile *f=TFile::Open(strfile.Data(),"READ");
  TPad *pickup=0x0;
  TCanvas *canvOut;
  TH1D *hSel=0x0;
  
  if(noCut==1 &&pickupPad && refit){
    TString strCanvName="cCanvasNoCutsSc";
    TString strHist="hNoCutsSc_PtBin";
    strHist+=ptbinUnique;
    TCanvas *c=(TCanvas*)f->Get(strCanvName.Data());
    TPad *p=(TPad*)c->GetPad(ptbinUnique+1);
    TH1D *h=(TH1D*)p->FindObject(strHist.Data());
    //    maxMassScFit=0.186;
      canvOut=new TCanvas("cSigmaout","cSigmaout",800,800);
      canvOut->cd();
      hSel=(TH1D*)h->Clone(Form("hNoCutsScClone_PtBin%d",ptbinUnique));
      Double_t signal,errsignal;
      TF1 *f=FitHisto(hSel,particle,-1,1,signal,errsignal);
      f->SetLineColor(kBlue);
      TF1 *fBack=(TF1*)f->Clone("fback");
      fBack->SetParameter(0,0);
      fBack->SetLineStyle(2);
      fBack->SetLineColor(kGray);
      fBack->SetRange(minMassScFit,maxMassScFit);
      fBack->Draw("same");
      TF1 *fSign=(TF1*)f->Clone("fsign");
      fSign->SetParameter(3,0);
      Printf("Fit done");	  
      WriteFitResultsAndAxisOnPad(hSel,f,(TPad*)(canvOut->cd()),signal,errsignal);
      
      TCanvas *cRes=new TCanvas("cRes","cRes",800,800);
      cRes->cd();
      TH1D *hRes=new TH1D(Form("hNoCutsScClone_PtBin%d",ptbinUnique),Form("hNoCutsScClone_PtBin%d",ptbinUnique),hSel->GetNbinsX(),hSel->GetBinLowEdge(1),hSel->GetBinLowEdge(hSel->GetNbinsX())+hSel->GetBinWidth(hSel->GetNbinsX()));
      //(TH1D*)hSel->Clone(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      //      TH1D *hRes=(TH1D*)hSel->Clone(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      Printf("Calculating residuals");
      for(Int_t kr=1;kr<=hRes->GetNbinsX();kr++){
	if(hRes->GetBinLowEdge(kr)>=minMassScFit && (hRes->GetBinLowEdge(kr)+hRes->GetBinWidth(kr))<maxMassScFit){
	  hRes->SetBinContent(kr,hSel->GetBinContent(kr)-fBack->Eval(hRes->GetBinCenter(kr)));//-fBack->Integral(hRes->GetBinLowEdge(kr),hRes->GetBinLowEdge(kr)+hRes->GetBinWidth(kr))/hRes->GetBinWidth(kr));
	  hRes->SetBinError(kr,hSel->GetBinError(kr));
	}
	else hRes->SetBinContent(kr,0);	
      }
      hRes->Draw();
      fSign->SetRange(minMassScFit,maxMassScFit);
      fSign->Draw("same");
      WriteFitResultsAndAxisOnPad(hRes,f,(TPad*)(cRes->cd()),signal,errsignal);
      return;

  }
  for(Int_t j=0;j<nptbins;j++){
    if(ptbinUnique>=0 && j!=ptbinUnique)continue;
    TString strCanvName="cCutScan";
    TString strHist="";
    if(strfile.Contains("CombinedLcSc_") || strfile.Contains("CombLcSc_")){
      if(particle==0){
	strCanvName.Append("Lc");
	strHist.Append("Lc");
      }
      if(particle==1){
	strCanvName.Append("Sc");
	strHist.Append("Sc");
      }
    }
 
    TCanvas *c=(TCanvas*)f->Get(Form("%s%d",strCanvName.Data(),j));
    c->Draw();
    if(pickupPad && refit){
      canvOut=new TCanvas("cSigmaout","cSigmaout",800,800);
      canvOut->cd();
      TPad *p=(TPad*)c->GetPad(var2fix+var1fix*5+1);
      TH1D *h=(TH1D*)p->FindObject(Form("hVarScan%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      hSel=(TH1D*)h->Clone(Form("hVarScanClone%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      Double_t signal,errsignal;
      TF1 *f=FitHisto(hSel,particle,-1,1,signal,errsignal);
      f->SetLineColor(kBlue);
      TF1 *fBack=(TF1*)f->Clone("fback");
      fBack->SetParameter(0,0);
      fBack->SetLineStyle(2);
      fBack->SetLineColor(kGray);
      fBack->SetRange(minMassScFit,maxMassScFit);
      fBack->Draw("same");
      TF1 *fSign=(TF1*)f->Clone("fsign");
      fSign->SetParameter(3,0);
      Printf("Fit done");	  
      WriteFitResultsAndAxisOnPad(hSel,f,(TPad*)(canvOut->cd()),signal,errsignal);
      
      TCanvas *cRes=new TCanvas("cRes","cRes",800,800);
      cRes->cd();
      TH1D *hRes=new TH1D(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix),Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix),hSel->GetNbinsX(),hSel->GetBinLowEdge(1),hSel->GetBinLowEdge(hSel->GetNbinsX())+hSel->GetBinWidth(hSel->GetNbinsX()));
      //(TH1D*)hSel->Clone(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      //      TH1D *hRes=(TH1D*)hSel->Clone(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      for(Int_t kr=1;kr<=hRes->GetNbinsX();kr++){
	if(hRes->GetBinLowEdge(kr)>=minMassScFit && (hRes->GetBinLowEdge(kr)+hRes->GetBinWidth(kr))<maxMassScFit){
	  hRes->SetBinContent(kr,hSel->GetBinContent(kr)-fBack->Eval(hRes->GetBinCenter(kr)));//-fBack->Integral(hRes->GetBinLowEdge(kr),hRes->GetBinLowEdge(kr)+hRes->GetBinWidth(kr))/hRes->GetBinWidth(kr));
	  hRes->SetBinError(kr,hSel->GetBinError(kr));
	}
	else hRes->SetBinContent(kr,0);
	
      }
      hRes->Draw();
      fSign->SetRange(minMassScFit,maxMassScFit);
      fSign->Draw("same");
      WriteFitResultsAndAxisOnPad(hRes,f,(TPad*)(cRes->cd()),signal,errsignal);

      return;
    }

    for(Int_t iv1=0;iv1<5;iv1++){
      if(var1fix>=0 && var1fix!=iv1)continue;
      for(Int_t iv2=0;iv2<5;iv2++){
	if(var2fix>=0 && var2fix!=iv2)continue;
	Printf("Getting pad %d",iv2+iv1*5+1);
	TPad *p=(TPad*)c->GetPad(iv2+iv1*5+1);
	Printf(Form("getting histo hVarScan%sPt%d_i%d_j%d",strHist.Data(),j,iv2,iv1));
	TH1D *h=(TH1D*)p->FindObject(Form("hVarScan%sPt%d_i%d_j%d",strHist.Data(),j,iv2,iv1));
	Double_t signal,errsignal;
	if(refit){
	  TF1 *f=FitHisto(h,particle,-1,1,signal,errsignal);
	  Printf("Fit done");	  
	  WriteFitResultsAndAxisOnPad(h,f,p,signal,errsignal);
	}
	if(pickupPad)pickup=p;
      }
    }
  }
  if(pickupPad && !refit){
    canvOut=new TCanvas("cSigmaout","cSigmaout",800,800);
    pickup->Draw();
  }
  
}


void StudySBfits(TH1D **hSB,Int_t nbins){
  
  Double_t minFS[4]={0.148,0.140,0.150,0.148};
  Double_t maxFS[4]={0.210,0.220,0.200,0.224};
  TString strFunc[6]={"Power","dstar","customLog","pol2","pol3","pol4"};
  Double_t defminMassScFit=minMassScFit;
  Double_t defmaxMassScFit=maxMassScFit;
  TCanvas ***c=new TCanvas**[6];
  for(Int_t istrfunc=0;istrfunc<6;istrfunc++){
    c[istrfunc]=new TCanvas*[3];
    for(Int_t irange=0;irange<4;irange++){
      minMassScFit=minFS[irange];
      maxMassScFit=maxFS[irange];
      c[istrfunc][irange]=new TCanvas(Form("cCanvasNoCutsScSBsumF%sRange%d_%d",strFunc[istrfunc].Data(),(Int_t)(minFS[irange]*1000),(Int_t)(maxFS[irange]*1000)),Form("cCanvasNoCutsScSBsumF%sRange%d_%d",strFunc[istrfunc].Data(),(Int_t)(minFS[irange]*1000),(Int_t)(maxFS[irange]*1000)),800,800);
      DivideCanvas(c[istrfunc][irange],nbins);
      for(Int_t ipt=0;ipt<nbins;ipt++){
	Double_t signalSc=0,errsignalSc=-1.;//,significanceSc=0;
	TPad *p=(TPad*)c[istrfunc][irange]->cd(ipt+1);
	TH1D *h=(TH1D*)hSB[ipt]->Clone(Form("hSBtoSubtractLplusRbin%d_F%sRange%d_%d",ipt,strFunc[istrfunc].Data(),(Int_t)(minFS[irange]*1000),(Int_t)(maxFS[irange]*1000)));
	h->Draw();
	TF1 *f;
	if(istrfunc==0)f=FitHisto(h,1,-1,1,signalSc,errsignalSc,"REMLI");
	else if(istrfunc==1)f=FitHisto(h,1,5,1,signalSc,errsignalSc,"REMLI");
	else if(istrfunc==2)f=FitHisto(h,1,6,1,signalSc,errsignalSc,"REMLI");
	else if(istrfunc==3)f=FitHisto(h,1,2,1,signalSc,errsignalSc,"REMLI");//h->Fit("pol2","REMLI","",minMassScFit,maxMassScFit);//f=FitHisto(h,1,5,1,signalSc,errsignalSc,"REMLI");
	else if(istrfunc==4)f=FitHisto(h,1,3,1,signalSc,errsignalSc,"REMLI");//h->Fit("pol3","REMLI","",minMassScFit,maxMassScFit);//f=FitHisto(h,1,5,1,signalSc,errsignalSc,"REMLI");
	else if(istrfunc==5)f=FitHisto(h,1,4,1,signalSc,errsignalSc,"REMLI");//h->Fit("pol4","REMLI","",minMassScFit,maxMassScFit);//f=FitHisto(h,1,5,1,signalSc,errsignalSc,"REMLI");
	if(f)f->SetName(Form("fScNoCutsSBLRfit_bin%d_F%sRange%d_%d",ipt,strFunc[istrfunc].Data(),(Int_t)(minFS[irange]*1000),(Int_t)(maxFS[irange]*1000)));
      }
      c[istrfunc][irange]->Draw();
    }
  }
  minMassScFit=defminMassScFit;
  maxMassScFit=defmaxMassScFit;
  
  return;
}


void AnalyseOutput(Int_t readMC=0,Int_t tryFitData=0,Int_t particle=0,Int_t useCuts=1,Bool_t lowField=kFALSE,Int_t varscan0=3/*2=Lxy,3=nLxy,4=cosThetaP,5=nDCA*/,Int_t varscan1=4,Int_t caseSelForUseCuts2=-1,TString strfile="AnalysisResults_16_18_QM19.root",TString strDir="PWG3_D2H_XicpKpi"){
  const Int_t nbins=6;
  const Int_t nscansVar0=5; 
  const Int_t nscansVar1=5; 
  const Int_t nvar=7;
  const Int_t nbinsTree=11;
  Double_t lcMassCutDefault=0.007;
  Double_t sigmaCMincosthetaStarDefCut=-0.8;
  Double_t ptbins[nbins+1]={1.,2.,3.,5.,8,16,24.};
  Double_t ptbinsTree[nbinsTree+1]={0.,1.,2.,3.,4.,6.,8,10.,12.,16,20.,24.};
  Double_t lowEdgeVar0[nbins][nscansVar0];
  Double_t lowEdgeVar1[nbins][nscansVar1];
  Double_t upEdgeVar0[nbins][nscansVar0];
  Double_t upEdgeVar1[nbins][nscansVar1];
  Double_t varLowStand[nvar]={-1,-1,0.010,4. ,-1.1,0.,-1};// used only when cuts are scanned ; variabls: pt, mass, lxy, nlxy, thetaPoint, NDCA, SpecialSel
  Double_t varUpStand[nvar]={ -1,-1,1.,   99.,1.1 ,2.,99};// used only when cuts are scanned
  TString strSel="SpecSelNo";


  switch (caseSelForUseCuts2){
  
  case 0:
    if(varscan0==3){
      varLowStand[2]=0.0100;    
    }    
    if(particle==0 || particle ==1){
      varLowStand[3]=2.5;
      varUpStand[5]=3.5;// value for good preliminary so far: 2.5
    }
    break;
  case 1:
    if(varscan0==3){
      varLowStand[2]=0.0100;    
    }    
    if(particle==0 || particle ==1){
      varLowStand[3]=2.5;
      varUpStand[5]=3.5;// value for good preliminary so far: 2.5
    }
    varLowStand[6]=2.001;
  case -1:
  default:
    // values used for preliminary (current) versions
    if(varscan1==5){
      varLowStand[4]=0.94;
    }
    if(varscan0==3){
      varLowStand[2]=0.0150;    
    }    
    if(particle==0 || particle ==1){
      varLowStand[3]=3.5;
      varUpStand[5]=2.5;// value for good preliminary so far: 2.5
    }
    break;
  }
  


  ////////////////////////////////////////////////////////
  const Int_t nScanLxy=3;
  const Int_t nScanNLxy=3;
  const Int_t nScanThPth=3;
  const Int_t nScanNDCA=3;
  const Int_t nScanCsThetStar=6;
  const Int_t nScanLcMass=3;

  Double_t scanGlobalLxyMin[nScanLxy]={0.01,0.02,0.035};
  Double_t scanGlobalLxyMax[nScanLxy]={99.,99.,99.};
  Double_t scanGlobalNLxyMin[nScanNLxy]={1.,3.5,6.};
  if(particle==0 || particle==1){
    scanGlobalNLxyMin[1]=2.0;
    scanGlobalNLxyMin[2]=3.5;
  }
  //Double_t scanGlobalNLxyMin[nScanNLxy]={4.5,5.5,7.};
  Double_t scanGlobalNLxyMax[nScanNLxy]={99.,99.,99.};
  Double_t scanGlobalThetaPtMin[nScanThPth]={0.85,0.9,0.95};
  Double_t scanGlobalThetaPtMax[nScanThPth]={1.,1.,1.};
  Double_t scanGlobalNDCAMin[nScanNDCA]={0.,0.,0.};
  Double_t scanGlobalNDCAMax[nScanNDCA]={1.,2.,4.};

  Double_t scanGlobalCosThStarMin[nScanCsThetStar]={-1.,-0.85,-0.4,-0.2,0.,0.2};
  Double_t scanGlobalCosThStarMax[nScanCsThetStar]={1.,1.,1.};
  Double_t scanGlobalLcMass[nScanLcMass]={0.030,0.015,0.007};

//   const Int_t nScanLxy=6;
//   const Int_t nScanNLxy=6;
//   const Int_t nScanThPth=6;
//   const Int_t nScanNDCA=6;

//   Double_t scanGlobalLxyMin[nScanLxy]={0.,0.1,0.15,0.2,0.25,0.35};
//   Double_t scanGlobalLxyMax[nScanLxy]={99.,99.,99.,99.,99.,99.};
//   Double_t scanGlobalNLxyMin[nScanNLxy]={1.,2.,3.,4.,5.,6.};
//   Double_t scanGlobalNLxyMax[nScanNLxy]={99.,99.,99.,99.,99.,99.};
//   Double_t scanGlobalThetaPtMin[nScanThPth]={0.8,0.85,0.9,0.92,0.94,0.98};
//   Double_t scanGlobalThetaPtMax[nScanThPth]={1.,1.,1.,1.,1.,1.};
//   Double_t scanGlobalNDCAMin[nScanNDCA]={0.,0.,0.,0.,0.,0.};
//   Double_t scanGlobalNDCAMax[nScanNDCA]={1.,1.5,2.,2.5,3.,4.};


  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
  Double_t scanNlxyLowEdgeLc[nbins][5]={
    {0.,1.,1.5,2.,2.5},
    {0.,1.,1.5,2.,2.5},
    {1.,1.5,2.,2.5,3.},
    {1.,2.,3.,4.,5.},
    {1.,2.,3.,4.,5.},
    {1.,2.,3.,4.,5.}
  };
  Double_t scanNlxyUpEdgeLc[nbins][5]={
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.}
  };
  Double_t scanLxyLowEdgeLc[nbins][5]={
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
  };
  Double_t scanLxyUpEdgeLc[nbins][5]={
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.}
  };
  Double_t scanThPointLowEdgeLc[nbins][5]={
    {0.8,0.85,0.9,0.92,0.95},
    {0.8,0.85,0.9,0.92,0.95},
    {-1,0.8,0.82,0.85,0.9},
    {0.8,0.85,0.9,0.92,0.95},
    {0.8,0.85,0.9,0.92,0.95},
    {0.8,0.85,0.9,0.92,0.95},
  };
  Double_t scanThPointUpEdgeLc[nbins][5]={
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.}
  };

  Double_t scanNDCALowEdgeLc[nbins][5]={
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.}
  };
  Double_t scanNDCAUpEdgeLc[nbins][5]={
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.}
  };

  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
  Double_t scanNlxyLowEdgeXic[nbins][5]={
    {1.,3.,5.,7.,8.},
    {1.,3.,5.,7.,8.},
    {1.,3.,5.,7.,8.},
    {1.,3.,5.,7.,8.},
    {1.,3.,5.,7.,8.},
    {1.,3.,5.,7.,8.}   
  };
  Double_t scanNlxyUpEdgeXic[nbins][5]={
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.}
  };
  Double_t scanLxyLowEdgeXic[nbins][5]={
    {0.0,0.01,0.02,0.035,0.05},
    {0.0,0.01,0.02,0.035,0.05},
    {0.0,0.01,0.02,0.035,0.05},
    {0.0,0.01,0.02,0.035,0.05},    
    {0.0,0.01,0.02,0.035,0.05},
    {0.0,0.01,0.02,0.035,0.05}
  };
  Double_t scanLxyUpEdgeXic[nbins][5]={
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.}
  };

  Double_t scanThPointLowEdgeXic[nbins][5]={
    {0.85,0.9,0.94,0.97,0.99},
    {0.85,0.9,0.94,0.97,0.99},
    {0.9,0.92,0.94,0.97,0.99},
    {0.9,0.92,0.94,0.97,0.99},
    {0.9,0.92,0.94,0.97,0.99},
    {0.9,0.92,0.94,0.97,0.99}
  };
  Double_t scanThPointUpEdgeXic[nbins][5]={
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.}
  };

 Double_t scanNDCALowEdgeXic[nbins][5]={
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.}
  };
  Double_t scanNDCAUpEdgeXic[nbins][5]={
    {0.5,1.,1.5,2.5,3.5},
    {0.5,1.,1.5,2.5,3.5},
    {0.5,1.,1.5,2.5,3.5},
    {0.5,1.,1.5,2.5,3.5},
    {0.5,1.,1.5,2.5,3.5},
    {0.5,1.,1.5,2.5,3.5}
  };




  for(Int_t k=0;k<nbins;k++){
    for(Int_t j=0;j<nscansVar0;j++){
      if(varscan0==2){
	if(particle==0 || particle==1){
	  lowEdgeVar0[k][j]=scanLxyLowEdgeLc[k][j];	
	  upEdgeVar0[k][j]=scanLxyUpEdgeLc[k][j];	
	}
	else if(particle==2){
	  lowEdgeVar0[k][j]=scanLxyLowEdgeXic[k][j];	
	  upEdgeVar0[k][j]=scanLxyUpEdgeXic[k][j];	
	}
      }
      if(varscan0==3){
	if(particle==0 || particle==1){
	  lowEdgeVar0[k][j]=scanNlxyLowEdgeLc[k][j];	
	  upEdgeVar0[k][j]=scanNlxyUpEdgeLc[k][j];	
	}
	else if(particle==2){
	  lowEdgeVar0[k][j]=scanNlxyLowEdgeXic[k][j];	
	  upEdgeVar0[k][j]=scanNlxyUpEdgeXic[k][j];	
	}
      }
    }
    for(Int_t j=0;j<nscansVar1;j++){
      if(varscan1==4){	
	if(particle==0 || particle==1){
	lowEdgeVar1[k][j]=scanThPointLowEdgeLc[k][j];	
	upEdgeVar1[k][j]=scanThPointUpEdgeLc[k][j];	
	}
	else if(particle==2){
	lowEdgeVar1[k][j]=scanThPointLowEdgeXic[k][j];	
	upEdgeVar1[k][j]=scanThPointUpEdgeXic[k][j];	
	}
      }
      if(varscan1==5){	
	if(particle==0 || particle==1){
	lowEdgeVar1[k][j]=scanNDCALowEdgeLc[k][j];	
	upEdgeVar1[k][j]=scanNDCAUpEdgeLc[k][j];	
	}
	else if(particle==2){
	lowEdgeVar1[k][j]=scanNDCALowEdgeXic[k][j];	
	upEdgeVar1[k][j]=scanNDCAUpEdgeXic[k][j];	
	}
      }
    }
  }
  
  
  TFile *f=TFile::Open(strfile.Data(),"READ");
  /*TDirectory *ddir=f->GetDirectory(strDir.Data());
  TList *lis=(TList*)ddir->Get("outputList");*/
  THnSparseF *hsparse;
  if(particle==0 || particle==2 )hsparse=(THnSparseF*)f->Get("fhSparseAnalysis");
  else if(particle==1)hsparse=(THnSparseF*)f->Get("fhSparseAnalysisSigma");
  TTree *treeVar;
  TH1D ***hVarSign=new TH1D**[17];
  TH1D ***hVarBack=new TH1D**[17];
  TCanvas **cVar=new TCanvas*[17];
  TCanvas *cPtTree;
  TH1D *hPtDistrTree;
  Float_t var[18];
  TString varNames[18]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","flagMC"};
  Double_t rangeL[18]={0.,-1,0,0,0,0,0,0,0,0,-0.1,-0.1,-0.1,-8,-8,-8,-0.1,-2};
  Double_t rangeU[18]={50,1,0.2,20,25,25,25,5,0.4,0.01,0.1,0.1,0.1,8,8,8,0.1,14};
  Double_t nCandSignal[nbinsTree];
  Double_t nCandBack[nbinsTree];


  if(readMC==2){
    treeVar=(TTree*)f->Get("T");
    for(Int_t j=0;j<18;j++){      
      treeVar->SetBranchAddress(varNames[j].Data(),&var[j]);
      if(j==0){
	hPtDistrTree=new TH1D("hPtDistrTree","hPtDistrTree",200,0,40.);
	cPtTree=new TCanvas("cPtTree","cPtTree",800,800);
      }
      else {
	cVar[j-1]=new TCanvas(varNames[j].Data(),varNames[j].Data(),800,800);
	cVar[j-1]->Divide(3,4);
	hVarSign[j-1]=new TH1D*[nbinsTree];
	hVarBack[j-1]=new TH1D*[nbinsTree];
	for(Int_t jpt=0;jpt<nbinsTree;jpt++){
	  hVarSign[j-1][jpt]=new TH1D(Form("Sign_%s_bin%d",varNames[j].Data(),jpt),Form("Sign_%s, %.1f<pt<%.1f",varNames[j].Data(),ptbinsTree[jpt],ptbinsTree[jpt+1]),200,rangeL[j],rangeU[j]);
	  hVarBack[j-1][jpt]=new TH1D(Form("Back_%s_bin%d",varNames[j].Data(),jpt),Form("Back_%s, %.1f<pt<%.1f",varNames[j].Data(),ptbinsTree[jpt],ptbinsTree[jpt+1]),200,rangeL[j],rangeU[j]);
	  nCandSignal[jpt]=0.;
	  nCandBack[jpt]=0.;
	}
      }
    }
    Printf("Tree and histos created");

    for(Int_t itr=0;itr<treeVar->GetEntries();itr++){
      treeVar->GetEntry(itr);
      Int_t binNumber=0;
      
      for(Int_t ib=0;ib<nbinsTree;ib++){
	if(var[0]>ptbinsTree[ib])binNumber=ib;
	else break;
      }
      for(Int_t j=0;j<18;j++){
	if(j==0)hPtDistrTree->Fill(var[0]);
	else {
	  if(var[17]>=1){
	    if(var[2]>=0.01&&var[1]>0.8&&var[7]<1.5){
	      hVarSign[j-1][binNumber]->Fill(var[j]);
	      nCandSignal[binNumber]++;
	    }
	  }
	  else{	    
	    if(var[2]>=0.01&&var[1]>0.8&&var[7]<1.5){
	      nCandBack[binNumber]++;
	      hVarBack[j-1][binNumber]->Fill(var[j]);
	    }
	  }
	}
      }
    }
    cPtTree->cd();
    hPtDistrTree->Draw();
    for(Int_t j=0;j<17;j++){
      for(Int_t ib=0;ib<nbinsTree;ib++){
	cVar[j]->cd(ib+1);
	hVarSign[j][ib]->Sumw2();
	hVarSign[j][ib]->Scale(1./nCandSignal[ib]);//hVarSign[j][ib]->Integral());
	hVarSign[j][ib]->SetLineColor(kBlack);
	hVarSign[j][ib]->SetLineWidth(2);
	hVarSign[j][ib]->Draw();
	hVarBack[j][ib]->Sumw2();
	hVarBack[j][ib]->Scale(1./nCandBack[ib]);//hVarBack[j][ib]->Integral());
	hVarBack[j][ib]->SetLineColor(kRed);
	hVarBack[j][ib]->SetLineWidth(2);
	hVarBack[j][ib]->Draw("sames");	
      }
    }
    return;
  }
  TAxis *axpt=hsparse->GetAxis(0); // axes: pt;mass;lxy;nLxy;cosThatPoint;normImpParXY; Xic mass: 2.468, Lc:2.286
  TAxis *axLxy=hsparse->GetAxis(2);
  TAxis *axNlxy=hsparse->GetAxis(3);
  TAxis *axPoint=hsparse->GetAxis(4);
  TAxis *axNDCA=hsparse->GetAxis(5);
  TAxis *axSel=hsparse->GetAxis(6);
  TAxis *axSpecialPIDs=hsparse->GetAxis(7);
  axSpecialPIDs->SetRangeUser(0,0);
  TAxis *axVar0=hsparse->GetAxis(varscan0);
  TAxis *axVar1=hsparse->GetAxis(varscan1);
  TString strvar[7]={"pt","mass","lxy","nLxy","cosThetaP","nDCA","SelSpec"};

  TAxis *axMassSigma=0x0;
  TAxis *axMassLcSigma=0x0;
  TAxis *axMassCosThetaStar=0x0;
  if(particle==1){
    //    axMassSigma=hsparse->GetAxis(6);
    if(version>0){
      axMassLcSigma=hsparse->GetAxis(8);
      axMassCosThetaStar=hsparse->GetAxis(9);
    }
    else {
      axMassLcSigma=hsparse->GetAxis(7);
      axMassCosThetaStar=hsparse->GetAxis(8);
    }
  }

  // DO IMMEDIATELY CANVAS W/O CUTS
  TCanvas *cCanvasNoCuts=new TCanvas("cCanvasNoCuts","cCanvasNoCuts",800,800);
  if(nbins<=4)cCanvasNoCuts->Divide(2,2);
  else if(nbins<=6)cCanvasNoCuts->Divide(2,3);
  else if(nbins<=9)cCanvasNoCuts->Divide(3,3);
  else if(nbins<=12)cCanvasNoCuts->Divide(4,3);
  else if(nbins<=16)cCanvasNoCuts->Divide(4,4);
  else cCanvasNoCuts->Divide(5,5);

  for(Int_t ipt=0;ipt<nbins;ipt++){
    cCanvasNoCuts->cd(ipt+1);
    axpt->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    TH1D *h=hsparse->Projection(1);
    h->SetName(Form("hNoCuts_PtBin%d",ipt));
    h->SetTitle(Form("hNoCuts %.1f<pt<%.1f GeV/c",ptbins[ipt],ptbins[ipt+1]));    
    h->Draw();;    
  }
  cCanvasNoCuts->Draw();


  TH1D **hmassData=new TH1D*[nbins];
  TH1D **hCosThetaStarSigma;
  if(particle==1)hCosThetaStarSigma=new TH1D*[nbins];
  TH1D **hmassMCacc=new TH1D*[nbins];
  
  TH1D *hYieldData=new TH1D("hYieldData","hYieldData",nbins,ptbins);
  TH1D *hYieldMCacc=new TH1D("hYieldMCacc","hYieldMCacc",nbins,ptbins);
  TH1D *hEffRecoGenAcc=new TH1D("hEffRecoGenAcc","hEffRecoGenAcc",nbins,ptbins);

  
  //   Double_t cutsMinLc[nbins][6]={/*boni* DEVI CAMBIARE LXY!!!/
  //     {2.1,ptbins[0],0.,-1.1,0.,-2.},
  //     {2.1,ptbins[1],2.5,0.81,0.,-2.},
  //     {2.1,ptbins[2],2.5,0.81,0.,-2.},
  //     {2.1,ptbins[3],2.5,0.81,0.,-1.},
  //     {2.1,ptbins[4],2.5,0.81,0.,-1.},
  //     {2.1,ptbins[5],2.5,0.81,0.,-1.}  
  //   };
  
  //   Double_t cutsMinLc[nbins][6]={/*boni? DEVI CAMBIARE LXY!!!*/
  //     {2.1,ptbins[0],0.,-1.1,0.,-2.},
  //     {2.1,ptbins[1],3.,0.81,0.,-2.},
  //     {2.1,ptbins[2],3.,0.81,0.,-2.},
  //     {2.1,ptbins[3],2.5,0.81,0.,-1.},
  //     {2.1,ptbins[4],2.5,0.81,0.,-1.},
  //     {2.1,ptbins[5],2.5,0.81,0.,-1.}  
  //   };
  
  //   Double_t cutsMinLc[nbins][7]={/*boni?*/
  //     {2.1,ptbins[0],0.005,2.,0.83,0.,-2.},
  //     {2.1,ptbins[1],0.005,2.,0.85,0.,-2.},
  //     {2.1,ptbins[2],0.005,3.,0.85,0.,-2.},
  //     {2.1,ptbins[3],0.0100,2.5,0.88,0.,-1.},
  //     {2.1,ptbins[4],0.0100,2.5,0.88,0.,-1.},
  //     {2.1,ptbins[5],0.0100,2.5,0.88,0.,-1.}  
  //   };
  
  
  //   Double_t cutsMaxLc[nbins][7]={
  //     {2.65,ptbins[1],99.,99,1.1, 1.1,99.},
  //     {2.65,ptbins[2],99.,99,1.1, 1.,99.},
//     {2.65,ptbins[3],99.,99,1.1, 1.,99.},
//     {2.65,ptbins[4],99.,99,1.1, 2.5,99.},
//     {2.65,ptbins[5],99.,99,1.1, 2.5,99.},
//     {2.65,ptbins[6],99.,99,1.1, 2.5,99.}  
//   };

//   Double_t cutsMinLc[nbins][7]={/*durezza intermedia*/
//     {2.1,ptbins[0],0.015,1.,0.8,0.,-2.},
//     {2.1,ptbins[1],0.02,1.5,0.8,0.,-2.},
//     {2.1,ptbins[2],0.02,1.5,0.8,0.,-2.},
//     {2.1,ptbins[3],0.0250,1.5,0.85,0.,-1.},
//     {2.1,ptbins[4],0.0250,1.5,0.9,0.,-1.},
//     {2.1,ptbins[5],0.0250,2.,0.95,0.,-1.}  
//   };

  Double_t cutsMinLc[nbins][7]={/*piu tight buoni per Sigma_C*/
    {2.1,ptbins[0],0.015,1.,0.8,0.,-2.},
    {2.1,ptbins[1],0.025,2.,0.8,0.,-2.},
    {2.1,ptbins[2],0.03,3.,0.8,0.,-2.},
    {2.1,ptbins[3],0.040,3.,0.85,0.,-1.},
    {2.1,ptbins[4],0.0350,3.,0.9,0.,-1.},
    {2.1,ptbins[5],0.0350,3.,0.95,0.,-1.}  
  };


  Double_t cutsMaxLc[nbins][7]={
    {2.65,ptbins[1],99.,99,1.1, 99.1,99.},
    {2.65,ptbins[2],99.,99,1.1, 99.,99.},
    {2.65,ptbins[3],99.,99,1.1, 99.,99.},
    {2.65,ptbins[4],99.,99,1.1, 99.5,99.},
    {2.65,ptbins[5],99.,99,1.1, 99.5,99.},
    {2.65,ptbins[6],99.,99,1.1, 99.5,99.}  
  };


 //  Double_t cutsMinXic[nbins][6]={/*boni*/ DEVI CAMBIARE LXY!! AGGIUNGI ASSE
//     {2.1,ptbins[0],1.5,0.93,0.,-2.},
//     {2.1,ptbins[1],3.5,0.93,0.,-2.},
//     {2.1,ptbins[2],3.5,0.91,0.,-2.},
//     {2.1,ptbins[3],2.5,0.9,0.,-1.},
//     {2.1,ptbins[4],3.5,0.92,0.,-1.},
//     {2.1,ptbins[5],4.5,0.92,0.,-1.}  
//   };
//   Double_t cutsMinXic[nbins][6]={ /*pure boni*/
//     {2.1,ptbins[0],1.5,0.93,0.,-2.},
//     {2.1,ptbins[1],3.5,0.93,0.,-2.},
//     {2.1,ptbins[2],3.5,0.95,0.,-2.},
//     {2.1,ptbins[3],2.5,0.95,0.,-1.},
//     {2.1,ptbins[4],3.5,0.95,0.,-1.},
//     {2.1,ptbins[5],4.5,0.95,0.,-1.}  
//     };



//   Double_t cutsMinXic[nbins][7]={
//     {2.1,ptbins[0],0.02,1.5,0.95,0.,-2.},
//     {2.1,ptbins[1],0.01,3.5,0.93,0.,-2.},/*    {2.1,ptbins[1],0.01,3.5,0.9,0.,-2.},*/
//     {2.1,ptbins[2],0.015,3.5,0.95,0.,-2.},
//     {2.1,ptbins[3],0.020,3.5,0.9,0.,-1.},
//     {2.1,ptbins[4],0.0200,3.5,0.9,0.,-1.},
//     {2.1,ptbins[5],0.0200,3.5,0.9,0.,-1.}  
//   };

//   Double_t cutsMaxXic[nbins][7]={
//     {2.65,ptbins[1],99,99,1.1, 2.5,99.},
//     {2.65,ptbins[2],99,99,1.1, 3.,99.},
//     {2.65,ptbins[3],99,99,1.1, 3.,99.},
//     {2.65,ptbins[4],99,99,1.1, 3.5,99.},
//     {2.65,ptbins[5],99,99,1.1, 4.5,99.},
//     {2.65,ptbins[6],99,99,1.1, 4.5,99.}  
//   };

  Double_t cutsMinXic[nbins][7]={
    {2.1,ptbins[0],0.025,5.,0.94,0.,-2.},
    {2.1,ptbins[1],0.025,5,0.94,0.,-2.},/*    {2.1,ptbins[1],0.01,3.5,0.9,0.,-2.},*/
    {2.1,ptbins[2],0.025,5,0.94,0.,-2.},
    {2.1,ptbins[3],0.0250,5,0.94,0.,-1.},
    {2.1,ptbins[4],0.0250,5,0.94,0.,-1.},
    {2.1,ptbins[5],0.0250,5,0.94,0.,-1.}  
  };

  Double_t cutsMaxXic[nbins][7]={
    {2.65,ptbins[1],99,99,1.1, 1.,99.},
    {2.65,ptbins[2],99,99,1.1, 1.,99.},
    {2.65,ptbins[3],99,99,1.1, 1.,99.},
    {2.65,ptbins[4],99,99,1.1, 1.,99.},
    {2.65,ptbins[5],99,99,1.1, 1.,99.},
    {2.65,ptbins[6],99,99,1.1, 1.5,99.}  
  };

  
  //   axLxy->SetRangeUser(3.,99.);
  //   axPoint->SetRangeUser(0.9,0.9999);
  //   axNDCA->SetRangeUser(0.,1.5);
  //   axSel->SetRangeUser(1.1,8.1);
  TCanvas *cCanvasCosThetaStar;
  if(particle==1){
    cCanvasCosThetaStar=new TCanvas("cCanvasCosThetaStar","cCanvasCosThetaStar",800,800);
    cCanvasCosThetaStar->Divide(3,3);
}
  TCanvas *cAllBins=new TCanvas("cAllBins","cAllBins",800,800);
  cAllBins->Divide(3,3);
  for(Int_t ipt=0;ipt<nbins;ipt++){
    cAllBins->cd(ipt+1);
    axpt->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    //     axLxy->SetRangeUser(2.5,99);
    //     axPoint->SetRangeUser(0.81,1.1);
    //     axNDCA->SetRangeUser(0.,3.5);
    //     axSel->SetRangeUser(1.1,99);
    if(useCuts){
      if(particle==0 || particle ==1){	
	axLxy->SetRangeUser(cutsMinLc[ipt][2]*1.0001,cutsMaxLc[ipt][2]*0.9999);
	axNlxy->SetRangeUser(cutsMinLc[ipt][3]*1.0001,cutsMaxLc[ipt][3]*0.9999);
	axPoint->SetRangeUser(cutsMinLc[ipt][4]*1.0001,cutsMaxLc[ipt][4]*0.9999);
	axNDCA->SetRangeUser(cutsMinLc[ipt][5]*1.0001,cutsMaxLc[ipt][5]*0.9999);
	axSel->SetRangeUser(cutsMinLc[ipt][6]*1.0001,cutsMaxLc[ipt][6]*0.9999);
	if(particle==1){
	  if(!lowField)axMassLcSigma->SetRangeUser(2.286-0.015,2.286+0.015);
	  else axMassLcSigma->SetRangeUser(2.286-0.025,2.286+0.025);
	  axMassCosThetaStar->SetRangeUser(-0.8,1);
	}
      }
      else if(particle==2){
	axLxy->SetRangeUser(cutsMinXic[ipt][2]*1.0001,cutsMaxXic[ipt][2]*0.9999);
	axNlxy->SetRangeUser(cutsMinXic[ipt][3]*1.0001,cutsMaxXic[ipt][3]*0.9999);
	axPoint->SetRangeUser(cutsMinXic[ipt][4]*1.0001,cutsMaxXic[ipt][4]*0.9999);
	axNDCA->SetRangeUser(cutsMinXic[ipt][5]*1.0001,cutsMaxXic[ipt][5]*0.9999);
	axSel->SetRangeUser(cutsMinXic[ipt][6]*1.0001,cutsMaxXic[ipt][6]*0.9999);
      }
    }
    else if(particle==1){
      axLxy->SetRangeUser(0.0150,99.);
      axPoint->SetRangeUser(0.082,1.);
      if(!lowField)axMassLcSigma->SetRangeUser(2.286-0.007,2.286+0.007);
      else axMassLcSigma->SetRangeUser(2.286-0.025,2.286+0.025);     
      axMassCosThetaStar->SetRangeUser(-0.8,1);
    }
    hmassData[ipt]=hsparse->Projection(1);
    hmassData[ipt]->SetName(Form("hmassData%d",ipt));
    hmassData[ipt]->Sumw2();
    if(particle==0)hmassData[ipt]->Rebin(2);
    hmassData[ipt]->Draw();

    if(useCuts&&particle==1){
      axLxy->SetRangeUser(0.,99.);
      axNlxy->SetRangeUser(0.,99.);
      axPoint->SetRangeUser(-1,1);
      axNDCA->SetRangeUser(0,99);
      axSel->SetRangeUser(-2,14);
      axMassCosThetaStar->SetRangeUser(-1.1,1.1);
      axMassLcSigma->SetRangeUser(-1,3);//2.286-0.010,2.286+0.010);
      hCosThetaStarSigma[ipt]=hsparse->Projection(8);
      hCosThetaStarSigma[ipt]->SetName(Form("hCosThetaStarSigma_%d",ipt));
      cCanvasCosThetaStar->cd(ipt+1);
      hCosThetaStarSigma[ipt]->Draw();      
    }
  }

  // Global desperate scan
  TString strpt;
  // Merge files
  if(useCuts==6){// Start global desperate scan
    TCanvas *cScan=new TCanvas("cScanSum","cScanSum",800,800);
    for(Int_t ipt=0;ipt<nbins;ipt++){
      strpt=Form("Pt%.0fto99",ptbins[ipt]);
      gROOT->ProcessLine(Form(".!mkdir %s",strpt.Data()));      
      for(Int_t ilxy=0;ilxy<nScanLxy;ilxy++){
	for(Int_t inlxy=0;inlxy<nScanNLxy;inlxy++){
	  for(Int_t iCosthPth=0;iCosthPth<nScanThPth;iCosthPth++){
	    for(Int_t indca=0;indca<nScanNDCA;indca++){
	      for(Int_t iLcMass=0;iLcMass<nScanLcMass;iLcMass++){
		if(iLcMass>0 && particle!=1)continue;		  
		for(Int_t iCsThetStar=0;iCsThetStar<nScanCsThetStar;iCsThetStar++){
		  if(iCsThetStar>0 && particle!=1)continue;
		  TH1D *hSum=0x0;
		  TString strScanSumName(Form("%s_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strpt.Data(),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data()));
		  if(particle==1){
		    if(scanGlobalCosThStarMin[iCsThetStar]<0){
		      strScanSumName=Form("%s_LcMass%.0f_CosThStMinMinus%.0f_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strpt.Data(),scanGlobalLcMass[iLcMass]*1.e+3,TMath::Abs(scanGlobalCosThStarMin[iCsThetStar]*100),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data());
		    }
		    else{
		      strScanSumName=Form("%s_LcMass%.0f_CosThStMinPlus%.0f_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strpt.Data(),scanGlobalLcMass[iLcMass]*1.e+3,TMath::Abs(scanGlobalCosThStarMin[iCsThetStar]*100),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data());
		    }
		  }
		  for(Int_t jpt=0;jpt<nbins;jpt++){
		    TString strptFileDiff=Form("Pt%.0fto%.0f",ptbins[jpt],ptbins[jpt+1]);		
		    TString strScanName(Form("%s_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strptFileDiff.Data(),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data()));
		    if(particle==1){
		      if(scanGlobalCosThStarMin[iCsThetStar]<0){
			strScanName=Form("%s_LcMass%.0f_CosThStMinMinus%.0f_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strptFileDiff.Data(),scanGlobalLcMass[iLcMass]*1.e+3,TMath::Abs(scanGlobalCosThStarMin[iCsThetStar]*100),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data());
		      }
		      else{
			strScanName=Form("%s_LcMass%.0f_CosThStMinPlus%.0f_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strptFileDiff.Data(),scanGlobalLcMass[iLcMass]*1.e+3,TMath::Abs(scanGlobalCosThStarMin[iCsThetStar]*100),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data());
		      }
		    }
		    Printf("Opening file: %s",Form("%s/%s.root.root",strptFileDiff.Data(),strScanName.Data()));
		    TFile *fcase=TFile::Open(Form("%s/%s.root.root",strptFileDiff.Data(),strScanName.Data()),"READ");
		    TCanvas *c=(TCanvas*)fcase->Get("cScan");
		    TH1D *h=(TH1D*)c->FindObject(Form("%s.root",strScanName.Data()));
		    Printf("Adding histo %p to hSum %p",h,hSum);
		    if(hSum)hSum->Add(h);
		    else {
		      hSum=new TH1D(*h);
		      hSum->SetName(strScanSumName.Data());
		      hSum->SetTitle(strScanSumName.Data());
		    }
		    fcase->Close();
		  }
		  cScan->cd();
		  hSum->Draw();
		  cScan->SaveAs(Form("%s/%s.root",strpt.Data(),strScanSumName.Data()));
		  cScan->SaveAs(Form("%s/%s.png",strpt.Data(),strScanSumName.Data()));
		  cScan->SaveAs(Form("%s/%s.eps",strpt.Data(),strScanSumName.Data()));		    
		  delete hSum;
		}
	      }
	    }
	  }
	}      
      }    
    }
    return;
  }
  
  if(useCuts==4 || useCuts==5){

    TCanvas *cScan=new TCanvas("cScan","cScan",800,800);

    for(Int_t ipt=0;ipt<nbins;ipt++){
      if(useCuts==4){
	strpt=Form("Pt%.0fto%.0f",ptbins[ipt],ptbins[ipt+1]);
	gROOT->ProcessLine(Form(".!mkdir %s",strpt.Data()));
	axpt->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
      }
      else if(useCuts==5){
	strpt=Form("Pt%.0fto99",ptbins[ipt]);
	gROOT->ProcessLine(Form(".!mkdir %s",strpt.Data()));
	axpt->SetRangeUser(ptbins[ipt]*1.001,999.);
      }

      for(Int_t iv=2;iv<7;iv++){
	Double_t temp=varLowStand[2];
	if((iv==2)&&(particle==0 || particle==1)&&(ipt<2)){
	  varLowStand[2]=0.;
	}
	hsparse->GetAxis(iv)->SetRangeUser(varLowStand[iv]*1.0001,varUpStand[iv]*0.999);    
	varLowStand[2]=temp;
      }
      if(particle==1){
	if(!lowField)axMassLcSigma->SetRangeUser(2.286-0.007,2.286+0.007);
	else axMassLcSigma->SetRangeUser(2.286-0.025,2.286+0.025);
	axMassCosThetaStar->SetRangeUser(-0.5,1);
      }
      
      for(Int_t ilxy=0;ilxy<nScanLxy;ilxy++){
	axLxy->SetRangeUser(scanGlobalLxyMin[ilxy]*1.0001,scanGlobalLxyMax[ilxy]*0.9999);
	for(Int_t inlxy=0;inlxy<nScanNLxy;inlxy++){
	  axNlxy->SetRangeUser(scanGlobalNLxyMin[inlxy]*1.0001,scanGlobalNLxyMax[inlxy]*0.9999);
	  for(Int_t iCosthPth=0;iCosthPth<nScanThPth;iCosthPth++){
	    axPoint->SetRangeUser(scanGlobalThetaPtMin[iCosthPth],scanGlobalThetaPtMax[iCosthPth]);
	    for(Int_t indca=0;indca<nScanNDCA;indca++){
	      axNDCA->SetRangeUser(scanGlobalNDCAMin[indca],scanGlobalNDCAMax[indca]);          
	      TString strScanName(Form("%s_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strpt.Data(),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data()));
	      for(Int_t iLcMass=0;iLcMass<nScanLcMass;iLcMass++){
		if(iLcMass>0 && particle!=1)continue;
		if(particle==1)axMassLcSigma->SetRangeUser(2.286-scanGlobalLcMass[iLcMass],2.286+scanGlobalLcMass[iLcMass]);
		for(Int_t iCsThetStar=0;iCsThetStar<nScanCsThetStar;iCsThetStar++){
		  if(iCsThetStar>0 && particle!=1)continue;
		  if(particle==1){
		    axMassCosThetaStar->SetRangeUser(scanGlobalCosThStarMin[iCsThetStar],scanGlobalCosThStarMax[iCsThetStar]);
		    if(scanGlobalCosThStarMin[iCsThetStar]<0){
		      strScanName=Form("%s_LcMass%.0f_CosThStMinMinus%.0f_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strpt.Data(),scanGlobalLcMass[iLcMass]*1.e+3,TMath::Abs(scanGlobalCosThStarMin[iCsThetStar]*100),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data());
		    }
		    else {
		      strScanName=Form("%s_LcMass%.0f_CosThStMinPlus%.0f_Lxy%.0f_nLxy%.1f_thPt%.2f_ndca%.1f_%s",strpt.Data(),scanGlobalLcMass[iLcMass]*1.e+3,TMath::Abs(scanGlobalCosThStarMin[iCsThetStar]*100),scanGlobalLxyMin[ilxy]*1.e+4,scanGlobalNLxyMin[inlxy],scanGlobalThetaPtMin[iCosthPth],scanGlobalNDCAMax[indca],strSel.Data());
		    }
		  }
		  
		  cScan->cd();
		  TH1D *hscan=hsparse->Projection(1);
		  hscan->SetName(strScanName.Data());
		  hscan->SetTitle(strScanName.Data());
		  hscan->Sumw2();
		  //	      hscan->Rebin(2);	      
		  hscan->Draw();
		  cScan->SaveAs(Form("%s/%s.root",strpt.Data(),strScanName.Data()));
		  cScan->SaveAs(Form("%s/%s.png",strpt.Data(),strScanName.Data()));
		  cScan->SaveAs(Form("%s/%s.eps",strpt.Data(),strScanName.Data()));
		  delete hscan;
		}
	      }
	    }
	  }
	}      
      }    
    }
  }
  
  // THIS IS THE PART WITH SCAN OF CUTS DEFINED BY THE ARRAYES varLowStand and varUpStand
  TCanvas **cCutScan;
  if(useCuts==2||useCuts==3){
    cCutScan=new TCanvas*[nbins];
    for(Int_t iv=2;iv<7;iv++){
      hsparse->GetAxis(iv)->SetRangeUser(varLowStand[iv]*1.0001,varUpStand[iv]*0.999);     
    }
    if(particle==1){

      if(!lowField)axMassLcSigma->SetRangeUser(2.286-lcMassCutDefault,2.286+lcMassCutDefault);
      else axMassLcSigma->SetRangeUser(2.286-0.025,2.286+0.025);
      axMassCosThetaStar->SetRangeUser(sigmaCMincosthetaStarDefCut,1);
    }
    for(Int_t jpt=0;jpt<nbins;jpt++){
      cCutScan[jpt]=new TCanvas(Form("cCutScan%d",jpt),Form("cCutScan%d",jpt),800,800);
      cCutScan[jpt]->Divide(nscansVar0,nscansVar1);     
      if(useCuts==2)axpt->SetRangeUser(ptbins[jpt]*1.001,ptbins[jpt+1]*0.999);
      else if(useCuts==3)axpt->SetRangeUser(ptbins[jpt]*1.001,999);
      cCutScan[jpt]->Draw();     
      for(Int_t ivar0=0;ivar0<nscansVar0;ivar0++){
	axVar0->SetRangeUser(lowEdgeVar0[jpt][ivar0]*1.0001,upEdgeVar0[jpt][ivar0]*0.9999);
	for(Int_t ivar1=0;ivar1<nscansVar1;ivar1++){
	  cCutScan[jpt]->cd(ivar1*nscansVar0+ivar0+1);
	  //	 Printf("Setting values %f,%f to axis %s",lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999,axVar1->GetName());
	  axVar1->SetRangeUser(lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999);
	  TH1D *h=hsparse->Projection(1);
	  h->Sumw2();
	  h->SetName(Form("hVarScanPt%d_i%d_j%d",jpt,ivar0,ivar1));
	  h->SetTitle(Form("Pt%d, %s: %.2f - %.2f, %s: %.2f-%.2f",jpt,strvar[varscan0].Data(),lowEdgeVar0[jpt][ivar0],upEdgeVar0[jpt][ivar0],strvar[varscan1].Data(),lowEdgeVar1[jpt][ivar1],upEdgeVar1[jpt][ivar1]));	 
	  h->Draw();
	  if(tryFitData==1){
	    Double_t signal,errsignal;
	    TF1 *f=FitHisto(h,particle,-1,1,signal,errsignal);
	    f->SetName(Form("f%d_%d_%d",jpt,ivar0,ivar1));
	  }
	}
      }
    }
  }
  
  if(readMC){
    if(tryFitData==2){
      for(Int_t ipt=0;ipt<nbins;ipt++){
	Double_t signal,errsignal;
	TF1 *f=FitHisto(hmassData[ipt],particle,-1,1,signal,errsignal);
	f->SetName(Form("f%d",ipt));
	// 	cAllBins->cd(ipt+1);
	// 	f->SetLineColor(kRed);
	// 	f->SetLineStyle(1);
	// 	Printf("F at 2.2: %f",f->Eval(2.2));
	// 	f->Draw("same");
	hYieldData->SetBinContent(ipt+1,signal);
	hYieldData->SetBinError(ipt+1,errsignal);
      }
    }
    TH2F *hMC=(TH2F*)f->Get("fhistMCSpectrumAccLc");
    TH1D *hGenLimAcc=hMC->ProjectionX("hGenLimAcc",2,2);
    TH1D *hGenAcc=hMC->ProjectionX("hGenAcc",3,3);
    TH1D *hRecoAcc=hMC->ProjectionX("hRecoCuts",8,8);
    TH1D *hRecoAccFilt=hMC->ProjectionX("hRecoCutsFiltering",13,13);
    TCanvas *cMC=new TCanvas("cMC","cMC",800,800);
    cMC->cd();
    hGenAcc->SetLineColor(kRed);
    hGenLimAcc->SetLineColor(kBlack);
    hGenLimAcc->Draw();
    hGenAcc->Draw("same");
  
    TAxis *axptMC=hGenLimAcc->GetXaxis();    
    for(Int_t ipt=0;ipt<nbins;ipt++){
      Int_t bl=axptMC->FindBin(ptbins[ipt]*1.001);
      Int_t bu=axptMC->FindBin(ptbins[ipt+1]*0.999);
      if(TMath::Abs(ptbins[ipt]-axptMC->GetBinLowEdge(bl))>1.e-3){
	Printf("Bin Mismatch (low edge %d), cannot calculate efficiency",ipt);
	return;
      }
      if(TMath::Abs(ptbins[ipt+1]-axptMC->GetBinUpEdge(bu))>1.e-6){
	Printf("Bin Mismatch (up edge %d), cannot calculate efficiency",ipt);
	return;
      }
      Double_t yieldGenAcc=hGenLimAcc->Integral(bl,bu);
      hYieldMCacc->SetBinContent(ipt+1,yieldGenAcc);
      hYieldMCacc->SetBinError(ipt+1,TMath::Sqrt(yieldGenAcc));
   
    }
    hEffRecoGenAcc->Divide(hYieldData,hYieldMCacc,1.,1.,"b");
    TCanvas *cEff=new TCanvas("cEff","cEff",800,800);
    cEff->cd();
    hEffRecoGenAcc->Draw();
    
    TCanvas *cCompareRecos=new TCanvas("cCompareRecos","cCompareRecos",800,800);
    cCompareRecos->cd();
    hRecoAcc->Draw();
    hRecoAcc->SetLineColor(kBlue);
    
    hRecoAccFilt->Draw("Sames");
    hRecoAccFilt->SetLineColor(kRed);
    hRecoAccFilt->SetLineStyle(2);

  }
  
  TCanvas *cPtIntCase1=new TCanvas("cPtIntCase1","cPtIntCase1",800,800);
  cPtIntCase1->cd();
  TH1D *hmass=0x0;
  for(Int_t j=0;j<nbins;j++){
    if(!hmass){
      hmass=new TH1D(*hmassData[j]);
      hmass->SetName(Form("hmassPtMore%.1f",ptbins[j]));
    }
    else{
      hmass->Add(hmassData[j]);
    }
  }
  hmass->Sumw2();
  hmass->Draw();  

  TCanvas *cPtIntCase2=new TCanvas("cPtIntCase2","cPtIntCase2",800,800);
  cPtIntCase2->cd();
  TH1D *hmass2=0x0;
  for(Int_t j=1;j<nbins;j++){
    if(!hmass2){
      hmass2=new TH1D(*hmassData[j]);
      hmass2->SetName(Form("hmass2PtMore%.1f",ptbins[j]));
    }
    else{
      hmass2->Add(hmassData[j]);
    }
  }
  hmass2->Sumw2();
  hmass2->Draw();  

  if(useCuts==2 || useCuts==3){
    TString strCosThetStarSigma=Form("%d",(Int_t)TMath::Abs(sigmaCMincosthetaStarDefCut*100));
    if(sigmaCMincosthetaStarDefCut<0)strCosThetStarSigma.Prepend("Minus");
    else strCosThetStarSigma.Prepend("Plus");
    TFile *foutScan=new TFile(Form("fScanData_LcMassCut%.0f_CosThetStar%s_varScan%d_%d_caseCuts%s.root",lcMassCutDefault*1000.,strCosThetStarSigma.Data(),varscan0,varscan1,caseSelForUseCuts2==-1 ? Form("Def"):Form("%d",caseSelForUseCuts2)),"RECREATE");
    foutScan->cd();
    for(Int_t jpt=0;jpt<nbins;jpt++){
      cCutScan[jpt]->Write();
    }
    cCanvasNoCuts->Write();
    foutScan->Close();
  }
  
}


void SetFilteringCuts(const Int_t nptbins,const Double_t *ptbins,Double_t cutsMin[nptbins][10],Double_t cutsMax[nptbins][10],Int_t *isCommonCut,Bool_t isForBackgroundStudy=kFALSE){

  // inv. mass [GeV]	pTK [GeV/c]	    pTP [GeV/c]	       d0K [cm]   lower limit!	       d0Pi [cm]  lower limit!	
  //   0.35                0.4                 0.5                       0                             0
  // dist12 (cm)	sigmavert (cm)	 dist prim-sec (cm)	pM=Max{pT1,pT2,pT3} (GeV/c)	cosThetaPoint	Sum d0^2 (cm^2)	    dca cut (cm)	cut on pTpion [GeV/c]
  // 0.01                 0.06                0.005                    0.7                           0                 0               0.05                    0.4


  // these are global variables
  version=1;
  minMassScFit=0.148;
  maxMassScFit=0.190;
  fixSigma=0.0012;
  fixSigmaBW=0.0007;
  widthSigmaBW=0.00189;

  // these 
  Double_t lcMassCutDefault=isForBackgroundStudy ? 999. : 0.012*0.999;
  Int_t itsRefit=1;
  Double_t sigmaCMincosthetaStarDefCut=-1.1;

//   const Int_t nbins=(useGlobalPt ? nbinsGlobal : 6);
//   Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,16,24.};
//   Double_t *ptbins=(useGlobalPt ? &ptbinsGlobal[0] : &ptbinsLocal[0]);//new Double_t[nbins]);

//   const Int_t nptbinsPrel=6;
//   Double_t ptbinsPrel[nptbinsPrel+1]={1.,2.,3.,5.,8,16,24.};
//   if(nptbinsPrel!=nptbins){
//     Printf("SetQM2019MassPlotCuts:  incompatible number of ptbins ");
//     return;
//   }
 
    for(Int_t k=0;k<nptbins;k++){
//       if(TMath::Abs(ptbins[k]-ptbinsPrel[k])>1.e-9 || TMath::Abs(ptbins[k+1]-ptbinsPrel[k+1])>1.e-9){
// 	Printf("SetQM2019MassPlotCuts: Incompatible bin edges");
// 	return;
//       }
      cutsMin[k][0]=-1;// pt, dummy value
      cutsMax[k][0]=99;// pt, dummy value

      cutsMin[k][1]=-1;// mass, dummy value
      cutsMax[k][1]=999;// mass, dummy value

      cutsMin[k][2]=-1; // lxy
      cutsMax[k][2]=999; // lxy

      cutsMin[k][3]=-1;//nlxy
      cutsMax[k][3]=99;//nlxy

      cutsMin[k][4]=-1.1;// thetaPoint (set to filtering cut)
      cutsMax[k][4]=1.1;// thetaPoint (set to filtering cut)

      cutsMin[k][5]=-1;// topomatic
      cutsMax[k][5]=99;// topomatic

      if(itsRefit){
	cutsMin[k][6]=-0.001;// ITS refit: 0 -> ITSrefit = kTRUE
	cutsMax[k][6]=0.001;// ITS refit: 0 -> ITSrefit = kTRUE
      }
      else {
	cutsMin[k][6]=-99;// ITS refit: 0 -> ITSrefit = kTRUE
	cutsMax[k][6]=99;// ITS refit: 0 -> ITSrefit = kTRUE
      }

      cutsMin[k][7]=-0.0001;// pID sel, 0= Bayes
      cutsMax[k][7]=0.0001;// pID sel, 0= Bayes

      cutsMin[k][8]=2.286-lcMassCutDefault;// lcMassCut, a range mass be set, dummy value
      cutsMax[k][8]=2.286+lcMassCutDefault;// lcMassCut, a range mass be set, dummy value

      cutsMin[k][9]=-1.1;
      cutsMax[k][9]=1.1;
    }

    isCommonCut[0]=0;
    isCommonCut[1]=0;
    isCommonCut[2]=0;
    isCommonCut[3]=0;
    isCommonCut[4]=0;
    isCommonCut[5]=0;
    isCommonCut[6]=1;
    isCommonCut[7]=1;
    isForBackgroundStudy ? isCommonCut[8]=0 : isCommonCut[8]=1;
    isCommonCut[9]=1;

    return;

}


void SetCuts(Int_t setType,const Int_t nptbins,const Double_t *ptbins,Double_t cutsMin[nptbins][10],Double_t cutsMax[nptbins][10],Int_t *isCommonCut,Bool_t isForBackgroundStudy=kFALSE){
  if(setType<2)return;
  // inv. mass [GeV]	pTK [GeV/c]	    pTP [GeV/c]	       d0K [cm]   lower limit!	       d0Pi [cm]  lower limit!	
  //   0.35                0.4                 0.5                       0                             0
  // dist12 (cm)	sigmavert (cm)	 dist prim-sec (cm)	pM=Max{pT1,pT2,pT3} (GeV/c)	cosThetaPoint	Sum d0^2 (cm^2)	    dca cut (cm)	cut on pTpion [GeV/c]
  // 0.01                 0.06                0.005                    0.7                           0                 0               0.05                    0.4


  // these are global variables
  version=1;
  minMassScFit=0.148;
  maxMassScFit=0.190;
  fixSigma=0.0012;
  fixSigmaBW=0.0007;
  widthSigmaBW=0.00189;

  // these 
  Double_t lcMassCutDefault=isForBackgroundStudy ? 999. : 0.012*0.999;
  Int_t itsRefit=1;
  Double_t sigmaCMincosthetaStarDefCut=-1.1;

//   const Int_t nbins=(useGlobalPt ? nbinsGlobal : 6);
//   Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,16,24.};
//   Double_t *ptbins=(useGlobalPt ? &ptbinsGlobal[0] : &ptbinsLocal[0]);//new Double_t[nbins]);

//   const Int_t nptbinsPrel=6;
//   Double_t ptbinsPrel[nptbinsPrel+1]={1.,2.,3.,5.,8,16,24.};
//   if(nptbinsPrel!=nptbins){
//     Printf("SetQM2019MassPlotCuts:  incompatible number of ptbins ");
//     return;
//   }
  if(setType==2){

    for(Int_t k=0;k<nptbins;k++){
      cutsMin[k][0]=-1;// pt, dummy value
      cutsMax[k][0]=99;// pt, dummy value

      cutsMin[k][1]=-1;// mass, dummy value
      cutsMax[k][1]=999;// mass, dummy value

      cutsMin[k][2]=0.0100; // lxy
      cutsMax[k][2]=999; // lxy

      cutsMin[k][3]=2.;//nlxy
      cutsMax[k][3]=99;//nlxy

      cutsMin[k][4]=-1.1;// thetaPoint (set to filtering cut)
      cutsMax[k][4]=1.1;// thetaPoint (set to filtering cut)

      cutsMin[k][5]=0;// topomatic
      cutsMax[k][5]=2.5;// topomatic

      if(itsRefit){
	cutsMin[k][6]=-0.001;// ITS refit: 0 -> ITSrefit = kTRUE
	cutsMax[k][6]=0.001;// ITS refit: 0 -> ITSrefit = kTRUE
      }
      else {
	cutsMin[k][6]=-99;// ITS refit: 0 -> ITSrefit = kTRUE
	cutsMax[k][6]=99;// ITS refit: 0 -> ITSrefit = kTRUE
      }

      cutsMin[k][7]=-0.0001;// pID sel, 0= Bayes
      cutsMax[k][7]=0.0001;// pID sel, 0= Bayes

      cutsMin[k][8]=2.286-lcMassCutDefault;// lcMassCut, a range mass be set, dummy value
      cutsMax[k][8]=2.286+lcMassCutDefault;// lcMassCut, a range mass be set, dummy value

      cutsMin[k][9]=sigmaCMincosthetaStarDefCut;
      cutsMax[k][9]=1.1;
    }

    cutsMin[2][3]=2.5;
    cutsMin[2][4]=0.85;
    
    cutsMin[3][3]=2.;
    cutsMin[3][4]=-1.1;
    
    cutsMin[4][2]=0.;// the filtering cut will be used on Lxy
    cutsMin[4][3]=0.;
    cutsMin[4][4]=-1.1;
    cutsMax[4][5]=99.;
    cutsMin[4][9]=-1.1;

    isCommonCut[0]=0;
    isCommonCut[1]=0;
    isCommonCut[2]=0;
    isCommonCut[3]=0;
    isCommonCut[4]=0;
    isCommonCut[5]=0;
    isCommonCut[6]=1;
    isCommonCut[7]=1;
    isForBackgroundStudy ? isCommonCut[8]=0 : isCommonCut[8]=1;
    isCommonCut[9]=1;
  }
  
  return;

}

void SetQM2019MassPlotCuts(const Int_t nptbins,const Double_t *ptbins,Double_t cutsMin[nptbins][10],Double_t cutsMax[nptbins][10],Int_t *isCommonCut,Bool_t isForBackgroundStudy=kFALSE){
   // AliPhysics: VO_ALICE@AliPhysics::vAN-20191004_ROOT6-1 
  // FILTERING CUTS: object in /alice/cern.ch/user/m/mfaggin/XiC/CutObjects/2019Sep19_withDCAtovtx_enlargedMass/cuts_ppRef_5TeV_Prod_Data_XicCutsClass_withDCAcut_enlargedMass.root or in local path: /Users/administrator/ALICE/CHARM/XsiC/2018June12/TestRunWithTask/pp13TeV/2019_Oct07_Trains2512_2515/QM_PRELIMINARYMASSPLOTS
  //  cut object name: XictopKpiProdCuts
  // 
  // 
  // wagon settings (AddTask): kFALSE,0,0.,0.,"","alien:///alice/cern.ch/user/m/mfaggin/XiC/CutObject_first_LEGOtrains/D0toKpiCuts.root","D0toKpiCuts","alien:///alice/cern.ch/user/m/mfaggin/XiC/CutObjects/2019Sep19_withDCAtovtx_enlargedMass/cuts_ppRef_5TeV_Prod_Data_XicCutsClass_withDCAcut_enlargedMass.root","XictopKpiProdCuts","AR"
  // inv. mass [GeV]	pTK [GeV/c]	    pTP [GeV/c]	       d0K [cm]   lower limit!	       d0Pi [cm]  lower limit!	
  //   0.35                0.4                 0.5                       0                             0
  // dist12 (cm)	sigmavert (cm)	 dist prim-sec (cm)	pM=Max{pT1,pT2,pT3} (GeV/c)	cosThetaPoint	Sum d0^2 (cm^2)	    dca cut (cm)	cut on pTpion [GeV/c]
  // 0.01                 0.06                0.005                    0.7                           0                 0               0.05                    0.4

// fCutsRD[0][0] = 0.35	
// fCutsRD[1][0] = 0.4	
// fCutsRD[2][0] = 0.5
// fCutsRD[3][0] = 0	
// fCutsRD[4][0] = 0	
// fCutsRD[5][0] = 0.01	
// fCutsRD[6][0] = 0.06	
// fCutsRD[7][0] = 0.005
// fCutsRD[8][0] = 0.7	
// fCutsRD[9][0] = 0
// fCutsRD[10][0] = 0
// fCutsRD[11][0] = 0.05
// fCutsRD[12][0] = 0.4	

// fUsePreselect=0 
// Detectors used for PID: TPC TOF 
// Minimum TPC PID clusters = 0
// Maximum momentum for using TPC PID = 999999.000000
// TOF Mismatch probablility cut = 0.010000

// TRACK CUTS: those defined in the cut object
//
//
// SOFT PION CUTS (as defined in the task, check with correct AliPhysics version):
//       
//
// EVENT CUTS: defined in the D0 cut object
//



  // these are global variables
  version=1;
  minMassScFit=0.148;
  maxMassScFit=0.190;
  fixSigma=0.0012;
  fixSigmaBW=0.0007;
  widthSigmaBW=0.00189;

  // these 
  Double_t lcMassCutDefault=0.012*0.999;
  if(isForBackgroundStudy)lcMassCutDefault=999.;
  Int_t itsRefit=1;
  Double_t sigmaCMincosthetaStarDefCut=-1.1;
  const Int_t nptbinsPrel=6;
  Double_t ptbinsPrel[nptbinsPrel+1]={1.,2.,3.,5.,8,16,24.};
  /*Double_t ptbinsPrel[nptbinsPrel+1]={1.,2.,3.,5.,8,12,24.};*/
  if(nptbinsPrel!=nptbins){
    Printf("SetQM2019MassPlotCuts:  incompatible number of ptbins ");
    return;
  }
 
    for(Int_t k=0;k<nptbins;k++){
      if(TMath::Abs(ptbins[k]-ptbinsPrel[k])>1.e-9 || TMath::Abs(ptbins[k+1]-ptbinsPrel[k+1])>1.e-9){
	Printf("SetQM2019MassPlotCuts: Incompatible bin edges");
	return;
      }
      cutsMin[k][0]=-1;// pt, dummy value
      cutsMax[k][0]=99;// pt, dummy value

      cutsMin[k][1]=-1;// mass, dummy value
      cutsMax[k][1]=999;// mass, dummy value

      cutsMin[k][2]=0.0100; // lxy
      cutsMax[k][2]=999; // lxy

      cutsMin[k][3]=2.;//nlxy
      cutsMax[k][3]=99;//nlxy

      cutsMin[k][4]=-1.1;// thetaPoint (set to filtering cut)
      cutsMax[k][4]=1.1;// thetaPoint (set to filtering cut)

      cutsMin[k][5]=0;// topomatic
      cutsMax[k][5]=2.5;// topomatic

      if(itsRefit){
	cutsMin[k][6]=-0.001;// ITS refit: 0 -> ITSrefit = kTRUE
	cutsMax[k][6]=0.001;// ITS refit: 0 -> ITSrefit = kTRUE
      }
      else {
	cutsMin[k][6]=-99;// ITS refit: 0 -> ITSrefit = kTRUE
	cutsMax[k][6]=99;// ITS refit: 0 -> ITSrefit = kTRUE
      }

      cutsMin[k][7]=-0.0001;// pID sel, 0= Bayes
      cutsMax[k][7]=0.0001;// pID sel, 0= Bayes

      cutsMin[k][8]=2.286-lcMassCutDefault;// lcMassCut, a range mass be set, dummy value
      cutsMax[k][8]=2.286+lcMassCutDefault;// lcMassCut, a range mass be set, dummy value

      cutsMin[k][9]=sigmaCMincosthetaStarDefCut;
      cutsMax[k][9]=1.1;
    }

    cutsMin[2][3]=2.5;
    cutsMin[2][4]=0.85;
    
    cutsMin[3][3]=2.;
    cutsMin[3][4]=-1.1;
    
    cutsMin[4][2]=0.;// the filtering cut will be used on Lxy
    cutsMin[4][3]=0.;
    cutsMin[4][4]=-1.1;
    cutsMax[4][5]=99.;
    cutsMin[4][9]=-1.1;

    isCommonCut[0]=0;
    isCommonCut[1]=0;
    isCommonCut[2]=0;
    isCommonCut[3]=0;
    isCommonCut[4]=0;
    isCommonCut[5]=0;
    isCommonCut[6]=1;
    isCommonCut[7]=1;
    isForBackgroundStudy ? isCommonCut[8]=0 : isCommonCut[8]=1;
    isCommonCut[9]=1;

    return;
}


void SetLcCutsFromSigmaC(const Int_t nptbins,const Double_t cutsScMin[nptbins][10],const Double_t cutsScMax[nptbins][10],const Int_t *isPtIndepCutSc,Double_t cutsLcMin[nptbins][8],Double_t cutsLcMax[nptbins][8],Int_t *isPtIndepCutLc,Int_t checkOrigin=-1){


  for(Int_t ip=0;ip<nptbins;ip++){
    for(Int_t ivar=0;ivar<8;ivar++){// pt(dummy), inv mass (dummy), lxy, nlxy, thetaPoint, topomatic, specialSelAxis (to be changed after), PID sel
      cutsLcMin[ip][ivar]=cutsScMin[ip][ivar];// pt, dummy value
      cutsLcMax[ip][ivar]=cutsScMax[ip][ivar];// pt, dummy value
    }    
    if(checkOrigin==-1){
      cutsLcMin[ip][6]=-1;
      cutsLcMax[ip][6]=99;
    }
    else if(checkOrigin==4){
      cutsLcMin[ip][6]=6.;
      cutsLcMax[ip][6]=7.;
    }
    else if(checkOrigin==5){
      cutsLcMin[ip][6]=8.;
      cutsLcMax[ip][6]=9.;
    }
  }
  
  for(Int_t ivar=0;ivar<8;ivar++){
    isPtIndepCutLc[ivar]=isPtIndepCutSc[ivar];
  }
  isPtIndepCutLc[6]=1;

  return; 
}

void CalculateCombinedEff(TH1D** hCombined, TH1D** hEffDirect, TH1D** hEffRes2, TH1D** hEffRes3, TH1D** hEffRes4, TString suff){

    (*hCombined) = (TH1D*) (*hEffDirect)->Clone();
    int nbins = (*hCombined)->GetNbinsX();

    // loop over the pT bins to compute the combined efficiency
    for(int i=1; i<=nbins; i++){
      // direct decay
      double centre_direct  = (*hEffDirect)->GetBinCenter(i);
      double content_direct = (*hEffDirect)->GetBinContent(i);
      double error_direct   = (*hEffDirect)->GetBinError(i);
      std::cout << "   [direct channel]         bin centre: " << centre_direct << "   bin content: " << content_direct << "+-" << error_direct << std::endl;
      // resonant decay (2)
      double centre_res_2  = (*hEffRes2)->GetBinCenter(i);
      double content_res_2 = (*hEffRes2)->GetBinContent(i);
      double error_res_2   = (*hEffRes2)->GetBinError(i);
      std::cout << "   [resonant channel (2)]   bin centre: " << centre_res_2 << "   bin content: " << content_res_2 << "+-" << error_res_2 << std::endl;
      // resonant decay (3)
      double centre_res_3  = (*hEffRes3)->GetBinCenter(i);
      double content_res_3 = (*hEffRes3)->GetBinContent(i);
      double error_res_3   = (*hEffRes3)->GetBinError(i);
      std::cout << "   [resonant channel (3)]   bin centre: " << centre_res_3 << "   bin content: " << content_res_3 << "+-" << error_res_3 << std::endl;
      // resonant decay (4)
      double centre_res_4  = (*hEffRes4)->GetBinCenter(i);
      double content_res_4 = (*hEffRes4)->GetBinContent(i);
      double error_res_4   = (*hEffRes4)->GetBinError(i);
      std::cout << "   [resonant channel (4)]   bin centre: " << centre_res_4 << "   bin content: " << content_res_4 << "+-" << error_res_4 << std::endl;

      //  check consistent binning
      if(abs(centre_direct-centre_res_2)>0.01){
        std::cout << "ERROR: wrong binning between direct and resonant (2) channels" << std::endl;
        return;
      }
      if(abs(centre_res_3-centre_res_2)>0.01){
        std::cout << "ERROR: wrong binning between resonant (2) and resonant (3) channels" << std::endl;
        return;
      }
      if(abs(centre_res_3-centre_res_4)>0.01){
        std::cout << "ERROR: wrong binning between resonant (3) and resonant (4) channels" << std::endl;
        return;
      }

      // do the calculation
      double content_combined = ( content_direct*br_direct + content_res_2*br_res_2 + content_res_3*br_res_3 + content_res_4*br_res_4 )/( br_direct + br_res_2 + br_res_3 + br_res_4 );
      double err_combined = sqrt( pow( (content_direct-content_combined)*error_direct,2 ) + pow(br_direct*err_br_direct,2)
                                 +pow( (content_res_2 -content_combined)*error_res_2,2 )  + pow(br_res_2 *err_br_res_2 ,2)
                                 +pow( (content_res_3 -content_combined)*error_res_3,2 )  + pow(br_res_3 *err_br_res_3 ,2)
                                 +pow( (content_res_4 -content_combined)*error_res_4,2 )  + pow(br_res_4 *err_br_res_4 ,2)
                                )/( br_direct + br_res_2 + br_res_3 + br_res_4 );

      // assign the values
      (*hCombined)->SetBinContent(i,content_combined);
      (*hCombined)->SetBinError(i,err_combined);
      std::cout << "  ---> combined efficiency: " << (*hCombined)->GetBinContent(i) << "+-" << (*hCombined)->GetBinError(i) << std::endl;
    }

    //
    // plot comparison
    //
    TCanvas* cancomp = new TCanvas(Form("cancomp_%s",suff.Data()),Form("eff comparison %s",suff.Data()),1000,1000);
    //cancomp->Divide(2,1);
    cancomp->cd(1);
    (*hEffDirect)->SetMarkerStyle(24);
    (*hEffRes2)->SetMarkerStyle(25);
    (*hEffRes3)->SetMarkerStyle(26);
    (*hEffRes4)->SetMarkerStyle(32);
    (*hCombined)->GetYaxis()->SetRangeUser(0,0.25);
    (*hCombined)->Draw();
    //(*hCombined)->GetYaxis()->UnZoom();
    (*hEffDirect)->Draw("same");
    (*hEffRes2)->Draw("same");
    (*hEffRes3)->Draw("same");
    (*hEffRes4)->Draw("same");
    gPad->SetTicks();
    gPad->SetGridy();
    TLegend* legeffcomp = new TLegend(0.25,0.5,0.5,0.85);
    legeffcomp->SetHeader(Form("efficiency %s", suff.Data()));
    legeffcomp->AddEntry((*hCombined),"combined");
    legeffcomp->AddEntry((*hEffDirect),"direct channel (1)");
    legeffcomp->AddEntry((*hEffRes2),"resonant channel (2)");
    legeffcomp->AddEntry((*hEffRes3),"resonant channel (3)");
    legeffcomp->AddEntry((*hEffRes4),"resonant channel (4)");
    legeffcomp->Draw();

}

void StudyEfficiency(int lowch, int upch, TString strfile="AnalysisResults.root",TString strSuff="AR",TString strDir="PWG3_D2H_XicpKpi",Bool_t useGlobalPt=kFALSE,Int_t cutSet=0
, bool useScpt=false  // mfaggin
){

  TDatime *t=new TDatime();
  TString tdate(Form("%d%d%d",t->GetDate(),t->GetMinute(),t->GetSecond()));

  const Int_t nbins=(useGlobalPt ? nbinsGlobal : 6);
  const Int_t nVarSc=10;
  const Int_t nVarLc=8;
  Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,16,24.};
  /*Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,12,24.};*/
  Int_t skipBinLocal[6]={1,1,0,0,0,1};
  Double_t *ptbins=(useGlobalPt ? &ptbinsGlobal[0] : &ptbinsLocal[0]);//new Double_t[nbins]);
  Int_t *skipBin=(useGlobalPt ? &skipBinGlobal[0] : &skipBinLocal[0]);//new Int_t[nbins]);

  // Lc
  TH1D *hLcMCRecoPID=new TH1D("hLcMCRecoPID","hLcMCRecoPID",nbins,ptbins);
  TH1D *hLcMCGenLimAcc=new TH1D("hLcMCGenLimAcc","hLcMCGenLimAcc",nbins,ptbins);
  TH1D *hLcMCGenAcc=new TH1D("hLcMCGenAcc","hLcMCGenAcc",nbins,ptbins);
  // Lc(<-Sc) or Sc
  TH1D *hLcScMCRecoPID  =new TH1D(useScpt?"hScMCRecoPID"  :"hLcScMCRecoPID"  ,useScpt?"hScMCRecoPID"  :"hLcScMCRecoPID"  ,nbins,ptbins);
  TH1D *hLcScMCGenLimAcc=new TH1D(useScpt?"hScMCGenLimAcc":"hLcScMCGenLimAcc",useScpt?"hScMCGenLimAcc":"hLcScMCGenLimAcc",nbins,ptbins);
  TH1D *hLcScMCGenAcc   =new TH1D(useScpt?"hScMCGenAcc"   :"hLcScMCGenAcc"   ,useScpt?"hScMCGenAcc"   :"hLcScMCGenAcc"   ,nbins,ptbins);
  // Sc
  //TH1D *hScMCRecoPID  =new TH1D("hScMCRecoPID"  ,"hScMCRecoPID",nbins,ptbins);
  //TH1D *hScMCGenLimAcc=new TH1D("hScMCGenLimAcc","hScMCGenLimAcc",nbins,ptbins);
  //TH1D *hScMCGenAcc   =new TH1D("hScMCGenAcc"   ,"hScMCGenAcc",nbins,ptbins);

  Double_t cutsMinSigmaC[nbins][nVarSc];
  Double_t cutsMaxSigmaC[nbins][nVarSc];
  Double_t cutsMinLc[nbins][nVarLc];
  Double_t cutsMaxLc[nbins][nVarLc];
  Int_t isPtCommonCutSigmaC[nVarSc];
  Int_t isPtCommonCutLc[nVarLc];
  
  TFile *f=TFile::Open(strfile.Data(),"READ");
  THnSparseF *hsparseLcAllAxes,*hsparseScAllAxes;
  THnSparseF *hsparseMCLcFromSc;
  THnF *hNdMCLc, *hNdMcSigmaC;  // new
  TH3F *h3dMCLc, *h3dMCSigmaC;  // new (second TH3F)
  if(strDir.Contains("-1")){
    hsparseLcAllAxes=(THnSparseF*)f->Get("fhSparseAnalysis");
    hsparseScAllAxes=(THnSparseF*)f->Get("fhSparseAnalysisSigma");
    //h3dMCLc=(TH3F*)f->Get("fhistMCSpectrumAccLc");
    hNdMCLc=(THnF*)f->Get("fhistMCSpectrumAccLc");
    hsparseMCLcFromSc=(THnSparseF*)f->Get("fhistMCSpectrumAccLcFromSc");
    hNdMcSigmaC=(THnF*)f->Get("fhistMCSpectrumAccSc");
  }
  else{
    TDirectory *ddir=f->GetDirectory(Form("%s%s",strDir.Data(),strSuff.Data()));
    TList *lis=(TList*)ddir->Get(Form("outputList%s",strSuff.Data()));
    hsparseLcAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysis");
    hsparseScAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysisSigma");
    //h3dMCLc=(TH3F*)lis->FindObject("fhistMCSpectrumAccLc");
    hNdMCLc=(THnF*)lis->FindObject("fhistMCSpectrumAccLc");
    hsparseMCLcFromSc=(THnSparseF*)lis->FindObject("fhistMCSpectrumAccLcFromSc");
    hNdMcSigmaC=(THnF*)lis->FindObject("fhistMCSpectrumAccSc");
  }
  TString strvar[7]={"pt","mass","lxy","nLxy","cosThetaP","nDCA","SelSpec"};
  
  //
  //  Apply the selection of the decay channel
  //
  //  Generation level
  //
  hNdMCLc->GetAxis(3)->SetRangeUser(lowch,upch);
  h3dMCLc = (TH3F*) hNdMCLc->Projection(0,1,2,"d");
  std::cout << "===> [Lc - generated] Set decay channel: " << lowch << "," << upch << std::endl;
  //
  hsparseMCLcFromSc->GetAxis(6)->SetRangeUser(lowch,upch);
  std::cout << "===> [Lc(<-Sc) - generated] Set decay channel: " << lowch << "," << upch << std::endl;
  //
  // SigmaC
  hNdMcSigmaC->GetAxis(3)->SetRangeUser(lowch,upch);  // info about the Lc decay channel
  h3dMCSigmaC = (TH3F*) hNdMcSigmaC->Projection(0,1,2,"d");
  std::cout << "===> [Sc - generated] Set decay channel for Lc: " << lowch << "," << upch << std::endl;

  //
  //  Reconstructed level
  //
  hsparseLcAllAxes->GetAxis(8)->SetRangeUser(lowch,upch);
  std::cout << "===> [Lc - reconstructed] Set decay channel: " << lowch << "," << upch << std::endl;
  //
  hsparseScAllAxes->GetAxis(13)->SetRangeUser(lowch,upch);
  std::cout << "===> [Lc(<-Sc) - reconstructed] Set decay channel: " << lowch << "," << upch << std::endl;

  //
  //  Select only prompt reconstructed for Sc (for Lc it is done through the SetLcCutsFromSigmaC function)
  //
  hsparseScAllAxes->GetAxis(11)->SetRangeUser(4,4);
  std::cout << "######### HELLO ########" << std::endl;

  //
  // Gen level
  //

  // Lc(<-Sc) or Sc
  TH1D *hLcFromScGenLimAcc, *hLcFromScGenAcc;
  if(!useScpt){ // Lc(<-Sc)
    std::cout << "---> Taking GenAcc and GenLimAcc stuff from Lc(<-Sc) sparse" << std::endl;
    TAxis *axStepMClcFromSc=(TAxis*)hsparseMCLcFromSc->GetAxis(1);
    TAxis *axOriginMClcFromSc=(TAxis*)hsparseMCLcFromSc->GetAxis(2);
    axStepMClcFromSc->SetRangeUser(kGenLimAcc,kGenLimAcc);
    axOriginMClcFromSc->SetRange(1,1);// prompt 

    hLcFromScGenLimAcc=hsparseMCLcFromSc->Projection(0);
    hLcFromScGenLimAcc->SetName(Form("hLcFromScGenLimAcc_%d_%d",lowch,upch));

    axStepMClcFromSc->SetRangeUser(kGenAcc,kGenAcc);
    hLcFromScGenAcc=hsparseMCLcFromSc->Projection(0);
    hLcFromScGenAcc->SetName(Form("hLcFromScGenAcc_%d_%d",lowch,upch));
  }
  else{ // Sc prompt (1,1 bin)
    std::cout << "---> Taking GenAcc and GenLimAcc stuff from Sc THnF (projected to a TH3F)" << std::endl;
    hLcFromScGenLimAcc=h3dMCSigmaC->ProjectionX(Form("hScMCGenLimAccYield_%d_%d",lowch,upch),h3dMCSigmaC->GetYaxis()->FindBin(kGenLimAcc),h3dMCSigmaC->GetYaxis()->FindBin(kGenLimAcc),1,1);
    hLcFromScGenAcc   =h3dMCSigmaC->ProjectionX(Form("hScMCGenAccYield_%d_%d",lowch,upch)   ,h3dMCSigmaC->GetYaxis()->FindBin(kGenAcc)   ,h3dMCSigmaC->GetYaxis()->FindBin(kGenAcc),1,1);
  }
  
  // Lc prompt (1,1 bin)
  TH1D *hLcGenLimAcc=h3dMCLc->ProjectionX(Form("hLcMCGenLimAccYield_%d_%d",lowch,upch),h3dMCLc->GetYaxis()->FindBin(kGenLimAcc),h3dMCLc->GetYaxis()->FindBin(kGenLimAcc),1,1);
  TH1D *hLcGenAcc=h3dMCLc->ProjectionX(Form("hLcMCGenAccYield_%d_%d",lowch,upch),h3dMCLc->GetYaxis()->FindBin(kGenAcc),h3dMCLc->GetYaxis()->FindBin(kGenAcc),1,1);


  for(Int_t jpt=0;jpt<nbins;jpt++){
    hLcScMCGenLimAcc->SetBinContent(hLcScMCGenLimAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcFromScGenLimAcc->Integral(hLcFromScGenLimAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcFromScGenLimAcc->GetXaxis()->FindBin(ptbins[jpt+1]*0.9999)));
    hLcMCGenLimAcc->SetBinContent(hLcMCGenLimAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcGenLimAcc->Integral(hLcGenLimAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcGenLimAcc->GetXaxis()->FindBin(ptbins[jpt+1]*0.9999)));

    hLcScMCGenAcc->SetBinContent(hLcScMCGenAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcFromScGenAcc->Integral(hLcFromScGenAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcFromScGenAcc->GetXaxis()->FindBin(ptbins[jpt+1]*0.9999)));
    hLcMCGenAcc->SetBinContent(hLcMCGenAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcGenAcc->Integral(hLcGenAcc->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLcGenAcc->GetXaxis()->FindBin(ptbins[jpt+1]*0.9999)));

  }

  // RECO PART
  // SETTING CUTS    
  if(cutSet==0){
    SetQM2019MassPlotCuts(nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC);
  }
  else if(cutSet==1){
    SetFilteringCuts(nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC);
  }
  else {
    SetCuts(cutSet,nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC);
  }
  SetLcCutsFromSigmaC(nbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC,cutsMinLc,cutsMaxLc,isPtCommonCutLc,4);
    
    // check cut values:
  for(Int_t ipt=0;ipt<nbins;ipt++){
    Printf(" Pt bin %d",ipt);
    Printf("Lc cuts:");
    for(Int_t ivarLc=0;ivarLc<nVarLc;ivarLc++){
      Printf(" %f < x < %f (is common-pt cut :%d)",cutsMinLc[ipt][ivarLc],cutsMaxLc[ipt][ivarLc],isPtCommonCutLc[ivarLc]);       
    }
    Printf("Sc cuts:");
    for(Int_t ivarSc=0;ivarSc<nVarSc;ivarSc++){
      Printf(" %f < x < %f (is common-pt cut :%d)",cutsMinSigmaC[ipt][ivarSc],cutsMaxSigmaC[ipt][ivarSc],isPtCommonCutSigmaC[ivarSc]);
    }
  }

  // REDUCED SPARSES TO SPEED UP PROJECTIONS
    Int_t naxisSp=hsparseScAllAxes->GetNdimensions();
    //if(naxisSp>=12){
    if(naxisSp>12){ // mfaggin
      TAxis *axRot=hsparseScAllAxes->GetAxis(12);
      TString strAxRotTitle(axRot->GetTitle());
      if(!(strAxRotTitle.EqualTo("isRotated"))){
	Printf("Rotation axis not at expected dimension: something wrong with code version? Axis 12 pointer %p, title is %s",axRot,strAxRotTitle.Data());
	//	return;
      }
      axRot->SetRange(1,1);
    }

  Int_t dimredSc[nVarSc];
  Int_t ndimredSc=0;
  for(Int_t j=0;j<nVarSc;j++){
    if(isPtCommonCutSigmaC[j]==1){
      TAxis *ax=hsparseScAllAxes->GetAxis(j);
      ax->SetRangeUser(cutsMinSigmaC[0][j]< 0 ? cutsMinSigmaC[0][j]*0.99999 : cutsMinSigmaC[0][j]*1.00001,cutsMaxSigmaC[0][j] < 0 ? cutsMaxSigmaC[0][j]*1.000001 : cutsMaxSigmaC[0][j]*0.999999);// the pt index should not matter!
    }
    else{
      dimredSc[ndimredSc]=j;
      if(j==0 && useScpt){  // mfaggin
        std::cout << "###" << std::endl;
        std::cout << "### [StudyEfficiency] BEWARE: putting pT of Sigma c in position 0 before reducing the sparse!" << std::endl;
        std::cout << "###" << std::endl;
        dimredSc[ndimredSc]=10; // pT of Sc
      }
      ndimredSc++;
    }
  }
  //     TAxis *axM=hsparseScAllAxes->GetAxis(8);
  //     Double_t   lcMassCutDefault=0.012;
  //     axM->SetRangeUser(2.286-lcMassCutDefault,2.286+lcMassCutDefault);
  
  THnSparseF *hsparseSc=(THnSparseF*)hsparseScAllAxes->Projection(ndimredSc,dimredSc,"A");
  hsparseSc->SetName(Form("hsparseSc_%d_%d",lowch,upch));
  TAxis *axptSc=hsparseSc->GetAxis(0);
  std::cout << "---> [StudyEfficiency] used pT (reco) for Sc analysis: " << axptSc->GetTitle() << std::endl;
  
  Int_t dimredLc[nVarLc];
  Int_t ndimredLc=0;
  for(Int_t j=0;j<nVarLc;j++){
    if(isPtCommonCutLc[j]==1){
      TAxis *ax=hsparseLcAllAxes->GetAxis(j);
      ax->SetRangeUser(cutsMinLc[0][j]< 0 ? cutsMinLc[0][j]*0.999 : cutsMinLc[0][j]*1.001,cutsMaxLc[0][j] < 0 ? cutsMaxLc[0][j]*1.0001 : cutsMaxLc[0][j]*0.999);// the pt index should not matter!
    }
    else{
      dimredLc[ndimredLc]=j;
      ndimredLc++;
    }
  }
  THnSparseF *hsparseLc=(THnSparseF*)hsparseLcAllAxes->Projection(ndimredLc,dimredLc,"A");
  hsparseLc->SetName(Form("hsparseLc_%d_%d",lowch,upch));
  TAxis *axptLc=hsparseLc->GetAxis(0);    


  TCanvas *cLcMassPlots=new TCanvas(Form("cLcMassPlots_%d_%d",lowch,upch),Form("cLcMassPlots_%d_%d",lowch,upch),800,800);
  DivideCanvas(cLcMassPlots,nbins);
  TCanvas *cScMassPlots=new TCanvas(Form("cScMassPlots_%d_%d",lowch,upch),Form("cScMassPlots_%d_%d",lowch,upch),800,800);
  DivideCanvas(cScMassPlots,nbins);
  
  // NOW START WORK
  for(Int_t jpt=0;jpt<nbins;jpt++){
    if(skipBin[jpt])continue;
    axptLc->SetRangeUser(ptbins[jpt]*1.0001,ptbins[jpt+1]*0.9999);
    axptSc->SetRangeUser(ptbins[jpt]*1.0001,ptbins[jpt+1]*0.9999);
    // apply cuts
    Int_t reduceLc=0;
    Int_t reduceSc=0;
    // Lc
    for(Int_t j=2;j<nVarLc;j++){// exclude pt and mass
      if(isPtCommonCutLc[j]==0){
	TAxis *ax=hsparseLc->GetAxis(j-reduceLc);	   
	ax->SetRangeUser(cutsMinLc[jpt][j]< 0 ? cutsMinLc[jpt][j]*0.999 : cutsMinLc[jpt][j]*1.001,cutsMaxLc[jpt][j] < 0 ? cutsMaxLc[jpt][j]*1.0001 : cutsMaxLc[jpt][j]*0.999);
      }
      else reduceLc++;
    }
    // Sc
    for(Int_t j=2;j<nVarSc;j++){// exclude pt and mass
      if(isPtCommonCutSigmaC[j]==0){
	TAxis *ax=hsparseSc->GetAxis(j-reduceSc);
	ax->SetRangeUser(cutsMinSigmaC[jpt][j]< 0 ? cutsMinSigmaC[jpt][j]*0.999 : cutsMinSigmaC[jpt][j]*1.001,cutsMaxSigmaC[jpt][j] < 0 ? cutsMaxSigmaC[jpt][j]*1.0001 : cutsMaxSigmaC[jpt][j]*0.999);
      }
      else reduceSc++;
    }
    
    TPad *p=(TPad*)cLcMassPlots->cd(jpt+1);
    TH1D *hLc=hsparseLc->Projection(1);
    hLc->Sumw2();
    hLc->SetName(Form("hMCLcPt%d",jpt));
    hLc->SetTitle(Form("Lc %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
    hLc->Draw();
    hLcMCRecoPID->SetBinContent(hLcMCRecoPID->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hLc->Integral());
    
    p=(TPad*)cScMassPlots->cd(jpt+1);
    TH1D *hSc=hsparseSc->Projection(1);
    hSc->Sumw2();
    hSc->SetName(Form("hMCScPt%d",jpt));
    hSc->SetTitle(Form("Sc %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
    hSc->Draw();
    hLcScMCRecoPID->SetBinContent(hLcScMCRecoPID->GetXaxis()->FindBin(ptbins[jpt]*1.0001),hSc->Integral());
    
  }

  TCanvas *cCompareRecoGenSpectra=new TCanvas(Form("cCompareRecoGenSpectra_%d_%d",lowch,upch),Form("cCompareRecoGenSpectra_%d_%d",lowch,upch),800,800);
  cCompareRecoGenSpectra->Divide(1,2);
  cCompareRecoGenSpectra->cd(1);
  gPad->SetLogy();
  gPad->SetTicks();
  hLcMCGenLimAcc->SetLineColor(kBlack);
  hLcMCGenLimAcc->Sumw2();
  hLcMCGenLimAcc->SetMarkerStyle(20);
  hLcMCGenLimAcc->SetMarkerColor(hLcMCGenLimAcc->GetLineColor());
  hLcMCGenLimAcc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hLcMCGenLimAcc->Draw();

  hLcMCGenAcc->SetLineColor(kBlue);
  hLcMCGenAcc->Sumw2();
  hLcMCGenAcc->SetMarkerStyle(20);
  hLcMCGenAcc->SetMarkerColor(hLcMCGenAcc->GetLineColor());
  hLcMCGenAcc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hLcMCGenAcc->Draw("sames");

  hLcMCRecoPID->SetLineColor(kRed);
  hLcMCRecoPID->Sumw2();
  hLcMCRecoPID->SetMarkerStyle(20);
  hLcMCRecoPID->SetMarkerColor(hLcMCRecoPID->GetLineColor());
  hLcMCRecoPID->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hLcMCRecoPID->Draw("sames");

  TLegend* legGenReco = new TLegend(0.5,0.6,0.75,0.8);
  legGenReco->AddEntry(hLcMCGenLimAcc, "GenLimAcc");
  legGenReco->AddEntry(hLcMCGenAcc   , "GenAcc");
  legGenReco->AddEntry(hLcMCRecoPID  , "RecoPID");
  legGenReco->Draw();


  cCompareRecoGenSpectra->cd(2);
  gPad->SetLogy();
  gPad->SetTicks();
  hLcScMCGenLimAcc->SetLineColor(kBlack);
  hLcScMCGenLimAcc->Sumw2();
  hLcScMCGenLimAcc->SetMarkerStyle(20);
  hLcScMCGenLimAcc->SetMarkerColor(hLcScMCGenLimAcc->GetLineColor());
  hLcScMCGenLimAcc->GetXaxis()->SetTitle(Form("#it{p}_{T}^{%s} (GeV/#it{c})",useScpt?"#Sigma_{c}":"#Lambda_{c}(<-#Sigma_{c})"));
  hLcScMCGenLimAcc->Draw();

  hLcScMCGenAcc->SetLineColor(kBlue);
  hLcScMCGenAcc->Sumw2();
  hLcScMCGenAcc->SetMarkerStyle(20);
  hLcScMCGenAcc->SetMarkerColor(hLcScMCGenAcc->GetLineColor());
  hLcScMCGenAcc->GetXaxis()->SetTitle(Form("#it{p}_{T}^{%s} (GeV/#it{c})",useScpt?"#Sigma_{c}":"#Lambda_{c}(<-#Sigma_{c})"));
  hLcScMCGenAcc->Draw("sames");

  hLcScMCRecoPID->SetLineColor(kRed);
  hLcScMCRecoPID->Sumw2();
  hLcScMCRecoPID->SetMarkerStyle(20);
  hLcScMCRecoPID->SetMarkerColor(hLcScMCRecoPID->GetLineColor());
  hLcScMCRecoPID->GetXaxis()->SetTitle(Form("#it{p}_{T}^{%s} (GeV/#it{c})",useScpt?"#Sigma_{c}":"#Lambda_{c}(<-#Sigma_{c})"));
  hLcScMCRecoPID->Draw("sames");

  legGenReco->Draw();

  TCanvas *cEfficiency=new TCanvas(Form("cEfficiency_%d_%d",lowch,upch),Form("cEfficiency_%d_%d",lowch,upch),800,800);
  //
  //  choose the right histiogram to treat, according to the decay channel
  //
  /*hEffRecoLc=(TH1D*)hLcMCRecoPID->Clone("hEfficiencyLc_RecoGenLimAcc");
  hEffRecoLc->Divide(hEffRecoLc,hLcMCGenLimAcc,1.,1.,"B");
  hEffRecoLc->SetLineColor(kBlack);
  hEffRecoLc->SetMarkerColor(hEffRecoLc->GetLineColor());
  hEffRecoLc->Draw();
  hEffRecoLcSc=(TH1D*)hLcScMCRecoPID->Clone("hEfficiencyLcSc_RecoGenLimAcc");
  hEffRecoLcSc->Divide(hEffRecoLcSc,hLcScMCGenLimAcc,1.,1.,"B");
  hEffRecoLcSc->SetLineColor(kRed);
  hEffRecoLcSc->Draw("same");*/
  // [...] create a pointer to the proper histogram and use it
  TH1D** heff_ptr_Lc   = nullptr;
  TH1D** heff_ptr_LcSc = nullptr;
  if(lowch==1 && upch==4){  // direct channel
    std::cout << "=========> [StudyEfficiency] efficiency computed without channel separation" << std::endl;
    heff_ptr_Lc   = &hEffRecoLc;
    heff_ptr_LcSc = &hEffRecoLcSc;
  }
  else if(lowch==1 && upch==1){ // resonant channel (2)
    std::cout << "=========> [StudyEfficiency] efficiency computed for the direct channel" << std::endl;
    heff_ptr_Lc   = &hEffRecoLc_direct;
    heff_ptr_LcSc = &hEffRecoLcSc_direct;
  }
  else if(lowch==2 && upch==2){
    std::cout << "=========> [StudyEfficiency] efficiency computed for the resonant channel (2)" << std::endl;
    heff_ptr_Lc   = &hEffRecoLc_res2;
    heff_ptr_LcSc = &hEffRecoLcSc_res2;
  }
  else if(lowch==3 && upch==3){
    std::cout << "=========> [StudyEfficiency] efficiency computed for the resonant channel (3)" << std::endl;
    heff_ptr_Lc   = &hEffRecoLc_res3;
    heff_ptr_LcSc = &hEffRecoLcSc_res3;
  }
  else if(lowch==4 && upch==4){
    std::cout << "=========> [StudyEfficiency] efficiency computed for the resonant channel (4)" << std::endl;
    heff_ptr_Lc   = &hEffRecoLc_res4;
    heff_ptr_LcSc = &hEffRecoLcSc_res4;
  }
  (*heff_ptr_Lc) = (TH1D*)hLcMCRecoPID->Clone("hEfficiencyLc_RecoGenLimAcc");
  (*heff_ptr_Lc)->Divide((*heff_ptr_Lc),hLcMCGenLimAcc,1.,1.,"B");
  (*heff_ptr_Lc)->SetLineColor(kBlack);
  (*heff_ptr_Lc)->SetMarkerColor((*heff_ptr_Lc)->GetLineColor());
  (*heff_ptr_Lc)->Draw();
  (*heff_ptr_LcSc) = (TH1D*)hLcScMCRecoPID->Clone(useScpt?"hEfficiencySc_RecoGenLimAcc":"hEfficiencyLcSc_RecoGenLimAcc");
  (*heff_ptr_LcSc)->Divide((*heff_ptr_LcSc),hLcScMCGenLimAcc,1.,1.,"B");
  (*heff_ptr_LcSc)->SetLineColor(kRed);
  (*heff_ptr_LcSc)->Draw("same");

  gPad->SetTicks();
  TLegend* legeff = new TLegend(0.3,0.4,0.65,0.5);
  legeff->AddEntry((*heff_ptr_Lc),(*heff_ptr_Lc)->GetName());
  legeff->AddEntry((*heff_ptr_LcSc),(*heff_ptr_LcSc)->GetName());
  legeff->Draw();


  TCanvas *cEfficiencyGenAccOverLimAcc=new TCanvas(Form("cEfficiencyGenAccOverLimAcc_%d_%d",lowch,upch),Form("cEfficiencyGenAccOverLimAcc_%d_%d",lowch,upch),800,800);
  TH1D *hGenAccOverLimAccLc=(TH1D*)hLcMCGenAcc->Clone("hEfficiencyLc_GenAccLimAcc");
  hGenAccOverLimAccLc->Divide(hGenAccOverLimAccLc,hLcMCGenLimAcc,1.,1.);
  hGenAccOverLimAccLc->SetLineColor(kBlack);
  hGenAccOverLimAccLc->SetMarkerColor(hGenAccOverLimAccLc->GetLineColor());
  hGenAccOverLimAccLc->Draw();

  TH1D *hGenAccOverLimAccLcSc=(TH1D*)hLcScMCGenAcc->Clone(useScpt?"hEfficiencySc_GenAccLimAcc":"hEfficiencyLcSc_GenAccLimAcc");
  hGenAccOverLimAccLcSc->Divide(hGenAccOverLimAccLcSc,hLcScMCGenLimAcc,1.,1.,"B");
  hGenAccOverLimAccLcSc->SetLineColor(kRed);
  hGenAccOverLimAccLcSc->SetMarkerColor(hGenAccOverLimAccLcSc->GetLineColor());
  hGenAccOverLimAccLcSc->Draw("same");

  TLegend* legGenAcc_over_GenLimAcc = new TLegend(0.3,0.4,0.65,0.5);
  legGenAcc_over_GenLimAcc->AddEntry(hGenAccOverLimAccLc,hGenAccOverLimAccLc->GetName());
  legGenAcc_over_GenLimAcc->AddEntry(hGenAccOverLimAccLcSc,hGenAccOverLimAccLcSc->GetName());
  legeff->Draw();

}

//
//  Compute efficiency in one of the two cases
//    1. with separation of decay channels
//    2. without separation of decay channels
//
void CalculateEfficiency(bool do_channel_separation, TString strfile="AnalysisResults.root",TString strSuff="AR",TString strDir="PWG3_D2H_XicpKpi",Bool_t useGlobalPt=kFALSE,Int_t cutSet=0,bool useScpt=false){

  std::cout << "=========================================================" << std::endl;
  std::cout << "========== CalculateEfficiency function called ==========" << std::endl;
  std::cout << "=========================================================" << std::endl;

  // no channel separation
  if(!do_channel_separation){
    std::cout << "======= Efficiency computation =======" << std::endl;
    std::cout << "===      NO channel separation     ===" << std::endl;
    StudyEfficiency(1,4,strfile,strSuff,strDir,useGlobalPt,cutSet,useScpt);
  }

  // with channel separation
  else{
    std::cout << "======= Efficiency computation =======" << std::endl;
    std::cout << "===         Direct channel         ===" << std::endl;
    StudyEfficiency(1,1,strfile,strSuff,strDir,useGlobalPt,cutSet,useScpt);
    std::cout << "======= Efficiency computation =======" << std::endl;
    std::cout << "===      Resonant channel (2)      ===" << std::endl;
    StudyEfficiency(2,2,strfile,strSuff,strDir,useGlobalPt,cutSet,useScpt);
    std::cout << "======= Efficiency computation =======" << std::endl;
    std::cout << "===      Resonant channel (3)      ===" << std::endl;
    StudyEfficiency(3,3,strfile,strSuff,strDir,useGlobalPt,cutSet,useScpt);
    std::cout << "======= Efficiency computation =======" << std::endl;
    std::cout << "===      Resonant channel (4)      ===" << std::endl;
    StudyEfficiency(4,4,strfile,strSuff,strDir,useGlobalPt,cutSet,useScpt);

    //
    //  Calculate now the combined efficiency as the weighted sum of the 
    //  efficiency of each decay channel, considered separately
    //
    std::cout                                                                << std::endl;
    std::cout << "=========================================================" << std::endl;
    std::cout << "             Compute combined efficiency        "          << std::endl;
    std::cout << "=========================================================" << std::endl;
    std::cout << "Combined efficiency computed as"                           << std::endl;
    std::cout <<                                                                std::endl;
    std::cout << "   eff = sum_i (eff_i*B_i)/sum_i(B_i)"                     << std::endl;
    std::cout <<                                                                std::endl;
    std::cout << "=========================================================" << std::endl;

    // Lc
    std::cout << " ---> calculation for Lc" << std::endl;
    CalculateCombinedEff( &hEffRecoLc, &hEffRecoLc_direct, &hEffRecoLc_res2, &hEffRecoLc_res3, &hEffRecoLc_res4, "Lc" );
    // LcSc
    std::cout << " ---> calculation for LcSc" << std::endl;
    CalculateCombinedEff( &hEffRecoLcSc, &hEffRecoLcSc_direct, &hEffRecoLcSc_res2, &hEffRecoLcSc_res3, &hEffRecoLcSc_res4, useScpt?"Sc":"LcSc" );
    
  }

  if(!hEffRecoLc || !hEffRecoLcSc){
    std::cout << "ERROR: no total efficiency computed" << std::endl;
    std::cout << "hEffRecoLc "  << hEffRecoLc   << std::endl;
    std::cout << (useScpt?"hEffRecoSc":"hEffRecoLcSc") << hEffRecoLcSc << std::endl;
    return;
  }

  return;
}

void StudyFit(std::vector<TString> vec_filenames={},TString strSuff="AR",TString strDir="PWG3_D2H_XicpKpi",Bool_t useGlobalPt=kFALSE,Bool_t produceQM2019Performance=kTRUE,Int_t cutSet=0
, bool useScpt=false  // mfaggin
){

    TDatime *t=new TDatime();
    TString tdate(Form("%d%d%d",t->GetDate(),t->GetMinute(),t->GetSecond()));
    const Int_t nbins=(useGlobalPt ? nbinsGlobal : 6);
    Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,16,24.};
    /*Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,12,24.};*/
    Int_t skipBinLocal[6]={1,1,0,0,0,1};
    const Int_t nVarSc=10;
    const Int_t nVarLc=8;
    Double_t *ptbins=(useGlobalPt ? &ptbinsGlobal[0] : &ptbinsLocal[0]);//new Double_t[nbins]);
    Int_t *skipBin=(useGlobalPt ? &skipBinGlobal[0] : &skipBinLocal[0]);//new Int_t[nbins]);

    hRawYieldLc=new TH1D("hRawYieldLc","hRawYieldLc",nbins,ptbins);
    hRawYieldLcSc=new TH1D(useScpt?"hRawYieldSc":"hRawYieldLcSc",useScpt?"hRawYieldSc":"hRawYieldLcSc",nbins,ptbins);

    Double_t cutsMinSigmaC[nbins][nVarSc];
    Double_t cutsMaxSigmaC[nbins][nVarSc];
    Double_t cutsMinLc[nbins][nVarLc];
    Double_t cutsMaxLc[nbins][nVarLc];
    Int_t isPtCommonCutSigmaC[nVarSc];
    Int_t isPtCommonCutLc[nVarLc];
    
//    std::vector<TString> vec_filenames = {
//      //
//      // 2701
//      //
//       "/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_1.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_2.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_3.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_4.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_5.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_6.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_7.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_8.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_9.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_10.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_11.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_12.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2701/AnalysisResults_2701_data_child_13.root"
//      //
//      //  2702
//      //
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/AnalysisResults_2702_child1.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/AnalysisResults_2702_child2.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/child3/merged_2702_child3.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/AnalysisResults_2702_child4.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/AnalysisResults_2702_child5.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/child6/merged_2702_child6.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/AnalysisResults_2702_child7.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/child8/merged_2702_child8.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/AnalysisResults_2702_child9.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2702/child10/merged_2702_child10.root"
//      //
//      //  2703
//      //
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2703/merged_2703_child1.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2703/AnalysisResults_2703_child2_hadPID.root"
//      //
//      //  2704
//      //
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_1.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_2.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_3.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_4_posMagField.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_4_negMagField.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_5.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_6_hadPID.root"
//      ,"/media/mattia/VirMachine/SigmaC_analysis/Data_ptSc_rotations/2704/AnalysisResults_2704_data_child_7.root"
//    };
    THnSparseF *hsparseLc, *hsparseSc;

    for(int ifile=0; ifile<vec_filenames.size(); ifile++){
          TString strfile = vec_filenames.at(ifile);
          std::cout << "[StudyFit] Looking at file " << strfile.Data() << "..." << std::endl;
          TFile *f=TFile::Open(strfile.Data(),"READ");
          THnSparseF *hsparseLcAllAxes,*hsparseScAllAxes;
          if(strDir.Contains("-1")){
            hsparseLcAllAxes=(THnSparseF*)f->Get("fhSparseAnalysis");
            hsparseScAllAxes=(THnSparseF*)f->Get("fhSparseAnalysisSigma");
            if(ifile==0)  normCount=(AliNormalizationCounter*)f->Get(Form("normalizationCounter%s",strSuff.Data()));
            else          normCount->Add((AliNormalizationCounter*)(f->Get(Form("normalizationCounter%s",strSuff.Data()))->Clone()));
          }
          else{
            TDirectory *ddir=f->GetDirectory(Form("%s%s",strDir.Data(),strSuff.Data()));
            TList *lis=(TList*)ddir->Get(Form("outputList%s",strSuff.Data()));
            hsparseLcAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysis");
            hsparseScAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysisSigma");
            if(ifile==0) normCount=(AliNormalizationCounter*)ddir->Get(Form("normalizationCounter%s",strSuff.Data()));
            else        normCount->Add((AliNormalizationCounter*)(ddir->Get(Form("normalizationCounter%s",strSuff.Data()))->Clone()));
          }

          //     TAxis *axptLcAllAxes=hsparseLcAllAxes->GetAxis(0); // axes: pt;mass;lxy;nLxy;cosThatPoint;normImpParXY; Xic mass: 2.468, Lc:2.286
          //     TAxis *axLxyLc=hsparseLcAllAxes->GetAxis(2);
          //     TAxis *axNlxyLc=hsparseLcAllAxes->GetAxis(3);
          //     TAxis *axPointLc=hsparseLcAllAxes->GetAxis(4);
          //     TAxis *axNDCALc=hsparseLcAllAxes->GetAxis(5);
          //     TAxis *axSelLc=hsparseLcAllAxes->GetAxis(6);
          //     axSelLc->SetRangeUser(-1,99);
          //     TAxis *axSpecialPIDsLc=hsparseLcAllAxes->GetAxis(7);
          //     axSpecialPIDsLc->SetRangeUser(pid,pid);


          //     TAxis *axptScAllAxes=hsparseScAllAxes->GetAxis(0); // axes: pt;mass;lxy;nLxy;cosThatPoint;normImpParXY; Xic mass: 2.468, Sc:2.286
          //     TAxis *axLxySc=hsparseScAllAxes->GetAxis(2);
          //     TAxis *axNlxySc=hsparseScAllAxes->GetAxis(3);
          //     TAxis *axPointSc=hsparseScAllAxes->GetAxis(4);
          //     TAxis *axNDCASc=hsparseScAllAxes->GetAxis(5);
          //     TAxis *axITSrefitSc=hsparseScAllAxes->GetAxis(6);
          //     TAxis *axSpecialPIDsSc=hsparseScAllAxes->GetAxis(7);
          //     axSpecialPIDsSc->SetRangeUser(pid,pid);

          TString strvar[7]={"pt","mass","lxy","nLxy","cosThetaP","nDCA","SelSpec"};

          //     TAxis *axMassSigma=0x0;
          //     TAxis *axMassLcSigma=0x0;
          //     TAxis *axMassCosThetaStar=0x0;

          //     if(version>0){
          //       axMassLcSigma=hsparseScAllAxes->GetAxis(8);
          //       axMassCosThetaStar=hsparseScAllAxes->GetAxis(9);
          //     }
          //     else {
          //       axMassLcSigma=hsparseScAllAxes->GetAxis(7);
          //       axMassCosThetaStar=hsparseScAllAxes->GetAxis(8);
          //     }


          // SETTING CUTS    
          if(cutSet==0){
            SetQM2019MassPlotCuts(nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC);
          }
          else if(cutSet==1){
            SetFilteringCuts(nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC);
          }
          else {
            SetCuts(cutSet,nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC);
          }
          SetLcCutsFromSigmaC(nbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC,cutsMinLc,cutsMaxLc,isPtCommonCutLc);

          // check cut values:
          for(Int_t ipt=0;ipt<nbins;ipt++){
            Printf(" Pt bin %d",ipt);
            Printf("Lc cuts:");
            for(Int_t ivarLc=0;ivarLc<nVarLc;ivarLc++){
      	Printf(" %f < x < %f (is common-pt cut :%d)",cutsMinLc[ipt][ivarLc],cutsMaxLc[ipt][ivarLc],isPtCommonCutLc[ivarLc]);       
            }
            Printf("Sc cuts:");
            for(Int_t ivarSc=0;ivarSc<nVarSc;ivarSc++){
      	Printf(" %f < x < %f (is common-pt cut :%d)",cutsMinSigmaC[ipt][ivarSc],cutsMaxSigmaC[ipt][ivarSc],isPtCommonCutSigmaC[ivarSc]);
            }

          }

          Int_t naxisSp=hsparseScAllAxes->GetNdimensions();
          if(naxisSp>=12){
            TAxis *axRot=hsparseScAllAxes->GetAxis(12);
            TString strAxRotTitle(axRot->GetTitle());
            if(!(strAxRotTitle.EqualTo("isRotated"))){
      	Printf("Rotation axis not at expected dimension: something wrong with code version? Axis 12 pointer %p, title is %s",axRot,strAxRotTitle.Data());
      	//	return;
            }
            axRot->SetRange(1,1);
          }
          // REDUCED SPARSES TO SPEED UP PROJECTIONS
          Int_t dimredSc[nVarSc];
          Int_t ndimredSc=0;
          for(Int_t j=0;j<nVarSc;j++){  // Lc(<-Sc) or Sc
            if(isPtCommonCutSigmaC[j]==1){
      	      TAxis *ax=hsparseScAllAxes->GetAxis(j);
      	      ax->SetRangeUser(cutsMinSigmaC[0][j]< 0 ? cutsMinSigmaC[0][j]*0.99999 : cutsMinSigmaC[0][j]*1.00001,cutsMaxSigmaC[0][j] < 0 ? cutsMaxSigmaC[0][j]*1.000001 : cutsMaxSigmaC[0][j]*0.999999);// the pt index should not matter!
            }
            else{
      	      dimredSc[ndimredSc]=j;
              if(j==0 && useScpt){  // mfaggin
                std::cout << "###" << std::endl;
                std::cout << "### [StudyFit] BEWARE: putting pT of Sigma c in position 0 before reducing the sparse!" << std::endl;
                std::cout << "###" << std::endl;
                dimredSc[ndimredSc]=10; // pT of Sc
              }
      	      ndimredSc++;
            }
          }
          //     TAxis *axM=hsparseScAllAxes->GetAxis(8);
          //     Double_t   lcMassCutDefault=0.012;
          //     axM->SetRangeUser(2.286-lcMassCutDefault,2.286+lcMassCutDefault);

          if(ifile==0)  hsparseSc=(THnSparseF*) (hsparseScAllAxes->Projection(ndimredSc,dimredSc,"A"))->Clone();
          else          hsparseSc->Add( (THnSparseF*) ((THnSparseF*)hsparseScAllAxes->Projection(ndimredSc,dimredSc,"A"))->Clone() );
          

          Int_t dimredLc[nVarLc];
          Int_t ndimredLc=0;
          for(Int_t j=0;j<nVarLc;j++){
            if(isPtCommonCutLc[j]==1){
      	      TAxis *ax=hsparseLcAllAxes->GetAxis(j);
      	      ax->SetRangeUser(cutsMinLc[0][j]< 0 ? cutsMinLc[0][j]*0.999 : cutsMinLc[0][j]*1.001,cutsMaxLc[0][j] < 0 ? cutsMaxLc[0][j]*1.0001 : cutsMaxLc[0][j]*0.999);// the pt index should not matter!
            }
            else{
      	      dimredLc[ndimredLc]=j;
      	      ndimredLc++;
            }
          }
          if(ifile==0)  hsparseLc=(THnSparseF*)(hsparseLcAllAxes->Projection(ndimredLc,dimredLc,"A"))->Clone();
          else          hsparseLc->Add( (THnSparseF*)(hsparseLcAllAxes->Projection(ndimredLc,dimredLc,"A"))->Clone() );
          std::cout << "... done!" << std::endl;
          delete f;
          delete hsparseScAllAxes;
          delete hsparseLcAllAxes;
    }
    hsparseSc->SetName("hsparseSc");
    TAxis *axptSc=hsparseSc->GetAxis(0);
    std::cout << "---> [StudyFit] used pT (reco) for Sc analysis: " << axptSc->GetTitle() << std::endl;
    hsparseLc->SetName("hsparseLc");
    TAxis *axptLc=hsparseLc->GetAxis(0); 
    
    
    TCanvas *cLcMassPlots=new TCanvas("cLcMassPlots","cLcMassPlots",800,800);
    DivideCanvas(cLcMassPlots,nbins);
    TCanvas *cScMassPlots=new TCanvas("cScMassPlots","cScMassPlots",800,800);
    DivideCanvas(cScMassPlots,nbins);
    TCanvas *cScRes=new TCanvas("cScRes","cScRes",800,800);
    DivideCanvas(cScRes,nbins);
    TCanvas *cPerfQM2019;
    if(produceQM2019Performance)cPerfQM2019=new TCanvas("cPerfQM2019","cPerfQM2019",1200,600);
    cPerfQM2019->Divide(2,1);

    // NOW START WORK
    for(Int_t jpt=0;jpt<nbins;jpt++){
      Int_t reduceLc=0;
      Int_t reduceSc=0;
      if(skipBin[jpt])continue;
      axptLc->SetRangeUser(ptbins[jpt]*1.0001,ptbins[jpt+1]*0.9999);
      axptSc->SetRangeUser(ptbins[jpt]*1.0001,ptbins[jpt+1]*0.9999);
      // apply cuts
      // Lc
      for(Int_t j=2;j<nVarLc;j++){// exclude pt and mass
	if(isPtCommonCutLc[j]==0){
	  TAxis *ax=hsparseLc->GetAxis(j-reduceLc);	   
	  ax->SetRangeUser(cutsMinLc[jpt][j]< 0 ? cutsMinLc[jpt][j]*0.999 : cutsMinLc[jpt][j]*1.001,cutsMaxLc[jpt][j] < 0 ? cutsMaxLc[jpt][j]*1.0001 : cutsMaxLc[jpt][j]*0.999);
	}
	else reduceLc++;
      }
      // Sc
      for(Int_t j=2;j<nVarSc;j++){// exclude pt and mass
	if(isPtCommonCutSigmaC[j]==0){
	  TAxis *ax=hsparseSc->GetAxis(j-reduceSc);
	  ax->SetRangeUser(cutsMinSigmaC[jpt][j]< 0 ? cutsMinSigmaC[jpt][j]*0.999 : cutsMinSigmaC[jpt][j]*1.001,cutsMaxSigmaC[jpt][j] < 0 ? cutsMaxSigmaC[jpt][j]*1.0001 : cutsMaxSigmaC[jpt][j]*0.999);
	}
	else reduceSc++;
      }
      
      
      
      TPad *p=(TPad*)cLcMassPlots->cd(jpt+1);
      TH1D *hLc=hsparseLc->Projection(1);
      hLc->Sumw2();
      hLc->SetName(Form("hLcPt%d",jpt));
      hLc->SetTitle(Form("Lc %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hLc->Draw();
      Double_t signalLc=-1,errsignalLc=-1.,significanceLc=0;
      Double_t signalSc=-1,errsignalSc=-1.;
      
      Printf("fitting Lc");
      TF1 *f=FitHisto(hLc,0,-1,1,signalLc,errsignalLc,"RLEM");
      f->SetName(Form("fLc%d",jpt));
      significanceLc=errsignalLc>1.e-6 ? signalLc/errsignalLc : 0;
      //	    hSignificanceLc[jpt]->SetBinContent(ivar0+1,nscansVar1-ivar1,significanceLc);
      hRawYieldLc->SetBinContent(jpt+1,signalLc);
      hRawYieldLc->SetBinError(jpt+1,errsignalLc);
      WriteFitResultsAndAxisOnPad(hLc,f,p,signalLc,errsignalLc);
      
      
      p=(TPad*)cScMassPlots->cd(jpt+1);
      //	 Printf("Setting values %f,%f to axis %s",lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999,axVar1->GetName());
      TH1D *hSc=hsparseSc->Projection(1);
      hSc->Sumw2();
      hSc->SetName(Form("hScPt%d",jpt));
      hSc->SetTitle(Form("Sc %.0f<#it{p}_{T}^{%s}<%.0f",ptbins[jpt],useScpt?"#Sigma_{c}":"#Lambda_{c}(<-#Sigma_{c})",ptbins[jpt+1]));	 
      hSc->Draw();
      if(significanceLc>3.8){
	Printf("fitting Sc");
	TF1 *f=FitHisto(hSc,1,3,2,signalSc,errsignalSc,"RLEM");
	f->SetLineColor(kBlue);
	f->SetLineWidth(3);
	hSc->GetFunction("fFit")->SetLineColor(kBlue);
	hSc->GetFunction("fFit")->SetLineWidth(3);
	//f->Draw("same");
	f->SetName(Form("fSc%d",jpt));
	//	    hSignificanceSc[jpt]->SetBinContent(ivar0+1,nscansVar1-ivar1,errsignalSc>1.e-6 ? signalSc/errsignalSc : 0);
	if(errsignalSc>1.e-6 && errsignalLc>1.e-6){
	  //	      hRatioSigmaCoverLc[jpt]->SetBinContent(ivar0+1,nscansVar1-ivar1,signalSc/signalLc);
	  WriteFitResultsAndAxisOnPad(hSc,f,p,signalSc,errsignalSc);
	}
	hRawYieldLcSc->SetBinContent(jpt+1,signalSc);
	hRawYieldLcSc->SetBinError(jpt+1,errsignalSc);


	TF1 *fBack=f->DrawCopy("same");
	fBack->SetName("fback");
	fBack->SetParameter(0,0);
	fBack->SetLineStyle(2);
	fBack->SetLineWidth(3);
	fBack->SetLineColor(kRed);
	fBack->SetRange(minMassScFit,maxMassScFit);
	//	fBack->Draw("same");
	TF1 *fSign=new TF1(*f);//->Clone("fsign");
	fSign->SetName("fsign");
	fSign->SetParameter(3,0);
	if(produceQM2019Performance&(!useScpt) && (jpt==2 || jpt==3)){
	  TPad *pdPerf=(TPad*)cPerfQM2019->cd(1+jpt-2);
	  pdPerf->SetTopMargin(0.04);
	  pdPerf->SetBottomMargin(0.12);
	  if(jpt==2){
	    pdPerf->SetLeftMargin(0.14);
	    pdPerf->SetRightMargin(0.02);
	  }
	  if(jpt==3){
	    pdPerf->SetLeftMargin(0.14);
	    pdPerf->SetRightMargin(0.02);
	  }
	  pdPerf->SetTickx();
	  pdPerf->SetTicky();
	  TH1D *hPerf=(TH1D*)hSc->DrawCopy();
	  hPerf->SetTitle("");
	  hPerf->SetYTitle(Form("Counts / %.0f MeV/#it{c}^{2}",hPerf->GetBinWidth(1)*1000.));
	  hPerf->GetXaxis()->SetLabelSize(0.04);
	  hPerf->GetYaxis()->SetLabelSize(0.04);
	  hPerf->SetMarkerStyle(20);
	  hPerf->SetMarkerSize(0.8);
	  hPerf->SetMarkerColor(kBlack);
	  hPerf->SetLineColor(kBlack);
	  hPerf->SetLineWidth(2);
	  fBack->Draw("same");
	  hSc->GetFunction("fFit")->Draw("same");
	  //	  f->Draw("same");
	 
	  Double_t fwhm=2.*(fSign->GetParameter(1)-fSign->GetX(fSign->GetMaximum()/2.,minMassScFit,fSign->GetParameter(1)));
	  Printf("Signal FWHM is: %f",fwhm);
	  WriteFitResultsAndAxisOnPad(hPerf,f,pdPerf,signalSc,errsignalSc,fwhm,kTRUE,Form("%.0f < #it{p}_{T}^{#Lambda_{c}^{+}} < %.0f GeV/#it{c}",ptbins[jpt],ptbins[jpt+1]));
	  hPerf->GetXaxis()->SetRangeUser(minMassScFit,maxMassScFit);//TMath::Max(0.138,minMassScFit-(minMassScFit-0.138)*0.01),maxMassScFit+(maxMassScFit-0.138)*0.01);
	  hPerf->GetYaxis()->SetRangeUser(0.,hPerf->GetMaximum()*1.35);
	}
	TH1D *hRes=new TH1D(Form("hScResiduals_PtBin%d",jpt),Form("hScResiduals_PtBin%d",jpt),hSc->GetNbinsX(),hSc->GetBinLowEdge(1),hSc->GetBinLowEdge(hSc->GetNbinsX())+hSc->GetBinWidth(hSc->GetNbinsX()));
      //(TH1D*)hSel->Clone(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
      //      TH1D *hRes=(TH1D*)hSel->Clone(Form("hVarScanRes%sPt%d_i%d_j%d",strHist.Data(),j,var2fix,var1fix));
	Printf("Calculating Sc residuals for bin %d",jpt);
	for(Int_t kr=1;kr<=hRes->GetNbinsX();kr++){
	  if(hRes->GetBinLowEdge(kr)>=minMassScFit && (hRes->GetBinLowEdge(kr)+hRes->GetBinWidth(kr))<maxMassScFit){
	    hRes->SetBinContent(kr,hSc->GetBinContent(kr)-fBack->Eval(hSc->GetBinCenter(kr)));//-fBack->Integral(hRes->GetBinLowEdge(kr),hRes->GetBinLowEdge(kr)+hRes->GetBinWidth(kr))/hRes->GetBinWidth(kr));
	    hRes->SetBinError(kr,hSc->GetBinError(kr));
	  }
	  else hRes->SetBinContent(kr,0);	
	}
	cScRes->cd(jpt+1);
	hRes->Draw();
	fSign->SetRange(minMassScFit,maxMassScFit);
	fSign->Draw("same");
	WriteFitResultsAndAxisOnPad(hRes,f,(TPad*)(cScRes->cd(jpt+1)),signalSc,errsignalSc);
      }
    }
}

void GetCrossSectionLcSc(std::vector<TString> vec_files,TString strFileMC,Int_t cutset=0,TString strSuffData="AR",TString strDirData="PWG3_D2H_XicpKpi",TString strSuffMC="AR",TString strDirMC="PWG3_D2H_XicpKpi"
, bool useGlobalpt = false  // mfaggin
, bool do_ch_separation = false // mfaggin
, bool useScpt = false  // mfaggin
){
  
  //StudyEfficiency(strFileMC,strSuffMC,strDirMC,kTRUE,cutset); // arossi
  //StudyFit(strFileData,strSuffData,strDirData,kTRUE,kTRUE,cutset);  // arossi
  
  //StudyEfficiency(strFileMC,strSuffMC,strDirMC,useGlobalpt,cutset);
  CalculateEfficiency(do_ch_separation,strFileMC,strSuffMC,strDirMC,useGlobalpt,cutset,useScpt);
  StudyFit(vec_files,strSuffData,strDirData,useGlobalpt,kTRUE,cutset,useScpt);

  Int_t skipBinLocal[6]={1,1,0,0,0,1};
  Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,16,24.};
  Double_t *ptbins=(useGlobalpt ? &ptbinsGlobal[0] : &ptbinsLocal[0]);//new Double_t[nbins]);
  Int_t *skipBin=(useGlobalpt ? &skipBinGlobal[0] : &skipBinLocal[0]);//new Int_t[nbins]);

  Double_t nev=normCount->GetNEventsForNorm();

  TH1D *hCrossSectionLc=(TH1D*)hRawYieldLc->Clone("hCrossSectionLc");
  hCrossSectionLc->Divide(hEffRecoLc);
  for(Int_t i=1;i<=hCrossSectionLc->GetNbinsX();i++){
    if(skipBin[i-1])hCrossSectionLc->SetBinContent(i,0);
    else {
      hCrossSectionLc->SetBinContent(i,hCrossSectionLc->GetBinContent(i)/2./(ptbins[i]-ptbins[i-1])/nev*57.8*1000.);//microbarn/GeV/c
      hCrossSectionLc->SetBinError(i,hCrossSectionLc->GetBinError(i)/2./(ptbins[i]-ptbins[i-1])/nev*57.8*1000.);//microbarn/GeV/c
    }
  }

  TH1D *hCrossSectionLcSc=(TH1D*)hRawYieldLcSc->Clone(useScpt?"hCrossSectionSc":"hCrossSectionLcSc");
  hCrossSectionLcSc->Divide(hEffRecoLcSc);
  for(Int_t i=1;i<=hCrossSectionLcSc->GetNbinsX();i++){
    if(skipBin[i-1])hCrossSectionLcSc->SetBinContent(i,0);
    else {
      hCrossSectionLcSc->SetBinContent(i,hCrossSectionLcSc->GetBinContent(i)/2./(ptbins[i]-ptbins[i-1])/nev*57.8*1000.);//microbarn/GeV/c
      hCrossSectionLcSc->SetBinError(i,hCrossSectionLcSc->GetBinError(i)/2./(ptbins[i]-ptbins[i-1])/nev*57.8*1000.);//microbarn/GeV/c
    }
  }

  TCanvas *cRawYield=new TCanvas("cRawYield","cRawYield",800,800);
  hRawYieldLc->SetLineColor(kBlack);
  hRawYieldLc->SetMarkerColor(hRawYieldLc->GetLineColor());
  hRawYieldLc->SetMarkerStyle(20);
  hRawYieldLc->GetYaxis()->SetRangeUser(10,100000);
  hRawYieldLc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hRawYieldLc->Draw();
  hRawYieldLcSc->SetLineColor(kRed);
  hRawYieldLcSc->SetMarkerColor(hRawYieldLcSc->GetLineColor());
  hRawYieldLcSc->SetMarkerStyle(20);
  hRawYieldLcSc->GetXaxis()->SetTitle(Form("#it{p}_{T}^{%s} (GeV/#it{c})",useScpt?"#Sigma_{c}":"#Lambda_{c}(<-#Sigma_{c})"));
  hRawYieldLcSc->Draw("Same");
  gPad->SetTicks();
  gPad->SetLogy();
  TLegend* legRawYield = new TLegend(0.5,0.5,0.8,0.8);
  legRawYield->AddEntry(hRawYieldLc,hRawYieldLc->GetName());
  legRawYield->AddEntry(hRawYieldLcSc,hRawYieldLcSc->GetName());
  legRawYield->Draw();

  TCanvas *cCrossSectionCompare=new TCanvas("cCrossSectionCompare","cCrossSectionCompare",800,800);
  cCrossSectionCompare->cd();
  hCrossSectionLc->SetLineColor(kBlack);
  hCrossSectionLc->SetMarkerColor(hCrossSectionLc->GetLineColor());
  hCrossSectionLc->SetMarkerStyle(20);
  hCrossSectionLc->SetYTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub#upoint#it{c}/GeV)");
  hCrossSectionLc->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCrossSectionLc->GetYaxis()->SetRangeUser(0.001,10);
  hCrossSectionLc->Draw();
  hCrossSectionLcSc->SetLineColor(kRed);
  hCrossSectionLcSc->SetMarkerColor(hCrossSectionLcSc->GetLineColor());
  hCrossSectionLcSc->SetMarkerStyle(20);
  hCrossSectionLcSc->SetYTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub#upoint#it{c}/GeV)");
  hCrossSectionLcSc->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCrossSectionLcSc->Draw("same");
  gPad->SetTicks();
  gPad->SetLogy();
  TLegend* legCrossSection = new TLegend(0.5,0.5,0.8,0.8);
  legCrossSection->SetHeader("NO branching ratio correction"); // mfaggin
  legCrossSection->AddEntry(hCrossSectionLc,hCrossSectionLc->GetName());
  legCrossSection->AddEntry(hCrossSectionLcSc,hCrossSectionLcSc->GetName());
  legCrossSection->Draw();
  
  TFile *fPrel=TFile::Open("/home/mattia/Documenti/PhD/D2H_Analysis/SigmaC/files_Cristina/HFPtSpectrum_LcpKpi_STD_pp_13TeV_NbNbx2_22Oct_corrected_noWeigNb.root");
  TH1D *hPrel=(TH1D*)fPrel->Get("histoSigmaCorr");
  hPrel->SetName("histoSigmaCorr_Lcprel");  // mfaggin
  hPrel->Scale(1.e-6);// go to microbarn from picobarn
  hPrel->SetLineColor(kBlue);
  hPrel->SetMarkerColor(hPrel->GetLineColor());
  hPrel->Draw("same");
  legCrossSection->AddEntry(hPrel,"#Lambda_{c} preliminary");

  //
  //  NB: it makes sense only with same binning!
  //
  TCanvas *cLcCompare=new TCanvas("cLcCompare","cLcCompare",800,800);
  cLcCompare->cd();
  TH1D *hCrossSectionRatioLcMineLcPrel=(TH1D*)hCrossSectionLc->Clone("hCrossSectionRatioLcMineLcPrel");  
  hCrossSectionRatioLcMineLcPrel->Divide(hPrel);
  hCrossSectionRatioLcMineLcPrel->Draw();
  //

  if(!useScpt){
    TCanvas *cRatiosLcOverSc=new TCanvas("cRatiosLcOverSc","cRatiosLcOverSc",800,800);
    cRatiosLcOverSc->cd();
    TH1D *hCrossSectionRatioLcSc=(TH1D*)hCrossSectionLcSc->Clone("hCrossSectionRatioLcSc");
    hCrossSectionRatioLcSc->Divide(hCrossSectionLc);
    hCrossSectionRatioLcSc->Draw();
  }

  //
  //  NB: it makes sense only with same binning!
  //
  TCanvas *cCompareToD0=new TCanvas("cCompareToD0","cCompareToD0",800,800);
  cCompareToD0->cd();
  
  TFile *fPrelD0=TFile::Open("/home/mattia/Documenti/PhD/D2H_Analysis/SigmaC/files_Cristina/HFPtSpectrum_D0_multInt_22OctnoWeig_Nb.root");
  TH1D *hPrelD0=(TH1D*)fPrelD0->Get("histoSigmaCorr");
  hPrelD0->SetName("histoSigmaCorr_D0");
  hPrelD0->Scale(1.e-6);// go to microbarn from picobarn
  hPrelD0->Scale(1/0.0389);
  hPrelD0->SetYTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub#upoint#it{c}/GeV)");
  hPrelD0->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  hPrelD0->Draw("same");
  TH1D *hCrossSectionLcBRcorr=(TH1D*)hCrossSectionLc->Clone("hCrossSectionLcBRcorr");
  hCrossSectionLcBRcorr->Scale(1./0.0635);
  hCrossSectionLcBRcorr->Draw("same");

  TH1D *hCrossSectionLcScBRcorr=(TH1D*)hCrossSectionLcSc->Clone("hCrossSectionLcScBRcorr");
  hCrossSectionLcScBRcorr->Scale(1./0.0635);
  hCrossSectionLcScBRcorr->Draw("same");

  TCanvas *cRatioToD0=new TCanvas("cRatioToD0","cRatioToD0",800,800);
  cRatioToD0->cd();
  TH1D *hCrossSectionLcBRcorrRatioToD0=(TH1D*)hCrossSectionLcBRcorr->Clone("hCrossSectionLcBRcorrRatioToD0");
  TH1D *hCrossSectionLcScBRcorrRatioToD0=(TH1D*)hCrossSectionLcScBRcorr->Clone("hCrossSectionLcScBRcorrRatioToD0");
  hCrossSectionLcBRcorrRatioToD0->Divide(hPrelD0);
  hCrossSectionLcScBRcorrRatioToD0->Divide(hPrelD0);
  TH1D *hCrossSectionLcScBRcorrRatioToD0ScaledIsospin=(TH1D*)hCrossSectionLcScBRcorrRatioToD0->Clone("hCrossSectionLcScBRcorrRatioToD0ScalIsospin");
  hCrossSectionLcScBRcorrRatioToD0ScaledIsospin->Scale(3./2.);
  hCrossSectionLcBRcorrRatioToD0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCrossSectionLcBRcorrRatioToD0->SetYTitle("ratio to D0");
  hCrossSectionLcBRcorrRatioToD0->Draw();
  hCrossSectionLcScBRcorrRatioToD0->Draw("same");
  hCrossSectionLcScBRcorrRatioToD0ScaledIsospin->SetLineStyle(2);
  hCrossSectionLcScBRcorrRatioToD0ScaledIsospin->Draw("same");

  /*TCanvas *cCompareRatioToModels=new TCanvas("cCompareRatioToModels","cCompareRatioToModels",800,800);
  cCompareRatioToModels->cd();

  TFile *f=TFile::Open("/Users/administrator/soft/PYTHIA8/RunOnGrid/PYTHIA8243/HFbaryon/13TeV/2019Oct26/SIGMACplots/LcFromScOverD0Ratio_Mode2VsMonash.root");
  TCanvas *cModels=(TCanvas*)f->Get("cLundPoster");
  cCompareRatioToModels->cd();
  TH2F *hDrawRatioStyle_copy=(TH2F*)cModels->FindObject("hDrawRatioStyle_copy");
  hDrawRatioStyle_copy->Draw();
  TH1D *h1DtempLcFromSigmaY05_Mon2013def=(TH1D*)cModels->FindObject("h1DtempLcFromSigmaY05_Mon2013def");
  h1DtempLcFromSigmaY05_Mon2013def->Draw("HIST C Same");
  TH1D *h1DtempLcY05_Mon2013def=(TH1D*)cModels->FindObject("h1DtempLcY05_Mon2013def");
  h1DtempLcY05_Mon2013def->Draw("HIST C same");
  TH1D *h1DtempLcFromSigmaY05_Md2_13tev=(TH1D*)cModels->FindObject("h1DtempLcFromSigmaY05_Md2_13tev");
  h1DtempLcFromSigmaY05_Md2_13tev->Draw("HIST C same");
  TH1D *h1DtempLcY05_Md2_13tev=(TH1D*)cModels->FindObject("h1DtempLcY05_Md2_13tev");
  h1DtempLcY05_Md2_13tev->Draw("HIST C same");
  hCrossSectionLcBRcorrRatioToD0->Draw("same");
  hCrossSectionLcScBRcorrRatioToD0ScaledIsospin->Draw("same");*/

  TH1D* hPrelLcOverD0=(TH1D*)hPrel->Clone("hLcOverD0WithPreliminaryInputs");
  hPrelLcOverD0->Scale(1./0.0635);
  hPrelLcOverD0->Divide(hPrelD0);
  hPrelLcOverD0->Draw("same");
  //

  /*TLegend *leg=(TLegend*)cModels->FindObject("TPave");
  leg->Draw("same");*/
  

}

void StudyRotationalBackground(TString strfile="AnalysisResults.root",TString strSuff="AR",TString strDir="PWG3_D2H_XicpKpi",Bool_t useGlobalPt=kFALSE,Bool_t produceQM2019Performance=kTRUE,Int_t cutSet=0
, bool useScpt = false  // mfaggin
){


    TDatime *t=new TDatime();
    TString tdate(Form("%d%d%d",t->GetDate(),t->GetMinute(),t->GetSecond()));
    const Int_t nbins=(useGlobalPt ? nbinsGlobal : 6);
    Double_t ptbinsLocal[7]={1.,2.,3.,5.,8,16,24.};
    Int_t skipBinLocal[6]={1,1,0,0,0,1};
    const Int_t nVarSc=10;
    //    const Int_t nVarLc=8;
    Double_t *ptbins=(useGlobalPt ? &ptbinsGlobal[0] : &ptbinsLocal[0]);//new Double_t[nbins]);
    Int_t *skipBin=(useGlobalPt ? &skipBinGlobal[0] : &skipBinLocal[0]);//new Int_t[nbins]);

    hRawYieldLc=new TH1D("hRawYieldLc","hRawYieldLc",nbins,ptbins);
    hRawYieldLcSc=new TH1D("hRawYieldLcSc","hRawYieldLcSc",nbins,ptbins);

    Double_t cutsMinSigmaC[nbins][nVarSc];
    Double_t cutsMaxSigmaC[nbins][nVarSc];
    //    Double_t cutsMinLc[nbins][nVarLc];
    //    Double_t cutsMaxLc[nbins][nVarLc];
    Int_t isPtCommonCutSigmaC[nVarSc];
    //    Int_t isPtCommonCutLc[nVarLc];
    
    TFile *f=TFile::Open(strfile.Data(),"READ");
    THnSparseF *hsparseScAllAxes;
    if(strDir.Contains("-1")){
      //      hsparseLcAllAxes=(THnSparseF*)f->Get("fhSparseAnalysis");
      hsparseScAllAxes=(THnSparseF*)f->Get("fhSparseAnalysisSigma");
    }
    else{
      TDirectory *ddir=f->GetDirectory(Form("%s%s",strDir.Data(),strSuff.Data()));
      TList *lis=(TList*)ddir->Get(Form("outputList%s",strSuff.Data()));
      //      hsparseLcAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysis");
      hsparseScAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysisSigma");
    }
    
    TString strvar[7]={"pt","mass","lxy","nLxy","cosThetaP","nDCA","SelSpec"};
    
    // SETTING CUTS    
    if(cutSet==0){
      SetQM2019MassPlotCuts(nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC,kTRUE);
    }
    else if(cutSet==1){
      SetFilteringCuts(nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC,kTRUE);
    }
    else {
      SetCuts(cutSet,nbins,ptbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC,kTRUE);
    }
    //    SetLcCutsFromSigmaC(nbins,cutsMinSigmaC,cutsMaxSigmaC,isPtCommonCutSigmaC,cutsMinLc,cutsMaxLc,isPtCommonCutLc);
    TCanvas *cScMassRatioLplusRSBRotOverPeak=new TCanvas("cScMassRatioLplusRSBRotOverPeak","cScMassRatioLplusRSBRotOverPeak",800,800);
    DivideCanvas(cScMassRatioLplusRSBRotOverPeak,nbins);
    
    TCanvas **cScMassPlots=new TCanvas*[nbins];
    for(Int_t k=0;k<nbins;k++){
      cScMassPlots[k]=new TCanvas(Form("cScMassPlots_%d",k),Form("cScMassPlots %.1f<pt<%.1f",ptbins[k],ptbins[k+1]),800,800);
      cScMassPlots[k]->Divide(2,4);
    }

    // check cut values:
    for(Int_t ipt=0;ipt<nbins;ipt++){
      Printf(" Pt bin %d",ipt);
//       Printf("Lc cuts:");
//       for(Int_t ivarLc=0;ivarLc<nVarLc;ivarLc++){
// 	Printf(" %f < x < %f (is common-pt cut :%d)",cutsMinLc[ipt][ivarLc],cutsMaxLc[ipt][ivarLc],isPtCommonCutLc[ivarLc]);       
//       }
      Printf("Sc cuts:");
      for(Int_t ivarSc=0;ivarSc<nVarSc;ivarSc++){
	Printf(" %f < x < %f (is common-pt cut :%d)",cutsMinSigmaC[ipt][ivarSc],cutsMaxSigmaC[ipt][ivarSc],isPtCommonCutSigmaC[ivarSc]);
      }
      
    }

    // REDUCED SPARSES TO SPEED UP PROJECTIONS
    Printf("REDUCING SPARSES");
    Int_t naxisSp=hsparseScAllAxes->GetNdimensions();
    if(naxisSp>=12){
      TAxis *axRot=hsparseScAllAxes->GetAxis(12);
      TString strAxRotTitle(axRot->GetTitle());
      if(!(strAxRotTitle.EqualTo("isRotated"))){
	Printf("Rotation axis not at expected dimension: something wrong with code version? Axis 12 pointer %p, title is %s",axRot,strAxRotTitle.Data());
	//	return;
      }
    }
    
    Int_t dimredSc[nVarSc+1];
    Int_t ndimredSc=0;
    for(Int_t j=0;j<nVarSc;j++){
      if(isPtCommonCutSigmaC[j]==1){
	      TAxis *ax=hsparseScAllAxes->GetAxis(j);
	      ax->SetRangeUser(cutsMinSigmaC[0][j]< 0 ? cutsMinSigmaC[0][j]*0.99999 : cutsMinSigmaC[0][j]*1.00001,cutsMaxSigmaC[0][j] < 0 ? cutsMaxSigmaC[0][j]*1.000001 : cutsMaxSigmaC[0][j]*0.999999);// the pt index should not matter!
      }
      else{
	      dimredSc[ndimredSc]=j;
        if(j==0 && useScpt){  // mfaggin
          std::cout << "###" << std::endl;
          std::cout << "### [StudyRotationalBackground] BEWARE: putting pT of Sigma c in position 0 before reducing the sparse!" << std::endl;
          std::cout << "###" << std::endl;
          dimredSc[ndimredSc]=10; // pT of Sc
        }
	      ndimredSc++;
      }
    }
    dimredSc[ndimredSc]=12;// adding rotation axis
    ndimredSc++;
    
    THnSparseF *hsparseSc=(THnSparseF*)hsparseScAllAxes->Projection(ndimredSc,dimredSc,"A");
    hsparseSc->SetName("hsparseSc");
    TAxis *axptSc=hsparseSc->GetAxis(0);
    std::cout << "---> [StudyRotationalBackground] used pT (reco) for Sc analysis: " << axptSc->GetTitle() << std::endl;

    //     Int_t dimredLc[nVarLc];
    //     Int_t ndimredLc=0;
    //     for(Int_t j=0;j<nVarLc;j++){
    //       if(isPtCommonCutLc[j]==1){
    // 	TAxis *ax=hsparseLcAllAxes->GetAxis(j);
    // 	ax->SetRangeUser(cutsMinLc[0][j]< 0 ? cutsMinLc[0][j]*0.999 : cutsMinLc[0][j]*1.001,cutsMaxLc[0][j] < 0 ? cutsMaxLc[0][j]*1.0001 : cutsMaxLc[0][j]*0.999);// the pt index should not matter!
    //       }
    //       else{
    // 	dimredLc[ndimredLc]=j;
    // 	ndimredLc++;
    //       }
    //     }
    //     THnSparseF *hsparseLc=(THnSparseF*)hsparseLcAllAxes->Projection(ndimredLc,dimredLc,"A");
    //     hsparseLc->SetName("hsparseLc");
    //     TAxis *axptLc=hsparseLc->GetAxis(0);    
    

    // NOW START WORK

    // COMPARE rotated and unrotated background from Lc SB
    // COMPARE rotated and unrotated background below Lc peak
    
    Printf("NOW START PT LOOP");
    for(Int_t jpt=0;jpt<nbins;jpt++){
      //      Int_t reduceLc=0;
      Int_t reduceSc=0;
      if(skipBin[jpt])continue;
      //      axptLc->SetRangeUser(ptbins[jpt]*1.0001,ptbins[jpt+1]*0.9999);
      axptSc->SetRangeUser(ptbins[jpt]*1.0001,ptbins[jpt+1]*0.9999);
      
      // apply cuts
//       // Lc
//       for(Int_t j=2;j<nVarLc;j++){// exclude pt and mass
// 	if(isPtCommonCutLc[j]==0){
// 	  TAxis *ax=hsparseLc->GetAxis(j-reduceLc);	   
// 	  ax->SetRangeUser(cutsMinLc[jpt][j]< 0 ? cutsMinLc[jpt][j]*0.999 : cutsMinLc[jpt][j]*1.001,cutsMaxLc[jpt][j] < 0 ? cutsMaxLc[jpt][j]*1.0001 : cutsMaxLc[jpt][j]*0.999);
// 	}
// 	else reduceLc++;
//       }
      // Sc
 
      TAxis *axRotRedSp,*axLcMassRedSp;
      for(Int_t j=2;j<nVarSc;j++){// exclude pt and mass
	TAxis *ax;
	if(isPtCommonCutSigmaC[j]==0){
	  Printf("Getting ax %d-%d",j,reduceSc);
	  ax=hsparseSc->GetAxis(j-reduceSc);
	  if(j==8){
	    TString strTitleAxMassLc=ax->GetTitle();
	    if(!strTitleAxMassLc.EqualTo("LcMass")){
	      Printf("Something wrong in reduced-sparse axis indices, Lc Mass axis not at the expected number");
	      return;
	    }
	    axLcMassRedSp=ax;      
	  }
	  else  {
	    ax->SetRangeUser(cutsMinSigmaC[jpt][j]< 0 ? cutsMinSigmaC[jpt][j]*0.999 : cutsMinSigmaC[jpt][j]*1.001,cutsMaxSigmaC[jpt][j] < 0 ? cutsMaxSigmaC[jpt][j]*1.0001 : cutsMaxSigmaC[jpt][j]*0.999);
	  }
	}
	else {
	  Printf("Axis %d was reduced",j);
	  reduceSc++;
	}
      }
    
      axRotRedSp=hsparseSc->GetAxis(ndimredSc-1);
      TString strTitleAxRotReduceSparseLc=axRotRedSp->GetTitle();
      if(!strTitleAxRotReduceSparseLc.EqualTo("isRotated")){
	Printf("Something wrong in reduced-sparse axis indices, rotation axis not at the expected number");
	//	      return;
      }
      

      Printf("PROJECTIONS");
      axLcMassRedSp->SetRange(0,axLcMassRedSp->FindBin(2.286-0.018)); // Left SB
      axRotRedSp->SetRange(-1,0); // select rotated background    
      TPad *p=(TPad*)cScMassPlots[jpt]->cd(1);
      //	 Printf("Setting values %f,%f to axis %s",lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999,axVar1->GetName());
      TH1D *hScRotLSBlc=hsparseSc->Projection(1);
      Printf("First projection done");
      hScRotLSBlc->Sumw2();
      hScRotLSBlc->SetName(Form("hScRotLSBlcPt%d",jpt));
      hScRotLSBlc->SetTitle(Form("Rot Sc LeftLcSB %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hScRotLSBlc->SetLineColor(kRed);
      hScRotLSBlc->SetMarkerColor(kRed);
      hScRotLSBlc->Draw();
      
      axRotRedSp->SetRange(1,1); 
      TH1D *hScLSBlc=hsparseSc->Projection(1);
      hScLSBlc->Sumw2();
      hScLSBlc->SetName(Form("hScLSBlcPt%d",jpt));
      hScLSBlc->SetTitle(Form(" Sc LeftLcSB %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hScLSBlc->Draw("same");

      p=(TPad*)cScMassPlots[jpt]->cd(2);
      TH1D *hRatioRotLSBoverRot=(TH1D*)hScRotLSBlc->Clone(Form("hRatioScRotLSBlcOverLSBlcPt%d",jpt));
      hRatioRotLSBoverRot->Divide(hScLSBlc);
      hRatioRotLSBoverRot->Draw();
      

      p=(TPad*)cScMassPlots[jpt]->cd(3);
      axRotRedSp->SetRange(-1,0); 
      axLcMassRedSp->SetRange(axLcMassRedSp->FindBin(2.286+0.018),axLcMassRedSp->GetNbins()+1); // Right SB
      TH1D *hScRotRSBlc=hsparseSc->Projection(1);
      hScRotRSBlc->Sumw2();
      hScRotRSBlc->SetName(Form("hScRotRSBlcPt%d",jpt));
      hScRotRSBlc->SetTitle(Form("Rot Sc RightLcSB %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hScRotRSBlc->SetLineColor(kBlue);
      hScRotRSBlc->SetMarkerColor(kBlue);
      hScRotRSBlc->Draw();

      axRotRedSp->SetRange(1,1); 
      TH1D *hScRSBlc=hsparseSc->Projection(1);
      hScRSBlc->Sumw2();
      hScRSBlc->SetName(Form("hScRSBlcPt%d",jpt));
      hScRSBlc->SetTitle(Form(" Sc RightLcSB %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hScRSBlc->Draw("same");

      p=(TPad*)cScMassPlots[jpt]->cd(4);
      TH1D *hRatioRotRSBoverRot=(TH1D*)hScRotRSBlc->Clone(Form("hRatioScRotRSBlcOverRSBlcPt%d",jpt));
      hRatioRotRSBoverRot->Divide(hScRSBlc);
      hRatioRotRSBoverRot->Draw();


      p=(TPad*)cScMassPlots[jpt]->cd(5);
      axRotRedSp->SetRange(-1,0); 
      axLcMassRedSp->SetRange(axLcMassRedSp->FindBin(2.286-0.012),axLcMassRedSp->FindBin(2.286+0.012)); // Peak region
      TH1D *hScRotPeaklc=hsparseSc->Projection(1);
      hScRotPeaklc->Sumw2();
      hScRotPeaklc->SetName(Form("hScRotPeaklcPt%d",jpt));
      hScRotPeaklc->SetTitle(Form("Rot Sc Lc Peak %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hScRotPeaklc->SetLineColor(kViolet);
      hScRotPeaklc->SetMarkerColor(kViolet);
      hScRotPeaklc->Draw();

      axRotRedSp->SetRange(1,1); 
      TH1D *hScPeaklc=hsparseSc->Projection(1);
      hScPeaklc->Sumw2();
      hScPeaklc->SetName(Form("hScPeaklcPt%d",jpt));
      hScPeaklc->SetTitle(Form(" Sc Lc Peak %.0f<pt<%.0f",ptbins[jpt],ptbins[jpt+1]));	 
      hScPeaklc->Draw("same");

      p=(TPad*)cScMassPlots[jpt]->cd(6);
      TH1D *hRatioRotPeakoverRot=(TH1D*)hScRotPeaklc->Clone(Form("hRatioScRotPeaklcOverPeaklcPt%d",jpt));
      hRatioRotPeakoverRot->Divide(hScPeaklc);
      hRatioRotPeakoverRot->Draw();


      p=(TPad*)cScMassPlots[jpt]->cd(7);
      TH1D *hRatioLSBoverRSB=(TH1D*)hScLSBlc->Clone(Form("hRatioScLSBlcOverRSBlcPt%d",jpt));
      hRatioLSBoverRSB->Divide(hScRSBlc);
      hRatioLSBoverRSB->Draw();

      TH1D *hRatioRotLSBoverRotRSB=(TH1D*)hScRotLSBlc->Clone(Form("hRatioRotScLSBlcOverRotRSBlcPt%d",jpt));
      hRatioRotLSBoverRotRSB->Divide(hScRotRSBlc);
      hRatioRotLSBoverRotRSB->Draw("same");

      TH1D *hRatioLplusRSBoverPeak=(TH1D*)hScLSBlc->Clone(Form("hRatioScLplusRSBlcOverPeakLcPt%d",jpt));
      hRatioLplusRSBoverPeak->Add(hScRSBlc);
      hRatioLplusRSBoverPeak->Divide(hScPeaklc);
      hRatioLplusRSBoverPeak->SetLineColor(kPink);
      hRatioLplusRSBoverPeak->SetMarkerColor(kPink);
      hRatioLplusRSBoverPeak->Draw("same");

      TH1D *hRatioRotLplusRSBoverRotPeak=(TH1D*)hScRotLSBlc->Clone(Form("hRatioRotScLplusRSBlcOverRotPeaklcPt%d",jpt));
      hRatioRotLplusRSBoverRotPeak->Add(hScRotRSBlc);
      TH1D *hRatioRotLplusRSBoverPeak=(TH1D*)hRatioRotLplusRSBoverRotPeak->Clone(Form("hRatioRotScLplusRSBlcOverPeaklcPt%d",jpt));
      hRatioRotLplusRSBoverRotPeak->Divide(hScRotPeaklc);
      hRatioRotLplusRSBoverRotPeak->SetLineColor(kViolet);
      hRatioRotLplusRSBoverRotPeak->SetMarkerColor(kViolet);
      hRatioRotLplusRSBoverRotPeak->Draw("same");


      p=(TPad*)cScMassPlots[jpt]->cd(8);
      TH1D *hDoubleRatioRotOverNotRotLSBoverRSB=(TH1D*)hRatioRotLSBoverRotRSB->Clone(Form("hDoubleRatioRotOverNotRotLSBoverRSBPt%d",jpt));
      hDoubleRatioRotOverNotRotLSBoverRSB->Divide(hRatioLSBoverRSB);
      hDoubleRatioRotOverNotRotLSBoverRSB->Draw();

      TH1D *hDoubleRatioRotOverNotRotLplusRSBoverPeak=(TH1D*)hRatioRotLplusRSBoverRotPeak->Clone(Form("hDoubleRatioRotOverNotRotLplusRSBoverPeakPt%d",jpt));
      hDoubleRatioRotOverNotRotLplusRSBoverPeak->Divide(hRatioLplusRSBoverPeak);
      hDoubleRatioRotOverNotRotLplusRSBoverPeak->Draw("same");


      cScMassRatioLplusRSBRotOverPeak->cd(jpt+1);
      hRatioRotLplusRSBoverPeak->Divide(hScPeaklc);
      hRatioRotLplusRSBoverPeak->Draw();
    }
    
    

}

void CombinedAnalyseSigmaCLambdaC(Int_t readMC=0,Int_t tryFitData=0,Int_t useCuts=2,Bool_t lowField=kFALSE,Int_t varscan0=3/*2=Lxy,3=nLxy,4=cosThetaP,5=nDCA*/,Int_t varscan1=4,Int_t caseSelForUseCuts2=-1,TString strfile="AnalysisResults.root",TString strSuff="AR",TString strDir="PWG3_D2H_XicpKpi"
, bool useScpt=false  // mfaggin
){

  // TO-DO LIST:
  //    affianca plot dei residui come sottrazione del fondo dal fit, ai casi con gli studi vs. i tagli
  //    aggiungi fit di sottrazione del fondo da Lc side-bands
  //    includi fit Lc side bands --> funzione di fit = LcSB + pol1
  //    Mix-event e Rotational back
  //    OTTIMIZZAZIONE TAGLI

    TDatime *t=new TDatime();
    TString tdate(Form("%d%d%d",t->GetDate(),t->GetMinute(),t->GetSecond()));

  const Int_t nbins=6;
  const Int_t nscansVar0=5; 
  const Int_t nscansVar1=5; 
  const Int_t nvar=7;
  Int_t pid=0; 
  //  const Int_t nbinsTree=11;
  Double_t lcMassCutDefault=0.012*0.999;
  Int_t itsRefit=1;
  Double_t sigmaCMincosthetaStarDefCut=-1.1;
  Double_t ptbins[nbins+1]={1.,2.,3.,5.,8,16,24.};
  //  Double_t ptbinsTree[nbinsTree+1]={0.,1.,2.,3.,4.,6.,8,10.,12.,16,20.,24.};
  Int_t skipBin[nbins]={1,1,0,0,0,1};
  Double_t lowEdgeVar0[nbins][nscansVar0];
  Double_t lowEdgeVar1[nbins][nscansVar1];
  Double_t upEdgeVar0[nbins][nscansVar0];
  Double_t upEdgeVar1[nbins][nscansVar1];
  Double_t varLowStand[nvar]={-1,-1,0.010,4. ,-1.1,0.,-1};// used only when cuts are scanned ; variabls: pt, mass, lxy, nlxy, thetaPoint, NDCA, SpecialSel
  Double_t varUpStand[nvar]={ -1,-1,1.,   99.,1.1 ,2.,99};// used only when cuts are scanned; original NDCA for "preliminary" was 2
  TString strSel="SpecSelNo";
  

  switch (caseSelForUseCuts2){
  
  case 0:
    if(varscan0==3){
      varLowStand[2]=0.0100;    
    }    
    varLowStand[3]=2.5;
    varUpStand[5]=3.5;// value for good preliminary so far: 2.5
    break;
  case 1:
    if(varscan0==3){
      varLowStand[2]=0.0100;    
    }  
    varLowStand[3]=2.5;
    varUpStand[5]=3.5;// value for good preliminary so far: 2.5
 
    varLowStand[6]=2.001;
  case -1:
  default:
    // values used for preliminary (current) versions
    if(varscan1==5){
      varLowStand[4]=0.94;
    }
    if(varscan0==3){
      varLowStand[2]=0.0100;  //  value for good preliminary so far: 0.0150 ; otherwise 100
    }    
    varLowStand[3]=3.5;
    varUpStand[5]=2.5;// value for good preliminary so far: 2.5 ; otherwise 99
    
    break;
  }
  
  
  
  ////////////////////////////////////////////////////////
  //   const Int_t nScanLxy=3;
  //   const Int_t nScanNLxy=3;
  //   const Int_t nScanThPth=3;
  //   const Int_t nScanNDCA=3;
  //   const Int_t nScanCsThetStar=6;
  //   const Int_t nScanLcMass=3;
  
  //   Double_t scanGlobalLxyMin[nScanLxy]={0.01,0.02,0.035};
  //   Double_t scanGlobalLxyMax[nScanLxy]={99.,99.,99.};
  //   Double_t scanGlobalNLxyMin[nScanNLxy]={1.,3.5,6.};
  //   if(particle==0 || particle==1){
  //     scanGlobalNLxyMin[1]=2.0;
  //     scanGlobalNLxyMin[2]=3.5;
  //   }
  //   //Double_t scanGlobalNLxyMin[nScanNLxy]={4.5,5.5,7.};
  //   Double_t scanGlobalNLxyMax[nScanNLxy]={99.,99.,99.};
  //   Double_t scanGlobalThetaPtMin[nScanThPth]={0.85,0.9,0.95};
  //   Double_t scanGlobalThetaPtMax[nScanThPth]={1.,1.,1.};
  //   Double_t scanGlobalNDCAMin[nScanNDCA]={0.,0.,0.};
  //   Double_t scanGlobalNDCAMax[nScanNDCA]={1.,2.,4.};
  
  //   Double_t scanGlobalCosThStarMin[nScanCsThetStar]={-1.,-0.85,-0.4,-0.2,0.,0.2};
  //   Double_t scanGlobalCosThStarMax[nScanCsThetStar]={1.,1.,1.};
  //   Double_t scanGlobalLcMass[nScanLcMass]={0.030,0.015,0.007};
  
  // //   const Int_t nScanLxy=6;
  // //   const Int_t nScanNLxy=6;
  // //   const Int_t nScanThPth=6;
  // //   const Int_t nScanNDCA=6;
  
  // //   Double_t scanGlobalLxyMin[nScanLxy]={0.,0.1,0.15,0.2,0.25,0.35};
  // //   Double_t scanGlobalLxyMax[nScanLxy]={99.,99.,99.,99.,99.,99.};
  // //   Double_t scanGlobalNLxyMin[nScanNLxy]={1.,2.,3.,4.,5.,6.};
  // //   Double_t scanGlobalNLxyMax[nScanNLxy]={99.,99.,99.,99.,99.,99.};
  // //   Double_t scanGlobalThetaPtMin[nScanThPth]={0.8,0.85,0.9,0.92,0.94,0.98};
  // //   Double_t scanGlobalThetaPtMax[nScanThPth]={1.,1.,1.,1.,1.,1.};
  // //   Double_t scanGlobalNDCAMin[nScanNDCA]={0.,0.,0.,0.,0.,0.};
  // //   Double_t scanGlobalNDCAMax[nScanNDCA]={1.,1.5,2.,2.5,3.,4.};
  
  
  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
  Double_t scanNlxyLowEdgeLc[nbins][5]={
    {0.,1.,1.5,2.,2.5},
    {0.,1.,1.5,2.,2.5},
    {0.,1.5,2.,2.5,3.},
    {0.,2.,3.,4.,5.},
    {0.,2.,3.,4.,5.},
    {0.,2.,3.,4.,5.}
  };
  Double_t scanNlxyUpEdgeLc[nbins][5]={
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.}
  };
  Double_t scanLxyLowEdgeLc[nbins][5]={
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
    {0.0,0.005,0.01,0.015,0.02},
  };
  Double_t scanLxyUpEdgeLc[nbins][5]={
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.},
    {99.,99.,99.,99.,99.}
  };
  Double_t scanThPointLowEdgeLc[nbins][5]={
    {-1,0.,0.8,0.9,0.95},
    {-1,0.,0.8,0.9,0.95},
    {-1,0.8,0.82,0.85,0.9},
    {-1.,0.85,0.9,0.92,0.95},
    {0.8,0.85,0.9,0.92,0.95},
    {0.8,0.85,0.9,0.92,0.95},
  };
  Double_t scanThPointUpEdgeLc[nbins][5]={
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.},
    {1.,1.,1.,1.,1.}
  };

  Double_t scanNDCALowEdgeLc[nbins][5]={
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.}
  };
  Double_t scanNDCAUpEdgeLc[nbins][5]={
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.},
    {1.,1.5,2.,2.5,3.}
  };

//   ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
//   Double_t scanNlxyLowEdgeXic[nbins][5]={
//     {1.,3.,5.,7.,8.},
//     {1.,3.,5.,7.,8.},
//     {1.,3.,5.,7.,8.},
//     {1.,3.,5.,7.,8.},
//     {1.,3.,5.,7.,8.},
//     {1.,3.,5.,7.,8.}   
//   };
//   Double_t scanNlxyUpEdgeXic[nbins][5]={
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.}
//   };
//   Double_t scanLxyLowEdgeXic[nbins][5]={
//     {0.0,0.01,0.02,0.035,0.05},
//     {0.0,0.01,0.02,0.035,0.05},
//     {0.0,0.01,0.02,0.035,0.05},
//     {0.0,0.01,0.02,0.035,0.05},    
//     {0.0,0.01,0.02,0.035,0.05},
//     {0.0,0.01,0.02,0.035,0.05}
//   };
//   Double_t scanLxyUpEdgeXic[nbins][5]={
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.},
//     {99.,99.,99.,99.,99.}
//   };

//   Double_t scanThPointLowEdgeXic[nbins][5]={
//     {0.85,0.9,0.94,0.97,0.99},
//     {0.85,0.9,0.94,0.97,0.99},
//     {0.9,0.92,0.94,0.97,0.99},
//     {0.9,0.92,0.94,0.97,0.99},
//     {0.9,0.92,0.94,0.97,0.99},
//     {0.9,0.92,0.94,0.97,0.99}
//   };
//   Double_t scanThPointUpEdgeXic[nbins][5]={
//     {1.,1.,1.,1.,1.},
//     {1.,1.,1.,1.,1.},
//     {1.,1.,1.,1.,1.},
//     {1.,1.,1.,1.,1.},
//     {1.,1.,1.,1.,1.},
//     {1.,1.,1.,1.,1.}
//   };

//  Double_t scanNDCALowEdgeXic[nbins][5]={
//     {0.,0.,0.,0.,0.},
//     {0.,0.,0.,0.,0.},
//     {0.,0.,0.,0.,0.},
//     {0.,0.,0.,0.,0.},
//     {0.,0.,0.,0.,0.},
//     {0.,0.,0.,0.,0.}
//   };
//   Double_t scanNDCAUpEdgeXic[nbins][5]={
//     {0.5,1.,1.5,2.5,3.5},
//     {0.5,1.,1.5,2.5,3.5},
//     {0.5,1.,1.5,2.5,3.5},
//     {0.5,1.,1.5,2.5,3.5},
//     {0.5,1.,1.5,2.5,3.5},
//     {0.5,1.,1.5,2.5,3.5}
//   };




  for(Int_t k=0;k<nbins;k++){
    for(Int_t j=0;j<nscansVar0;j++){
      if(varscan0==2){
	lowEdgeVar0[k][j]=scanLxyLowEdgeLc[k][j];	
	upEdgeVar0[k][j]=scanLxyUpEdgeLc[k][j];	
      }
      if(varscan0==3){	
	lowEdgeVar0[k][j]=scanNlxyLowEdgeLc[k][j];	
	upEdgeVar0[k][j]=scanNlxyUpEdgeLc[k][j];	
      }
    }
    for(Int_t j=0;j<nscansVar1;j++){
      if(varscan1==4){		
	lowEdgeVar1[k][j]=scanThPointLowEdgeLc[k][j];	
	upEdgeVar1[k][j]=scanThPointUpEdgeLc[k][j];		
      }
      if(varscan1==5){	
	lowEdgeVar1[k][j]=scanNDCALowEdgeLc[k][j];	
	upEdgeVar1[k][j]=scanNDCAUpEdgeLc[k][j];	
      }
    }
  }
  
  
  TFile *f=TFile::Open(strfile.Data(),"READ");
  THnSparseF *hsparseLcAllAxes,*hsparseScAllAxes;
  if(strDir.Contains("-1")){
    hsparseLcAllAxes=(THnSparseF*)f->Get("fhSparseAnalysis");
    hsparseScAllAxes=(THnSparseF*)f->Get("fhSparseAnalysisSigma");
  }
  else{
    TDirectory *ddir=f->GetDirectory(Form("%s%s",strDir.Data(),strSuff.Data()));
    TList *lis=(TList*)ddir->Get(Form("outputList%s",strSuff.Data()));
    hsparseLcAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysis");
    hsparseScAllAxes=(THnSparseF*)lis->FindObject("fhSparseAnalysisSigma");
  }
  
  //  TH1D ***hVarSign=new TH1D**[17];
  //  TH1D ***hVarBack=new TH1D**[17];
  //  TCanvas **cVar=new TCanvas*[17];
  //  TCanvas *cPtTree;
  //  TH1D *hPtDistrTree;
  //  Float_t var[18];
  //  TString varNames[18]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","flagMC"};
  //  Double_t rangeL[18]={0.,-1,0,0,0,0,0,0,0,0,-0.1,-0.1,-0.1,-8,-8,-8,-0.1,-2};
  //  Double_t rangeU[18]={50,1,0.2,20,25,25,25,5,0.4,0.01,0.1,0.1,0.1,8,8,8,0.1,14};
  //  Double_t nCandSignal[nbinsTree];
  //  Double_t nCandBack[nbinsTree];

  TAxis *axptLcAllAxes=hsparseLcAllAxes->GetAxis(0); // axes: pt;mass;lxy;nLxy;cosThatPoint;normImpParXY; Xic mass: 2.468, Lc:2.286
  TAxis *axLxyLc=hsparseLcAllAxes->GetAxis(2);
  TAxis *axNlxyLc=hsparseLcAllAxes->GetAxis(3);
  TAxis *axPointLc=hsparseLcAllAxes->GetAxis(4);
  TAxis *axNDCALc=hsparseLcAllAxes->GetAxis(5);
  TAxis *axSelLc=hsparseLcAllAxes->GetAxis(6);
  axSelLc->SetRangeUser(-1,99);
  //  TAxis *axVar0Lc=hsparseLcAllAxes->GetAxis(varscan0);
  //  TAxis *axVar1Lc=hsparseLcAllAxes->GetAxis(varscan1);
  TAxis *axSpecialPIDsLc=hsparseLcAllAxes->GetAxis(7);
  axSpecialPIDsLc->SetRangeUser(pid,pid);

  TAxis *axptScAllAxes=hsparseScAllAxes->GetAxis(0); // axes: pt;mass;lxy;nLxy;cosThatPoint;normImpParXY; Xic mass: 2.468, Sc:2.286
  TAxis *axLxySc=hsparseScAllAxes->GetAxis(2);
  TAxis *axNlxySc=hsparseScAllAxes->GetAxis(3);
  TAxis *axPointSc=hsparseScAllAxes->GetAxis(4);
  TAxis *axNDCASc=hsparseScAllAxes->GetAxis(5);
  TAxis *axITSrefitSc=hsparseScAllAxes->GetAxis(6);
  //  TAxis *axVar0Sc=hsparseScAllAxes->GetAxis(varscan0);
  //  TAxis *axVar1Sc=hsparseScAllAxes->GetAxis(varscan1);
  TAxis *axSpecialPIDsSc=hsparseScAllAxes->GetAxis(7);
  axSpecialPIDsSc->SetRangeUser(pid,pid);


  TString strvar[7]={"pt","mass","lxy","nLxy","cosThetaP","nDCA","SelSpec"};

  TAxis *axMassSigma=0x0;
  TAxis *axMassLcSigma=0x0;
  TAxis *axMassCosThetaStar=0x0;
  
  //    axMassSigma=hsparse->GetAxis(6);
  if(version>0){
    axMassLcSigma=hsparseScAllAxes->GetAxis(8);
    axMassCosThetaStar=hsparseScAllAxes->GetAxis(9);
  }
  else {
    axMassLcSigma=hsparseScAllAxes->GetAxis(7);
    axMassCosThetaStar=hsparseScAllAxes->GetAxis(8);
  }
  
  
  // SIDE-BANDS STUDY
  if(lowField)Printf("######## LEFT AND RIGHT SIDEBANDS MAY CONTAIN SIGNALS ########");
  //  axMassLcSigma->SetRangeUser(2.286-0.030,2.286-0.024);
  //  Printf("Bins are: %d, %d",axMassLcSigma->FindBin(2.286-0.024),axMassLcSigma->FindBin(2.286+0.024));
  axMassLcSigma->SetRange(0,axMassLcSigma->FindBin(2.286-0.018));
  TH2D *hSc2DFromSparseLSB=hsparseScAllAxes->Projection(1,0,"A");
  hSc2DFromSparseLSB->SetName("hSc2DFromSparseLSB");
  //  axMassLcSigma->SetRangeUser(2.286+0.024,2.286+0.030);
  axMassLcSigma->SetRange(axMassLcSigma->FindBin(2.286+0.018),axMassLcSigma->GetNbins()+1);
  TH2D *hSc2DFromSparseRSB=hsparseScAllAxes->Projection(1,0,"A");
  hSc2DFromSparseRSB->SetName("hSc2DFromSparseRSB");

  TH1D **hSBtoSubtract=new TH1D*[nbins];
  TCanvas *cCanvasNoCutsScSBratio=new TCanvas("cCanvasNoCutsScSBratio","cCanvasNoCutsScSBratio",800,800);
  DivideCanvas(cCanvasNoCutsScSBratio,nbins);
  
  for(Int_t ipt=0;ipt<nbins;ipt++){
    TPad *p=(TPad*)cCanvasNoCutsScSBratio->cd(ipt+1);
    //    axptScAllAxes->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    //    TH1D *h=hsparseScAllAxes->Projection(1);
    hSc2DFromSparseLSB->GetXaxis()->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    TH1D *hLSB=hSc2DFromSparseLSB->ProjectionY();
    hLSB->SetName(Form("hLSB_ptbin%d",ipt));
    hLSB->Sumw2();
    hSBtoSubtract[ipt]=(TH1D*)hLSB;
    hSBtoSubtract[ipt]->Sumw2();
    hSc2DFromSparseRSB->GetXaxis()->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    TH1D *hRSB=hSc2DFromSparseRSB->ProjectionY();
    hRSB->SetName(Form("hRSB_ptbin%d",ipt));
    hRSB->Sumw2();
    hSBtoSubtract[ipt]->Add(hRSB);
    hRSB->Divide(hLSB);
    hRSB->SetName(Form("hSBratioNoCutsScSBratio_PtBin%d",ipt));
    hRSB->SetTitle(Form("hSBratioNoCutsScSBratio %.1f<pt<%.1f GeV/c",ptbins[ipt],ptbins[ipt+1]));    
    //    hRSB->Sumw2();
    hRSB->Draw();    
  }
  cCanvasNoCutsScSBratio->Draw();
  Printf("Now studying SB");
  StudySBfits(hSBtoSubtract,nbins);
  
  return;

  // DO IMMEDIATELY CANVAS W/O CUTS
  TCanvas *cCanvasNoCutsLc=new TCanvas("cCanvasNoCutsLc","cCanvasNoCutsLc",800,800);
  DivideCanvas(cCanvasNoCutsLc,nbins);

  // reduce sparses to 2D to speed up
  TH2D *hLc2DFromSparse=hsparseLcAllAxes->Projection(1,0,"A");
  hLc2DFromSparse->SetName("hLc2DFromSparse");

  if(!lowField)axMassLcSigma->SetRangeUser(2.286-lcMassCutDefault,2.286+lcMassCutDefault);
  else axMassLcSigma->SetRangeUser(2.286-0.025,2.286+0.025);
  TH2D *hSc2DFromSparse=hsparseScAllAxes->Projection(1,0,"A");
  hSc2DFromSparse->SetName("hSc2DFromSparse");
  //   TCanvas *cCheck=new TCanvas("cCheck","cCheck",800,800);
  //   cCheck->cd();
  //   hSc2DFromSparse->Draw("colz");

  // NOW FITTING< STARTING FROM Lc
  for(Int_t ipt=0;ipt<nbins;ipt++){
    TPad *p=(TPad*)cCanvasNoCutsLc->cd(ipt+1);
    //    axptLcAllAxes->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    //     TH1D *h=hsparseLcAllAxes->Projection(1);
    hLc2DFromSparse->GetXaxis()->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    TH1D *h=hLc2DFromSparse->ProjectionY();
    h->SetName(Form("hNoCutsLc_PtBin%d",ipt));
    h->SetTitle(Form("hNoCutsLc %.1f<pt<%.1f GeV/c",ptbins[ipt],ptbins[ipt+1]));    
    h->Sumw2();
    h->Draw();;    
    Double_t signalLc=-1,errsignalLc=-1.;//,significanceSc=0;
    TF1 *f=FitHisto(h,0,-1,1,signalLc,errsignalLc,"REMI");
    f->SetName(Form("fLcNoCuts%d",ipt));
    WriteFitResultsAndAxisOnPad(h,f,p,signalLc,errsignalLc);
  }
  cCanvasNoCutsLc->Draw();
  
  TCanvas *cCanvasNoCutsSc=new TCanvas("cCanvasNoCutsSc","cCanvasNoCutsSc",800,800);
  DivideCanvas(cCanvasNoCutsSc,nbins);

  TCanvas *cCanvasNoCutsScSubtr=new TCanvas("cCanvasNoCutsScSubtr","cCanvasNoCutsScSubtr",800,800);
  DivideCanvas(cCanvasNoCutsScSubtr,nbins);
  
  // SIGMA C FITS
  for(Int_t ipt=0;ipt<nbins;ipt++){
    TPad *p=(TPad*)cCanvasNoCutsSc->cd(ipt+1);
    //    axptScAllAxes->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    //    TH1D *h=hsparseScAllAxes->Projection(1);
    hSc2DFromSparse->GetXaxis()->SetRangeUser(ptbins[ipt]*1.001,ptbins[ipt+1]*0.999);
    TH1D *h=hSc2DFromSparse->ProjectionY();
    h->SetName(Form("hNoCutsSc_PtBin%d",ipt));
    h->SetTitle(Form("hNoCutsSc %.1f<pt<%.1f GeV/c",ptbins[ipt],ptbins[ipt+1]));    
    h->Sumw2();
    h->Draw();   
    TH1D *hSub=(TH1D*)h->Clone(Form("hNoCutsScSBSub_PtBin%d",ipt)); 
    Double_t signalSc=-1,errsignalSc=-1.;//,significanceSc=0;
    TF1 *f=FitHisto(h,1,3,2,signalSc,errsignalSc,"RLEMI");
    f->SetName(Form("fScNoCuts%d",ipt));
    WriteFitResultsAndAxisOnPad(h,f,p,signalSc,errsignalSc);
    TF1 *fSign=f->DrawCopy("same");
    fSign->SetName(Form("fScNoCutsNoBack%d",ipt));
    fSign->SetParameter(3,0);
    fSign->SetRange(minMassScFit,maxMassScFit);
    fSign->SetLineColor(kBlue);
    fSign->Draw("same");
    
    TPad *p2=(TPad*)cCanvasNoCutsScSubtr->cd(ipt+1);
    if(readMC!=2)hSub->Add(hSBtoSubtract[ipt],-1.*h->Integral(h->FindFixBin(0.175),h->FindFixBin(0.220))/hSBtoSubtract[ipt]->Integral(hSBtoSubtract[ipt]->FindFixBin(0.175),hSBtoSubtract[ipt]->FindFixBin(0.220)));
    signalSc=-1,errsignalSc=-1.;//,significanceSc=0;
    TF1 *fsub=FitHisto(hSub,1,-3,2,signalSc,errsignalSc,"REMI");
    fsub->SetName(Form("fScNoCutsSub%d",ipt));
    WriteFitResultsAndAxisOnPad(hSub,fsub,p2,signalSc,-errsignalSc);
    //    h->Rebin(4);
  }

  cCanvasNoCutsSc->Draw();

  //
  // THIS IS THE PART WITH SCAN OF CUTS DEFINED BY THE ARRAYES varLowStand and varUpStand
  TCanvas **cCutScanLc,**cCutScanSc;
  TCanvas **cCutScanLcSignif,**cCutScanScSignif;
  TCanvas **cCutScanLcSoverB,**cCutScanScSoverB;
  TCanvas **cCutScanLcRawYield,**cCutScanScRawYield;
  TCanvas **cCutScanScOverLcRawYield;
  TH2D **hSignificanceLc;
  TH2D **hSignificanceSc;
  TH2D **hRatioSigmaCoverLc;

  if(useCuts==2||useCuts==3){
    cCutScanLc=new TCanvas*[nbins];
    cCutScanSc=new TCanvas*[nbins];
    hSignificanceLc=new TH2D*[nbins];    
    hSignificanceLc=new TH2D*[nbins];    
    hSignificanceSc=new TH2D*[nbins];    
    
    cCutScanLcSignif=new TCanvas*[nbins];
    cCutScanScSignif=new TCanvas*[nbins];
    //     cCutScanLcSoverB=new TCanvas*[nbins];
    //     cCutScanScSoverB=new TCanvas*[nbins];
    //cCutScanLcRawYield=new TCanvas*[nbins];
    //    cCutScanScRawYield=new TCanvas*[nbins];
    cCutScanScOverLcRawYield=new TCanvas*[nbins];
    
    for(Int_t iv=2;iv<6;iv++){
      hsparseLcAllAxes->GetAxis(iv)->SetRangeUser(varLowStand[iv]*1.0001,varUpStand[iv]*0.999);     
      hsparseScAllAxes->GetAxis(iv)->SetRangeUser(varLowStand[iv]*1.0001,varUpStand[iv]*0.999);     
    }
    if(!lowField)axMassLcSigma->SetRangeUser(2.286-lcMassCutDefault,2.286+lcMassCutDefault);
    else axMassLcSigma->SetRangeUser(2.286-0.025,2.286+0.025);
    axMassCosThetaStar->SetRangeUser(sigmaCMincosthetaStarDefCut,1);
    if(itsRefit){
      axITSrefitSc->SetRange(1,1);
    }


    // REDUCE SPARSE DIMENSION TO SPEED UP
    const Int_t ndimred=4; Int_t dimred[ndimred]={0,1,varscan0,varscan1};
    if(useScpt){
      std::cout << "###" << std::endl;
      std::cout << "### [CombinedAnalyseSigmaCLambdaC] BEWARE: putting pT of Sigma c in position 0 before reducing the sparse!" << std::endl;
      std::cout << "###" << std::endl;
      dimred[0]=10; // pT of Sc
    }
    THnSparseF *hsparseLc=(THnSparseF*)hsparseLcAllAxes->Projection(ndimred,dimred,"A");
    hsparseLc->SetName("hsparseLc");
    THnSparseF *hsparseSc=(THnSparseF*)hsparseScAllAxes->Projection(ndimred,dimred,"A");
    hsparseSc->SetName("hsparseSc");

    TAxis *axVar0Lc=hsparseLc->GetAxis(2);
    TAxis *axVar1Lc=hsparseLc->GetAxis(3);

  TAxis *axVar0Sc=hsparseSc->GetAxis(2);
  TAxis *axVar1Sc=hsparseSc->GetAxis(3);
  TAxis *axptLc=hsparseLc->GetAxis(0); 
  TAxis *axptSc=hsparseSc->GetAxis(0);
  std::cout << "---> [CombinedAnalyseSigmaCLambdaC] used pT (reco) for Sc analysis: " << axptSc->GetTitle() << std::endl;
    //
    for(Int_t jpt=0;jpt<nbins;jpt++){
      if(skipBin[jpt])continue;
      cCutScanLc[jpt]=new TCanvas(Form("cCutScanLc%d",jpt),Form("cCutScanLc%d",jpt),800,800);
      cCutScanLc[jpt]->Divide(nscansVar0,nscansVar1);     
      cCutScanSc[jpt]=new TCanvas(Form("cCutScanSc%d",jpt),Form("cCutScanSc%d",jpt),800,800);
      cCutScanSc[jpt]->Divide(nscansVar0,nscansVar1);     
      cCutScanLcSignif[jpt]=new TCanvas(Form("cCutScanLcSignif%d",jpt),Form("cCutScanLcSignif%d",jpt),800,800);
      cCutScanScSignif[jpt]=new TCanvas(Form("cCutScanScSignif%d",jpt),Form("cCutScanScSignif%d",jpt),800,800);
      cCutScanScOverLcRawYield[jpt]=new TCanvas(Form("cCutScanScOverLcRawYield%d",jpt),Form("cCutScanScOverLcRawYield%d",jpt),800,800);

      hSignificanceLc[jpt]=new TH2D(Form("hSignificanceLc_PtBin%d",jpt),Form("hSignificanceLc_%d",jpt),nscansVar0,-0.5,nscansVar0-0.5,nscansVar1,-0.5,nscansVar1-0.5);
      hSignificanceSc[jpt]=new TH2D(Form("hSignificanceSc_PtBin%d",jpt),Form("hSignificanceSc_%d",jpt),nscansVar0,-0.5,nscansVar0-0.5,nscansVar1,-0.5,nscansVar1-0.5);
      hRatioSigmaCoverLc[jpt]=new TH2D(Form("hRatioSigmaCoverLc_PtBin%d",jpt),Form("hRatioSigmaCoverLc_%d",jpt),nscansVar0,-0.5,nscansVar0-0.5,nscansVar1,-0.5,nscansVar1-0.5);
      if(useCuts==2){
	axptLc->SetRangeUser(ptbins[jpt]*1.001,ptbins[jpt+1]*0.999);
	axptSc->SetRangeUser(ptbins[jpt]*1.001,ptbins[jpt+1]*0.999);
      }
      else if(useCuts==3){
	axptLc->SetRangeUser(ptbins[jpt]*1.001,999);
	axptSc->SetRangeUser(ptbins[jpt]*1.001,999);
      }
      cCutScanLc[jpt]->Draw();     
      cCutScanSc[jpt]->Draw();     
      for(Int_t ivar0=0;ivar0<nscansVar0;ivar0++){
	axVar0Lc->SetRangeUser(lowEdgeVar0[jpt][ivar0]*1.0001,upEdgeVar0[jpt][ivar0]*0.9999);
	axVar0Sc->SetRangeUser(lowEdgeVar0[jpt][ivar0]*1.0001,upEdgeVar0[jpt][ivar0]*0.9999);
	for(Int_t ivar1=0;ivar1<nscansVar1;ivar1++){
	  axVar1Lc->SetRangeUser(lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999);
	  axVar1Sc->SetRangeUser(lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999);

	  TPad *p=(TPad*)cCutScanLc[jpt]->cd(ivar1*nscansVar0+ivar0+1);
	  //	 Printf("Setting values %f,%f to axis %s",lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999,axVar1->GetName());
	  TH1D *hLc=hsparseLc->Projection(1);
	  hLc->Sumw2();
	  hLc->SetName(Form("hVarScanLcPt%d_i%d_j%d",jpt,ivar0,ivar1));
	  hLc->SetTitle(Form("Lc Pt%d, %s: %.2f - %.2f, %s: %.2f-%.2f",jpt,strvar[varscan0].Data(),lowEdgeVar0[jpt][ivar0],upEdgeVar0[jpt][ivar0],strvar[varscan1].Data(),lowEdgeVar1[jpt][ivar1],upEdgeVar1[jpt][ivar1]));	 
	  hLc->Draw();
	  Double_t signalLc=-1,errsignalLc=-1.,significanceLc=0;
	  Double_t signalSc=-1,errsignalSc=-1.;
	  if(tryFitData==1){
	    Printf("fitting Lc");
	    TF1 *f=FitHisto(hLc,0,-1,1,signalLc,errsignalLc);
	    f->SetName(Form("fLc%d_%d_%d",jpt,ivar0,ivar1));
	    significanceLc=errsignalLc>1.e-6 ? signalLc/errsignalLc : 0;
	    hSignificanceLc[jpt]->SetBinContent(ivar0+1,nscansVar1-ivar1,significanceLc);
	    WriteFitResultsAndAxisOnPad(hLc,f,p,signalLc,errsignalLc);
	  }

	  p=(TPad*)cCutScanSc[jpt]->cd(ivar1*nscansVar0+ivar0+1);
	  //	 Printf("Setting values %f,%f to axis %s",lowEdgeVar1[jpt][ivar1]*1.0001,upEdgeVar1[jpt][ivar1]*0.9999,axVar1->GetName());
	  TH1D *hSc=hsparseSc->Projection(1);
	  hSc->Sumw2();
	  hSc->SetName(Form("hVarScanScPt%d_i%d_j%d",jpt,ivar0,ivar1));
	  hSc->SetTitle(Form("Sc Pt%d, %s: %.2f - %.2f, %s: %.2f-%.2f",jpt,strvar[varscan0].Data(),lowEdgeVar0[jpt][ivar0],upEdgeVar0[jpt][ivar0],strvar[varscan1].Data(),lowEdgeVar1[jpt][ivar1],upEdgeVar1[jpt][ivar1]));	 
	  hSc->Draw();
	  if(tryFitData==1 && significanceLc>3.8){
	    Printf("fitting Sc");
	    TF1 *f=FitHisto(hSc,1,3,2,signalSc,errsignalSc);
	    f->SetName(Form("fSc%d_%d_%d",jpt,ivar0,ivar1));
	    hSignificanceSc[jpt]->SetBinContent(ivar0+1,nscansVar1-ivar1,errsignalSc>1.e-6 ? signalSc/errsignalSc : 0);
	    if(errsignalSc>1.e-6 && errsignalLc>1.e-6){
	      hRatioSigmaCoverLc[jpt]->SetBinContent(ivar0+1,nscansVar1-ivar1,signalSc/signalLc);
	      WriteFitResultsAndAxisOnPad(hSc,f,p,signalSc,errsignalSc);
	    }
	  }
	}
      }
      cCutScanScOverLcRawYield[jpt]->cd();
      hRatioSigmaCoverLc[jpt]->Draw("colz text");
      
      cCutScanLcSignif[jpt]->cd();
      hSignificanceLc[jpt]->Draw("colz text");
      
      cCutScanScSignif[jpt]->cd();
      hSignificanceSc[jpt]->Draw("colz text");
    }
  }



  TString strCosThetStarSigma=Form("%d",(Int_t)TMath::Abs(sigmaCMincosthetaStarDefCut*100));
    if(sigmaCMincosthetaStarDefCut<0)strCosThetStarSigma.Prepend("Minus");
    else strCosThetStarSigma.Prepend("Plus");
    TFile *foutScan=new TFile(Form("fScanCombLcSc_LcM%.0f_CsThSt%s_varScan%d_%d_caseCuts%s_itsref%d_pid%d_FitRang%d_%d_%s%s.root",lcMassCutDefault*1000.,strCosThetStarSigma.Data(),varscan0,varscan1,caseSelForUseCuts2==-1 ? Form("Def"):Form("%d",caseSelForUseCuts2),itsRefit,pid,(Int_t)(minMassScFit*1000),(Int_t)(maxMassScFit*1000),tdate.Data(),fixSigma>0 ? Form("fxSg%d",(Int_t)(fixSigma*10000)) : ""),"RECREATE");
    foutScan->cd();
    for(Int_t jpt=0;jpt<nbins;jpt++){
      if(skipBin[jpt])continue;
      cCutScanLc[jpt]->Write();
      cCutScanSc[jpt]->Write();
      cCutScanScOverLcRawYield[jpt]->Write();
      hSignificanceLc[jpt]->Write();
      hSignificanceSc[jpt]->Write();
    }
    cCanvasNoCutsLc->Write();
    cCanvasNoCutsSc->Write();
    foutScan->Close();
    // now save a snapshot of this macro
    gROOT->ProcessLine(Form(".!cp /Users/administrator/ALICE/CHARM/XsiC/2018June12/TestRunWithTask/AnalyseOutput.C AnalyseOutput_%s.C",tdate.Data()));
    
}
