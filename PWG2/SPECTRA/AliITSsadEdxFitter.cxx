/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: $ */

///////////////////////////////////////////////////////////////////////
// Class with the fits algorithms to be used in the identified       //
// spectra analysis using the ITS in stand-alone mode                //
// Author: E.Biolcati,biolcati@to.infn.it                            //
///////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TLatex.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLine.h>
#include <TH2F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <Rtypes.h>
#include "AliITSsadEdxFitter.h"


ClassImp(AliITSsadEdxFitter)
//______________________________________________________________________
AliITSsadEdxFitter::AliITSsadEdxFitter():TObject(){
  // standard constructor
  for(Int_t i=0; i<5; i++)  fFitpar[i] = 0.;
  for(Int_t i=0; i<5; i++)  fFitparErr[i] = 0.;
  SetRangeStep1();
  SetRangeStep2();
  SetRangeStep3();
  SetRangeFinalStep();
  SetLimitsOnSigmaPion();
  SetLimitsOnSigmaKaon();
  SetLimitsOnSigmaProton();
};

//________________________________________________________
Double_t AliITSsadEdxFitter::CalcSigma(Int_t code,Float_t x,Bool_t mc){
  // calculate the sigma 12-ott-2010  
  Double_t p[2]={0.};
  if(mc){
    if(code==211){
      p[0] = -1.20337e-04;
      p[1] = 1.13060e-01;
    }    
    else if(code==321 && x>0.15){
      p[0] = -2.39631e-03;
      p[1] = 1.15723e-01;
    }    
    else if(code==2212 && x>0.3){
      p[0] = -8.34576e-03;
      p[1] = 1.34237e-01;
    }    
    else return -1;
  }
  else{
    if(code==211){
      p[0] = -6.55200e-05;
      p[1] = 1.26657e-01;
    } 
    else if(code==321 && x>0.15){
      p[0] = -6.22639e-04;
      p[1] = 1.43289e-01;
    }    
    else if(code==2212 && x>0.3){
      p[0] = -2.13243e-03;
      p[1] = 1.68614e-01;
    } 
    else return -1;
  }
  return p[0]/(x*x)*TMath::Log(x)+p[1];
}

//_______________________________________________________
Int_t AliITSsadEdxFitter::CalcMean(Int_t code, Float_t x, Float_t mean0, Float_t &mean1, Float_t &mean2){
  // calculate the mean 12-ott-2010  
  Double_t p1[4]={0.};
  Double_t p2[4]={0.};
  if(code==211){
    p1[0] = 1.77049;
    p1[1] = -2.65469;
    p2[0] = 0.890856;
    p2[1] = -0.276719;
    mean1 = mean0 + p1[0]+ p1[1]*x + p1[2]*x*x + p1[3]*x*x*x;
    mean2 = mean1 + p2[0]+ p2[1]*x + p2[2]*x*x + p2[3]*x*x*x;
  }
  else if(code==321){
    p2[0] = 1.57664;
    p2[1] = -6.88635;
    p2[2] = 18.702;
    p2[3] = -16.3385;
    mean1 = 0.;
    mean2 = mean1 + p2[0]+ p2[1]*x + p2[2]*x*x + p2[3]*x*x*x;
  }
  else if(code==2212){
    p1[0] = 4.24861; 
    p1[1] = -19.178;
    p1[2] = 31.5947;
    p1[3] = -18.178;
    mean1 = mean0 + p1[0]+ p1[1]*x + p1[2]*x*x + p1[3]*x*x*x;
    mean2 = 0.; 
  }
  else return -1;
  return 0;
}

//________________________________________________________
Bool_t AliITSsadEdxFitter::IsGoodBin(Int_t bin,Int_t code){
  //method to select the bins used for the analysis only
  Bool_t retvalue=kTRUE;
  TLine *l[2];
  l[0]=new TLine(-2.1,0,2.,100.);
  l[1]=new TLine(-1.9,120,2.,0.);
  for(Int_t j=0;j<2;j++){
    l[j]->SetLineColor(4);
    l[j]->SetLineWidth(5);
    if(code==211 && (bin>14 || bin<1)){
      l[j]->Draw("same");	
      retvalue=kFALSE;
    }
    if(code==321 && (bin>12 || bin<5)){
      l[j]->Draw("same");
      retvalue=kFALSE;
    }
    if(code==2212 && (bin<8 || bin>16)){
      l[j]->Draw("same");	
      retvalue=kFALSE;
    }
  }
  return retvalue;
}

//________________________________________________________
Double_t SingleGausTail(const Double_t *x, const Double_t *par){
  //single gaussian with exponential tail
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  Double_t xx = x[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t xNormSquare1 = (xx - mean1) * (xx - mean1);
  Double_t sigmaSquare1 = sigma1 * sigma1;
  Double_t xdiv1 = mean1 + par[3] * sigma1;
  Double_t y1=0.0;
  if(xx < xdiv1) y1 = par[0]/(s2pi*par[2]) * TMath::Exp(-0.5 * xNormSquare1 / sigmaSquare1);
  else y1 = TMath::Exp(-par[4]*(xx-xdiv1)) * par[0] / (s2pi*par[2]) * TMath::Exp(-0.5*(par[3]*par[3]));
  return y1;
}

//________________________________________________________
Double_t DoubleGausTail(const Double_t *x, const Double_t *par){
  //sum of two gaussians with exponential tail
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  Double_t xx = x[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t xNormSquare1 = (xx - mean1) * (xx - mean1);
  Double_t sigmaSquare1 = sigma1 * sigma1;
  Double_t xdiv1 = mean1 + par[3] * sigma1;
  Double_t y1=0.0;
  if(xx < xdiv1) y1 = par[0]/(s2pi*par[2]) * TMath::Exp(-0.5 * xNormSquare1 / sigmaSquare1);
  else y1 = TMath::Exp(-par[4]*(xx-xdiv1)) * par[0] / (s2pi*par[2]) * TMath::Exp(-0.5*(par[3]*par[3]));

  Double_t mean2 = par[6];
  Double_t sigma2 = par[7];
  Double_t xNormSquare2 = (xx - mean2) * (xx - mean2);
  Double_t sigmaSquare2 = sigma2 * sigma2;
  Double_t xdiv2 = mean2 + par[8] * sigma2;
  Double_t y2=0.0;
  if(xx < xdiv2) y2 = par[5]/(s2pi*par[7]) * TMath::Exp(-0.5 * xNormSquare2 / sigmaSquare2);
  else y2 = TMath::Exp(-par[9]*(xx-xdiv2)) * par[5] / (s2pi*par[7]) * TMath::Exp(-0.5*(par[8]*par[8]));
  return y1+y2;
}

//________________________________________________________
Double_t FinalGausTail(const Double_t *x, const Double_t *par){
  //sum of three gaussians with exponential tail
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  Double_t xx = x[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t xNormSquare1 = (xx - mean1) * (xx - mean1);
  Double_t sigmaSquare1 = sigma1 * sigma1;
  Double_t xdiv1 = mean1 + par[3] * sigma1;
  Double_t y1=0.0;
  if(xx < xdiv1) y1 = par[0]/(s2pi*par[2]) * TMath::Exp(-0.5 * xNormSquare1 / sigmaSquare1);
  else y1 = TMath::Exp(-par[4]*(xx-xdiv1)) * par[0] / (s2pi*par[2]) * TMath::Exp(-0.5*(par[3]*par[3]));

  Double_t mean2 = par[6];
  Double_t sigma2 = par[7];
  Double_t xNormSquare2 = (xx - mean2) * (xx - mean2);
  Double_t sigmaSquare2 = sigma2 * sigma2;
  Double_t xdiv2 = mean2 + par[8] * sigma2;
  Double_t y2=0.0;
  if(xx < xdiv2) y2 = par[5]/(s2pi*par[7]) * TMath::Exp(-0.5 * xNormSquare2 / sigmaSquare2);
  else y2 = TMath::Exp(-par[9]*(xx-xdiv2)) * par[5] / (s2pi*par[7]) * TMath::Exp(-0.5*(par[8]*par[8]));

  Double_t mean3 = par[11];
  Double_t sigma3 = par[12];
  Double_t xNormSquare3 = (xx - mean3) * (xx - mean3);
  Double_t sigmaSquare3 = sigma3 * sigma3;
  Double_t xdiv3 = mean3 + par[13] * sigma3;
  Double_t y3=0.0;
  if(xx < xdiv3) y3 = par[10]/(s2pi*par[12]) * TMath::Exp(-0.5 * xNormSquare3 / sigmaSquare3);
  else y3 = TMath::Exp(-par[14]*(xx-xdiv3)) * par[10] / (s2pi*par[12]) * TMath::Exp(-0.5*(par[13]*par[13]));
  return y1+y2+y3;
}

//______________________________________________________________________
void AliITSsadEdxFitter::CalcResidual(TH1F *h,TF1 *fun,TGraph *gres) const{
  //code to calculate the difference fit function and histogram point (residual)
  //to use as goodness test for the fit
  Int_t ipt=0;
  Double_t x=0.,yhis=0.,yfun=0.;
  for(Int_t i=0;i<h->GetNbinsX();i++){
    x=(h->GetBinLowEdge(i+2)+h->GetBinLowEdge(i+1))/2;
    yfun=fun->Eval(x);
    yhis=h->GetBinContent(i+1);
    if(yhis>0. && yfun>0.) {
      gres->SetPoint(ipt,x,(yhis-yfun)/yhis);
      ipt++;
    }
  }
  return;
}


//________________________________________________________
Double_t SingleGausStep(const Double_t *x, const Double_t *par){
  //single normalized gaussian
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  Double_t xx = x[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t xNorm1Square = (xx - mean1) * (xx - mean1);
  Double_t sigma1Square = sigma1 * sigma1;
  Double_t step1 = par[0]/(s2pi*par[2]) * TMath::Exp(-0.5 * xNorm1Square / sigma1Square);
  return step1;
}

//________________________________________________________
Double_t DoubleGausStep(const Double_t *x, const Double_t *par){
  //sum of two normalized gaussians
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  Double_t xx = x[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t xNorm1Square = (xx - mean1) * (xx - mean1);
  Double_t sigma1Square = sigma1 * sigma1;
  Double_t step1 = par[0]/(s2pi*par[2]) * TMath::Exp( - 0.5 * xNorm1Square / sigma1Square );
  Double_t mean2 = par[4];
  Double_t sigma2 = par[5];
  Double_t xNorm2Square = (xx - mean2) * (xx - mean2);
  Double_t sigma2Square = sigma2 * sigma2;
  Double_t step2 = par[3]/(s2pi*par[5]) * TMath::Exp( - 0.5 * xNorm2Square / sigma2Square );
  return step1+step2;
}

//________________________________________________________
Double_t FinalGausStep(const Double_t *x, const Double_t *par){
  //sum of three normalized gaussians
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  Double_t xx = x[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t xNorm1Square = (xx - mean1) * (xx - mean1);
  Double_t sigma1Square = sigma1 * sigma1;
  Double_t step1 = par[0]/(s2pi*par[2]) * TMath::Exp( - 0.5 * xNorm1Square / sigma1Square );
  Double_t mean2 = par[4];
  Double_t sigma2 = par[5];
  Double_t xNorm2Square = (xx - mean2) * (xx - mean2);
  Double_t sigma2Square = sigma2 * sigma2;
  Double_t step2 = par[3]/(s2pi*par[5]) * TMath::Exp( - 0.5 * xNorm2Square / sigma2Square );
  Double_t mean3 = par[7];
  Double_t sigma3 = par[8];
  Double_t xNorm3Square = (xx - mean3) * (xx - mean3);
  Double_t sigma3Square = sigma3 * sigma3;
  Double_t step3 = par[6]/(s2pi*par[8]) * TMath::Exp( - 0.5 * xNorm3Square / sigma3Square);
  return step1+step2+step3;
}

//______________________________________________________________________
Double_t AliITSsadEdxFitter::GausPlusTail(const Double_t x, const Double_t mean, Double_t rms, Double_t c, Double_t slope, Double_t cut ) const{
  //gaussian with an exponential tail on the right side
  Double_t factor=1.0/(TMath::Sqrt(2.0*TMath::Pi()));
  Double_t returnvalue=0.0;
  Double_t n=0.5*(1.0+TMath::Erf(cut/TMath::Sqrt(2.0)))+TMath::Exp(-cut*cut*0.5)*factor/(TMath::Abs(rms)*slope);
  if (x<mean+cut*rms) returnvalue=TMath::Exp(-1.0*(x-mean)*(x-mean)/(2.0*rms*rms))*factor/TMath::Abs(rms);
  else returnvalue=TMath::Exp(slope*(mean+cut*rms-x))*TMath::Exp(-cut*cut*0.5)*factor/TMath::Abs(rms);
  return c*returnvalue/n;
}


//______________________________________________________________________
Double_t AliITSsadEdxFitter::GausOnBackground(const Double_t* x, const Double_t *par) const {
  //gaussian with a background parametrisation  
  //cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<" "<<par[5]<<" "<<x[0]<< endl;
  Double_t returnvalue=0.0;
  Double_t factor=1.0/(TMath::Sqrt(2.0*TMath::Pi()));
  if(par[6]<0.0) returnvalue+=TMath::Exp(-1.0*(x[0]-par[0])*(x[0]-par[0])/(2.0*par[1]*par[1]))*par[2]* factor/TMath::Abs(par[1]);
  else returnvalue+=GausPlusTail(x[0], par[0], par[1],par[2], par[6], 1.2);
  returnvalue+=par[3]*TMath::Exp((par[5]-x[0])*par[4]);
  return returnvalue;
}

//______________________________________________________________________
void AliITSsadEdxFitter::DrawFitFunction(TF1 *fun) const {
  //code to draw the final fit function and the single gaussians used
  //to extract the yields for the 3 species
  TF1 *fdraw1=new TF1("fdraw1",SingleGausStep,-3.5,3.5,3);
  TF1 *fdraw2=new TF1("fdraw2",SingleGausStep,-3.5,3.5,3);
  TF1 *fdraw3=new TF1("fdraw3",SingleGausStep,-3.5,3.5,3);
  fdraw1->SetParameter(0,fun->GetParameter(0));
  fdraw1->SetParameter(1,fun->GetParameter(1));
  fdraw1->SetParameter(2,fun->GetParameter(2));
  fdraw2->SetParameter(0,fun->GetParameter(3));
  fdraw2->SetParameter(1,fun->GetParameter(4));
  fdraw2->SetParameter(2,fun->GetParameter(5));
  fdraw3->SetParameter(0,fun->GetParameter(6));
  fdraw3->SetParameter(1,fun->GetParameter(7));
  fdraw3->SetParameter(2,fun->GetParameter(8));

  fdraw1->SetLineColor(6);//color code
  fdraw2->SetLineColor(2);
  fdraw3->SetLineColor(4);  

  fdraw1->SetLineStyle(2);//dot lines
  fdraw2->SetLineStyle(2);
  fdraw3->SetLineStyle(2);

  fun->SetLineWidth(3);
  fdraw1->DrawCopy("same");
  fdraw2->DrawCopy("same");
  fdraw3->DrawCopy("same");
  fun->DrawCopy("same");

  TLatex *ltx[3];
  for(Int_t j=0;j<3;j++){
    ltx[0]=new TLatex(0.13,0.25,"pions");
    ltx[1]=new TLatex(0.13,0.20,"kaons");
    ltx[2]=new TLatex(0.13,0.15,"protons");
    ltx[0]->SetTextColor(6);
    ltx[1]->SetTextColor(2);
    ltx[2]->SetTextColor(4);
    ltx[j]->SetNDC();
    ltx[j]->Draw();
  }
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::DoFit(TH1F *h, Int_t bin, Int_t signedcode, Bool_t mc, TGraph *gres){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: all free, second step: pion gaussian fixed, third step: kaon gaussian fixed
  // final step: refit all using the parameters and tollerance limits (+-20%)

  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  Double_t initialParametersStepTot[9];

  Int_t code=TMath::Abs(signedcode);
  //************************ drawing and label *******
  Double_t xbins[23]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  Double_t pt=(xbins[bin+1]+xbins[bin])/2;
  h->SetTitle(Form("p_{t}=[%1.2f,%1.2f], code=%d",xbins[bin],xbins[bin+1],signedcode));
  h->GetXaxis()->SetTitle("[ln dE/dx]_{meas} - [ln dE/dx(i)]_{calc}");
  h->GetYaxis()->SetTitle("counts");
  h->Draw("e");
  h->SetFillColor(11);
  Int_t xmax=-1,ymax=-1,zmax=-1;
  if(!IsGoodBin(bin,code)) return;
  h->GetMaximumBin(xmax,ymax,zmax);
  TString modfit = "M0R+";

  printf("\n$$$$$$$$$$$$$$$$$$$$$$$ BIN %d - hypothesis %d $$$$$$$$$$$$$$$$$$$$$$$$$$$\n",bin,code);
  Double_t ampl = h->GetMaximum()/(h->GetRMS()*s2pi);
  Double_t mean = h->GetBinLowEdge(xmax); //expected mean values
  Float_t expKaonMean=0., expProtonMean=0.;
  Int_t calcmean=CalcMean(code,pt,mean,expKaonMean,expProtonMean);
  if(calcmean<0) return;
  printf("mean -> %f %f %f\n",mean,expKaonMean,expProtonMean);
  printf("integration range -> (%1.2f,%1.2f) (%1.2f,%1.2f) (%1.2f,%1.2f)\n",fRangeStep1[0],fRangeStep1[1],fRangeStep2[0],fRangeStep2[1],fRangeStep3[0],fRangeStep3[1]);

  Float_t expPionSig = CalcSigma(211,pt,mc); //expected sigma values
  Float_t expKaonSig = CalcSigma(321,pt,mc);
  Float_t expProtonSig = CalcSigma(2212,pt,mc);
  printf("sigma -> %f %f %f\n",expPionSig,expKaonSig,expProtonSig);
  printf("sigma range -> (%1.2f,%1.2f) (%1.2f,%1.2f) (%1.2f,%1.2f)\n",fLimitsOnSigmaPion[0],fLimitsOnSigmaPion[1],fLimitsOnSigmaKaon[0],fLimitsOnSigmaKaon[1],fLimitsOnSigmaProton[0],fLimitsOnSigmaProton[1]);

  printf("___________________________________________________________________\n");
  printf("First Step: pions\n");
  fstep1 = new TF1("step1",SingleGausStep,-3.5,3.5,3);
  fstep1->SetParameter(0,ampl);//initial 
  fstep1->SetParameter(1,mean);
  fstep1->SetParLimits(0,0,ampl*1.2);
  fstep1->SetParLimits(1,mean+fRangeStep1[0],mean+fRangeStep1[1]);
  fstep1->SetParLimits(2,expPionSig*fLimitsOnSigmaPion[0],expPionSig*fLimitsOnSigmaPion[1]);
  //fstep1->FixParameter(2,expPionSig);
  if(expPionSig>0) h->Fit(fstep1,modfit.Data(),"",mean+fRangeStep1[0],mean+fRangeStep1[1]);//first fit
  else for(Int_t npar=0;npar<3;npar++) fstep1->FixParameter(npar,0.00001);

  printf("\n___________________________________________________________________\n");
  printf("Second Step: kaons\n"); 
  fstep2 = new TF1("fstep2",DoubleGausStep,-3.5,3.5,6);
  fstep2->FixParameter(0,fstep1->GetParameter(0));//fixed
  fstep2->FixParameter(1,fstep1->GetParameter(1));
  fstep2->FixParameter(2,fstep1->GetParameter(2));
  fstep2->SetParameter(3,fstep1->GetParameter(0)/8.);//initial
  fstep2->SetParameter(4,expKaonMean);
  fstep2->SetParLimits(3,fstep1->GetParameter(0)/100.,fstep1->GetParameter(0));//limits
  fstep2->SetParLimits(4,expKaonMean+fRangeStep2[0],expKaonMean+fRangeStep2[1]);
  fstep2->SetParLimits(5,expKaonSig*fLimitsOnSigmaKaon[0],expKaonSig*fLimitsOnSigmaKaon[1]);
  //fstep2->FixParameter(5,expKaonSig);

  if(expKaonSig>0) h->Fit(fstep2,modfit.Data(),"",expKaonMean+fRangeStep2[0],expKaonMean+fRangeStep2[1]);
  else for(Int_t npar=3;npar<6;npar++) fstep2->FixParameter(npar,0.00001);


  TLine *l[3];
  l[0] = new TLine(expKaonMean,0,expKaonMean,10000);
  l[1] = new TLine(expKaonMean+fRangeStep2[0],0,expKaonMean+fRangeStep2[0],10000);
  l[2] = new TLine(expKaonMean+fRangeStep2[1],0,expKaonMean+fRangeStep2[1],10000);
  for(Int_t dp=0;dp<3;dp++) {
    l[dp]->Draw("same");
    l[dp]->SetLineColor(4);
    l[dp]->SetLineWidth(4);
  }


  printf("\n____________________________________________________________________\n");
  printf("Third Step: protons \n");
  fstep3= new TF1("fstep3",FinalGausStep,-3.5,3.5,9);
  fstep3->FixParameter(0,fstep1->GetParameter(0));//fixed
  fstep3->FixParameter(1,fstep1->GetParameter(1));
  fstep3->FixParameter(2,fstep1->GetParameter(2));
  fstep3->FixParameter(3,fstep2->GetParameter(3));
  fstep3->FixParameter(4,fstep2->GetParameter(4));
  fstep3->FixParameter(5,fstep2->GetParameter(5));
  fstep3->SetParameter(6,fstep2->GetParameter(0)/16);//initial
  fstep3->SetParameter(7,expProtonMean);
  fstep3->SetParLimits(6,fstep2->GetParameter(0)/300,fstep2->GetParameter(0));//limits
  fstep3->SetParLimits(7,expProtonMean+fRangeStep3[0],expProtonMean+fRangeStep3[1]);
  fstep3->SetParLimits(8,expProtonSig*fLimitsOnSigmaProton[0],expProtonSig*fLimitsOnSigmaProton[1]);
  //fstep3->FixParameter(8,expProtonSig);
  if(expProtonSig>0) h->Fit(fstep3,modfit.Data(),"",expProtonMean+fRangeStep3[0],expProtonMean+fRangeStep3[1]);
  else for(Int_t npar=6;npar<9;npar++) fstep3->FixParameter(npar,0.00001);

  printf("\n_____________________________________________________________________\n");
  printf("Final Step: refit all \n"); 
  fstepTot = new TF1("funztot",FinalGausStep,-3.5,3.5,9);
  fstepTot->SetLineColor(1);
  initialParametersStepTot[0] = fstep1->GetParameter(0);//first gaussian
  initialParametersStepTot[1] = fstep1->GetParameter(1);
  initialParametersStepTot[2] = fstep1->GetParameter(2);

  initialParametersStepTot[3] = fstep2->GetParameter(3);//second gaussian
  initialParametersStepTot[4] = fstep2->GetParameter(4);
  initialParametersStepTot[5] = fstep2->GetParameter(5);

  initialParametersStepTot[6] = fstep3->GetParameter(6);//third gaussian
  initialParametersStepTot[7] = fstep3->GetParameter(7);
  initialParametersStepTot[8] = fstep3->GetParameter(8);  

  fstepTot->SetParameters(initialParametersStepTot);//initial parameter

  fstepTot->SetParLimits(0,initialParametersStepTot[0]*0.9,initialParametersStepTot[0]*1.1);//tolerance limit
  fstepTot->SetParLimits(1,initialParametersStepTot[1]*0.9,initialParametersStepTot[1]*1.1);
  fstepTot->SetParLimits(2,initialParametersStepTot[2]*0.9,initialParametersStepTot[2]*1.1);
  //fstepTot->FixParameter(2,initialParametersStepTot[2]);
  fstepTot->SetParLimits(3,initialParametersStepTot[3]*0.9,initialParametersStepTot[3]*1.1);
  fstepTot->SetParLimits(4,initialParametersStepTot[4]*0.9,initialParametersStepTot[4]*1.1);
  fstepTot->SetParLimits(5,initialParametersStepTot[5]*0.9,initialParametersStepTot[5]*1.1);
  //fstepTot->FixParameter(5,initialParametersStepTot[5]);
  fstepTot->SetParLimits(6,initialParametersStepTot[6]*0.9,initialParametersStepTot[6]*1.1);
  fstepTot->SetParLimits(7,initialParametersStepTot[7]*0.9,initialParametersStepTot[7]*1.1);
  fstepTot->SetParLimits(8,initialParametersStepTot[8]*0.9,initialParametersStepTot[8]*1.1);
  //fstepTot->FixParameter(8,initialParametersStepTot[8]);

  h->Fit(fstepTot,modfit.Data(),"",fRangeStep4[0],fRangeStep4[1]);//refit all


  //************************************* storing parameter to calculate the yields *******
  Int_t chpa=0;
  if(code==321) chpa=3;
  if(code==2212) chpa=6;
  for(Int_t j=0;j<3;j++) {
    fFitpar[j] = fstepTot->GetParameter(j+chpa);
    fFitparErr[j] = fstepTot->GetParError(j+chpa);
  }

  DrawFitFunction(fstepTot);
  CalcResidual(h,fstepTot,gres);
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::DoFitProton(TH1F *h, Int_t bin, Int_t signedcode, Bool_t mc, TGraph *gres){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: all free, second step: pion gaussian fixed, third step: kaon gaussian fixed
  // final step: refit all using the parameters and tollerance limits (+-20%)
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  Double_t initialParametersStepTot[9];

  Int_t code=TMath::Abs(signedcode);
  //************************ drawing and label *******
  Double_t xbins[23]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  Double_t pt=(xbins[bin+1]+xbins[bin])/2;
  h->SetTitle(Form("p_{t}=[%1.2f,%1.2f], code=%d",xbins[bin],xbins[bin+1],signedcode));
  h->GetXaxis()->SetTitle("[ln dE/dx]_{meas} - [ln dE/dx(i)]_{calc}");
  h->GetYaxis()->SetTitle("counts");
  h->Draw("e");
  h->SetFillColor(11);
  Int_t xmax=-1,ymax=-1,zmax=-1;
  h->GetMaximumBin(xmax,ymax,zmax);
  if(!IsGoodBin(bin,code)) return;

  printf("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ BIN %d - hypothesis %d $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n",bin,code);	
  Double_t ampl = h->GetMaximum()/(h->GetRMS()*s2pi);
  Double_t mean = h->GetBinCenter(xmax); //expected mean values
  Float_t expKaonMean=0., expProtonMean=0.;
  Int_t calcmean=CalcMean(code,pt,mean,expKaonMean,expProtonMean);
  if(calcmean<0) return;
  printf("mean -> %f %f %f\n",mean,expKaonMean,expProtonMean);
  printf("integration range -> (%1.2f,%1.2f) (%1.2f,%1.2f) (%1.2f,%1.2f)\n",fRangeStep1[0],fRangeStep1[1],fRangeStep2[0],fRangeStep2[1],fRangeStep3[0],fRangeStep3[1]);

  Float_t expPionSig = CalcSigma(211,pt,mc); //expected sigma values
  Float_t expKaonSig = CalcSigma(321,pt,mc);
  Float_t expProtonSig = CalcSigma(2212,pt,mc);
  printf("sigma -> %f %f %f\n",expPionSig,expKaonSig,expProtonSig);
  printf("sigma range -> (%1.2f,%1.2f) (%1.2f,%1.2f) (%1.2f,%1.2f)\n",fLimitsOnSigmaPion[0],fLimitsOnSigmaPion[1],fLimitsOnSigmaKaon[0],fLimitsOnSigmaKaon[1],fLimitsOnSigmaProton[0],fLimitsOnSigmaProton[1]);

  printf("___________________________________________________________________\n");
  printf("First Step: pions\n\n");
  fstep1 = new TF1("step1",SingleGausStep,-3.5,3.5,3);
  fstep1->SetParameter(0,ampl);//initial 
  fstep1->SetParameter(1,mean);
  fstep1->SetParLimits(0,0,ampl*1.2);
  fstep1->SetParLimits(1,mean+fRangeStep1[0],mean+fRangeStep1[1]);
  fstep1->SetParLimits(2,expPionSig*fLimitsOnSigmaPion[0],expPionSig*fLimitsOnSigmaPion[1]);
  //fstep1->FixParameter(2,expPionSig);
  if(expPionSig>0)  h->Fit(fstep1,"M0R+","",mean+fRangeStep1[0],mean+fRangeStep1[1]);//first fit
  else for(Int_t npar=0;npar<3;npar++) fstep1->FixParameter(npar,0.00001);

  printf("\n___________________________________________________________________\n");
  printf("Second Step: protons\n\n"); 
  fstep2 = new TF1("step2",SingleGausStep,-3.5,3.5,3);
  fstep2->SetParameter(0,fstep1->GetParameter(0)/8.);//initial
  fstep2->SetParameter(1,expProtonMean);
  fstep2->SetParLimits(0,fstep1->GetParameter(0)/100.,fstep1->GetParameter(0));//limits
  fstep2->SetParLimits(1,-3.5,3.5);
  fstep2->SetParLimits(2,expProtonSig*fLimitsOnSigmaProton[0],expProtonSig*fLimitsOnSigmaProton[1]);
  //fstep2->FixParameter(2,expProtonSig);
  if(expProtonSig>0) h->Fit(fstep2,"M0R+","",expProtonMean+fRangeStep3[0],expProtonMean+fRangeStep3[1]);
  else for(Int_t npar=0;npar<3;npar++) fstep2->FixParameter(npar,0.00001);

  printf("\n____________________________________________________________________\n");
  printf("Third Step: kaons \n\n");
  fstep3= new TF1("fstep3",FinalGausStep,-3.5,3.5,9);
  fstep3->FixParameter(0,fstep1->GetParameter(0));//fixed
  fstep3->FixParameter(1,fstep1->GetParameter(1));
  fstep3->FixParameter(2,fstep1->GetParameter(2));
  fstep3->SetParameter(3,fstep1->GetParameter(0)/10);//initial
  fstep3->SetParameter(4,expKaonMean);
  fstep3->FixParameter(6,fstep2->GetParameter(0));
  fstep3->FixParameter(7,fstep2->GetParameter(1));
  fstep3->FixParameter(8,fstep2->GetParameter(2));
  fstep3->SetParLimits(3,fstep1->GetParameter(0)/100.,fstep1->GetParameter(0));//limits
  fstep3->SetParLimits(4,-3.5,3.5);
  fstep3->SetParLimits(5,expKaonSig*fLimitsOnSigmaKaon[0],expKaonSig*fLimitsOnSigmaKaon[1]);
  //fstep3->FixParameter(5,expKaonSig);

  TLine *l[3];
  l[0] = new TLine(expProtonMean,0,expProtonMean,10000);
  l[1] = new TLine(expProtonMean+fRangeStep3[0],0,expProtonMean+fRangeStep3[0],10000);
  l[2] = new TLine(expProtonMean+fRangeStep3[1],0,expProtonMean+fRangeStep3[1],10000);
  for(Int_t dp=0;dp<3;dp++) {
    l[dp]->Draw("same");
    l[dp]->SetLineColor(2);
    l[dp]->SetLineWidth(4);
  }

  if(expKaonSig>0) h->Fit(fstep3,"M0R+","",expKaonMean+fRangeStep2[0],expKaonMean+fRangeStep2[1]);
  else for(Int_t npar=3;npar<6;npar++) fstep3->FixParameter(npar,0.00001);

  printf("\n_____________________________________________________________________\n");
  printf("Final Step: refit all \n\n"); 
  fstepTot = new TF1("funztot",FinalGausStep,-3.5,3.5,9);
  fstepTot->SetLineColor(1);
  initialParametersStepTot[0] = fstep1->GetParameter(0);//first gaussian
  initialParametersStepTot[1] = fstep1->GetParameter(1);
  initialParametersStepTot[2] = fstep1->GetParameter(2);

  initialParametersStepTot[3] = fstep3->GetParameter(3);//second gaussian
  initialParametersStepTot[4] = fstep3->GetParameter(4);
  initialParametersStepTot[5] = fstep3->GetParameter(5);

  initialParametersStepTot[6] = fstep2->GetParameter(0);//third gaussian
  initialParametersStepTot[7] = fstep2->GetParameter(1);
  initialParametersStepTot[8] = fstep2->GetParameter(2);  

  fstepTot->SetParameters(initialParametersStepTot);//initial parameter

  fstepTot->SetParLimits(0,initialParametersStepTot[0]*0.9,initialParametersStepTot[0]*1.1);//tolerance limit
  fstepTot->SetParLimits(1,initialParametersStepTot[1]*0.9,initialParametersStepTot[1]*1.1);
  fstepTot->SetParLimits(2,initialParametersStepTot[2]*0.9,initialParametersStepTot[2]*1.1);
  //fstepTot->FixParameter(2,initialParametersStepTot[2]);
  fstepTot->SetParLimits(3,initialParametersStepTot[3]*0.9,initialParametersStepTot[3]*1.1);
  fstepTot->SetParLimits(4,initialParametersStepTot[4]*0.9,initialParametersStepTot[4]*1.1);
  fstepTot->SetParLimits(5,initialParametersStepTot[5]*0.9,initialParametersStepTot[5]*1.1);
  //fstepTot->FixParameter(5,initialParametersStepTot[5]);
  fstepTot->SetParLimits(6,initialParametersStepTot[6]*0.9,initialParametersStepTot[6]*1.1);
  fstepTot->SetParLimits(7,initialParametersStepTot[7]*0.9,initialParametersStepTot[7]*1.1);
  fstepTot->SetParLimits(8,initialParametersStepTot[8]*0.9,initialParametersStepTot[8]*1.1);
  //fstepTot->FixParameter(8,initialParametersStepTot[8]);

  h->Fit(fstepTot,"M0R+","",fRangeStep4[0],fRangeStep4[1]);//refit all


  //************************************* storing parameter to calculate the yields *******
  Int_t chpa=0;
  if(code==321) chpa=3;
  if(code==2212) chpa=6;
  for(Int_t j=0;j<3;j++) {
    fFitpar[j] = fstepTot->GetParameter(j+chpa);
    fFitparErr[j] = fstepTot->GetParError(j+chpa);
  }

  DrawFitFunction(fstepTot);
  CalcResidual(h,fstepTot,gres);
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::DoFitTail(TH1F *h, Int_t bin, Int_t code){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: all free, second step: pion gaussian fixed, third step: kaon gaussian fixed
  // final step: refit all using the parameters and tollerance limits (+-20%)
  // WARNING: exponential tail added in the right of the Gaussian shape
  Double_t s2pi=TMath::Sqrt(2*TMath::Pi());
  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  Double_t initialParametersStepTot[15];


  //************************ drawing and label *******
  Double_t xbins[23]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  h->SetTitle(Form("p_{t}=[%1.2f,%1.2f], code=%d",xbins[bin],xbins[bin+1],code));
  h->GetXaxis()->SetTitle("[ln dE/dx]_{meas} - [ln dE/dx(i)]_{calc}");
  h->GetYaxis()->SetTitle("counts");
  h->Draw("");
  h->SetFillColor(11);
  Int_t xmax=-1,ymax=-1,zmax=-1;
  Int_t hmax=h->GetMaximumBin(xmax,ymax,zmax);
  hmax++;
  if(!IsGoodBin(bin,code)) return;
  Double_t ampl = h->GetMaximum()/(h->GetRMS()*s2pi);
  Double_t mean = h->GetBinLowEdge(xmax);

  printf("\n___________________________________________________________________\n First Step: pions\n\n");
  fstep1 = new TF1("step1",SingleGausTail,-3.5,3.5,5);
  fstep1->SetParameter(0,ampl);//initial 
  fstep1->SetParameter(1,mean);
  fstep1->SetParameter(3,1.2);
  fstep1->SetParameter(4,10.);

  fstep1->SetParLimits(0,0,ampl*1.2);
  fstep1->SetParLimits(1,-3.5,3.5);
  fstep1->SetParLimits(2,0.1,0.25);
  fstep1->SetParLimits(4,5.,20.);
  if(bin<8) fstep1->SetParLimits(4,13.,25.);

  h->Fit(fstep1,"M0R+","",mean-0.45,mean+0.45);//first fit

  printf("\n___________________________________________________________________\n Second Step: kaons\n\n"); 
  fstep2 = new TF1("fstep2",DoubleGausTail,-3.5,3.5,10);
  fstep2->FixParameter(0,fstep1->GetParameter(0));//fixed
  fstep2->FixParameter(1,fstep1->GetParameter(1));
  fstep2->FixParameter(2,fstep1->GetParameter(2));
  fstep2->FixParameter(3,fstep1->GetParameter(3));
  fstep2->FixParameter(4,fstep1->GetParameter(4));

  fstep2->SetParameter(5,fstep1->GetParameter(0)/8);//initial
  //fstep2->SetParameter(6,CalcP(code,322,bin));
  fstep2->SetParameter(7,0.145);
  fstep2->FixParameter(8,1.2);
  fstep2->SetParameter(9,13.);

  fstep2->SetParLimits(5,fstep1->GetParameter(0)/100,fstep1->GetParameter(0));//limits
  fstep2->SetParLimits(6,-3.5,3.5);
  fstep2->SetParLimits(7,0.12,0.2);
  fstep2->SetParLimits(9,9.,20.);
  if(bin<9) fstep2->SetParLimits(9,13.,25.);

  //h->Fit(fstep2,"M0R+","",CalcP(code,321,bin)-0.3,CalcP(code,321,bin)+0.3);//second fit
  if(bin<6 || bin>12) for(Int_t npar=5;npar<10;npar++) fstep2->FixParameter(npar,-0.0000000001);

  printf("\n____________________________________________________________________\n Third Step: protons \n\n");
  fstep3= new TF1("fstep3",FinalGausTail,-3.5,3.5,15);
  fstep3->FixParameter(0,fstep1->GetParameter(0));//fixed
  fstep3->FixParameter(1,fstep1->GetParameter(1));
  fstep3->FixParameter(2,fstep1->GetParameter(2));
  fstep3->FixParameter(3,fstep1->GetParameter(3));
  fstep3->FixParameter(4,fstep1->GetParameter(4));
  fstep3->FixParameter(5,fstep2->GetParameter(5));
  fstep3->FixParameter(6,fstep2->GetParameter(6));
  fstep3->FixParameter(7,fstep2->GetParameter(7));
  fstep3->FixParameter(8,fstep2->GetParameter(8));
  fstep3->FixParameter(9,fstep2->GetParameter(9));

  fstep3->SetParameter(10,fstep2->GetParameter(0)/8);//initial
  //fstep3->SetParameter(11,CalcP(code,2212,bin));
  fstep3->SetParameter(12,0.145);
  fstep3->FixParameter(13,1.2);
  fstep3->SetParameter(14,10.);

  fstep3->SetParLimits(10,fstep2->GetParameter(0)/100,fstep2->GetParameter(0));//limits
  fstep3->SetParLimits(11,-3.5,3.5);
  fstep3->SetParLimits(12,0.12,0.2);
  fstep3->SetParLimits(14,11.,25.);

  //h->Fit(fstep3,"M0R+","",CalcP(code,2212,bin)-0.3,CalcP(code,2212,bin)+0.3);//third fit

  printf("\n_____________________________________________________________________\n Final Step: refit all \n\n"); 
  fstepTot = new TF1("funztot",FinalGausTail,-3.5,3.5,15);
  fstepTot->SetLineColor(1);
  initialParametersStepTot[0] = fstep1->GetParameter(0);//first gaussian
  initialParametersStepTot[1] = fstep1->GetParameter(1);
  initialParametersStepTot[2] = fstep1->GetParameter(2);
  initialParametersStepTot[3] = fstep1->GetParameter(3);
  initialParametersStepTot[4] = fstep1->GetParameter(4);

  initialParametersStepTot[5] = fstep2->GetParameter(5);//second gaussian
  initialParametersStepTot[6] = fstep2->GetParameter(6);
  initialParametersStepTot[7] = fstep2->GetParameter(7);
  initialParametersStepTot[8] = fstep2->GetParameter(8);
  initialParametersStepTot[9] = fstep2->GetParameter(9);

  initialParametersStepTot[10] = fstep3->GetParameter(10);//third gaussian
  initialParametersStepTot[11] = fstep3->GetParameter(11);
  initialParametersStepTot[12] = fstep3->GetParameter(12);  
  initialParametersStepTot[13] = fstep3->GetParameter(13);  
  initialParametersStepTot[14] = fstep3->GetParameter(14);  

  fstepTot->SetParameters(initialParametersStepTot);//initial parameter


  fstepTot->SetParLimits(0,initialParametersStepTot[0]*0.9,initialParametersStepTot[0]*1.1);//tollerance limit
  fstepTot->SetParLimits(1,initialParametersStepTot[1]*0.9,initialParametersStepTot[1]*1.1);
  fstepTot->SetParLimits(2,initialParametersStepTot[2]*0.9,initialParametersStepTot[2]*1.1);
  fstepTot->SetParLimits(3,initialParametersStepTot[3]*0.9,initialParametersStepTot[3]*1.1);
  fstepTot->SetParLimits(4,initialParametersStepTot[4]*0.9,initialParametersStepTot[4]*1.1);
  fstepTot->SetParLimits(5,initialParametersStepTot[5]*0.9,initialParametersStepTot[5]*1.1);
  fstepTot->SetParLimits(6,initialParametersStepTot[6]*0.9,initialParametersStepTot[6]*1.1);
  fstepTot->SetParLimits(7,initialParametersStepTot[7]*0.9,initialParametersStepTot[7]*1.1);
  fstepTot->SetParLimits(8,initialParametersStepTot[8]*0.9,initialParametersStepTot[8]*1.1);
  fstepTot->SetParLimits(9,initialParametersStepTot[9]*0.9,initialParametersStepTot[9]*1.1);
  fstepTot->SetParLimits(10,initialParametersStepTot[10]*0.9,initialParametersStepTot[10]*1.1);
  fstepTot->SetParLimits(11,initialParametersStepTot[11]*0.9,initialParametersStepTot[11]*1.1);
  fstepTot->SetParLimits(12,initialParametersStepTot[12]*0.9,initialParametersStepTot[12]*1.1);
  fstepTot->SetParLimits(13,initialParametersStepTot[13]*0.9,initialParametersStepTot[13]*1.1);
  fstepTot->SetParLimits(14,initialParametersStepTot[14]*0.9,initialParametersStepTot[14]*1.1);

  if(bin<9) for(Int_t npar=10;npar<15;npar++) fstepTot->FixParameter(npar,-0.00000000001);

  h->Fit(fstepTot,"M0R+","",-3.5,3.5); //refit all


  //************************************* storing parameter to calculate the yields *******
  Int_t chpa=0;
  if(code==321) chpa=5;
  if(code==2212) chpa=10;
  for(Int_t j=0;j<5;j++) {
    fFitpar[j] = fstepTot->GetParameter(j+chpa);
    fFitparErr[j] = fstepTot->GetParError(j+chpa);
  }

  DrawFitFunction(fstepTot);
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::FillHisto(TH1F *hsps, Int_t bin, Float_t binsize, Int_t code){
  // fill the spectra histo calculating the yield
  // first bin has to be 1
  Double_t yield = 0;
  Double_t err = 0;
  Double_t ptbin = hsps->GetBinLowEdge(bin+1) - hsps->GetBinLowEdge(bin); 

  if(IsGoodBin(bin-1,code)) {
    yield = fFitpar[0] / ptbin / binsize;
    err = fFitparErr[0] / ptbin / binsize;
  }

  hsps->SetBinContent(bin,yield);
  hsps->SetBinError(bin,err);
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::FillHistoMC(TH1F *hsps, Int_t bin, Int_t code, TH1F *h){
  // fill the spectra histo calculating the yield (using the MC truth)
  // first bin has to be 1
  Double_t yield = 0;
  Double_t erryield=0;
  Double_t ptbin = hsps->GetBinLowEdge(bin+1) - hsps->GetBinLowEdge(bin);

  if(IsGoodBin(bin-1,code)){
    yield = h->GetEntries() / ptbin;
    erryield=TMath::Sqrt(h->GetEntries()) / ptbin;
  }
  hsps->SetBinContent(bin,yield);
  hsps->SetBinError(bin,erryield);
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::GetFitPar(Double_t *fitpar, Double_t *fitparerr) const {
  // getter of the fit parameters and the relative errors
  for(Int_t i=0;i<5;i++) {
    fitpar[i] = fFitpar[i];
    fitparerr[i] = fFitparErr[i];
  }
  return;
}


//________________________________________________________
void AliITSsadEdxFitter::PrintAll() const{
  printf("Range 1 = %f %f\n",fRangeStep1[0],fRangeStep1[1]);
  printf("Range 2 = %f %f\n",fRangeStep2[0],fRangeStep2[1]);
  printf("Range 3 = %f %f\n",fRangeStep3[0],fRangeStep3[1]);
  printf("Range F = %f %f\n",fRangeStep4[0],fRangeStep4[1]);
  printf(" Sigma1 = %f %f\n",fLimitsOnSigmaPion[0],fLimitsOnSigmaPion[1]);
  printf(" Sigma2 = %f %f\n",fLimitsOnSigmaKaon[0],fLimitsOnSigmaKaon[1]);
  printf(" Sigma3 = %f %f\n",fLimitsOnSigmaProton[0],fLimitsOnSigmaProton[1]);
}
