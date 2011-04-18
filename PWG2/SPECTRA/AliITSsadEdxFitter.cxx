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

/* $Id$ */

///////////////////////////////////////////////////////////////////////
// Class with the fits algorithms to be used in the identified       //
// spectra analysis using the ITS in stand-alone mode                //
// Author: E.Biolcati, biolcati@to.infn.it                           //
//         F.Prino, prino@to.infn.it                                 //
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
AliITSsadEdxFitter::AliITSsadEdxFitter():TObject(),
  fExpPionMean(0.),
  fExpKaonMean(0.),
  fExpProtonMean(0.),
  fExpPionMeanRange(0.),
  fExpKaonMeanRange(0.),
  fExpProtonMeanRange(0.),
  fExpPionSigma(0.),
  fExpKaonSigma(0.),
  fExpProtonSigma(0.),
  fExpPionAmpl(0.),
  fIsMC(kFALSE),
  fOptInit(0),
  fITSpid(0)
{
  // standard constructor
  for(Int_t i=0; i<9; i++)  fFitPars[i] = 0.;
  for(Int_t i=0; i<9; i++)  fFitParErrors[i] = 0.;
  for(Int_t i=0; i<3; i++)  fFitpar[i] = 0.;
  for(Int_t i=0; i<3; i++)  fFitparErr[i] = 0.;
  SetRangeStep1();
  SetRangeStep2();
  SetRangeStep3();
  SetRangeFinalStep();
  SetLimitsOnSigmaPion();
  SetLimitsOnSigmaKaon();
  SetLimitsOnSigmaProton();
  SetBinsUsedPion();
  SetBinsUsedKaon();
  SetBinsUsedProton();
  fITSpid=new AliITSPIDResponse(kFALSE);
};
//______________________________________________________________________
AliITSsadEdxFitter::AliITSsadEdxFitter(Bool_t isMC):TObject(),
  fExpPionMean(0.),
  fExpKaonMean(0.),
  fExpProtonMean(0.),
  fExpPionMeanRange(0.),
  fExpKaonMeanRange(0.),
  fExpProtonMeanRange(0.),
  fExpPionSigma(0.),
  fExpKaonSigma(0.),
  fExpProtonSigma(0.),
  fExpPionAmpl(0.),
  fIsMC(isMC), 
  fOptInit(0),
  fITSpid(0)
{
  // standard constructor
  for(Int_t i=0; i<9; i++)  fFitPars[i] = 0.;
  for(Int_t i=0; i<9; i++)  fFitParErrors[i] = 0.;
  for(Int_t i=0; i<3; i++)  fFitpar[i] = 0.;
  for(Int_t i=0; i<3; i++)  fFitparErr[i] = 0.;
  SetRangeStep1();
  SetRangeStep2();
  SetRangeStep3();
  SetRangeFinalStep();
  SetLimitsOnSigmaPion();
  SetLimitsOnSigmaKaon();
  SetLimitsOnSigmaProton();
  SetBinsUsedPion();
  SetBinsUsedKaon();
  SetBinsUsedProton();
  fITSpid=new AliITSPIDResponse(isMC);
};

//________________________________________________________
Double_t AliITSsadEdxFitter::CalcSigma(Int_t code,Float_t x) const {
  // calculate the sigma 12-ott-2010  
  Double_t p[2]={0.};
  Double_t xMinKaon=0.15; //minimum pt value to consider the kaon peak
  Double_t xMinProton=0.3;//minimum pt value to consider the proton peak
  if(fIsMC){
    if(code==211){
      p[0] = -1.20337e-04;
      p[1] = 1.13060e-01;
    }    
    else if(code==321 && x>xMinKaon){
      p[0] = -2.39631e-03;
      p[1] = 1.15723e-01;
    }    
    else if(code==2212 && x>xMinProton){
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
    else if(code==321 && x>xMinKaon){
      p[0] = -6.22639e-04;
      p[1] = 1.43289e-01;
    }    
    else if(code==2212 && x>xMinProton){
      p[0] = -2.13243e-03;
      p[1] = 1.68614e-01;
    } 
    else return -1;
  }
  return p[0]/(x*x)*TMath::Log(x)+p[1];
}

//_______________________________________________________
Int_t AliITSsadEdxFitter::InitializeMeanPosParam(Int_t code, Float_t x){
  // Set initial values for peak positions from ad-hoc parameterizations (Melinda, 5 april 2010)
  
  Double_t p1[5]={0.};
  Double_t p2[5]={0.};
  Double_t p3[5]={0.};
  Double_t ep1[5]={0.};
  Double_t ep2[5]={0.};
  Double_t ep3[5]={0.};
  if(!fIsMC){ // data parameters
    if(code==211){

      p1[0]=-4.21001e-04; p1[1]=-3.41133e-02; p1[2]=0.; p1[3]=0.; p1[4]=0.;
      ep1[0]= 3.62105e-06; ep1[1]= 1.52077e-04; ep1[2]= 0.; ep1[3]= 0.; ep1[4]= 0.;

      p2[0] = 8.04038e+00; p2[1] = 1.51139e-01; p2[2] = 8.89984e-02; p2[3]= 0.; p2[4]= 0.;
      ep2[0]= 6.47330e-02; ep2[1]= 3.09501e-03; ep2[2]= 6.03642e-04; ep2[3]= 0.; ep2[4]= 0.;

      p3[0] = 7.41939e+00; p3[1] =-1.31031e+00; p3[2] = 9.44093e-01; p3[3]= 0.; p3[4]= 0.;
      ep3[0]= 1.29862e-01; ep3[1]= 1.37370e-02; ep3[2]= 4.74995e-03; ep3[3]= 0.; ep3[4]= 0.;

      fExpPionMean=p1[0]/(x*x)*TMath::Log(x)+p1[1];
      fExpPionMeanRange=TMath::Sqrt(TMath::Power(1/(x*x)*TMath::Log(x),2.)*ep1[0]*ep1[0]+ep1[1]*ep1[1] );
      fExpKaonMean=p2[0]*TMath::Landau(x,p2[1],p2[2],kFALSE);
      fExpKaonMeanRange=TMath::Sqrt(ep2[0]*ep2[0]+(x*x)*ep2[1]*ep2[1]+(x*x*x*x)*ep2[2]*ep2[2]);
      fExpProtonMean=-9999.;
      fExpProtonMeanRange=0.;
      if(x>0.2){
	fExpProtonMean=p3[0]*TMath::Exp(-0.5*TMath::Power((x-p3[1])/p3[2],2.)) ;
	fExpProtonMeanRange=TMath::Sqrt(ep3[0]*ep3[0]+(x*x)*ep3[1]*ep3[1]+(x*x*x*x)*ep3[2]*ep3[2]);
      }

    }else if(code==321){

      p1[0] = 8.35645e-01; p1[1] = 3.61638e+00; p1[2] =-2.51068e-01; p1[3]= 0.; p1[4]= 0.;
      ep1[0]= 6.63101e-04; ep1[1]= 2.33140e-02; ep1[2]= 1.93276e-03; ep1[3]= 0.; ep1[4]= 0.;

      p2[0] = -4.75573e-02; p2[1] = -7.05764e-03; p2[2] =  7.96719e+00; p2[3]=  -7.48333e-02; p2[4]= 0.;
      ep2[0]= 1.70759e-02; ep2[1]= 2.95064e-04; ep2[2]= 2.70493e+00; ep2[3]= 1.14466e-02; ep2[4]= 0.;

      p3[0] = 4.20188e+00; p3[1] = 3.95855e-01; p3[2] = 2.01159e-01; p3[3]=  0.; p3[4]=  0.;
      ep3[0]= 2.18608e-02; ep3[1]= 4.15293e-03; ep3[2]= 8.91846e-03; ep3[3]= 0.; ep3[4]= 0.;

      fExpPionMean= p1[0]*TMath::Log(x)*TMath::Exp(-x*x*p1[1])+p1[2];
      fExpPionMeanRange=TMath::Sqrt(ep1[0]*ep1[0]+(x*x)*ep1[1]*ep1[1]+(x*x*x*x)*ep1[2]*ep1[2]);
      fExpKaonMean = p2[0]*TMath::Log(x)+p2[1]*TMath::Exp(-x*x*p2[2])/(x*x)+p2[3];
      fExpKaonMeanRange=TMath::Sqrt(ep2[0]*ep2[0]+(x*x)*ep2[1]*ep2[1]+(x*x*x*x)*ep2[2]*ep2[2]);
      fExpProtonMean = p3[0]*TMath::Landau(x, p3[1], p3[2],kFALSE);
      fExpProtonMeanRange=TMath::Sqrt(ep3[0]*ep3[0]+(x*x)*ep3[1]*ep3[1]+(x*x*x*x)*ep3[2]*ep3[2]);

    }else if(code==2212){

      p1[0] =1.48277e+00; p1[1] =1.14050e+00; p1[2] =1.01751e+02; p1[3]= -2.52768e-01; p1[4]= 0.;
      ep1[0]= 1.85491e-03; ep1[1]= 4.23200e-02; ep1[2]= 2.66855e+00; ep1[3]= 1.82240e-03; ep1[4]= 0.;
      
      p2[0]= 3.80541e-01 ; p2[1]=-3.79158e-02; p2[2]=-4.81653e-01; p2[3]= 0.; p2[4]= 0.;      
      ep2[0]= 1.43643e-01; ep2[1]= 2.47156e-02; ep2[2]= 1.06836e-01; ep2[3]= 0.; ep2[4]= 0.;
      
      p3[0] = 1.81168e-01; p3[1] = -6.69364e-01; p3[2] =  2.10685e+01; p3[3]=   5.14868e-02; p3[4]=   0.;
      ep3[0]= 4.61718e-02; ep3[1]= 8.21949e-02; ep3[2]= 4.46134e+00; ep3[3]= 2.69751e-02; ep3[4]= 0.;

      fExpPionMean= p1[0]*TMath::Log(x)+p1[1]*TMath::Exp(-x*x*p1[2])+p1[3];
      fExpPionMeanRange=TMath::Sqrt(ep1[0]*ep1[0]+(x*x)*ep1[1]*ep1[1]+(x*x*x*x)*ep1[2]*ep1[2]);
      fExpKaonMean = p2[0]*TMath::Log(x)+p2[1]/(x)+p2[2];
      fExpKaonMeanRange=TMath::Sqrt(ep2[0]*ep2[0]+(x*x)*ep2[1]*ep2[1]+(x*x*x*x)*ep2[2]*ep2[2]+(x*x*x*x*x*x)*ep2[3]*ep2[3]+(x*x*x*x*x*x*x*x)*ep2[4]*ep2[4]);
      fExpProtonMean= p3[0]*TMath::Log(x)+p3[1]*TMath::Exp(-x*x*p3[2])+p3[3];
      fExpProtonMeanRange=TMath::Sqrt(ep3[0]*ep3[0]+(x*x)*ep3[1]*ep3[1]+(x*x*x*x)*ep3[2]*ep3[2]+(x*x*x*x*x*x)*ep3[3]*ep3[3]+(x*x*x*x*x*x*x*x)*ep3[4]*ep3[4]);
    }else{
      printf("ERROR: Wrong particle code %d\n",code);
      return -1;
    }
  }else{ // MC parameters
    if(code==211){

      p1[0]= -3.28744e-04; p1[1]= -1.78123e-02; p1[2]= 0.; p1[3]= 0.; p1[4]= 0.;
      ep1[0]= 2.03102e-06; ep1[1]= 7.65859e-05; ep1[2]= 0.; ep1[3]= 0.; ep1[4]= 0.;    

      p2[0] = 8.84991e+00 ; p2[1] = 1.32096e-01 ; p2[2] = 9.42394e-02 ; p2[3]=0.; p2[4]=0.;
      ep2[0]= 4.78123e-02; ep2[1]= 1.99037e-03; ep2[2]= 2.85136e-04; ep2[3]= 0.; ep2[4]= 0.;

      p3[0] =  6.25025e+00; p3[1] = -1.08315e+00; p3[2] =  8.86474e-01; p3[3]= 0.; p3[4]= 0.;
      ep3[0]= 1.42379e+00; ep3[1]= 2.34807e-01; ep3[2]= 6.78015e-02; ep3[3]= 0.; ep3[4]= 0.;

      fExpPionMean= p1[0]/(x*x)*TMath::Log(x)+p1[1];
      fExpPionMeanRange=TMath::Sqrt(TMath::Power(1/(x*x)*TMath::Log(x),2.)*ep1[0]*ep1[0]+ep1[1]*ep1[1] );
      fExpKaonMean = p2[0]*TMath::Landau(x,p2[1],p2[2],kFALSE);
      fExpKaonMeanRange=TMath::Sqrt(ep2[0]*ep2[0]+(x*x)*ep2[1]*ep2[1]+(x*x*x*x)*ep2[2]*ep2[2]);
      fExpProtonMean=-9999.;
      fExpProtonMeanRange=0.;
      if(x>0.2){
	fExpProtonMean = p3[0]*TMath::Exp(-0.5*TMath::Power((x-p3[1])/p3[2],2.));
	fExpProtonMeanRange = TMath::Sqrt(ep3[0]*ep3[0]+(x*x)*ep3[1]*ep3[1]+(x*x*x*x)*ep3[2]*ep3[2]);
      }

    }else if(code==321){

      p1[0] = 9.19679e-01; p1[1] = 5.07077e-03; p1[2] = 5.48262e-05; p1[3]=  8.65674e-02; p1[4]= 0.;
      ep1[0]= 4.03325e-04; ep1[1]= 8.85421e-06; ep1[2]= 1.77437e+00; ep1[3]= 9.00492e-03; ep1[4]= 0.;

      p2[0] = -2.05127e-01; p2[1] = -2.05303e-03; p2[2] =  3.10090e-06; p2[3]=  -2.19936e-01; p2[4]= 0.;
      ep2[0]= 3.41607e-03; ep2[1]= 1.03661e-04; ep2[2]= 9.50776e-01; ep2[3]= 3.41148e-03; ep2[4]= 0.;

      p3[0] = 4.68022e+00; p3[1] = 3.02178e-01; p3[2] = 2.47185e-01; p3[3]= 0.; p3[4]= 0.;
      ep3[0]= 1.63926e-02; ep3[1]= 6.42072e-03; ep3[2]= 7.24079e-03; ep3[3]= 0.; ep3[4]= 0.;

      fExpPionMean=p1[0]*TMath::Log(x)+p1[1]*TMath::Exp(-x*x*p1[2])/(x*x)+p1[3];
      fExpPionMeanRange=TMath::Sqrt(ep1[0]*ep1[0]+(x*x)*ep1[1]*ep1[1]+(x*x*x*x)*ep1[2]*ep1[2]);
      fExpKaonMean=p2[0]*TMath::Log(x)+p2[1]*TMath::Exp(-x*x*p2[2])/(x*x+p2[3]);
      fExpKaonMeanRange=TMath::Sqrt(ep2[0]*ep2[0]+(x*x)*ep2[1]*ep2[1]+(x*x*x*x)*ep2[2]*ep2[2]);
      fExpProtonMean= p3[0]*TMath::Landau(x, p3[1], p3[2]);
      fExpProtonMeanRange=TMath::Sqrt(ep3[0]*ep3[0]+(x*x)*ep3[1]*ep3[1]+(x*x*x*x)*ep3[2]*ep3[2]);

    }else if(code==2212){

      p1[0] = 1.05826e+00; p1[1]=8.59088e-01; p1[2]=9.39446e+01; p1[3]= -4.86736e-01; p1[4]=  0.;
      ep1[0]= 1.07813e-03; ep1[1]= 1.15536e-02; ep1[2]= 9.74843e-01; ep1[3]= 1.24205e-03; ep1[4]= 0.;

      p2[0]=4.91075e-01; p2[1]=6.27855e-01; p2[2] = 1.09965e+01; p2[3] =-4.24656e-01; p2[4] = 0.;
      ep2[0]=3.00310e-02; ep2[1]= 4.42650e-02; ep2[2]= 3.40324e-01; ep2[3]= 2.46448e-02; ep2[4]= 0. ;

      p3[0] = 6.75789e-01; p3[1] = 1.23129e+00; p3[2] = 3.46294e+00; p3[3]= -2.48336e-02; p3[4]= 0.;
      ep3[0]= 2.90342e-02; ep3[1]= 4.34524e-02; ep3[2]= 7.18743e-02; ep3[3]= 1.08694e-02; ep3[4]= 0.;

      fExpPionMean= p1[0]*TMath::Log(x)+p1[1]*TMath::Exp(-x*x*p1[2])+p1[3];
      fExpPionMeanRange=TMath::Sqrt(ep1[0]*ep1[0]+(x*x)*ep1[1]*ep1[1]+(x*x*x*x)*ep1[2]*ep1[2]);
      fExpKaonMean = p2[0]*TMath::Log(x)+p2[1]*TMath::Exp(-x*x*p2[2])+p2[3];
      fExpKaonMeanRange=TMath::Sqrt(ep2[0]*ep2[0]+(x*x)*ep2[1]*ep2[1]+(x*x*x*x)*ep2[2]*ep2[2]+(x*x*x*x*x*x)*ep2[3]*ep2[3]+(x*x*x*x*x*x*x*x)*ep2[4]*ep2[4]);
      fExpProtonMean=  p3[0]*TMath::Log(x)+p3[1]*TMath::Exp(-x*x*p3[2])+p3[3];
      fExpProtonMeanRange=TMath::Sqrt(ep3[0]*ep3[0]+(x*x)*ep3[1]*ep3[1]+(x*x*x*x)*ep3[2]*ep3[2]+(x*x*x*x*x*x)*ep3[3]*ep3[3]+(x*x*x*x*x*x*x*x)*ep3[4]*ep3[4]);
    }else{
      printf("ERROR: Wrong particle code %d\n",code);
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________
Int_t AliITSsadEdxFitter::InitializeMeanPosBB(Int_t code, Float_t pt){
  // computes peak positions from parametrized BB

  Double_t massp=AliPID::ParticleMass(AliPID::kProton);
  Double_t massk=AliPID::ParticleMass(AliPID::kKaon);
  Double_t masspi=AliPID::ParticleMass(AliPID::kPion);
  Bool_t isSA=kTRUE;
  Float_t sinthmed=0.8878;
  Float_t p=pt/sinthmed;
  Double_t bethep=fITSpid->Bethe(p,massp,isSA);
  Double_t bethek=fITSpid->Bethe(p,massk,isSA);
  Double_t bethepi=fITSpid->Bethe(p,masspi,isSA);
  Double_t betheref=bethepi;
  if(TMath::Abs(code)==321) betheref=bethek;
  else if(TMath::Abs(code)==2212) betheref=bethep;
  fExpPionMean=TMath::Log(bethepi)-TMath::Log(betheref);
  fExpKaonMean=TMath::Log(bethek)-TMath::Log(betheref);
  fExpProtonMean=TMath::Log(bethep)-TMath::Log(betheref);
  return 0;
}
//________________________________________________________
Bool_t AliITSsadEdxFitter::IsGoodBin(Int_t bin,Int_t code){
  //method to select the bins used for the analysis only
  Bool_t retvalue=kTRUE;
  TLine *l[2]; //for the cross 
  l[0]=new TLine(-2.1,0,2.,100.);
  l[1]=new TLine(-1.9,120,2.,0.);
  for(Int_t j=0;j<2;j++){
    l[j]->SetLineColor(4);
    l[j]->SetLineWidth(5);
  }

  if(code==211 && (bin<fBinsUsedPion[0] || bin>fBinsUsedPion[1])){
    for(Int_t j=0;j<2;j++) l[j]->Draw("same");	
    retvalue=kFALSE;
  }
  if(code==321 && (bin<fBinsUsedKaon[0] || bin>fBinsUsedKaon[1])){
    for(Int_t j=0;j<2;j++) l[j]->Draw("same");	
    retvalue=kFALSE;
  }
  if(code==2212 && (bin<fBinsUsedProton[0] || bin>fBinsUsedProton[1])){
    for(Int_t j=0;j<2;j++) l[j]->Draw("same");	
    retvalue=kFALSE;
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

//______________________________________________________________________
void AliITSsadEdxFitter::Initialize(TH1F* h, Int_t code,Int_t bin){
  //code to get the expected values to use for fitting

  Double_t xbins[23]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  Double_t pt=(xbins[bin+1]+xbins[bin])/2;

  //draw and label
  h->SetTitle(Form("p_{t}=[%1.2f,%1.2f], code=%d",xbins[bin],xbins[bin+1],code));
  h->GetXaxis()->SetTitle("[ln dE/dx]_{meas} - [ln dE/dx(i)]_{calc}");
  h->GetYaxis()->SetTitle("counts");
  h->Draw("e");
  h->SetFillColor(11);
	
  //expected values
  fExpPionAmpl = h->GetMaximum()/(h->GetRMS()*TMath::Sqrt(2*TMath::Pi()));
  if(fOptInit==0) InitializeMeanPosBB(code, pt);
  else if(fOptInit==1) InitializeMeanPosParam(code, pt);
  fExpPionSigma = CalcSigma(211,pt); //expected sigma values
  fExpKaonSigma = CalcSigma(321,pt);
  fExpProtonSigma = CalcSigma(2212,pt);

}


//________________________________________________________
void AliITSsadEdxFitter::DoFit(TH1F *h, Int_t bin, Int_t signedcode, TGraph *gres){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: all free, second step: pion gaussian fixed, third step: kaon gaussian fixed
  // final step: refit all using the parameters and tollerance limits (+-20%)

  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  TString modfit = "M0R+";
  Int_t code=TMath::Abs(signedcode);
  Initialize(h,code,bin);
  PrintInitialValues();
  if(!IsGoodBin(bin,code)) return;

  printf("___________________________________________________________________ First Step: pions\n");
  fstep1 = new TF1("step1",SingleGausStep,fRangeStep4[0],fRangeStep4[1],3);
  fstep1->SetParameter(0,fExpPionAmpl);       //initial ampl pion
  fstep1->SetParameter(1,fExpPionMean);       //initial mean pion
  fstep1->SetParameter(2,fExpPionSigma); //initial sigma pion
  fstep1->SetParLimits(0,0.,fExpPionAmpl*1.2);    //limits ampl pion
  fstep1->SetParLimits(1,fRangeStep4[0],fRangeStep4[1]);                                     //limits mean pion (dummy)
  fstep1->SetParLimits(2,fExpPionSigma*fLimitsOnSigmaPion[0],fExpPionSigma*fLimitsOnSigmaPion[1]); //limits sigma pion

  if(fExpPionSigma>0) h->Fit(fstep1,modfit.Data(),"",fExpPionMean+fRangeStep1[0],fExpPionMean+fRangeStep1[1]);//first fit
  else for(Int_t npar=0;npar<3;npar++) fstep1->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Second Step: kaons\n");
  fstep2 = new TF1("fstep2",DoubleGausStep,fRangeStep4[0],fRangeStep4[1],6);
  fstep2->FixParameter(0,fstep1->GetParameter(0)); //fixed ampl pion
  fstep2->FixParameter(1,fstep1->GetParameter(1)); //fixed mean pion
  fstep2->FixParameter(2,fstep1->GetParameter(2)); //fixed sigma pion
  fstep2->SetParameter(3,fstep1->GetParameter(0)/8.); //initial ampl kaon
  fstep2->SetParameter(4,fExpKaonMean);                //initial mean kaon
  fstep2->SetParameter(3,fExpKaonSigma);                 //initial sigma kaon
  fstep2->SetParLimits(3,0.,fstep1->GetParameter(0));                                         //limits ampl kaon
  fstep2->SetParLimits(4,fstep1->GetParameter(1),fRangeStep4[1]);                             //limits mean kaon 
  fstep2->SetParLimits(5,fExpKaonSigma*fLimitsOnSigmaKaon[0],fExpKaonSigma*fLimitsOnSigmaKaon[1]);  //limits sigma kaon

  if(fExpKaonSigma>0) h->Fit(fstep2,modfit.Data(),"",fExpKaonMean+fRangeStep2[0],fExpKaonMean+fRangeStep2[1]);//second fit
  else for(Int_t npar=3;npar<6;npar++) fstep2->FixParameter(npar,0.00001);


  printf("___________________________________________________________________ Third Step: protons\n");
  fstep3= new TF1("fstep3",FinalGausStep,fRangeStep4[0],fRangeStep4[1],9);
  fstep3->FixParameter(0,fstep1->GetParameter(0)); //fixed ampl pion
  fstep3->FixParameter(1,fstep1->GetParameter(1)); //fixed mean pion
  fstep3->FixParameter(2,fstep1->GetParameter(2)); //fixed sigma pion
  fstep3->FixParameter(3,fstep2->GetParameter(3)); //fixed ampl kaon
  fstep3->FixParameter(4,fstep2->GetParameter(4)); //fixed mean kaon
  fstep3->FixParameter(5,fstep2->GetParameter(5)); //fidex sigma kaon
  fstep3->SetParameter(6,fstep2->GetParameter(0)/16.); //initial ampl proton
  fstep3->SetParameter(7,fExpProtonMean);               //initial mean proton
  fstep3->SetParameter(8,fExpProtonSigma);                //initial sigma proton
  fstep3->SetParLimits(6,0.,fstep2->GetParameter(0));                                                //limits ampl proton
  fstep3->SetParLimits(7,fstep2->GetParameter(4),fRangeStep4[1]);                                    //limits mean proton
  fstep3->SetParLimits(8,fExpProtonSigma*fLimitsOnSigmaProton[0],fExpProtonSigma*fLimitsOnSigmaProton[1]); //limits sigma proton

  if(fExpProtonSigma>0) h->Fit(fstep3,modfit.Data(),"",fExpProtonMean+fRangeStep3[0],fExpProtonMean+fRangeStep3[1]);//third fit
  else for(Int_t npar=6;npar<9;npar++) fstep3->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Final Step: refit all\n");
  fstepTot = new TF1("funztot",FinalGausStep,fRangeStep4[0],fRangeStep4[1],9);
  fstepTot->SetLineColor(1);

  Double_t initialParametersStepTot[9];
  initialParametersStepTot[0] = fstep1->GetParameter(0);//first gaussian
  initialParametersStepTot[1] = fstep1->GetParameter(1);
  initialParametersStepTot[2] = fstep1->GetParameter(2);

  initialParametersStepTot[3] = fstep2->GetParameter(3);//second gaussian
  initialParametersStepTot[4] = fstep2->GetParameter(4);
  initialParametersStepTot[5] = fstep2->GetParameter(5);

  initialParametersStepTot[6] = fstep3->GetParameter(6);//third gaussian
  initialParametersStepTot[7] = fstep3->GetParameter(7);
  initialParametersStepTot[8] = fstep3->GetParameter(8);  

  fstepTot->SetParameters(initialParametersStepTot); //initial parameters
  fstepTot->SetParLimits(0,initialParametersStepTot[0]*0.9,initialParametersStepTot[0]*1.1); //tolerance limits ampl pion
  fstepTot->FixParameter(1,initialParametersStepTot[1]); //fixed mean pion
  fstepTot->FixParameter(2,initialParametersStepTot[2]); //fixed sigma pion
  fstepTot->SetParLimits(3,initialParametersStepTot[3]*0.9,initialParametersStepTot[3]*1.1); //tolerance limits ampl kaon
  fstepTot->FixParameter(4,initialParametersStepTot[4]); //fixed mean kaon
  fstepTot->FixParameter(5,initialParametersStepTot[5]); //fixed sigma kaon
  fstepTot->SetParLimits(6,initialParametersStepTot[6]*0.9,initialParametersStepTot[6]*1.1); //tolerance limits ampl proton
  fstepTot->FixParameter(7,initialParametersStepTot[7]); //fixed mean proton
  fstepTot->FixParameter(8,initialParametersStepTot[8]); //fixed sigma proton

  h->Fit(fstepTot,modfit.Data(),"",fRangeStep4[0],fRangeStep4[1]);//refit all

  for(Int_t j=0;j<9;j++) {
    fFitPars[j] = fstepTot->GetParameter(j);
    fFitParErrors[j] = fstepTot->GetParError(j);
  }

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
void AliITSsadEdxFitter::DoFitProton(TH1F *h, Int_t bin, Int_t signedcode, TGraph *gres){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: pion peak, second step: proton peak, third step: kaon peak
  // final step: refit all using the parameters
  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  TString modfit = "M0R+";
  Int_t code=TMath::Abs(signedcode);
  Initialize(h,code,bin);
  PrintInitialValues();
  if(!IsGoodBin(bin,code)) return;

  printf("___________________________________________________________________ First Step: pion\n");
  fstep1 = new TF1("step1",SingleGausStep,fRangeStep4[0],fRangeStep4[1],3);
  fstep1->SetParameter(0,fExpPionAmpl);       //initial ampl pion
  fstep1->SetParameter(1,fExpPionMean);       //initial mean pion
  fstep1->SetParameter(2,fExpPionSigma); //initial sigma pion
  fstep1->SetParLimits(0,0,fExpPionAmpl*1.2);                                                          //limits ampl pion
  fstep1->SetParLimits(1,fRangeStep4[0],fRangeStep4[1]);                                       //limits mean pion (dummy)
  fstep1->SetParLimits(2,fExpPionSigma*fLimitsOnSigmaPion[0],fExpPionSigma*fLimitsOnSigmaPion[1]);   //limits sigma pion

  if(fExpPionSigma>0)  h->Fit(fstep1,modfit,"",fExpPionMean+fRangeStep1[0],fExpPionMean+fRangeStep1[1]);//first fit
  else for(Int_t npar=0;npar<3;npar++) fstep1->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Second Step: proton\n");
  fstep2 = new TF1("step2",SingleGausStep,fRangeStep4[0],fRangeStep4[1],3);
  fstep2->SetParameter(0,fstep1->GetParameter(0)/16.);//initial ampl proton
  fstep2->SetParameter(1,fExpProtonMean);              //initial mean proton
  fstep2->SetParameter(2,fExpProtonSigma);               //initial sigma proton
  fstep2->SetParLimits(0,0.,fstep1->GetParameter(0));                                                //limits ampl proton
  fstep2->SetParLimits(1,fstep1->GetParameter(1),fRangeStep4[1]);                                    //limits mean proton
  fstep2->SetParLimits(2,fExpProtonSigma*fLimitsOnSigmaProton[0],fExpProtonSigma*fLimitsOnSigmaProton[1]); //limits sigma proton

  if(fExpProtonSigma>0) h->Fit(fstep2,modfit,"",fExpProtonMean+fRangeStep3[0],fExpProtonMean+fRangeStep3[1]);//second fit
  else for(Int_t npar=0;npar<3;npar++) fstep2->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Third Step: kaon\n");
  fstep3= new TF1("fstep3",FinalGausStep,fRangeStep4[0],fRangeStep4[1],9);
  fstep3->FixParameter(0,fstep1->GetParameter(0)); //fixed ampl pion
  fstep3->FixParameter(1,fstep1->GetParameter(1)); //fixed mean pion
  fstep3->FixParameter(2,fstep1->GetParameter(2)); //fixed sigma pion
  fstep3->FixParameter(6,fstep2->GetParameter(0)); //fixed ampl proton
  fstep3->FixParameter(7,fstep2->GetParameter(1)); //fixed mean proton
  fstep3->FixParameter(8,fstep2->GetParameter(2)); //fixed sigma proton
  fstep3->SetParameter(3,fstep1->GetParameter(0)/8.); //initial ampl kaon
  fstep3->SetParameter(4,fExpKaonMean);                //initial mean kaon
  fstep3->SetParameter(5,fExpKaonSigma);                 //initial sigma kaon
  fstep3->SetParLimits(3,fstep2->GetParameter(0),fstep1->GetParameter(0));                   //limits ampl kaon
  fstep3->SetParLimits(4,fstep1->GetParameter(1),fstep2->GetParameter(1));                   //limits mean kaon
  fstep3->SetParLimits(5,fExpKaonSigma*fLimitsOnSigmaKaon[0],fExpKaonSigma*fLimitsOnSigmaKaon[1]); //limits sigma kaon
  if(fExpKaonSigma>0) h->Fit(fstep3,modfit,"",fExpKaonMean+fRangeStep2[0],fExpKaonMean+fRangeStep2[1]);//third fit
  else for(Int_t npar=3;npar<6;npar++) fstep3->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Final Step: refit all\n");
  fstepTot = new TF1("funztot",FinalGausStep,fRangeStep4[0],fRangeStep4[1],9);
  fstepTot->SetLineColor(1);

  Double_t initialParametersStepTot[9];
  initialParametersStepTot[0] = fstep1->GetParameter(0);//first gaussian
  initialParametersStepTot[1] = fstep1->GetParameter(1);
  initialParametersStepTot[2] = fstep1->GetParameter(2);

  initialParametersStepTot[3] = fstep3->GetParameter(3);//second gaussian
  initialParametersStepTot[4] = fstep3->GetParameter(4);
  initialParametersStepTot[5] = fstep3->GetParameter(5);

  initialParametersStepTot[6] = fstep2->GetParameter(0);//third gaussian
  initialParametersStepTot[7] = fstep2->GetParameter(1);
  initialParametersStepTot[8] = fstep2->GetParameter(2);  

  fstepTot->SetParameters(initialParametersStepTot); //initial parameters
  fstepTot->SetParLimits(0,initialParametersStepTot[0]*0.9,initialParametersStepTot[0]*1.1); //tolerance limits ampl pion
  fstepTot->FixParameter(1,initialParametersStepTot[1]); //fixed mean pion
  fstepTot->FixParameter(2,initialParametersStepTot[2]); //fixed sigma pion
  fstepTot->SetParLimits(3,initialParametersStepTot[3]*0.9,initialParametersStepTot[3]*1.1); //tolerance limits ampl kaon
  fstepTot->FixParameter(4,initialParametersStepTot[4]); //fixed mean kaon
  fstepTot->FixParameter(5,initialParametersStepTot[5]); //fixed sigma kaon
  fstepTot->SetParLimits(6,initialParametersStepTot[6]*0.9,initialParametersStepTot[6]*1.1); //tolerance limits ampl proton
  fstepTot->FixParameter(7,initialParametersStepTot[7]); //fixed mean proton
  fstepTot->FixParameter(8,initialParametersStepTot[8]); //fixed sigma proton

  h->Fit(fstepTot,modfit,"",fRangeStep4[0],fRangeStep4[1]);//refit all

  for(Int_t j=0;j<9;j++) {
    fFitPars[j] = fstepTot->GetParameter(j);
    fFitParErrors[j] = fstepTot->GetParError(j);
  }
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
void AliITSsadEdxFitter::DoFitProtonFirst(TH1F *h, Int_t bin, Int_t signedcode, TGraph *gres){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: proton peak, second step: pion peak, third step: kaon peak
  // final step: refit all using the parameters
  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  TString modfit = "M0R+";
  Int_t code=TMath::Abs(signedcode);
  Initialize(h,code,bin);
  PrintInitialValues();
  if(!IsGoodBin(bin,code)) return;

  printf("___________________________________________________________________ First Step: proton\n");
  fstep1 = new TF1("step1",SingleGausStep,fRangeStep4[0],fRangeStep4[1],3);
  fstep1->SetParameter(0,fExpPionAmpl/16.);       //initial ampl proton`
  fstep1->SetParameter(1,fExpProtonMean);  //initial mean proton
  fstep1->SetParameter(2,fExpProtonSigma);   //initial sigma proton
  fstep1->SetParLimits(0,0,fExpPionAmpl);                                                                    //limits ampl proton
  fstep1->SetParLimits(1,fExpPionMean,fRangeStep4[1]);                                                       //limits mean proton (dummy)
  fstep1->SetParLimits(2,fExpProtonSigma*fLimitsOnSigmaProton[0],fExpProtonSigma*fLimitsOnSigmaProton[1]); //limits sigma proton

  if(fExpProtonSigma>0)  h->Fit(fstep1,modfit,"",fExpProtonMean+fRangeStep3[0],fExpProtonMean+fRangeStep3[1]);//first fit
  else for(Int_t npar=0;npar<3;npar++) fstep1->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Second Step: pion\n");
  fstep2 = new TF1("step2",DoubleGausStep,fRangeStep4[0],fRangeStep4[1],6);
  fstep2->FixParameter(0,fstep1->GetParameter(0)); //fixed ampl proton 
  fstep2->FixParameter(1,fstep1->GetParameter(1)); //fixed mean proton
  fstep2->FixParameter(2,fstep1->GetParameter(2)); //fixed sigma proton
  fstep2->SetParameter(3,fExpPionAmpl);             //initial ampl pion
  fstep2->SetParameter(4,fExpPionMean);             //initial mean pion
  fstep2->SetParameter(5,fExpPionSigma);       //initial sigma pion
  fstep2->SetParLimits(3,0.,1.2*fExpPionAmpl);                                                               //limits ampl pion
  fstep2->SetParLimits(4,fRangeStep4[0],fstep1->GetParameter(1));                                //limits mean pion
  fstep2->SetParLimits(5,fExpPionSigma*fLimitsOnSigmaPion[0],fExpPionSigma*fLimitsOnSigmaPion[1]);     //limits sigma pion

  if(fExpPionSigma>0) h->Fit(fstep2,modfit,"",fExpPionMean+fRangeStep1[0],fExpPionMean+fRangeStep1[1]);//second fit
  else for(Int_t npar=0;npar<3;npar++) fstep2->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Third Step: kaon\n");
  fstep3= new TF1("fstep3",FinalGausStep,fRangeStep4[0],fRangeStep4[1],9);
  fstep3->FixParameter(0,fstep1->GetParameter(0)); //fixed ampl proton
  fstep3->FixParameter(1,fstep1->GetParameter(1)); //fixed mean proton
  fstep3->FixParameter(2,fstep1->GetParameter(2)); //fixed sigma proton
  fstep3->FixParameter(3,fstep2->GetParameter(3)); //fixed ampl pion
  fstep3->FixParameter(4,fstep2->GetParameter(4)); //fixed mean pion
  fstep3->FixParameter(5,fstep2->GetParameter(5)); //fixed sigma pion
  fstep3->SetParameter(6,fstep2->GetParameter(0)/8.); //initial ampl kaon
  fstep3->SetParameter(7,fExpKaonMean);                //initial mean kaon
  fstep3->SetParameter(8,fExpKaonSigma);                 //initial sigma kaon
  fstep3->SetParLimits(6,fstep1->GetParameter(0),fstep2->GetParameter(3));                   //limits ampl kaon
  fstep3->SetParLimits(7,fstep2->GetParameter(4),fstep1->GetParameter(1));                   //limits mean kaon
  fstep3->SetParLimits(8,fExpKaonSigma*fLimitsOnSigmaKaon[0],fExpKaonSigma*fLimitsOnSigmaKaon[1]); //limits sigma kaon
  if(fExpKaonSigma>0) h->Fit(fstep3,modfit,"",fExpKaonMean+fRangeStep2[0],fExpKaonMean+fRangeStep2[1]);//third fit
  else for(Int_t npar=3;npar<6;npar++) fstep3->FixParameter(npar,0.00001);

  printf("___________________________________________________________________ Final Step: refit all\n");
  fstepTot = new TF1("funztot",FinalGausStep,fRangeStep4[0],fRangeStep4[1],9);
  fstepTot->SetLineColor(1);

  Double_t initialParametersStepTot[9];
  initialParametersStepTot[0] = fstep1->GetParameter(0);//first gaussian
  initialParametersStepTot[1] = fstep1->GetParameter(1);
  initialParametersStepTot[2] = fstep1->GetParameter(2);

  initialParametersStepTot[3] = fstep3->GetParameter(6);//second gaussian
  initialParametersStepTot[4] = fstep3->GetParameter(7);
  initialParametersStepTot[5] = fstep3->GetParameter(8);

  initialParametersStepTot[6] = fstep2->GetParameter(3);//third gaussian
  initialParametersStepTot[7] = fstep2->GetParameter(4);
  initialParametersStepTot[8] = fstep2->GetParameter(5);  

  fstepTot->SetParameters(initialParametersStepTot); //initial parameters
  fstepTot->SetParLimits(0,initialParametersStepTot[0]*0.9,initialParametersStepTot[0]*1.1); //tolerance limits ampl proton
  fstepTot->FixParameter(1,initialParametersStepTot[1]); //fixed mean pion
  fstepTot->FixParameter(2,initialParametersStepTot[2]); //fixed sigma pion
  fstepTot->SetParLimits(3,initialParametersStepTot[3]*0.9,initialParametersStepTot[3]*1.1); //tolerance limits ampl kaon
  fstepTot->FixParameter(4,initialParametersStepTot[4]); //fixed mean kaon
  fstepTot->FixParameter(5,initialParametersStepTot[5]); //fixed sigma kaon
  fstepTot->SetParLimits(6,initialParametersStepTot[6]*0.9,initialParametersStepTot[6]*1.1); //tolerance limits ampl pion
  fstepTot->FixParameter(7,initialParametersStepTot[7]); //fixed mean proton
  fstepTot->FixParameter(8,initialParametersStepTot[8]); //fixed sigma proton

  h->Fit(fstepTot,modfit,"",fRangeStep4[0],fRangeStep4[1]);//refit all

  for(Int_t j=0;j<9;j++) {
    fFitPars[j] = fstepTot->GetParameter(j);
    fFitParErrors[j] = fstepTot->GetParError(j);
  }
  //************************************* storing parameter to calculate the yields *******
  Int_t chpa=0;
  if(code==321) chpa=3;
  if(code==211) chpa=6;
  for(Int_t j=0;j<3;j++) {
    fFitpar[j] = fstepTot->GetParameter(j+chpa);
    fFitparErr[j] = fstepTot->GetParError(j+chpa);
  }

  DrawFitFunction(fstepTot);
  CalcResidual(h,fstepTot,gres);
  return;
}


//________________________________________________________
void AliITSsadEdxFitter::DoFitOnePeak(TH1F *h, Int_t bin, Int_t signedcode){
  // single-gaussian fit to log(dedx)-log(dedxBB) histogram
  TF1 *fstep1;
  TString modfit = "M0R+";
  Int_t code=TMath::Abs(signedcode);
  Initialize(h,code,bin);
  PrintInitialValues();
  if(!IsGoodBin(bin,code)) return;

  printf("___________________________________________________________________ Single Step\n");
  fstep1 = new TF1("step2",SingleGausStep,fRangeStep4[0],fRangeStep4[1],3);
  fstep1->SetParameter(0,fExpPionAmpl/16.);                   //initial ampl 
  fstep1->SetParameter(1,fExpProtonMean);              //initial mean 
  fstep1->SetParameter(2,fExpProtonSigma);               //initial sigma 
  fstep1->SetParLimits(0,0.,fExpPionAmpl);                                                                   //limits ampl proton
  fstep1->SetParLimits(1,fExpPionMean,fRangeStep4[1]);                                                       //limits mean proton
  //fstep1->SetParLimits(2,fExpProtonSigma*fLimitsOnSigmaProton[0],fExpProtonSigma*fLimitsOnSigmaProton[1]); //limits sigma proton

  if(fExpProtonSigma>0) h->Fit(fstep1,modfit,"",fExpProtonMean+fRangeStep3[0],fExpProtonMean+fRangeStep3[1]);//fit
  else for(Int_t npar=0;npar<3;npar++) fstep1->FixParameter(npar,0.00001);

  fstep1->SetLineColor(1);
  fstep1->Draw("same");

  fFitpar[0] = fstep1->GetParameter(0);
  fFitparErr[0] = fstep1->GetParError(0);
  return;
}

//________________________________________________________
void AliITSsadEdxFitter::DoFitTail(TH1F *h, Int_t bin, Int_t signedcode){
  // 3-gaussian fit to log(dedx)-log(dedxBB) histogram
  // pt bin from 0 to 20, code={211,321,2212} 
  // first step: all free, second step: pion gaussian fixed, third step: kaon gaussian fixed
  // final step: refit all using the parameters and tollerance limits (+-20%)
  // WARNING: exponential tail added in the right of the Gaussian shape
  Int_t code=TMath::Abs(signedcode);
  if(!IsGoodBin(bin,code)) return;

  TF1 *fstep1, *fstep2, *fstep3, *fstepTot;
  TString modfit = "M0R+";
  Initialize(h,code,bin);
  PrintInitialValues();

  printf("\n___________________________________________________________________\n First Step: pions\n\n");
  fstep1 = new TF1("step1",SingleGausTail,-3.5,3.5,5);
  fstep1->SetParameter(0,fExpPionAmpl);//initial 
  fstep1->SetParameter(1,fExpPionMean);
  fstep1->SetParameter(3,1.2);
  fstep1->SetParameter(4,10.);

  fstep1->SetParLimits(0,0,fExpPionAmpl*1.2);
  fstep1->SetParLimits(1,-3.5,3.5);
  fstep1->SetParLimits(2,0.1,0.25);
  fstep1->SetParLimits(4,5.,20.);
  if(bin<8) fstep1->SetParLimits(4,13.,25.);

  h->Fit(fstep1,modfit,"",fExpPionMean-0.45,fExpPionMean+0.45);//first fit

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

  printf("\n_____________________________________________________________________\n Final Step: refit all \n\n"); 
  fstepTot = new TF1("funztot",FinalGausTail,-3.5,3.5,15);
  fstepTot->SetLineColor(1);

  Double_t initialParametersStepTot[15];
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
  h->Fit(fstepTot,modfit,"",-3.5,3.5); //refit all

  for(Int_t j=0;j<9;j++) {
    fFitPars[j] = fstepTot->GetParameter(j);
    fFitParErrors[j] = fstepTot->GetParError(j);
  }

  //************************************* storing parameter to calculate the yields *******
  Int_t chpa=0;
  if(code==321) chpa=5;
  if(code==2212) chpa=10;
  for(Int_t j=0;j<3;j++) {
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
  for(Int_t i=0;i<3;i++) {
    fitpar[i] = fFitpar[i];
    fitparerr[i] = fFitparErr[i];
  }
  return;
}


//________________________________________________________
void AliITSsadEdxFitter::PrintAll() const{
  // print parameters

  printf("Range 1 = %f %f\n",fRangeStep1[0],fRangeStep1[1]);
  printf("Range 2 = %f %f\n",fRangeStep2[0],fRangeStep2[1]);
  printf("Range 3 = %f %f\n",fRangeStep3[0],fRangeStep3[1]);
  printf("Range F = %f %f\n",fRangeStep4[0],fRangeStep4[1]);
  printf(" Sigma1 = %f %f\n",fLimitsOnSigmaPion[0],fLimitsOnSigmaPion[1]);
  printf(" Sigma2 = %f %f\n",fLimitsOnSigmaKaon[0],fLimitsOnSigmaKaon[1]);
  printf(" Sigma3 = %f %f\n",fLimitsOnSigmaProton[0],fLimitsOnSigmaProton[1]);
}
//______________________________________________________________________
void AliITSsadEdxFitter::PrintInitialValues() const{
  // print initial values  

  printf("mean values        -> %f %f %f\n",fExpPionMean,fExpKaonMean,fExpProtonMean);
  printf("integration ranges -> (%1.2f,%1.2f) (%1.2f,%1.2f) (%1.2f,%1.2f)\n",
	 fRangeStep1[0],fRangeStep1[1],fRangeStep2[0],fRangeStep2[1],fRangeStep3[0],fRangeStep3[1]);
  printf("sigma values -> %f %f %f\n",fExpPionSigma,fExpKaonSigma,fExpProtonSigma);
  printf("sigma ranges -> (%1.2f,%1.2f) (%1.2f,%1.2f) (%1.2f,%1.2f)\n",
	 fLimitsOnSigmaPion[0],fLimitsOnSigmaPion[1],
	 fLimitsOnSigmaKaon[0],fLimitsOnSigmaKaon[1],
	 fLimitsOnSigmaProton[0],fLimitsOnSigmaProton[1]);
}
