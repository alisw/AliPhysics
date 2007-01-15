/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////
//Class for PID in the ITS                             //
//                                                     //
//                                                     //
/////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TF1.h>
#include <TMath.h>

#include "AliITSPidParItem.h"

ClassImp(AliITSPidParItem)
//____________________________________________________________________
AliITSPidParItem::AliITSPidParItem():
fPCenter(0),
fPWidth(0)
{
  // default constructor
  for(Int_t i=0;i<39;i++){
  fBuff[i]=0;
  }
}//____________________________________________________________________
AliITSPidParItem::AliITSPidParItem(Float_t center,Float_t width,Double_t *buff):
fPCenter(center),
fPWidth(width)
{
  // standard constructor
  for (Int_t i=0;i<39;i++) fBuff[i]=buff[i];

}

//____________________________________________________________________
void AliITSPidParItem::GetParameters(Double_t *buff) const{
  //get all the parameters
  for (Int_t i=0;i<39;i++) buff[i]=fBuff[i];

}
//____________________________________________________________________
void AliITSPidParItem::GetProtonPar(Double_t *buffp) const{
  //get the protons' parameters (Width Landau, Most Probable, Area, Width Gaussian, Chi2 fit, NDF fit, Integral fit)
  buffp[0]=fBuff[0];
  buffp[1]=fBuff[1];
  buffp[2]=fBuff[2];
  buffp[3]=fBuff[3];
  buffp[4]=fBuff[8];
  buffp[5]=fBuff[9];
  buffp[6]=fBuff[10];
}
//____________________________________________________________________
void AliITSPidParItem::GetKaonPar(Double_t *buffk) const{
  //get the kaons' parameters (Width Landau, Most Probable, Area, Width Gaussian, Chi2 fit, NDF fit, Integral fit)
  buffk[0]=fBuff[13];
  buffk[1]=fBuff[14];
  buffk[2]=fBuff[15];
  buffk[3]=fBuff[16];
  buffk[4]=fBuff[21];
  buffk[5]=fBuff[22];
  buffk[6]=fBuff[23];
}
//____________________________________________________________________
void AliITSPidParItem::GetPionPar(Double_t *buffpi) const{
  //get the pions' parameters (Width Landau, Most Probable, Area, Width Gaussian, Chi2 fit, NDF fit, Integral fit)
  buffpi[0]=fBuff[26];
  buffpi[1]=fBuff[27];
  buffpi[2]=fBuff[28];
  buffpi[3]=fBuff[29];
  buffpi[4]=fBuff[34];
  buffpi[5]=fBuff[35];
  buffpi[6]=fBuff[36];
}
//____________________________________________________________________
void AliITSPidParItem::GetPar0(Double_t *buff0) const{
  //Width Landau for protons, kaons, pions.
  buff0[0]=fBuff[0];
  buff0[1]=fBuff[13];
  buff0[2]=fBuff[26];
}
//____________________________________________________________________
void AliITSPidParItem::GetPar1(Double_t *buff1) const{
  //Most Probable for protons, kaons, pions.
  buff1[0]=fBuff[1];
  buff1[1]=fBuff[14];
  buff1[2]=fBuff[27];
}
//____________________________________________________________________
void AliITSPidParItem::GetPar2(Double_t *buff2) const{
  //Area for protons, kaons, pions.
  buff2[0]=fBuff[2];
  buff2[1]=fBuff[15];
  buff2[2]=fBuff[28];
}
//____________________________________________________________________
void AliITSPidParItem::GetPar3(Double_t *buff3) const{
  //Width Gaussian for protons, kaons, pions.
  buff3[0]=fBuff[3];
  buff3[1]=fBuff[16];
  buff3[2]=fBuff[29];
}
//____________________________________________________________________
void AliITSPidParItem::GetChisquare(Double_t *buffchi) const{
  //Chi2 of the fit for protons, kaons, pions.
  buffchi[0]=fBuff[8];
  buffchi[1]=fBuff[21];
  buffchi[2]=fBuff[34];
}
//____________________________________________________________________
void AliITSPidParItem::GetNDF(Double_t *buffndf) const{
  //NDF of the fit for protons, kaons, pions.
  buffndf[0]=fBuff[9];
  buffndf[1]=fBuff[22];
  buffndf[2]=fBuff[35];
}
//____________________________________________________________________
void AliITSPidParItem::GetProParFun(Double_t *pfun) const{
  //some Protons parameters: Width Landau, Most Probable, Area, Width Gaussian, Integral fit
  pfun[0]=fBuff[0];
  pfun[1]=fBuff[1];
  pfun[2]=fBuff[2];
  pfun[3]=fBuff[3];
  pfun[4]=fBuff[10];
}
//____________________________________________________________________
void AliITSPidParItem::GetKaoParFun(Double_t *kfun) const{
  //some Kaons parameters: Width Landau, Most Probable, Area, Width Gaussian, Integral fit
  kfun[0]=fBuff[13];
  kfun[1]=fBuff[14];
  kfun[2]=fBuff[15];
  kfun[3]=fBuff[16];
  kfun[4]=fBuff[23];
}
//____________________________________________________________________
void AliITSPidParItem::GetPiParFun(Double_t *pifun) const{
  //some Pions parameters: Width Landau, Most Probable, Area, Width Gaussian, Integral fit
  pifun[0]=fBuff[26];
  pifun[1]=fBuff[27];
  pifun[2]=fBuff[28];
  pifun[3]=fBuff[29];
  pifun[4]=fBuff[36];
}
//____________________________________________________________________
void AliITSPidParItem::GetRangeLim(Double_t *range) const{
  //Range limits for the response functions
  range[0]=fBuff[11];//proton low
  range[1]=fBuff[12];//proton high
  range[2]=fBuff[24];//kaon low
  range[3]=fBuff[25];//kaon high
  range[4]=fBuff[37];//pion low
  range[5]=fBuff[38];//pion high
}
//____________________________________________________________________
void AliITSPidParItem::GetProtonParErr(Double_t *bufferp)const{
  //errors on the protons' parameters
  bufferp[0]=fBuff[4];
  bufferp[1]=fBuff[5];
  bufferp[2]=fBuff[6];
  bufferp[3]=fBuff[7];
}
//____________________________________________________________________
void AliITSPidParItem::GetKaonParErr(Double_t *bufferk)const{
  //errors on the kaons' parameters
  bufferk[0]=fBuff[17];
  bufferk[1]=fBuff[18];
  bufferk[2]=fBuff[19];
  bufferk[3]=fBuff[20];
}
//____________________________________________________________________
void AliITSPidParItem::GetPionParErr(Double_t *bufferpi)const{
  //errors on the pions' parameters
  bufferpi[0]=fBuff[30];
  bufferpi[1]=fBuff[31];
  bufferpi[2]=fBuff[32];
  bufferpi[3]=fBuff[33];
}
//____________________________________________________________________
void AliITSPidParItem::PrintParameters() const {
  // Prints the data members of this class
cout<<"==============================***************"<<endl;
cout<<"Momentum (GeV/c) - center of the bin - "<<fPCenter<<endl;
cout<<" Width of momentum bin (GeV/c) "<<fPWidth<<endl;
for (Int_t i=0;i<39;i++) cout<<"Parameter"<<i<<" = "<<fBuff[i]<<endl;

}

//_______________________________________________________________________
TF1* AliITSPidParItem::CookFunIts(TString namefun,Double_t *par,Double_t rangei,Double_t rangef,TString comment){
  //
  TF1 *fun;
    if (par[4]!=0) {
      fun=new TF1(comment,Langaufun2,rangei,rangef,5);
      fun->SetParameters(par);
     
    }
    else {fun=new TF1(namefun,"0");}
    return fun;
}

//_________________________________________________________________________
Double_t AliITSPidParItem::Langaufun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}
//_______________________________________________________________________
Double_t AliITSPidParItem::Langaufun2(Double_t *x, Double_t *par){
  //
  return 1/par[4]*Langaufun(x,par);
}
