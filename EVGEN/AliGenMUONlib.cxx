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

// Library class for particle pt and y distributions used for 
// muon spectrometer simulations.
// To be used with AliGenParam.
// The following particle typed can be simulated:
// pi, K, phi, omega, eta, J/Psi, Upsilon, charm and beauty mesons. 
//
// andreas.morsch@cern.ch
//

#include "TMath.h"
#include "TRandom.h"
#include "TDatabasePDG.h"

#include "AliGenMUONlib.h"

ClassImp(AliGenMUONlib)
//
//  Pions
Double_t AliGenMUONlib::PtPion(const Double_t *px, const Double_t* /*dummy*/)
{
//
//     PT-PARAMETERIZATION CDF, PRL 61(88) 1819
//     POWER LAW FOR PT > 500 MEV
//     MT SCALING BELOW (T=160 MEV)
//
  const Double_t kp0 = 1.3;
  const Double_t kxn = 8.28;
  const Double_t kxlim=0.5;
  const Double_t kt=0.160;
  const Double_t kxmpi=0.139;
  const Double_t kb=1.;
  Double_t y, y1, xmpi2, ynorm, a;
  Double_t x=*px;
  //
  y1=TMath::Power(kp0/(kp0+kxlim),kxn);
  xmpi2=kxmpi*kxmpi;
  ynorm=kb*(TMath::Exp(-sqrt(kxlim*kxlim+xmpi2)/kt));
  a=ynorm/y1;
  if (x > kxlim)
    y=a*TMath::Power(kp0/(kp0+x),kxn);
  else
    y=kb*TMath::Exp(-sqrt(x*x+xmpi2)/kt);
  return y*x;
}
//
// y-distribution
//
Double_t AliGenMUONlib::YPion( const Double_t *py, const Double_t */*dummy*/)
{
// Pion y
  Double_t y=TMath::Abs(*py);
/*
  const Double_t ka    = 7000.;
  const Double_t kdy   = 4.;
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
*/
  return 1.16526e+04+y*-3.79886e+03+y*y*4.31130e+02;
  
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPion(TRandom *ran)
{
// Pion composition 
    if (ran->Rndm() < 0.5) {
	return  211;
    } else {
	return -211;
    }
}

//____________________________________________________________
//
// Mt-scaling

Double_t AliGenMUONlib::PtScal(Double_t pt, Int_t np)
{
  //    SCALING EN MASSE PAR RAPPORT A PTPI
  //    MASS PI,K,ETA,RHO,OMEGA,ETA',PHI
  const Double_t khm[10] = {.13957,.493,.5488,.769,.7826,.958,1.02,0,0,0};
  //     VALUE MESON/PI AT 5 GEV
  const Double_t kfmax[10]={1.,0.3,0.55,1.0,1.0,1.0,1.0,0,0,0};
  np--;
  Double_t f5=TMath::Power(((sqrt(100.018215)+2.)/(sqrt(100.+khm[np]*khm[np])+2.0)),12.3);
  Double_t fmax2=f5/kfmax[np];
  // PIONS
  Double_t ptpion=100.*PtPion(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
				 (sqrt(pt*pt+khm[np]*khm[np])+2.0)),12.3)/ fmax2;
  return fmtscal*ptpion;
}
//
// kaon
//
//                pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtKaon( const Double_t *px, const Double_t */*dummy*/)
{
// Kaon pT
  return PtScal(*px,2);
}

// y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YKaon( const Double_t *py, const Double_t */*dummy*/)
{
// Kaon y
  Double_t y=TMath::Abs(*py);
/*
  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
*/

  return 1.16526e+04+y*-3.79886e+03+y*y*4.31130e+02;
}

//                 particle composition
//
Int_t AliGenMUONlib::IpKaon(TRandom *ran)
{
// Kaon composition
    if (ran->Rndm() < 0.5) {
	return  321;
    } else {
	return -321;
    }
}

//                    J/Psi 
//
//
//                pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtJpsiPPdummy(Double_t x, Double_t energy)
{
// J/Psi pT
// pp
// from the fit of RHIC, CDF & LHC data, see arXiv:1103.2394
//
  const Double_t kpt0 = 1.04*TMath::Power(energy,0.101);
  const Double_t kxn  = 3.9;
  //
  Double_t pass1 = 1.+0.363*(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiPP7000(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
// pp 7 TeV
//
  return PtJpsiPPdummy(*px,7000);
}

Double_t AliGenMUONlib::PtJpsiPP2760(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
// pp 2.76 TeV
//
  return PtJpsiPPdummy(*px,2760);
}

Double_t AliGenMUONlib::PtJpsiPP8800(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
// pp 8.8 TeV
//
  return PtJpsiPPdummy(*px,8800);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760ShFdummy(Double_t x, Int_t n)
{
// J/Psi shadowing factor vs pT for PbPb min. bias and 11 centr. bins (in 2.5<y<4)
//
// PbPb 2.76 TeV, for EKS98, minimum bias shadowing factor = 0.66 in 4pi
// S.Grigoryan, details presented at the PWG3-Muon meeting (05.10.2011)
// https://indico.cern.ch/conferenceDisplay.py?confId=157367
//
  const Double_t f1[12] = {1, 1.128, 1.097, 1.037, 0.937, 0.821, 0.693, 0.558,
			   0.428, 0.317, 0.231, 0.156};
  const Double_t f2[12] = {1, 1.313, 1.202, 1.039, 0.814, 0.593, 0.391, 0.224,
			   0.106, 0.041, 0.013, 0.002};
  const Double_t c1[7] = {1.6077e+00, 7.6300e-02,-7.1880e-03, 3.4067e-04,
			  -9.2776e-06,1.5138e-07, 1.4652e-09}; 
  const Double_t c2[7] = {6.2047e-01, 5.7653e-02,-4.1414e-03, 1.0301e-04, 
			  9.6205e-07,-7.4098e-08, 5.0946e-09}; 
  Double_t y1, y2;
  Int_t j;
  y1 = c1[j = 6]; y2 = c2[6];
  while (j > 0) {y1 = y1 * x + c1[--j]; y2 = y2 * x + c2[j];}
  
  y1 /= 1.+c1[6]*TMath::Power(x,6);
  y2 /= 1.+c2[6]*TMath::Power(x,6);
  //  
  y1 = 1 + (y1-2)*f1[n] + (y2+1-y1)*f2[n];
  if(y1<0) y1=0;
  return y1;
}

Double_t AliGenMUONlib::PtJpsiPbPb2760(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, minimum bias 0-100 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 0) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c1(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 1st centrality bin 0-5 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 1) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c2(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 2nd centrality bin 5-10 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 2) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c3(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 3rd centrality bin 10-20 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 3) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c4(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 4th centrality bin 20-30 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 4) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c5(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 5th centrality bin 30-40 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 5) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c6(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 6th centrality bin 40-50 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 6) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c7(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 7th centrality bin 50-60 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 7) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c8(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 8th centrality bin 60-70 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 8) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c9(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 9th centrality bin 70-80 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 9) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c10(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 10th centrality bin 80-90 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 10) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbPb2760c11(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// PbPb 2.76 TeV, 11th centrality bin 90-100 %
//
  return PtJpsiPbPb2760ShFdummy(*px, 11) * PtJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPPb8800ShFdummy(Double_t x, Int_t n)
{
// J/Psi shadowing factor vs pT for pPb min. bias and 4 centr. bins (in 2.5<y<4)
//
// pPb 8.8 TeV, for EKS98, minimum bias shadowing factor = 0.78 in 4pi
//
  const Double_t f[5] = {1, 1.33, 1.05, 0.67, 0.23};
  const Double_t c[7] = {6.4922e-01, 6.4857e-03, 4.7268e-03,-9.5075e-04, 
			 8.4075e-05,-4.2006e-06, 4.9970e-07};
  Double_t y;
  Int_t j;
  y = c[j = 6];
  while (j > 0) y  = y * x + c[--j];
  y /= 1 + c[6]*TMath::Power(x,6);
  //  
  return 1 + (y-1)*f[n];
}

Double_t AliGenMUONlib::PtJpsiPPb8800(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// pPb 8.8 TeV, minimum bias 0-100 %
//
  return PtJpsiPPb8800ShFdummy(*px, 0) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPPb8800c1(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// pPb 8.8 TeV, 1st centrality bin 0-20 %
//
  return PtJpsiPPb8800ShFdummy(*px, 1) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPPb8800c2(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// pPb 8.8 TeV, 2nd centrality bin 20-40 %
//
  return PtJpsiPPb8800ShFdummy(*px, 2) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPPb8800c3(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// pPb 8.8 TeV, 3rd centrality bin 40-60 %
//
  return PtJpsiPPb8800ShFdummy(*px, 3) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPPb8800c4(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// pPb 8.8 TeV, 4th centrality bin 60-100 %
//
  return PtJpsiPPb8800ShFdummy(*px, 4) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbP8800ShFdummy(Double_t x, Int_t n)
{
// J/Psi shadowing factor vs pT for Pbp min. bias and 4 centr. bins (in 2.5<y<4)
//
// Pbp 8.8 TeV, for EKS98, minimum bias shadowing factor = 0.78 in 4pi
//
  const Double_t f[5] = {1, 1.33, 1.05, 0.67, 0.23};
  const Double_t c[7] = {8.7562e-01, 2.1944e-02, 7.8509e-03,-1.3979e-03, 
			 3.8513e-05, 4.2008e-06, 1.7088e-06};
  Double_t y;
  Int_t j;
  y = c[j = 6];
  while (j > 0) y  = y * x + c[--j];
  y /= 1 + c[6]*TMath::Power(x,6);
  //  
  return 1 + (y-1)*f[n];
}

Double_t AliGenMUONlib::PtJpsiPbP8800(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// Pbp 8.8 TeV, minimum bias 0-100 %
//
  return PtJpsiPbP8800ShFdummy(*px, 0) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbP8800c1(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// Pbp 8.8 TeV, 1st centrality bin 0-20 %
//
  return PtJpsiPbP8800ShFdummy(*px, 1) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbP8800c2(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// Pbp 8.8 TeV, 2nd centrality bin 20-40 %
//
  return PtJpsiPbP8800ShFdummy(*px, 2) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbP8800c3(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// Pbp 8.8 TeV, 3rd centrality bin 40-60 %
//
  return PtJpsiPbP8800ShFdummy(*px, 3) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsiPbP8800c4(const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
// Pbp 8.8 TeV, 4th centrality bin 60-100 %
//
  return PtJpsiPbP8800ShFdummy(*px, 4) * PtJpsiPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtJpsi( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaled( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// PbPb 5.5 TeV
// scaled from CDF data at 2 TeV
// see S.Grigoryan, PWG3 Meeting, 27th Oct 2008

  const Double_t kpt0 = 5.100;
  const Double_t kxn  = 4.102;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 14 TeV
// scaled from CDF data at 2 TeV

  const Double_t kpt0 = 5.630;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP10( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 10 TeV
// scaled from CDF data at 2 TeV

  const Double_t kpt0 = 5.334;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP9( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 8.8 TeV
// scaled from CDF data at 2 TeV
//
  const Double_t kpt0 = 5.245;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP7( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 7 TeV
// scaled from CDF data at 2 TeV

  const Double_t kpt0 = 5.072;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP4( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 3.94 TeV
// scaled from CDF data at 2 TeV
//
  const Double_t kpt0 = 4.647;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP3( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 2.76 TeV
// scaled from CDF data at 1.9 TeV
//
  const Double_t kpt0 = 4.435;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPP2( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
//
// pp 1.96 TeV
// fit of the CDF data at 1.96 TeV
//
  const Double_t kpt0 = 4.233;
  const Double_t kxn  = 4.071;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPPb9( const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
//
// pPb 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.79
//
  Double_t c[5] = {6.42774e-01, 1.86168e-02, -6.77296e-04, 8.93512e-06, 1.31586e-07};
  Double_t x=*px;
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  //  
  Double_t d = 1.+c[4]*TMath::Power(x,4);
  return y/d * AliGenMUONlib::PtJpsiCDFscaledPP9(px,dummy);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPbP9( const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
//
// Pbp 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.79
//
  Double_t c[5] = {8.58557e-01, 5.39791e-02, -4.75180e-03, 2.49463e-04, 5.52396e-05};
  Double_t x=*px;
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  //  
  Double_t d = 1.+c[4]*TMath::Power(x,4);
  return y/d * AliGenMUONlib::PtJpsiCDFscaledPP9(px,dummy);
}

Double_t AliGenMUONlib::PtJpsiCDFscaledPbPb4( const Double_t *px, const Double_t *dummy)
{
// J/Psi pT
//
// PbPb 3.94 TeV, for EKS98 with minimum bias shadowing factor 0.66
//
  Double_t c[5] = {6.01022e-01, 4.70988e-02, -2.27917e-03, 3.09885e-05, 1.31955e-06};
  Double_t x=*px;
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  //  
  Double_t d = 1.+c[4]*TMath::Power(x,4);
  return y/d * AliGenMUONlib::PtJpsiCDFscaledPP4(px,dummy);
}

Double_t AliGenMUONlib::PtJpsiFlat( const Double_t */*px*/, const Double_t */*dummy*/ )
{
  return 1.;
}

Double_t AliGenMUONlib::PtJpsiPbPb( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT spectrum
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Float_t x = px[0];
    Float_t c[8] = {
	-2.13098e+00, 9.46552e+00, -5.06799e+00, 1.27260e+00, 
	-1.83806e-01, 1.55853e-02, -7.23241e-04, 1.42105e-05
    };
    
    Double_t y;
    if (x < 10.) {
	Int_t j;
	y = c[j = 7];
	while (j > 0) y  = y * x +c[--j];
	y = x * TMath::Exp(y);
    } else {
	y = 0.;
    }
    return y;
}

Double_t AliGenMUONlib::PtJpsiBPbPb( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT spectrum
// B -> J/Psi X
    Double_t x0 =   4.0384;
    Double_t  n =   3.0288;
    
    Double_t x = px[0];
    Double_t y = x / TMath::Power((1. + (x/x0)*(x/x0)), n);
    
    return y;
}


Double_t AliGenMUONlib::PtJpsiPP( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT spectrum
//
// R. Vogt 2002
// pp 14 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Float_t x = px[0];
    Float_t c[4] = {8.47471e+00, -1.93567e+00, 1.50271e-01, -5.51212e-03};
 
    Double_t y;
    if (x < 10.) {
	Int_t j;
	y = c[j = 3];
	while (j > 0) y  = y * x +c[--j];
	y = x * TMath::Exp(y);
    } else {
	y = 0.;
    }
    return y;
}

//
//               y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YJpsiPPdummy(Double_t x, Double_t energy)
{
// J/Psi y
// pp
// from the fit of RHIC + LHC data, see arXiv:1103.2394
//
    x = x/TMath::Log(energy/3.097);
    x = x*x;
    Double_t y = TMath::Exp(-x/0.4/0.4/2);
    if(x > 1) y=0;
    return y;
}

Double_t AliGenMUONlib::YJpsiPPpoly(Double_t x, Double_t energy)
{
// J/Psi y
// pp
// from the fit of RHIC + LHC data, see arXiv:1103.2394
//
    x = x/TMath::Log(energy/3.097);
    x = x*x;
    Double_t y = 1 - 6.9*x*x;
    if(y < 0) y=0;
    return y;
}

Double_t AliGenMUONlib::YJpsiPP7000(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 7 TeV
//
  return YJpsiPPdummy(*px, 7000);
}

Double_t AliGenMUONlib::YJpsiPP2760(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 2.76 TeV
//
  return YJpsiPPdummy(*px, 2760);
}

Double_t AliGenMUONlib::YJpsiPP8800(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 8.8 TeV
//
  return YJpsiPPdummy(*px, 8800);
}

Double_t AliGenMUONlib::YJpsiPPpoly7000(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 7 TeV
//
  return YJpsiPPpoly(*px, 7000);
}

Double_t AliGenMUONlib::YJpsiPPpoly2760(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 2.76 TeV
//
  return YJpsiPPpoly(*px, 2760);
}

Double_t AliGenMUONlib::YJpsiPbPb2760ShFdummy(Double_t x, Int_t n)
{
// J/Psi shadowing factor vs y for PbPb min. bias and 11 centr. bins
//
// PbPb 2.76 TeV, for EKS98, minimum bias shadowing factor = 0.66 in 4pi
//
  const Double_t f1[12] = {1, 1.128, 1.097, 1.037, 0.937, 0.821, 0.693, 0.558,
			   0.428, 0.317, 0.231, 0.156};
  const Double_t f2[12] = {1, 1.313, 1.202, 1.039, 0.814, 0.593, 0.391, 0.224,
			   0.106, 0.041, 0.013, 0.002};
  const Double_t c1[5] = {1.5591e+00, 7.5827e-03, 2.0676e-03,-1.1717e-04, 1.5237e-06}; 
  const Double_t c2[5] = {6.0861e-01, 4.8854e-03, 1.3685e-03,-7.9182e-05, 1.0475e-06}; 

  x = x*x;
  Double_t y1, y2;
  Int_t j;
  y1 = c1[j = 4]; y2 = c2[4];
  while (j > 0) {y1 = y1 * x + c1[--j]; y2 = y2 * x + c2[j];}
  
  y1 = 1 + (y1-2)*f1[n] + (y2+1-y1)*f2[n];
  if(y1<0) y1=0;
  return y1;
}

Double_t AliGenMUONlib::YJpsiPbPb2760(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, minimum bias 0-100 %
//
  return YJpsiPbPb2760ShFdummy(*px, 0) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c1(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 1st centrality bin 0-5 %
//
  return YJpsiPbPb2760ShFdummy(*px, 1) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c2(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 2nd centrality bin 5-10 %
//
  return YJpsiPbPb2760ShFdummy(*px, 2) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c3(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 3rd centrality bin 10-20 %
//
  return YJpsiPbPb2760ShFdummy(*px, 3) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c4(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 4th centrality bin 20-30 %
//
  return YJpsiPbPb2760ShFdummy(*px, 4) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c5(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 5th centrality bin 30-40 %
//
  return YJpsiPbPb2760ShFdummy(*px, 5) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c6(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 6th centrality bin 40-50 %
//
  return YJpsiPbPb2760ShFdummy(*px, 6) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c7(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 7th centrality bin 50-60 %
//
  return YJpsiPbPb2760ShFdummy(*px, 7) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c8(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 8th centrality bin 60-70 %
//
  return YJpsiPbPb2760ShFdummy(*px, 8) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c9(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 9th centrality bin 70-80 %
//
  return YJpsiPbPb2760ShFdummy(*px, 9) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c10(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 10th centrality bin 80-90 %
//
  return YJpsiPbPb2760ShFdummy(*px, 10) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPbPb2760c11(const Double_t *px, const Double_t *dummy)
{
// J/Psi y
// PbPb 2.76 TeV, 11th centrality bin 90-100 %
//
  return YJpsiPbPb2760ShFdummy(*px, 11) * YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiPP8800dummy(Double_t px)
{
    return AliGenMUONlib::YJpsiPP8800(&px, (Double_t*) 0);
}

Double_t AliGenMUONlib::YJpsiPPb8800ShFdummy(Double_t x, Int_t n)
{
// J/Psi shadowing factor vs y for pPb min. bias and 4 centr. bins
//
// pPb 8.8 TeV, for EKS98, minimum bias shadowing factor = 0.78 in 4pi
//
    const Double_t f[5] = {1, 1.33, 1.05, 0.67, 0.23};
    const Double_t c[7] = {7.4372e-01, 2.3299e-02, 2.8678e-03, 1.9595e-03, 
			   3.2849e-04,-4.0547e-05,-7.9732e-06}; 
    Double_t y;
    Int_t j;
    y = c[j = 6];
    while (j > 0) y = y * x + c[--j];
    if(y<0) y=0;
    //
    return 1 +(y-1)*f[n];
}

Double_t AliGenMUONlib::YJpsiPPb8800(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pPb 8.8 TeV, minimum bias 0-100 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 0) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPPb8800c1(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pPb 8.8 TeV, 1st centrality bin 0-20 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 1) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPPb8800c2(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pPb 8.8 TeV, 2nd centrality bin 20-40 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 2) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPPb8800c3(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pPb 8.8 TeV, 3rd centrality bin 40-60 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 3) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPPb8800c4(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pPb 8.8 TeV, 4th centrality bin 60-100 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 4) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPbP8800(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// Pbp 8.8 TeV, minimum bias 0-100 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 0) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPbP8800c1(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// Pbp 8.8 TeV, 1st centrality bin 0-20 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 1) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPbP8800c2(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// Pbp 8.8 TeV, 2nd centrality bin 20-40 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 2) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPbP8800c3(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// Pbp 8.8 TeV, 3rd centrality bin 40-60 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 3) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsiPbP8800c4(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// Pbp 8.8 TeV, 4th centrality bin 60-100 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YJpsiPPb8800ShFdummy(x, 4) * YJpsiPP8800dummy(x);
}

Double_t AliGenMUONlib::YJpsi(const Double_t *py, const Double_t */*dummy*/)
{
// J/psi y
  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}

Double_t AliGenMUONlib::YJpsiFlat( const Double_t */*py*/, const Double_t */*dummy*/ )
{
  return 1.;
}


Double_t AliGenMUONlib::YJpsiPbPb( const Double_t *px, const Double_t */*dummy*/)
{

//
// J/Psi y
//
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Double_t c[5] = {-6.03425e+02, 4.98257e+02, -1.38794e+02, 1.62209e+01, -6.85955e-01};
    Double_t x = TMath::Abs(px[0]);
    Double_t y;
    
    if (x < 4.) {
	y = 31.754;
    } else if (x < 6) {
	Int_t j;
	y = c[j = 4];
	while (j > 0) y  = y * x + c[--j];
    } else {
	y =0.;
    }
    
    return y;
}

Double_t AliGenMUONlib::YJpsiCDFscaled( const Double_t *px, const Double_t* dummy)
{
    // J/Psi y 
    return AliGenMUONlib::YJpsiPbPb(px, dummy);
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP( const Double_t *px, const Double_t* dummy)
{
    // J/Psi y 
    return AliGenMUONlib::YJpsiPP(px, dummy);
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP10( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
//
// pp 10 TeV
// scaled from YJpsiPP(14 TeV) using 10 TeV / 14 TeV ratio of y-spectra in LO pQCD. 
// see S.Grigoryan, PWG3 Meeting, 27th Oct 2008
//

    Double_t c[5] = {2.46681e+01, 8.91486e+01, -3.21227e+01, 3.63075e+00, -1.32047e-01};

    Double_t x = TMath::Abs(px[0]);
    Double_t y;

    if (x < 3.2) {
        y = 98.523 - 1.3664 * x * x;
    } else if (x < 7.5) {
        Int_t j;
        y = c[j = 4];
        while (j > 0) y  = y * x + c[--j];
    } else {
        y =0.;
    }

    if(y<0) y=0;

    return y;
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP9( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
//
// pp 8.8 TeV
// rescaling of YJpsiPP(14 TeV) using 8.8 TeV / 14 TeV ratio of y-spectra in LO QCD 
//
    Double_t c[5] = {3.33882e+02, -1.30980e+02, 2.59082e+01, -3.08935e+00, 1.56375e-01};
    Double_t x = TMath::Abs(px[0]);
    Double_t y;

    if (x < 3.7) {
        y = 99.236 - 1.5498 * x * x;
    } else if (x < 7.4) {
	Int_t j;
	y = c[j = 4];
	while (j > 0) y  = y * x + c[--j];
    } else {
	y =0.;
    }

    if(y<0) y=0;

    return y;
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP9dummy(Double_t px)
{
    return AliGenMUONlib::YJpsiCDFscaledPP9(&px, (Double_t*) 0);
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP7( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
//
// pp 7 TeV
// scaled from YJpsiPP(14 TeV) using 7 TeV / 14 TeV ratio of y-spectra in LO pQCD. 
//

    Double_t c[5] = {6.71181e+02, -3.69240e+02, 8.89644e+01, -1.04937e+01, 4.80959e-01};

    Double_t x = TMath::Abs(px[0]);
    Double_t y;

    if (x < 4.0) {
        y = 100.78 - 1.8353 * x * x;
    } else if (x < 7.3) {
        Int_t j;
        y = c[j = 4];
        while (j > 0) y  = y * x + c[--j];
    } else {
        y =0.;
    }

    if(y<0) y=0;

    return y;
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP4( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
//
// pp 3.94 TeV
// rescaling of YJpsiPP(14 TeV) using 3.94 TeV / 14 TeV ratio of y-spectra in LO QCD 
//
    Double_t c[5] = {4.00785e+02, -1.41159e+01, -3.28599e+01, 5.53048e+00, -2.45151e-01};
    Double_t x = TMath::Abs(px[0]);
    Double_t y;

    if (x < 5.5) {
        y = 107.389 - 2.7454 * x * x;
    } else if (x < 7.0) {
	Int_t j;
	y = c[j = 4];
	while (j > 0) y  = y * x + c[--j];
    } else {
	y =0.;
    }

    if(y<0) y=0;

    return y;
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP3( const Double_t *px, const Double_t *dummy)
{
// J/Psi y 
    return AliGenMUONlib::YJpsiPP2760(px, dummy);
}

Double_t AliGenMUONlib::YJpsiCDFscaledPP2( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 1.96 TeV
//
  return YJpsiPPdummy(*px, 1960);
}

Double_t AliGenMUONlib::YJpsiPP( const Double_t *px, const Double_t */*dummy*/)
{

//
// J/Psi y
//
//
// R. Vogt 2002
// pp 14  TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//

    Double_t c[5] = {1.38532e+00, 1.00596e+02, -3.46378e+01, 3.94172e+00, -1.48319e-01};
    Double_t x = TMath::Abs(px[0]);
    Double_t y;
    
    if (x < 2.5) {
	y = 96.455 - 0.8483 * x * x;
    } else if (x < 7.9) {
	Int_t j;
	y = c[j = 4];
	while (j > 0) y  = y * x + c[--j];
    } else {
	y =0.;
    }
    
    return y;
}

Double_t AliGenMUONlib::YJpsiCDFscaledPPb9( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
//
// pPb 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.79
//
    Double_t c[7] = {7.52296e-01, 2.49917e-02, 3.36500e-03, 1.91187e-03, 2.92154e-04,
		     -4.16509e-05,-7.62709e-06}; 
    Double_t y;
    Double_t x = px[0] + 0.47;              // rapidity shift
    Int_t j;
    y = c[j = 6];
    while (j > 0) y = y * x + c[--j];
    if(y<0) y=0;

    return y * AliGenMUONlib::YJpsiCDFscaledPP9dummy(x);
}

Double_t AliGenMUONlib::YJpsiCDFscaledPbP9( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
//
// Pbp 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.79
//
    Double_t c[7] = {7.52296e-01, 2.49917e-02, 3.36500e-03, 1.91187e-03, 2.92154e-04,
		     -4.16509e-05,-7.62709e-06}; 
    Double_t y;
    Double_t x = -px[0] + 0.47;              // rapidity shift
    Int_t j;
    y = c[j = 6];
    while (j > 0) y = y * x + c[--j];
    if(y<0) y=0;

    return y * AliGenMUONlib::YJpsiCDFscaledPP9dummy(x);
}

Double_t AliGenMUONlib::YJpsiCDFscaledPbPb4( const Double_t *px, const Double_t *dummy)
{
// J/Psi y
//
// PbPb 3.94 TeV, for EKS98 with minimum bias shadowing factor 0.66
//
    Double_t c[4] = {5.95228e-01, 9.45069e-03, 2.44710e-04, -1.32894e-05}; 
    Double_t x = px[0]*px[0];
    Double_t y;
    Int_t j;
    y = c[j = 3];
    while (j > 0) y  = y * x + c[--j];
    if(y<0) y=0;

    return y * AliGenMUONlib::YJpsiCDFscaledPP4(px,dummy);
}

Double_t AliGenMUONlib::YJpsiBPbPb( const Double_t *px, const Double_t */*dummy*/)
{

//
// J/Psi from B->J/Psi X
//
//
    

    Double_t c[7] = {7.37025e-02, 0., -2.94487e-03, 0., 6.07953e-06, 0., 5.39219e-07};
    
    Double_t x = TMath::Abs(px[0]);
    Double_t y;
    
    if (x > 6.) {
	y = 0.;
    } else {
	Int_t j;
	y = c[j = 6];
	while (j > 0) y  = y * x + c[--j];
    } 
    
    return y;
}



//                 particle composition
//
Int_t AliGenMUONlib::IpJpsi(TRandom *)
{
// J/Psi composition
    return 443;
}
Int_t AliGenMUONlib::IpPsiP(TRandom *)
{
// Psi prime composition
    return 100443;
}
Int_t AliGenMUONlib::IpJpsiFamily(TRandom *)
{
// J/Psi composition
  Int_t ip;
  Float_t r = gRandom->Rndm();
  if (r < 0.98) {
    ip = 443;
  } else {
    ip = 100443;
  }
  return ip;
}



//                      Upsilon
//
//
//                  pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtUpsilonPPdummy(Double_t x, Double_t energy)
{
// Upsilon pT
// pp
// from the fit of CDF & LHC data, like for J/Psi in arXiv:1103.2394
//
  const Double_t kpt0 = 1.96*TMath::Power(energy,0.095);
  const Double_t kxn  = 3.4;
  //
  Double_t pass1 = 1.+0.471*(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonPP7000(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
// pp 7 TeV
//
  return PtUpsilonPPdummy(*px,7000);
}

Double_t AliGenMUONlib::PtUpsilonPP2760(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
// pp 2.76 TeV
//
  return PtUpsilonPPdummy(*px,2760);
}

Double_t AliGenMUONlib::PtUpsilonPP8800(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
// pp 8.8 TeV
//
  return PtUpsilonPPdummy(*px,8800);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760ShFdummy(Double_t x, Int_t n)
{
// Usilon shadowing factor vs pT for PbPb min. bias and 11 centr. bins (in 2.5<y<4)
//
// PbPb 2.76 TeV, for EKS98, minimum bias shadowing factor = 0.87 in 4pi
//
  const Double_t f1[12] = {1, 1.128, 1.097, 1.037, 0.937, 0.821, 0.693, 0.558,
			   0.428, 0.317, 0.231, 0.156};
  const Double_t f2[12] = {1, 1.313, 1.202, 1.039, 0.814, 0.593, 0.391, 0.224,
			   0.106, 0.041, 0.013, 0.002};
  const Double_t c1[7] = {1.9089e+00, 1.2969e-03, 8.9786e-05,-5.3062e-06,
			  -1.0046e-06,6.1446e-08, 1.0885e-09};
  const Double_t c2[7] = {8.8423e-01,-8.7488e-05, 5.9857e-04,-5.7959e-05, 
			  2.0059e-06,-2.7343e-08, 6.6053e-10};
  Double_t y1, y2;
  Int_t j;
  y1 = c1[j = 6]; y2 = c2[6];
  while (j > 0) {y1 = y1 * x + c1[--j]; y2 = y2 * x + c2[j];}
  
  y1 /= 1.+c1[6]*TMath::Power(x,6);
  y2 /= 1.+c2[6]*TMath::Power(x,6);
  //  
  y1 = 1 + (y1-2)*f1[n] + (y2+1-y1)*f2[n];
  if(y1<0) y1=0;
  return y1;
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, minimum bias 0-100 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 0) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c1(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 1st centrality bin 0-5 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 1) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c2(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 2nd centrality bin 5-10 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 2) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c3(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 3rd centrality bin 10-20 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 3) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c4(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 4th centrality bin 20-30 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 4) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c5(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 5th centrality bin 30-40 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 5) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c6(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 6th centrality bin 40-50 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 6) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c7(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 7th centrality bin 50-60 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 7) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c8(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 8th centrality bin 60-70 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 8) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c9(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 9th centrality bin 70-80 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 9) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c10(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 10th centrality bin 80-90 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 10) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbPb2760c11(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// PbPb 2.76 TeV, 11th centrality bin 90-100 %
//
  return PtUpsilonPbPb2760ShFdummy(*px, 11) * PtUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPPb8800ShFdummy(Double_t x, Int_t n)
{
// Upsilon shadowing factor vs pT for pPb min. bias and 4 centr. bins (in 2.5<y<4)
//
// pPb 8.8 TeV, for EKS98, minimum bias shadowing factor = 0.89 in 4pi
//
  const Double_t f[5] = {1, 1.33, 1.05, 0.67, 0.23};
  const Double_t c[5] = {7.6561e-01, 1.1360e-04, 4.9596e-04,-3.0287e-05, 3.7555e-06};
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  y /= 1 + c[4]*TMath::Power(x,4);
  //  
  return 1 + (y-1)*f[n];
}

Double_t AliGenMUONlib::PtUpsilonPPb8800(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// pPb 8.8 TeV, minimum bias 0-100 %
//
  return PtUpsilonPPb8800ShFdummy(*px, 0) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPPb8800c1(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// pPb 8.8 TeV, 1st centrality bin 0-20 %
//
  return PtUpsilonPPb8800ShFdummy(*px, 1) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPPb8800c2(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// pPb 8.8 TeV, 2nd centrality bin 20-40 %
//
  return PtUpsilonPPb8800ShFdummy(*px, 2) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPPb8800c3(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// pPb 8.8 TeV, 3rd centrality bin 40-60 %
//
  return PtUpsilonPPb8800ShFdummy(*px, 3) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPPb8800c4(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// pPb 8.8 TeV, 4th centrality bin 60-100 %
//
  return PtUpsilonPPb8800ShFdummy(*px, 4) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbP8800ShFdummy(Double_t x, Int_t n)
{
// Upsilon shadowing factor vs pT for Pbp min. bias and 4 centr. bins (in 2.5<y<4)
//
// Pbp 8.8 TeV, for EKS98, minimum bias shadowing factor = 0.89 in 4pi
//
  const Double_t f[5] = {1, 1.33, 1.05, 0.67, 0.23};
  const Double_t c[5] = {1.0975, 3.1905e-03,-2.0477e-04, 8.5270e-06, 2.5343e-06};
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  y /= 1 + c[4]*TMath::Power(x,4);
  //  
  return 1 + (y-1)*f[n];
}

Double_t AliGenMUONlib::PtUpsilonPbP8800(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// Pbp 8.8 TeV, minimum bias 0-100 %
//
  return PtUpsilonPbP8800ShFdummy(*px, 0) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbP8800c1(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// Pbp 8.8 TeV, 1st centrality bin 0-20 %
//
  return PtUpsilonPbP8800ShFdummy(*px, 1) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbP8800c2(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// Pbp 8.8 TeV, 2nd centrality bin 20-40 %
//
  return PtUpsilonPbP8800ShFdummy(*px, 2) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbP8800c3(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// Pbp 8.8 TeV, 3rd centrality bin 40-60 %
//
  return PtUpsilonPbP8800ShFdummy(*px, 3) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilonPbP8800c4(const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
// Pbp 8.8 TeV, 4th centrality bin 60-100 %
//
  return PtUpsilonPbP8800ShFdummy(*px, 4) * PtUpsilonPP8800(px, dummy);
}

Double_t AliGenMUONlib::PtUpsilon( const Double_t *px, const Double_t */*dummy*/ )
{
// Upsilon pT
  const Double_t kpt0 = 5.3;
  const Double_t kxn  = 2.5;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaled( const Double_t *px, const Double_t */*dummy*/ )
{
// Upsilon pT
  const Double_t kpt0 = 7.753;
  const Double_t kxn  = 3.042;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPP( const Double_t *px, const Double_t */*dummy*/ )
{
// Upsilon pT
//
// pp 14 TeV
//
// scaled from CDF data at 2 TeV

  const Double_t kpt0 = 8.610;
  const Double_t kxn  = 3.051;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPP10( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
//
// pp 10 TeV
//
// scaled from CDF data at 2 TeV

  const Double_t kpt0 = 8.235;
  const Double_t kxn  = 3.051;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPP9( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
//
// pp 8.8 TeV
// scaled from CDF data at 2 TeV
//
  const Double_t kpt0 = 8.048;
  const Double_t kxn  = 3.051;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPP7( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
//
// pp 7 TeV
//
// scaled from CDF data at 2 TeV

  const Double_t kpt0 = 7.817;
  const Double_t kxn  = 3.051;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPP4( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon pT
//
// pp 3.94 TeV
// scaled from CDF data at 2 TeV
//
  const Double_t kpt0 = 7.189;
  const Double_t kxn  = 3.051;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPPb9( const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
//
// pPb 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.89
//
  Double_t c[5] = {7.64952e-01, 1.12501e-04, 4.96038e-04, -3.03198e-05, 3.74035e-06};
  Double_t x=*px;
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  //  
  Double_t d = 1.+c[4]*TMath::Power(x,4);
  return y/d * AliGenMUONlib::PtUpsilonCDFscaledPP9(px,dummy);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPbP9( const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
//
// Pbp 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.89
//
  Double_t c[5] = {1.09881e+00, 3.08329e-03, -2.00356e-04, 8.28991e-06, 2.52576e-06};
  Double_t x=*px;
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  //  
  Double_t d = 1.+c[4]*TMath::Power(x,4);
  return y/d * AliGenMUONlib::PtUpsilonCDFscaledPP9(px,dummy);
}

Double_t AliGenMUONlib::PtUpsilonCDFscaledPbPb4( const Double_t *px, const Double_t *dummy)
{
// Upsilon pT
//
// PbPb 3.94 TeV, for EKS98 with minimum bias shadowing factor 0.85
//
  Double_t c[5] = {8.65872e-01, 2.05465e-03, 2.56063e-04, -1.65598e-05, 2.29209e-06};
  Double_t x=*px;
  Double_t y;
  Int_t j;
  y = c[j = 4];
  while (j > 0) y  = y * x + c[--j];
  //  
  Double_t d = 1.+c[4]*TMath::Power(x,4);
  return y/d * AliGenMUONlib::PtUpsilonCDFscaledPP4(px,dummy);
}

Double_t AliGenMUONlib::PtUpsilonFlat( const Double_t */*px*/, const Double_t */*dummy*/ )
{
  return 1.;
}

Double_t AliGenMUONlib::PtUpsilonPbPb( const Double_t *px, const Double_t */*dummy*/)
{
//
// Upsilon pT
//
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Float_t x = px[0];
    Double_t c[8] = {
	-1.03488e+01, 1.28065e+01, -6.60500e+00, 1.66140e+00,       
	-2.34293e-01, 1.86925e-02, -7.80708e-04, 1.30610e-05
    };
    Double_t y;
    if (x < 10.) {
	Int_t j;
	y = c[j = 7];
	while (j > 0) y  = y * x +c[--j];
	y = x * TMath::Exp(y);
    } else {
	y = 0.;
    }
    return y;
}

Double_t AliGenMUONlib::PtUpsilonPP( const Double_t *px, const Double_t */*dummy*/)
{
//
// Upsilon pT
//
//
// R. Vogt 2002
// pp 14 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Float_t x = px[0];
    Double_t c[8] = {-7.93955e+00, 1.06306e+01, -5.21392e+00, 1.19703e+00,   
		     -1.45718e-01, 8.95151e-03, -2.04806e-04, -1.13053e-06};
    
    Double_t y;
    if (x < 10.) {
	Int_t j;
	y = c[j = 7];
	while (j > 0) y  = y * x +c[--j];
	y = x * TMath::Exp(y);
    } else {
	y = 0.;
    }
    return y;
}

//
//                    y-distribution
//
//____________________________________________________________
Double_t AliGenMUONlib::YUpsilonPPdummy(Double_t x, Double_t energy)
{
// Upsilon y
// pp
// from the fit of CDF & LHC data, like for J/Psi in arXiv:1103.2394
//
    x = x/TMath::Log(energy/9.46);
    x = x*x;
    Double_t y = TMath::Exp(-x/0.4/0.4/2);
    if(x > 1) y=0;
    return y;
}

Double_t AliGenMUONlib::YUpsilonPPpoly(Double_t x, Double_t energy)
{
// Upsilon y
// pp
// from the fit of CDF & LHC data, like for J/Psi in arXiv:1103.2394
//
    x = x/TMath::Log(energy/9.46);
    x = x*x;
    Double_t y = 1 - 6.9*x*x;
    if(y < 0) y=0;
    return y;
}

Double_t AliGenMUONlib::YUpsilonPP7000(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pp 7 TeV
//
  return YUpsilonPPdummy(*px, 7000);
}

Double_t AliGenMUONlib::YUpsilonPP2760(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pp 2.76 TeV
//
  return YUpsilonPPdummy(*px, 2760);
}

Double_t AliGenMUONlib::YUpsilonPP8800(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pp 8.8 TeV
//
  return YUpsilonPPdummy(*px, 8800);
}

Double_t AliGenMUONlib::YUpsilonPPpoly7000(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pp 7 TeV
//
  return YUpsilonPPpoly(*px, 7000);
}

Double_t AliGenMUONlib::YUpsilonPPpoly2760(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pp 2.76 TeV
//
  return YUpsilonPPpoly(*px, 2760);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760ShFdummy(Double_t x, Int_t n)
{
// Upsilon shadowing factor vs y for PbPb min. bias and 11 centr. bins
//
// PbPb 2.76 TeV, for EKS98, minimum bias shadowing factor = 0.87 in 4pi
//
  const Double_t f1[12] = {1, 1.128, 1.097, 1.037, 0.937, 0.821, 0.693, 0.558,
			   0.428, 0.317, 0.231, 0.156};
  const Double_t f2[12] = {1, 1.313, 1.202, 1.039, 0.814, 0.593, 0.391, 0.224,
			   0.106, 0.041, 0.013, 0.002};
  const Double_t c1[5] = {1.8547e+00, 1.6717e-02,-2.1285e-04,-9.7662e-05, 2.5768e-06};
  const Double_t c2[5] = {8.6029e-01, 1.1742e-02,-2.7944e-04,-6.7973e-05, 1.8838e-06}; 

  x = x*x;
  Double_t y1, y2;
  Int_t j;
  y1 = c1[j = 4]; y2 = c2[4];
  while (j > 0) {y1 = y1 * x + c1[--j]; y2 = y2 * x + c2[j];}
  
  y1 = 1 + (y1-2)*f1[n] + (y2+1-y1)*f2[n];
  if(y1<0) y1=0;
  return y1;
}

Double_t AliGenMUONlib::YUpsilonPbPb2760(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, minimum bias 0-100 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 0) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c1(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 1st centrality bin 0-5 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 1) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c2(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 2nd centrality bin 5-10 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 2) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c3(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 3rd centrality bin 10-20 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 3) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c4(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 4th centrality bin 20-30 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 4) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c5(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 5th centrality bin 30-40 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 5) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c6(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 6th centrality bin 40-50 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 6) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c7(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 7th centrality bin 50-60 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 7) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c8(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 8th centrality bin 60-70 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 8) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c9(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 9th centrality bin 70-80 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 9) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c10(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 10th centrality bin 80-90 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 10) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPbPb2760c11(const Double_t *px, const Double_t *dummy)
{
// Upsilon y
// PbPb 2.76 TeV, 11th centrality bin 90-100 %
//
  return YUpsilonPbPb2760ShFdummy(*px, 11) * YUpsilonPP2760(px, dummy);
}

Double_t AliGenMUONlib::YUpsilonPP8800dummy(Double_t px)
{
    return AliGenMUONlib::YUpsilonPP8800(&px, (Double_t*) 0);
}

Double_t AliGenMUONlib::YUpsilonPPb8800ShFdummy(Double_t x, Int_t n)
{
// Upsilon shadowing factor vs y for pPb min. bias and 4 centr. bins
//
// pPb 8.8 TeV, for EKS98, minimum bias shadowing factor = 0.89 in 4pi
//
    const Double_t f[5] = {1, 1.33, 1.05, 0.67, 0.23};
    const Double_t c[7] = {8.6581e-01, 4.6111e-02, 7.6911e-03, 8.7313e-04,
			   -1.4700e-04,-5.0975e-05,-3.5718e-06}; 
    Double_t y;
    Int_t j;
    y = c[j = 6];
    while (j > 0) y = y * x + c[--j];
    if(y<0) y=0;
    //
    return 1 +(y-1)*f[n];
}

Double_t AliGenMUONlib::YUpsilonPPb8800(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pPb 8.8 TeV, minimum bias 0-100 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 0) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPPb8800c1(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pPb 8.8 TeV, 1st centrality bin 0-20 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 1) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPPb8800c2(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pPb 8.8 TeV, 2nd centrality bin 20-40 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 2) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPPb8800c3(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pPb 8.8 TeV, 3rd centrality bin 40-60 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 3) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPPb8800c4(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// pPb 8.8 TeV, 4th centrality bin 60-100 %
//
  Double_t x = px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 4) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPbP8800(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// Pbp 8.8 TeV, minimum bias 0-100 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 0) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPbP8800c1(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// Pbp 8.8 TeV, 1st centrality bin 0-20 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 1) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPbP8800c2(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// Pbp 8.8 TeV, 2nd centrality bin 20-40 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 2) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPbP8800c3(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// Pbp 8.8 TeV, 3rd centrality bin 40-60 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 3) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilonPbP8800c4(const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
// Pbp 8.8 TeV, 4th centrality bin 60-100 %
//
  Double_t x = -px[0] + 0.47;              // rapidity shift
  return YUpsilonPPb8800ShFdummy(x, 4) * YUpsilonPP8800dummy(x);
}

Double_t AliGenMUONlib::YUpsilon(const Double_t *py, const Double_t */*dummy*/)
{
// Upsilon y
  const Double_t ky0 = 3.;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yu=kb;
  else
    yu=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yu;
}


Double_t AliGenMUONlib::YUpsilonPbPb( const Double_t *px, const Double_t */*dummy*/)
{

//
// Upsilon y
//
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//

    Double_t c[7] = {3.40036e-01, -3.98882e-07, -4.48398e-03, 8.46411e-08, -6.10854e-04,
		     -2.99753e-09, 1.28895e-05};
    Double_t x = TMath::Abs(px[0]);
    if (x > 5.55) return 0.;
    Int_t j;
    Double_t y = c[j = 6];
    while (j > 0) y  = y * x +c[--j];
    return y;
}

Double_t AliGenMUONlib::YUpsilonCDFscaled( const Double_t *px, const Double_t *dummy)
{
    // Upsilon y
    return AliGenMUONlib::YUpsilonPbPb(px, dummy);
    
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPP( const Double_t *px, const Double_t *dummy)
{
    // Upsilon y
    return AliGenMUONlib::YUpsilonPP(px, dummy);
    
}

Double_t AliGenMUONlib::YUpsilonFlat( const Double_t */*px*/, const Double_t */*dummy*/)
{
    // Upsilon y
    return 1.;
    
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPP10( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
//
// pp 10 TeV
// scaled from YUpsilonPP(14 TeV) using 10 TeV / 14 TeV ratio of y-spectra in LO pQCD. 
// see S.Grigoryan, PWG3 Meeting, 27th Oct 2008
//
    Double_t c[4] = {1., -2.17877e-02, -6.52830e-04, 1.40578e-05};
    Double_t x = TMath::Abs(px[0]);
    if (x > 6.1) return 0.;
    Int_t j;
    Double_t y = c[j = 3];
    while (j > 0) y  = y * x*x +c[--j];
    return y;
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPP9( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
//
// pp 8.8 TeV
// rescaling of YUpsilonPP(14 TeV) using 8.8 TeV / 14 TeV ratio of y-spectra in LO QCD 
//
    Double_t c[4] = {1., -2.37621e-02, -6.29610e-04, 1.47976e-05};
    Double_t x = TMath::Abs(px[0]);
    if (x > 6.1) return 0.;
    Int_t j;
    Double_t y = c[j = 3];
    while (j > 0) y  = y * x*x +c[--j];
    return y;
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPP9dummy(Double_t px)
{
    return AliGenMUONlib::YUpsilonCDFscaledPP9(&px, (Double_t*) 0);
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPP7( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
//
// pp 7 TeV
// scaled from YUpsilonPP(14 TeV) using 7 TeV / 14 TeV ratio of y-spectra in LO pQCD. 
//
    Double_t c[4] = {1., -2.61009e-02, -6.83937e-04, 1.78451e-05};
    Double_t x = TMath::Abs(px[0]);
    if (x > 6.0) return 0.;
    Int_t j;
    Double_t y = c[j = 3];
    while (j > 0) y  = y * x*x +c[--j];
    return y;
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPP4( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
//
// pp 3.94 TeV
// rescaling of YUpsilonPP(14 TeV) using 3.94 TeV / 14 TeV ratio of y-spectra in LO QCD
//
    Double_t c[4] = {1., -3.91924e-02, -4.26184e-04, 2.10914e-05};
    Double_t x = TMath::Abs(px[0]);
    if (x > 5.7) return 0.;
    Int_t j;
    Double_t y = c[j = 3];
    while (j > 0) y  = y * x*x +c[--j];

    return y;
}

Double_t AliGenMUONlib::YUpsilonPP( const Double_t *px, const Double_t */*dummy*/)
{

//
// Upsilon y
//
//
// R. Vogt 2002
// p p  14. TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Double_t c[7] = {8.91936e-01, -6.46645e-07, -1.52774e-02, 4.28677e-08, -7.01517e-04, 
		     -6.20539e-10, 1.29943e-05};
    Double_t x = TMath::Abs(px[0]);
    if (x > 6.2) return 0.;
    Int_t j;
    Double_t y = c[j = 6];
    while (j > 0) y  = y * x +c[--j];
    return y;
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPPb9( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
//
// pPb 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.89
//
    Double_t c[7] = {8.71829e-01, 4.77467e-02, 8.09671e-03, 6.45294e-04, -2.15730e-04,
		     -4.67538e-05,-2.11683e-06}; 
    Double_t y;
    Double_t x = px[0] + 0.47;              // rapidity shift
    Int_t j;
    y = c[j = 6];
    while (j > 0) y = y * x + c[--j];
    if(y<0) y=0;

    return y * AliGenMUONlib::YUpsilonCDFscaledPP9dummy(x);
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPbP9( const Double_t *px, const Double_t */*dummy*/)
{
// Upsilon y
//
// Pbp 8.8 TeV, for EKS98 with minimum bias shadowing factor 0.89
//
    Double_t c[7] = {8.71829e-01, 4.77467e-02, 8.09671e-03, 6.45294e-04, -2.15730e-04,
		     -4.67538e-05,-2.11683e-06}; 
    Double_t y;
    Double_t x = -px[0] + 0.47;              // rapidity shift
    Int_t j;
    y = c[j = 6];
    while (j > 0) y = y * x + c[--j];
    if(y<0) y=0;

    return y * AliGenMUONlib::YUpsilonCDFscaledPP9dummy(x);
}

Double_t AliGenMUONlib::YUpsilonCDFscaledPbPb4( const Double_t *px, const Double_t *dummy)
{
// Upsilon y
//
// PbPb 3.94 TeV, for EKS98 with minimum bias shadowing factor 0.85
//
    Double_t c[4] = {8.27837e-01, 1.70115e-02, -1.26046e-03, 1.52091e-05}; 
    Double_t x = px[0]*px[0];
    Double_t y;
    Int_t j;
    y = c[j = 3];
    while (j > 0) y  = y * x + c[--j];
    if(y<0) y=0;

    return y * AliGenMUONlib::YUpsilonCDFscaledPP4(px,dummy);
}


//                 particle composition
//
Int_t AliGenMUONlib::IpUpsilon(TRandom *)
{
// y composition
    return 553;
}
Int_t AliGenMUONlib::IpUpsilonP(TRandom *)
{
// y composition
    return 100553;
}
Int_t AliGenMUONlib::IpUpsilonPP(TRandom *)
{
// y composition
    return 200553;
}
Int_t AliGenMUONlib::IpUpsilonFamily(TRandom *)
{
// y composition
  Int_t ip;
  Float_t r = gRandom->Rndm();
  
  if (r < 0.712) {
    ip = 553;
  } else if (r < 0.896) {
    ip = 100553;
  } else {
    ip = 200553;
  }
  return ip;
}


//
//                        Phi
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtPhi( const Double_t *px, const Double_t */*dummy*/)
{
// Phi pT
  return PtScal(*px,7);
}
//    y-distribution
Double_t AliGenMUONlib::YPhi( const Double_t *px, const Double_t */*dummy*/)
{
// Phi y
    Double_t *dum=0;
    return YJpsi(px,dum);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPhi(TRandom *)
{
// Phi composition
    return 333;
}

//
//                        omega
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtOmega( const Double_t *px, const Double_t */*dummy*/)
{
// Omega pT
  return PtScal(*px,5);
}
//    y-distribution
Double_t AliGenMUONlib::YOmega( const Double_t *px, const Double_t */*dummy*/)
{
// Omega y
    Double_t *dum=0;
    return YJpsi(px,dum);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpOmega(TRandom *)
{
// Omega composition
    return 223;
}


//
//                        Eta
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtEta( const Double_t *px, const Double_t */*dummy*/)
{
// Eta pT
  return PtScal(*px,3);
}
//    y-distribution
Double_t AliGenMUONlib::YEta( const Double_t *px, const Double_t */*dummy*/)
{
// Eta y
    Double_t *dum=0;
    return YJpsi(px,dum);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpEta(TRandom *)
{
// Eta composition
    return 221;
}

//
//                        Charm
//
//
//                    pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtCharm( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
  const Double_t kpt0 = 2.25;
  const Double_t kxn  = 3.17;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtCharmCentral( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
  const Double_t kpt0 = 2.12;
  const Double_t kxn  = 2.78;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtCharmF0M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// FiMjSkPP define theoretical uncertainties around F0M0S0PP as follows:
// PtCharmFiMjSkPP = PtCharmF0M0S0PP * (dN(i,j,k)/dpt / dN(0,0,0)/dpt)_MNR
//       i=0,1,2;  j=0,1,2;  k=0,1,...,6
// dN(i,j,k)/dpt - spectra obtained by A.Dainese (hep-ph/0601164, p.88; 
// http://www-zeus.desy.de/~corradi/benchmarks) from NLO pQCD (MNR)
// calculations for the following inputs: 
// Peterson fragmentation function (F) with \epsilon_c = 0.02, 0.002 & 0.11 
// for i=0,1 & 2 respectively; quark mass (M) of 1.5, 1.3 & 1.7 GeV 
// for j=0,1 & 2 respectively; 
// factorisation \mu_F = a*mt and renormalisation \mu_R = b*mt scales (S) 
// with a/b = 1/1, 1/0.5, 0.5/1, 0.5/0.5, 1/2, 2/1 & 2/2 
// for k = 0, 1, 2, 3, 4, 5 & 6 respectively; CTEQ6.1 PDF set 
// (PDF uncertainty not considered since is small, see hep-ph/0601164, p.89).
// June 2008, Smbat.Grigoryan@cern.ch

// Charm pT
// Pythia6.214 (kCharmppMNRwmi, PDF = CTEQ5L, quark mass = 1.2 GeV, PtHard > 2.76 GeV/c)
// for pp collisions at 14 TeV with one c-cbar pair per event.
// Corresponding NLO total cross section is 5.68 mb


  const Double_t kpt0 = 2.2930;
  const Double_t kxn  = 3.1196;
  Double_t c[3]={-5.2180e-01,1.8753e-01,2.8669e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF1M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 6.06 mb
  const Double_t kpt0 = 2.8669;
  const Double_t kxn  = 3.1044;
  Double_t c[3]={-4.6714e-01,1.5005e-01,4.5003e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF2M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 6.06 mb
  const Double_t kpt0 = 1.8361;
  const Double_t kxn  = 3.2966;
  Double_t c[3]={-6.1550e-01,2.6498e-01,1.0728e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M1S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 7.69 mb
  const Double_t kpt0 = 2.1280;
  const Double_t kxn  = 3.1397;
  Double_t c[3]={-5.4021e-01,2.0944e-01,2.5211e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M2S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 4.81 mb
  const Double_t kpt0 = 2.4579;
  const Double_t kxn  = 3.1095;
  Double_t c[3]={-5.1497e-01,1.7532e-01,3.2429e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M0S1PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 14.09 mb
  const Double_t kpt0 = 2.1272;
  const Double_t kxn  = 3.1904;
  Double_t c[3]={-4.6088e-01,2.1918e-01,2.3055e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M0S2PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 1.52 mb
  const Double_t kpt0 = 2.8159;
  const Double_t kxn  = 3.0857;
  Double_t c[3]={-6.4691e-01,2.0289e-01,2.4922e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M0S3PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 3.67 mb
  const Double_t kpt0 = 2.7297;
  const Double_t kxn  = 3.3019;
  Double_t c[3]={-6.2216e-01,1.9031e-01,1.5341e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M0S4PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 3.38 mb
  const Double_t kpt0 = 2.3894;
  const Double_t kxn  = 3.1075;
  Double_t c[3]={-4.9742e-01,1.7032e-01,2.5994e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M0S5PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 10.37 mb
  const Double_t kpt0 = 2.0187;
  const Double_t kxn  = 3.3011;
  Double_t c[3]={-3.9869e-01,2.9248e-01,1.1763e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}
Double_t AliGenMUONlib::PtCharmF0M0S6PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm pT
// Corresponding NLO total cross section is 7.22 mb
  const Double_t kpt0 = 2.1089;
  const Double_t kxn  = 3.1848;
  Double_t c[3]={-4.6275e-01,1.8114e-01,2.1363e-02};
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn)*(1.+c[0]*x+c[1]*x*x)/(1.+c[2]*x*x);
}

//                  y-distribution
Double_t AliGenMUONlib::YCharm( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y :: Carrer & Dainese : ALICE-INT-2003-019 v.3 (hep-ph/0311225) 
// Pythia tuned to reproduce the distribution given by the HVQMNR program based on NLO calculations (pQCD)
// shadowing + kt broadening 

    Double_t x=px[0];
    Double_t c[2]={-2.42985e-03,-2.31001e-04};
    Double_t y=1+(c[0]*TMath::Power(x,2))+(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>8) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// FiMjSkPP define theoretical uncertainties around F0M0S0PP as follows:
// YCharmFiMjSkPP = YCharmF0M0S0PP * (dN(i,j,k)/dy / dN(0,0,0)/dy)_MNR
//       i=0,1,2;  j=0,1,2;  k=0,1,...,6
// dN(i,j,k)/dy - spectra obtained by A.Dainese (hep-ph/0601164, p.88; 
// http://www-zeus.desy.de/~corradi/benchmarks) from NLO pQCD (MNR) 
// calculations for the following inputs: 
// Peterson fragmentation function (F) with \epsilon_c = 0.02, 0.002 & 0.11 
// for i=0,1 & 2 respectively; quark mass (M) of 1.5, 1.3 & 1.7 GeV 
// for j=0,1 & 2 respectively; 
// factorisation \mu_F = a*mt and renormalisation \mu_R = b*mt scales (S) 
// with a/b = 1/1,1/0.5, 0.5/1, 0.5/0.5, 1/2, 2/1 & 2/2 for 
// k = 0, 1, 2, 3, 4, 5 & 6 respectively; CTEQ6.1 PDF set
// (PDF uncertainty not considered since is small, see hep-ph/0601164, p.89).
// June 2008, Smbat.Grigoryan@cern.ch

// Charm y
// Pythia6.214 (kCharmppMNRwmi, PDF = CTEQ5L, quark mass = 1.2 GeV, PtHard > 2.76 GeV/c)
// for pp collisions at 14 TeV with one c-cbar pair per event.
// Corresponding NLO total cross section is 5.68 mb

    Double_t x=px[0];
    Double_t c[2]={7.0909e-03,6.1967e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF1M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 6.06 mb
    Double_t x=px[0];
    Double_t c[2]={6.9707e-03,6.0971e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF2M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 6.06 mb
    Double_t x=px[0];
    Double_t c[2]={7.1687e-03,6.5303e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M1S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 7.69 mb
    Double_t x=px[0];
    Double_t c[2]={5.9090e-03,7.1854e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M2S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 4.81 mb
    Double_t x=px[0];
    Double_t c[2]={8.0882e-03,5.5872e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S1PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 14.09 mb
    Double_t x=px[0];
    Double_t c[2]={7.2520e-03,6.2691e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S2PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 1.52 mb
    Double_t x=px[0];
    Double_t c[2]={1.1040e-04,1.4498e-04};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S3PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 3.67 mb
    Double_t x=px[0];
    Double_t c[2]={-3.1328e-03,1.8270e-04};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S4PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 3.38 mb
    Double_t x=px[0];
    Double_t c[2]={7.0865e-03,6.2532e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S5PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 10.37 mb
    Double_t x=px[0];
    Double_t c[2]={7.7070e-03,5.3533e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}
Double_t AliGenMUONlib::YCharmF0M0S6PP( const Double_t *px, const Double_t */*dummy*/)
{
// Charm y
// Corresponding NLO total cross section is 7.22 mb
    Double_t x=px[0];
    Double_t c[2]={7.9195e-03,5.3823e-05};
    Double_t y=1-(c[0]*TMath::Power(x,2))-(c[1]*TMath::Power(x,4));
    Double_t ycharm;
    
    if (TMath::Abs(x)>9) {
      ycharm=0.;
    }
    else {
      ycharm=TMath::Power(y,3);
    }
    
    return ycharm;
}


Int_t AliGenMUONlib::IpCharm(TRandom *ran)
{  
// Charm composition
    Float_t random;
    Int_t ip;
//    411,421,431,4122
    random = ran->Rndm();
//  Taux de production Carrer & Dainese : ALICE-INT-2003-019 v.3  
//  >>>>> cf. tab 4 p 11
  
    if (random < 0.30) {                       
        ip=421;
    } else if (random < 0.60) {
        ip=-421;
    } else if (random < 0.70) {
        ip=411;
    } else if (random < 0.80) {
        ip=-411;
    } else if (random < 0.86) {
        ip=431;
    } else if (random < 0.92) {
        ip=-431;	
    } else if (random < 0.96) {
        ip=4122;
    } else {
        ip=-4122;
    }
    
    return ip;
}

//
//                        Beauty
//
//
//                    pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtBeauty( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
  const Double_t kpt0 = 6.53;
  const Double_t kxn  = 3.59;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtBeautyCentral( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
  const Double_t kpt0 = 6.14;
  const Double_t kxn  = 2.93;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// FiMjSkPP define theoretical uncertainties around F0M0S0PP as follows:
// PtBeautyFiMjSkPP = PtBeautyF0M0S0PP * (dN(i,j,k)/dpt / dN(0,0,0)/dpt)_MNR
//       i=0,1,2;  j=0,1,2;  k=0,1,...,6
// dN(i,j,k)/dpt - spectra obtained by A.Dainese (hep-ph/0601164, p.88; 
// http://www-zeus.desy.de/~corradi/benchmarks) from NLO pQCD (MNR) 
// calculations for the following inputs: 
// Peterson fragmentation function (F) with \epsilon_b = 0.001, 0.0002 & 0.004 
// for i=0,1 & 2 respectively; quark mass (M) of 4.75, 4.5 & 5.0 GeV 
// for j=0,1 & 2 respectively; 
// factorisation \mu_F = a*mt and renormalisation \mu_R = b*mt scales (S) 
// with a/b = 1/1, 1/0.5, 0.5/1, 0.5/0.5, 1/2, 2/1 & 2/2 for 
// k = 0, 1, 2, 3, 4, 5 & 6 respectively; CTEQ6.1 PDF set
// (PDF uncertainty not considered since is small, see hep-ph/0601164, p.89).
// June 2008, Smbat.Grigoryan@cern.ch

// Beauty pT
// Pythia6.214 (kBeautyppMNRwmi, PDF = CTEQ5L, quark mass = 4.75 GeV, PtHard > 2.76 GeV/c)
// for pp collisions at 14 TeV with one b-bbar pair per event.
// Corresponding NLO total cross section is 0.494 mb

  const Double_t kpt0 = 8.0575;
  const Double_t kxn  = 3.1921;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF1M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.445 mb
  const Double_t kpt0 = 8.6239;
  const Double_t kxn  = 3.2911;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF2M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.445 mb
  const Double_t kpt0 = 7.3367;
  const Double_t kxn  = 3.0692;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M1S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.518 mb
  const Double_t kpt0 = 7.6409;
  const Double_t kxn  = 3.1364;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M2S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.384 mb
  const Double_t kpt0 = 8.4948;
  const Double_t kxn  = 3.2546;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S1PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.648 mb
  const Double_t kpt0 = 7.6631;
  const Double_t kxn  = 3.1621;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S2PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.294 mb
  const Double_t kpt0 = 8.7245;
  const Double_t kxn  = 3.2213;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S3PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.475 mb
  const Double_t kpt0 = 8.5296;
  const Double_t kxn  = 3.2187;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S4PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.324 mb
  const Double_t kpt0 = 7.9440;
  const Double_t kxn  = 3.1614;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S5PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.536 mb
  const Double_t kpt0 = 8.2408;
  const Double_t kxn  = 3.3029;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtBeautyF0M0S6PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty pT
// Corresponding NLO total cross section is 0.420 mb
  const Double_t kpt0 = 7.8041;
  const Double_t kxn  = 3.2094;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

//                     y-distribution
Double_t AliGenMUONlib::YBeauty( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y :: Carrer & Dainese : ALICE-INT-2003-019 v.3 (hep-ph/0311225) 
// Pythia tuned to reproduce the distribution given by the HVQMNR program based on NLO calculations (pQCD)
// shadowing + kt broadening 

    Double_t x=px[0];
    Double_t c[2]={-1.27590e-02,-2.42731e-04};
    Double_t y=1+c[0]*TMath::Power(x,2)+c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// FiMjSkPP define theoretical uncertainties around F0M0S0PP as follows:
// YBeautyFiMjSkPP = YBeautyF0M0S0PP * (dN(i,j,k)/dy / dN(0,0,0)/dy)_MNR
//       i=0,1,2;  j=0,1,2;  k=0,1,...,6
// dN(i,j,k)/dy - spectra obtained by A.Dainese (hep-ph/0601164, p.88; 
// http://www-zeus.desy.de/~corradi/benchmarks) from NLO pQCD (MNR) 
// calculations for the following inputs: 
// Peterson fragmentation function (F) with \epsilon_b = 0.001, 0.0002 & 0.004 
// for i=0,1 & 2 respectively; quark mass (M) of 4.75, 4.5 & 5.0 GeV 
// for j=0,1 & 2 respectively; 
// factorisation \mu_F = a*mt and renormalisation \mu_R = b*mt scales (S) 
// with a/b = 1/1, 1/0.5, 0.5/1, 0.5/0.5, 1/2, 2/1 & 2/2 
// for k = 0, 1, 2, 3, 4, 5 & 6 respectively; CTEQ6.1 PDF set 
// (PDF uncertainty not considered since is small, see hep-ph/0601164, p.89).
// June 2008, Smbat.Grigoryan@cern.ch

// Beauty y
// Pythia6.214 (kBeautyppMNRwmi, PDF = CTEQ5L, quark mass = 4.75 GeV, PtHard > 2.76 GeV/c)
// for pp collisions at 14 TeV with one b-bbar pair per event.
// Corresponding NLO total cross section is 0.494 mb


    Double_t x=px[0];
    Double_t c[2]={1.2350e-02,9.2667e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF1M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.445 mb
    Double_t x=px[0];
    Double_t c[2]={1.2292e-02,9.1847e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF2M0S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.445 mb
    Double_t x=px[0];
    Double_t c[2]={1.2436e-02,9.3709e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M1S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.518 mb
    Double_t x=px[0];
    Double_t c[2]={1.1714e-02,1.0068e-04};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M2S0PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.384 mb
    Double_t x=px[0];
    Double_t c[2]={1.2944e-02,8.5500e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S1PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.648 mb
    Double_t x=px[0];
    Double_t c[2]={1.2455e-02,9.2713e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S2PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.294 mb
    Double_t x=px[0];
    Double_t c[2]={1.0897e-02,1.1878e-04};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S3PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.475 mb
    Double_t x=px[0];
    Double_t c[2]={1.0912e-02,1.1858e-04};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S4PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.324 mb
    Double_t x=px[0];
    Double_t c[2]={1.2378e-02,9.2490e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S5PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.536 mb
    Double_t x=px[0];
    Double_t c[2]={1.2886e-02,8.2912e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}
Double_t AliGenMUONlib::YBeautyF0M0S6PP( const Double_t *px, const Double_t */*dummy*/)
{
// Beauty y
// Corresponding NLO total cross section is 0.420 mb
    Double_t x=px[0];
    Double_t c[2]={1.3106e-02,8.0115e-05};
    Double_t y=1-c[0]*TMath::Power(x,2)-c[1]*TMath::Power(x,4);
    Double_t ybeauty;
    
    if (TMath::Abs(x)>7.6) {
      ybeauty=0.;
    }
    else {
      ybeauty=TMath::Power(y,3);
    }
    
    return ybeauty;
}

Int_t AliGenMUONlib::IpBeauty(TRandom *ran)
{  
// Beauty Composition
    Float_t random;
    Int_t ip;
    random = ran->Rndm(); 
    
//  Taux de production Carrer & Dainese : ALICE-INT-2003-019 v.3  
//  >>>>> cf. tab 4 p 11
    
 if (random < 0.20) {                       
        ip=511;
    } else if (random < 0.40) {
        ip=-511;
    } else if (random < 0.605) {
        ip=521;
    } else if (random < 0.81) {
        ip=-521;
    } else if (random < 0.87) {
        ip=531;
    } else if (random < 0.93) {
        ip=-531;	
    } else if (random < 0.965) {
        ip=5122;
    } else {
        ip=-5122;
    }
    
 return ip;
}


typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);
GenFunc AliGenMUONlib::GetPt(Int_t param,  const char* tname) const
{
// Return pointer to pT parameterisation
    TString sname = TString(tname);
    GenFunc func;
    switch (param) 
    {
    case kPhi:
	func=PtPhi;
	break;
    case kOmega:
	func=PtOmega;
	break;
    case kEta:
	func=PtEta;
	break;
    case kJpsiFamily:
    case kPsiP:
    case kChic1:
    case kChic2:
    case kJpsi:
	if (sname == "Vogt" || sname == "Vogt PbPb") {
	    func=PtJpsiPbPb;
	} else if (sname == "Vogt pp") {
	    func=PtJpsiPP;
	} else if (sname == "pp 7") {
	    func=PtJpsiPP7000;
	} else if (sname == "pp 2.76") {
	    func=PtJpsiPP2760;
	} else if (sname == "PbPb 2.76") {
	    func=PtJpsiPbPb2760;
	} else if (sname == "PbPb 2.76c1") {
	    func=PtJpsiPbPb2760c1;
	} else if (sname == "PbPb 2.76c2") {
	    func=PtJpsiPbPb2760c2;
	} else if (sname == "PbPb 2.76c3") {
	    func=PtJpsiPbPb2760c3;
	} else if (sname == "PbPb 2.76c4") {
	    func=PtJpsiPbPb2760c4;
	} else if (sname == "PbPb 2.76c5") {
	    func=PtJpsiPbPb2760c5;
	} else if (sname == "PbPb 2.76c6") {
	    func=PtJpsiPbPb2760c6;
	} else if (sname == "PbPb 2.76c7") {
	    func=PtJpsiPbPb2760c7;
	} else if (sname == "PbPb 2.76c8") {
	    func=PtJpsiPbPb2760c8;
	} else if (sname == "PbPb 2.76c9") {
	    func=PtJpsiPbPb2760c9;
	} else if (sname == "PbPb 2.76c10") {
	    func=PtJpsiPbPb2760c10;
	} else if (sname == "PbPb 2.76c11") {
	    func=PtJpsiPbPb2760c11;
	} else if (sname == "pp 7 poly") {
	    func=PtJpsiPP7000;
	} else if (sname == "pp 2.76 poly") {
	    func=PtJpsiPP2760;
	} else if (sname == "pp 8.8") {
	    func=PtJpsiPP8800;
	} else if (sname == "pPb 8.8") {
	    func=PtJpsiPPb8800;
	} else if (sname == "pPb 8.8c1") {
	    func=PtJpsiPPb8800c1;
	} else if (sname == "pPb 8.8c2") {
	    func=PtJpsiPPb8800c2;
	} else if (sname == "pPb 8.8c3") {
	    func=PtJpsiPPb8800c3;
	} else if (sname == "pPb 8.8c4") {
	    func=PtJpsiPPb8800c4;
	} else if (sname == "Pbp 8.8") {
	    func=PtJpsiPbP8800;
	} else if (sname == "Pbp 8.8c1") {
	    func=PtJpsiPbP8800c1;
	} else if (sname == "Pbp 8.8c2") {
	    func=PtJpsiPbP8800c2;
	} else if (sname == "Pbp 8.8c3") {
	    func=PtJpsiPbP8800c3;
	} else if (sname == "Pbp 8.8c4") {
	    func=PtJpsiPbP8800c4;
	} else if (sname == "CDF scaled") {
	    func=PtJpsiCDFscaled;
	} else if (sname == "CDF pp") {
	    func=PtJpsiCDFscaledPP;
	} else if (sname == "CDF pp 10") {
	    func=PtJpsiCDFscaledPP10;
	} else if (sname == "CDF pp 8.8") {
	    func=PtJpsiCDFscaledPP9;
	} else if (sname == "CDF pp 7" || sname == "CDF pp 7 flat y") {
	    func=PtJpsiCDFscaledPP7;
	} else if (sname == "CDF pp 3.94") {
	    func=PtJpsiCDFscaledPP4;
	} else if (sname == "CDF pp 2.76") {
	    func=PtJpsiCDFscaledPP3;
	} else if (sname == "CDF pp 1.9") {
	    func=PtJpsiCDFscaledPP2;
	} else if (sname == "CDF pPb 8.8") {
	    func=PtJpsiCDFscaledPPb9;
	} else if (sname == "CDF Pbp 8.8") {
	    func=PtJpsiCDFscaledPbP9;
	} else if (sname == "CDF PbPb 3.94") {
	    func=PtJpsiCDFscaledPbPb4;
	} else if (sname == "Flat" || sname == "CDF pp 7 flat pt") {
	    func=PtJpsiFlat;
	} else {
	    func=PtJpsi;
	}
	break;
    case kJpsiFromB:
	func = PtJpsiBPbPb;
	break;
    case kUpsilonFamily:
    case kUpsilonP:
    case kUpsilonPP:
    case kUpsilon:
	if (sname == "Vogt" || sname == "Vogt PbPb") {
	    func=PtUpsilonPbPb;
	} else if (sname == "Vogt pp") {
	    func=PtUpsilonPP;
	} else if (sname == "pp 7") {
	    func=PtUpsilonPP7000;
	} else if (sname == "pp 2.76") {
	    func=PtUpsilonPP2760;
	} else if (sname == "PbPb 2.76") {
	    func=PtUpsilonPbPb2760;
	} else if (sname == "PbPb 2.76c1") {
	    func=PtUpsilonPbPb2760c1;
	} else if (sname == "PbPb 2.76c2") {
	    func=PtUpsilonPbPb2760c2;
	} else if (sname == "PbPb 2.76c3") {
	    func=PtUpsilonPbPb2760c3;
	} else if (sname == "PbPb 2.76c4") {
	    func=PtUpsilonPbPb2760c4;
	} else if (sname == "PbPb 2.76c5") {
	    func=PtUpsilonPbPb2760c5;
	} else if (sname == "PbPb 2.76c6") {
	    func=PtUpsilonPbPb2760c6;
	} else if (sname == "PbPb 2.76c7") {
	    func=PtUpsilonPbPb2760c7;
	} else if (sname == "PbPb 2.76c8") {
	    func=PtUpsilonPbPb2760c8;
	} else if (sname == "PbPb 2.76c9") {
	    func=PtUpsilonPbPb2760c9;
	} else if (sname == "PbPb 2.76c10") {
	    func=PtUpsilonPbPb2760c10;
	} else if (sname == "PbPb 2.76c11") {
	    func=PtUpsilonPbPb2760c11;
	} else if (sname == "pp 7 poly") {
	    func=PtUpsilonPP7000;
	} else if (sname == "pp 2.76 poly") {
	    func=PtUpsilonPP2760;
	} else if (sname == "pp 8.8") {
	    func=PtUpsilonPP8800;
	} else if (sname == "pPb 8.8") {
	    func=PtUpsilonPPb8800;
	} else if (sname == "pPb 8.8c1") {
	    func=PtUpsilonPPb8800c1;
	} else if (sname == "pPb 8.8c2") {
	    func=PtUpsilonPPb8800c2;
	} else if (sname == "pPb 8.8c3") {
	    func=PtUpsilonPPb8800c3;
	} else if (sname == "pPb 8.8c4") {
	    func=PtUpsilonPPb8800c4;
	} else if (sname == "Pbp 8.8") {
	    func=PtUpsilonPbP8800;
	} else if (sname == "Pbp 8.8c1") {
	    func=PtUpsilonPbP8800c1;
	} else if (sname == "Pbp 8.8c2") {
	    func=PtUpsilonPbP8800c2;
	} else if (sname == "Pbp 8.8c3") {
	    func=PtUpsilonPbP8800c3;
	} else if (sname == "Pbp 8.8c4") {
	    func=PtUpsilonPbP8800c4;
	} else if (sname == "CDF scaled") {
	    func=PtUpsilonCDFscaled;
	} else if (sname == "CDF pp") {
	    func=PtUpsilonCDFscaledPP;
	} else if (sname == "CDF pp 10") {
	    func=PtUpsilonCDFscaledPP10;
	} else if (sname == "CDF pp 8.8") {
	    func=PtUpsilonCDFscaledPP9;
	} else if (sname == "CDF pp 7") {
	    func=PtUpsilonCDFscaledPP7;
	} else if (sname == "CDF pp 3.94") {
	    func=PtUpsilonCDFscaledPP4;
	} else if (sname == "CDF pPb 8.8") {
	    func=PtUpsilonCDFscaledPPb9;
	} else if (sname == "CDF Pbp 8.8") {
	    func=PtUpsilonCDFscaledPbP9;
	} else if (sname == "CDF PbPb 3.94") {
	    func=PtUpsilonCDFscaledPbPb4;
	} else if (sname == "Flat") {
	    func=PtUpsilonFlat;
	} else {
	    func=PtUpsilon;
	}
	break;  
    case kCharm:
	if (sname == "F0M0S0 pp") {
	    func=PtCharmF0M0S0PP;
	} else if (sname == "F1M0S0 pp") {
	    func=PtCharmF1M0S0PP;
	} else if (sname == "F2M0S0 pp") {
	    func=PtCharmF2M0S0PP;
	} else if (sname == "F0M1S0 pp") {
	    func=PtCharmF0M1S0PP;
	} else if (sname == "F0M2S0 pp") {
	    func=PtCharmF0M2S0PP;
	} else if (sname == "F0M0S1 pp") {
	    func=PtCharmF0M0S1PP;
	} else if (sname == "F0M0S2 pp") {
	    func=PtCharmF0M0S2PP;
	} else if (sname == "F0M0S3 pp") {
	    func=PtCharmF0M0S3PP;
	} else if (sname == "F0M0S4 pp") {
	    func=PtCharmF0M0S4PP;
	} else if (sname == "F0M0S5 pp") {
	    func=PtCharmF0M0S5PP;
	} else if (sname == "F0M0S6 pp") {
	    func=PtCharmF0M0S6PP;
	} else if (sname == "central") {
	    func=PtCharmCentral;
	} else {
	    func=PtCharm;
	}
	break;
    case kBeauty:
	if (sname == "F0M0S0 pp") {
	    func=PtBeautyF0M0S0PP;
	} else if (sname == "F1M0S0 pp") {
	    func=PtBeautyF1M0S0PP;
	} else if (sname == "F2M0S0 pp") {
	    func=PtBeautyF2M0S0PP;
	} else if (sname == "F0M1S0 pp") {
	    func=PtBeautyF0M1S0PP;
	} else if (sname == "F0M2S0 pp") {
	    func=PtBeautyF0M2S0PP;
	} else if (sname == "F0M0S1 pp") {
	    func=PtBeautyF0M0S1PP;
	} else if (sname == "F0M0S2 pp") {
	    func=PtBeautyF0M0S2PP;
	} else if (sname == "F0M0S3 pp") {
	    func=PtBeautyF0M0S3PP;
	} else if (sname == "F0M0S4 pp") {
	    func=PtBeautyF0M0S4PP;
	} else if (sname == "F0M0S5 pp") {
	    func=PtBeautyF0M0S5PP;
	} else if (sname == "F0M0S6 pp") {
	    func=PtBeautyF0M0S6PP;
	} else if (sname == "central") {
	    func=PtBeautyCentral;
	} else {
	    func=PtBeauty;
	}
	break;
    case kPion:
        if (sname == "2010 Pos PP") {
            func=PtPionPos2010PP;
        } else if (sname == "2010 Neg PP") {
            func=PtPionNeg2010PP;
        } else {
	    func=PtPion;
        }
	break;
    case kKaon:
        if (sname == "2010 Pos PP") {
            func=PtKaonPos2010PP;
        } else if (sname == "2010 Neg PP") {
            func=PtKaonNeg2010PP;
        } else {
	    func=PtKaon;
        }
	break;
    case kChic0:
	func=PtChic0;
	break;
    case kChic:
	func=PtChic;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetPt> unknown parametrisation\n");
    }
    return func;
}

GenFunc AliGenMUONlib::GetY(Int_t param, const char* tname) const
{
  //    
  // Return pointer to y- parameterisation
  //
    TString sname = TString(tname);
    GenFunc func;
    switch (param) 
    {
    case kPhi:
	func=YPhi;
	break;
    case kEta:
	func=YEta;
	break;
    case kOmega:
	func=YOmega;
	break;
    case kJpsiFamily:
    case kPsiP:
    case kChic1:
    case kChic2:
    case kJpsi:
	if (sname == "Vogt" || sname == "Vogt PbPb") {
	    func=YJpsiPbPb;
	} else if (sname == "Vogt pp"){
	    func=YJpsiPP;
	} else if (sname == "pp 7") {
	    func=YJpsiPP7000;
	} else if (sname == "pp 2.76") {
	    func=YJpsiPP2760;
	} else if (sname == "PbPb 2.76") {
	    func=YJpsiPbPb2760;
	} else if (sname == "PbPb 2.76c1") {
	    func=YJpsiPbPb2760c1;
	} else if (sname == "PbPb 2.76c2") {
	    func=YJpsiPbPb2760c2;
	} else if (sname == "PbPb 2.76c3") {
	    func=YJpsiPbPb2760c3;
	} else if (sname == "PbPb 2.76c4") {
	    func=YJpsiPbPb2760c4;
	} else if (sname == "PbPb 2.76c5") {
	    func=YJpsiPbPb2760c5;
	} else if (sname == "PbPb 2.76c6") {
	    func=YJpsiPbPb2760c6;
	} else if (sname == "PbPb 2.76c7") {
	    func=YJpsiPbPb2760c7;
	} else if (sname == "PbPb 2.76c8") {
	    func=YJpsiPbPb2760c8;
	} else if (sname == "PbPb 2.76c9") {
	    func=YJpsiPbPb2760c9;
	} else if (sname == "PbPb 2.76c10") {
	    func=YJpsiPbPb2760c10;
	} else if (sname == "PbPb 2.76c11") {
	    func=YJpsiPbPb2760c11;
	} else if (sname == "pp 7 poly") {
	    func=YJpsiPPpoly7000;
	} else if (sname == "pp 2.76 poly") {
	    func=YJpsiPPpoly2760;
	} else if (sname == "pp 8.8") {
	    func=YJpsiPP8800;
	} else if (sname == "pPb 8.8") {
	    func=YJpsiPPb8800;
	} else if (sname == "pPb 8.8c1") {
	    func=YJpsiPPb8800c1;
	} else if (sname == "pPb 8.8c2") {
	    func=YJpsiPPb8800c2;
	} else if (sname == "pPb 8.8c3") {
	    func=YJpsiPPb8800c3;
	} else if (sname == "pPb 8.8c4") {
	    func=YJpsiPPb8800c4;
	} else if (sname == "Pbp 8.8") {
	    func=YJpsiPbP8800;
	} else if (sname == "Pbp 8.8c1") {
	    func=YJpsiPbP8800c1;
	} else if (sname == "Pbp 8.8c2") {
	    func=YJpsiPbP8800c2;
	} else if (sname == "Pbp 8.8c3") {
	    func=YJpsiPbP8800c3;
	} else if (sname == "Pbp 8.8c4") {
	    func=YJpsiPbP8800c4;
	} else if (sname == "CDF scaled") {
	    func=YJpsiCDFscaled;
	} else if (sname == "CDF pp") {
	    func=YJpsiCDFscaledPP;
	} else if (sname == "CDF pp 10") {
	    func=YJpsiCDFscaledPP10;
	} else if (sname == "CDF pp 8.8") {
	    func=YJpsiCDFscaledPP9;
	} else if (sname == "CDF pp 7" || sname == "CDF pp 7 flat pt") {
	    func=YJpsiCDFscaledPP7;
	} else if (sname == "CDF pp 3.94") {
	    func=YJpsiCDFscaledPP4;
	} else if (sname == "CDF pp 2.76") {
	    func=YJpsiCDFscaledPP3;
	} else if (sname == "CDF pp 1.9") {
	    func=YJpsiCDFscaledPP2;
	} else if (sname == "CDF pPb 8.8") {
	    func=YJpsiCDFscaledPPb9;
	} else if (sname == "CDF Pbp 8.8") {
	    func=YJpsiCDFscaledPbP9;
	} else if (sname == "CDF PbPb 3.94") {
	    func=YJpsiCDFscaledPbPb4;
	} else if (sname == "Flat" || sname == "CDF pp 7 flat y") {
	    func=YJpsiFlat;
	} else {
	    func=YJpsi;
	}
	break;
    case kJpsiFromB:
	func = YJpsiBPbPb;
	break;
    case kUpsilonFamily:
    case kUpsilonP:
    case kUpsilonPP:
    case kUpsilon:
	if (sname == "Vogt" || sname == "Vogt PbPb") {
	    func=YUpsilonPbPb;
	} else if (sname == "Vogt pp") {
	    func = YUpsilonPP;
	} else if (sname == "pp 7") {
	    func=YUpsilonPP7000;
	} else if (sname == "pp 2.76") {
	    func=YUpsilonPP2760;
	} else if (sname == "PbPb 2.76") {
	    func=YUpsilonPbPb2760;
	} else if (sname == "PbPb 2.76c1") {
	    func=YUpsilonPbPb2760c1;
	} else if (sname == "PbPb 2.76c2") {
	    func=YUpsilonPbPb2760c2;
	} else if (sname == "PbPb 2.76c3") {
	    func=YUpsilonPbPb2760c3;
	} else if (sname == "PbPb 2.76c4") {
	    func=YUpsilonPbPb2760c4;
	} else if (sname == "PbPb 2.76c5") {
	    func=YUpsilonPbPb2760c5;
	} else if (sname == "PbPb 2.76c6") {
	    func=YUpsilonPbPb2760c6;
	} else if (sname == "PbPb 2.76c7") {
	    func=YUpsilonPbPb2760c7;
	} else if (sname == "PbPb 2.76c8") {
	    func=YUpsilonPbPb2760c8;
	} else if (sname == "PbPb 2.76c9") {
	    func=YUpsilonPbPb2760c9;
	} else if (sname == "PbPb 2.76c10") {
	    func=YUpsilonPbPb2760c10;
	} else if (sname == "PbPb 2.76c11") {
	    func=YUpsilonPbPb2760c11;
	} else if (sname == "pp 7 poly") {
	    func=YUpsilonPPpoly7000;
	} else if (sname == "pp 2.76 poly") {
	    func=YUpsilonPPpoly2760;
	} else if (sname == "pp 8.8") {
	    func=YUpsilonPP8800;
	} else if (sname == "pPb 8.8") {
	    func=YUpsilonPPb8800;
	} else if (sname == "pPb 8.8c1") {
	    func=YUpsilonPPb8800c1;
	} else if (sname == "pPb 8.8c2") {
	    func=YUpsilonPPb8800c2;
	} else if (sname == "pPb 8.8c3") {
	    func=YUpsilonPPb8800c3;
	} else if (sname == "pPb 8.8c4") {
	    func=YUpsilonPPb8800c4;
	} else if (sname == "Pbp 8.8") {
	    func=YUpsilonPbP8800;
	} else if (sname == "Pbp 8.8c1") {
	    func=YUpsilonPbP8800c1;
	} else if (sname == "Pbp 8.8c2") {
	    func=YUpsilonPbP8800c2;
	} else if (sname == "Pbp 8.8c3") {
	    func=YUpsilonPbP8800c3;
	} else if (sname == "Pbp 8.8c4") {
	    func=YUpsilonPbP8800c4;
	} else if (sname == "CDF scaled") {
	    func=YUpsilonCDFscaled;
	} else if (sname == "CDF pp") {
	    func=YUpsilonCDFscaledPP;
	} else if (sname == "CDF pp 10") {
	    func=YUpsilonCDFscaledPP10;
	} else if (sname == "CDF pp 8.8") {
	    func=YUpsilonCDFscaledPP9;
	} else if (sname == "CDF pp 7") {
	    func=YUpsilonCDFscaledPP7;
	} else if (sname == "CDF pp 3.94") {
	    func=YUpsilonCDFscaledPP4;
	} else if (sname == "CDF pPb 8.8") {
	    func=YUpsilonCDFscaledPPb9;
	} else if (sname == "CDF Pbp 8.8") {
	    func=YUpsilonCDFscaledPbP9;
	} else if (sname == "CDF PbPb 3.94") {
	    func=YUpsilonCDFscaledPbPb4;
	} else if (sname == "Flat") {
	    func=YUpsilonFlat;
	} else {
	    func=YUpsilon;
	}
	break;
    case kCharm:
	if (sname == "F0M0S0 pp") {
	    func=YCharmF0M0S0PP;
	} else if (sname == "F1M0S0 pp") {
	    func=YCharmF1M0S0PP;
	} else if (sname == "F2M0S0 pp") {
	    func=YCharmF2M0S0PP;
	} else if (sname == "F0M1S0 pp") {
	    func=YCharmF0M1S0PP;
	} else if (sname == "F0M2S0 pp") {
	    func=YCharmF0M2S0PP;
	} else if (sname == "F0M0S1 pp") {
	    func=YCharmF0M0S1PP;
	} else if (sname == "F0M0S2 pp") {
	    func=YCharmF0M0S2PP;
	} else if (sname == "F0M0S3 pp") {
	    func=YCharmF0M0S3PP;
	} else if (sname == "F0M0S4 pp") {
	    func=YCharmF0M0S4PP;
	} else if (sname == "F0M0S5 pp") {
	    func=YCharmF0M0S5PP;
	} else if (sname == "F0M0S6 pp") {
	    func=YCharmF0M0S6PP;
	} else {
	    func=YCharm;
	}
	break;
    case kBeauty:
	if (sname == "F0M0S0 pp") {
	    func=YBeautyF0M0S0PP;
	} else if (sname == "F1M0S0 pp") {
	    func=YBeautyF1M0S0PP;
	} else if (sname == "F2M0S0 pp") {
	    func=YBeautyF2M0S0PP;
	} else if (sname == "F0M1S0 pp") {
	    func=YBeautyF0M1S0PP;
	} else if (sname == "F0M2S0 pp") {
	    func=YBeautyF0M2S0PP;
	} else if (sname == "F0M0S1 pp") {
	    func=YBeautyF0M0S1PP;
	} else if (sname == "F0M0S2 pp") {
	    func=YBeautyF0M0S2PP;
	} else if (sname == "F0M0S3 pp") {
	    func=YBeautyF0M0S3PP;
	} else if (sname == "F0M0S4 pp") {
	    func=YBeautyF0M0S4PP;
	} else if (sname == "F0M0S5 pp") {
	    func=YBeautyF0M0S5PP;
	} else if (sname == "F0M0S6 pp") {
	    func=YBeautyF0M0S6PP;
	} else {
	    func=YBeauty;
	}
	break;
    case kPion:
        if (sname == "2010 Pos PP") {
            func=YKaonPion2010PP;
        } else if (sname == "2010 Neg PP") {
            func=YKaonPion2010PP;
        } else {
	    func=YPion;
        }
	break;
    case kKaon:
        if (sname == "2010 Pos PP") {
            func=YKaonPion2010PP;
        } else if (sname == "2010 Neg PP") {
            func=YKaonPion2010PP;
        } else {
	    func=YKaon;
        }
	break;
    case kChic0:
	func=YChic0;
	break;
    case kChic:
	func=YChic;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetY> unknown parametrisation\n");
    }
    return func;
}

//
//                    Chi
//
//
//                pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtChic0( const Double_t *px, const Double_t */*dummy*/)
{
// Chi_c1 pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtChic1( const Double_t *px, const Double_t */*dummy*/)
{
// Chi_c1 pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtChic2( const Double_t *px, const Double_t */*dummy*/)
{
// Chi_c2 pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
Double_t AliGenMUONlib::PtChic( const Double_t *px, const Double_t */*dummy*/)
{
// Chi_c family pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

//
//               y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YChic0(const Double_t *py, const Double_t */*dummy*/)
{
// Chi-1c y
  const Double_t ky0 = 4.;
 const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}

Double_t AliGenMUONlib::YChic1(const Double_t *py, const Double_t */*dummy*/)
{
// Chi-1c y
  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}

Double_t AliGenMUONlib::YChic2(const Double_t *py, const Double_t */*dummy*/)
{
// Chi-2c y
  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}

Double_t AliGenMUONlib::YChic(const Double_t *py, const Double_t */*dummy*/)
{
// Chi_c family y
  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}

//                 particle composition
//
Int_t AliGenMUONlib::IpChic0(TRandom *)
{
// Chi composition
    return 10441;
}
//
Int_t AliGenMUONlib::IpChic1(TRandom *)
{
// Chi composition
    return 20443;
}
Int_t AliGenMUONlib::IpChic2(TRandom *)
{
// Chi_c2 prime composition
    return 445;
}
Int_t AliGenMUONlib::IpChic(TRandom *)
{
// Chi composition
  Int_t ip;
  Float_t r = gRandom->Rndm();
  if (r < 0.001) {
    ip = 10441;
  } else if( r < 0.377 ) {
    ip = 20443;
  } else {
    ip = 445;
  }
  return ip;
}


//_____________________________________________________________

typedef Int_t (*GenFuncIp) (TRandom *);
GenFuncIp AliGenMUONlib::GetIp(Int_t param,  const char* tname) const
{
// Return pointer to particle type parameterisation
    TString sname = TString(tname);
    GenFuncIp func;
    switch (param) 
    {
    case kPhi:
	func=IpPhi;
	break;
    case kEta:
	func=IpEta;
	break;
    case kOmega:
	func=IpOmega;
	break;
    case kJpsiFamily:
      	func=IpJpsiFamily;
	break;
    case kPsiP:
      	func=IpPsiP;
	break;
    case kJpsi:
    case kJpsiFromB:
	func=IpJpsi;
	break;
    case kUpsilon:
	func=IpUpsilon;
	break;
    case kUpsilonFamily:
      func=IpUpsilonFamily;
      break;
    case kUpsilonP:
	func=IpUpsilonP;
	break;
    case kUpsilonPP:
	func=IpUpsilonPP;
	break;
    case kCharm:
	func=IpCharm;
	break;
    case kBeauty:
	func=IpBeauty;
	break;
    case kPion:
        if (sname == "2010 Pos PP") {
            func=IpPionPos;
        } else if (sname == "2010 Neg PP") {
            func=IpPionNeg;
        } else {
	    func=IpPion;
        }
	break;
    case kKaon:
        if (sname == "2010 Pos PP") {
            func=IpKaonPos;
        } else if (sname == "2010 Neg PP") {
            func=IpKaonNeg;
        } else {
	    func=IpKaon;
        }
	break;
    case kChic0:
	func=IpChic0;
	break;
    case kChic1:
	func=IpChic1;
	break;
    case kChic2:
	func=IpChic2;
	break;
    case kChic:
        func=IpChic;
        break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetIp> unknown parametrisation\n");
    }
    return func;
}



Float_t AliGenMUONlib::Interpolate(Float_t x, Float_t* y, Float_t x0, 
				   Float_t dx,
				   Int_t n, Int_t no)
{
//
// Neville's alorithm for interpolation
//
// x:  x-value
// y:  Input array
// x0: minimum x 
// dx: step size
//  n: number of data points
// no: order of polynom 
//
    Float_t*  c = new Float_t[n];
    Float_t*  d = new Float_t[n];
    Int_t m, i;
    for (i = 0; i < n; i++) {
	c[i] = y[i];
	d[i] = y[i];
    }
    
    Int_t   ns  = int((x - x0)/dx);
    
    Float_t y1  = y[ns];
    ns--;    
    for (m = 0; m < no; m++) {
	for (i = 0; i < n-m; i++) {	
	    Float_t ho = x0 + Float_t(i) * dx - x;
	    Float_t hp = x0 + Float_t(i+m+1) * dx - x;
	    Float_t w  = c[i+1] - d[i];
	    Float_t den = ho-hp;
	    den = w/den;
	    d[i] = hp * den;
	    c[i] = ho * den;
	}
	Float_t dy;
	
	if (2*ns < (n-m-1)) {
	    dy  = c[ns+1];
	} else {
	    dy  = d[ns--];
	}
	y1 += dy;}
    delete[] c;
    delete[] d;

    return y1;
}

//=============================================================================
Double_t AliGenMUONlib::PtPionPos2010PP(const Double_t *px, const Double_t* /*dummy*/)
{
// Pos pion
  const Double_t par[3] = {2.27501, 0.116141, 5.59591};
  Double_t pt = px[0];
  Double_t m0 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t mt = TMath::Sqrt(m0*m0 + pt*pt);
  Double_t nc = par[1]*par[2];
  Double_t t1 = (par[2]-1.)/nc/(nc/(par[2]-2.)+m0);
  Double_t t2 = TMath::Power(1.+(mt-m0)/nc, -1.*par[2]);
  Double_t fn = par[0] * pt * t1 * t2;
  return fn;
}

//=============================================================================
Double_t AliGenMUONlib::PtPionNeg2010PP(const Double_t *px, const Double_t* /*dummy*/)
{
// Neg pion
  const Double_t par[3] = {2.25188, 0.12176, 5.91166};
  Double_t pt = px[0];
  Double_t m0 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t mt = TMath::Sqrt(m0*m0 + pt*pt);
  Double_t nc = par[1]*par[2];
  Double_t t1 = (par[2]-1.)/nc/(nc/(par[2]-2.)+m0);
  Double_t t2 = TMath::Power(1.+(mt-m0)/nc, -1.*par[2]);
  Double_t fn = par[0] * pt * t1 * t2;
  return fn;
}

//=============================================================================
Double_t AliGenMUONlib::PtKaonPos2010PP(const Double_t *px, const Double_t* /*dummy*/)
{
// Pos kaons
  const Double_t par[3] = {0.279386, 0.195466, 6.59587};
  Double_t pt = px[0];
  Double_t m0 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t mt = TMath::Sqrt(m0*m0 + pt*pt);
  Double_t nc = par[1]*par[2];
  Double_t t1 = (par[2]-1.)/nc/(nc/(par[2]-2.)+m0);
  Double_t t2 = TMath::Power(1.+(mt-m0)/nc, -1.*par[2]);
  Double_t fn = par[0] * pt * t1 * t2;
  return fn;
}

//=============================================================================
Double_t AliGenMUONlib::PtKaonNeg2010PP(const Double_t *px, const Double_t* /*dummy*/)
{
// Neg kaons
  const Double_t par[3] = {0.278927, 0.189049, 6.43006};
  Double_t pt = px[0];
  Double_t m0 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t mt = TMath::Sqrt(m0*m0 + pt*pt);
  Double_t nc = par[1]*par[2];
  Double_t t1 = (par[2]-1.)/nc/(nc/(par[2]-2.)+m0);
  Double_t t2 = TMath::Power(1.+(mt-m0)/nc, -1.*par[2]);
  Double_t fn = par[0] * pt * t1 * t2;
  return fn;
}

//=============================================================================
Double_t AliGenMUONlib::YKaonPion2010PP(const Double_t *px, const Double_t* /*dummy*/)
{
// pions and kaons
  Double_t y = px[0];
  Double_t sigma = 2.35;
  Double_t kernal = y/2./sigma;
  Double_t fxn = TMath::Exp(-1.*kernal*kernal);
  return fxn;
}

//=============================================================================
Int_t AliGenMUONlib::IpPionPos(TRandom *)
{
// Pos pions
    return 211;
}

//=============================================================================
Int_t AliGenMUONlib::IpPionNeg(TRandom *)
{
// Neg pions
    return -211;
}

//=============================================================================
Int_t AliGenMUONlib::IpKaonPos(TRandom *)
{
// pos Kaons
    return 321;
}

//=============================================================================
Int_t AliGenMUONlib::IpKaonNeg(TRandom *)
{
// neg Kaons 
    return -321;
}
