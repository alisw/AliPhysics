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

//-------------------------------------------------------------------------
// author: Sergey Kiselev, ITEP, Moscow
// e-mail: Sergey.Kiselev@cern.ch
// tel.: 007 495 129 95 45
//-------------------------------------------------------------------------
// Generator of direct thermal photons for the reaction A+B, sqrt(S)
// main assumptions:
// 1+1 Bjorken scaling hydrodinamics.
// 1st order phase transition
// QGP + Mixed (QGP+HHG) + HHG (Hot Hadron Gas) phases, 
// an ideal massless parton gas and ideal massless HHG 
// see 
// F.D.Steffen, nucl-th/9909035
// F.D.Steffen and M.H.Thoma, Phys.Lett. B510, 98 (2001)
// T.Peitzmann and M.H.Thoma, Phys.Rep., 364, 175 (2002) 
//
// photon rates for QGP: Phys.Rep., 364, 175 (2002), section 2.1.1
//
// photon rates for HHG
// prates for i rho --> pi gamma, pi pi --> rho gamma and rho --> pi pi gamma:
// Song and Fai, Phys.Rev. C58, 1689 (1998)
// rates for omega --> pi gamma: Phys.Rev. D44, 2774 (1991)
//
// input parameters:
//       fAProjectile, fATarget - number of nucleons in a nucleus A and B
//       fMinImpactParam - minimal impct parameter, fm
//       fMaxImpactParam - maximal impct parameter, fm
//       fEnergyCMS - sqrt(S) per nucleon pair, AGeV
//       fTau0 - initial proper time, fm/c
//       fT0 - initial temperature, GeV
//       fTc - critical temperature, GeV
//       fTf - freeze-out temperature, GeV
//       fGhhg - effective number of degrees of freedom in HHG
//
//       fYMin - minimal rapidity of photons 
//       fYMax - maximal rapidity of photons
//              in [fYMin,fYMax] uniform distribution of gamma is assumed
//       fPtMin - minimal p_t value of gamma, GeV/c
//       fPtMax - maximal p_t value of gamma, GeV/c
//-------------------------------------------------------------------------
// comparison with SPS and RHIC data, prediction for LHC can be found in
// arXiv:0811.2634 [nucl-th]
//-------------------------------------------------------------------------

//Begin_Html
/*
<img src="picts/AliGeneratorClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TArrayF.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>

#include "AliConst.h"
#include "AliGenEventHeader.h"
#include "AliGenThermalPhotons.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliGenThermalPhotons)

// -----------------------------------------------------------------------------------------------------
static Double_t rateQGP(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - tau (fm), proper time
// x[1] - yprime, space rapidity
// par[0] - p_T (GeV), photon transverse momentum 
// par[1] - y, photon rapidity in the c.m.s. A+A
// par[2] - tau0 (fm), initial proper time 
// par[3] - T_0 (GeV), initial temperature
// par[4] - T_c (GeV), critical temperature
// par[5] - iProcQGP, process number, 0, 1, 2
//
// output:
// tau EdR/d^3p = tau EdN/d^4xd^3p (fm fm^{-4}GeV^{-2}) 
//---------------------------------------------------
  Double_t tau=x[0],yprime=x[1];
  Double_t pT=par[0],y=par[1],tau0=par[2],t0=par[3],tc=par[4];
  Int_t iProcQGP=Int_t(par[5]), nFl=3;

  Double_t e=pT*TMath::CosH(yprime-y),t=t0*TMath::Power(tau0/tau,1./3.);

  const Double_t alpha=1./137.;
// factor to convert from GeV^2 to fm^{-4}GeV^{-2}: (1/hc)**4=(1/0.197)**4=659.921
  const Double_t factor=659.921;
  Double_t alphaS=3.*TMath::TwoPi()/((33.-2.*nFl)*TMath::Log(8.*t/tc));
  const Double_t abc[3]={0.0338, 0.0281, 0.0135} ; // a, b, c for nFf=3
  Double_t rate=1.;

  switch (iProcQGP) {

    case 0:
// 1-loop
      rate=factor*abc[iProcQGP]*alpha*alphaS*TMath::Exp(-e/t)*t*t*TMath::Log(0.2317*e/alphaS/t);
    break ;

    case 1:
// 2-loop: bremsstrahlung
      rate=factor*abc[iProcQGP]*alpha*alphaS*TMath::Exp(-e/t)*t*t;
    break ;

    case 2:
// 2-loop: annihilation with scattering
      rate=factor*abc[iProcQGP]*alpha*alphaS*TMath::Exp(-e/t)*t*e;
    break ;
  
    default:
      printf("NO iProcQGP=%i \n",iProcQGP);
  }

  return tau*rate;
}   

// -----------------------------------------------------------------------------------------------------
static Double_t fromQGP(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - p_T (GeV), photon p_T
// par[0] - tau0 (fm), initial proper time 
// par[1] - T_0 (GeV), initial temperature
// par[2] - tauCQGP (fm), end of QGP
// par[3] - yNucl, rapidity of projectile nucleus
// par[4] - T_c (GeV), critical temperature
// par[5] - y, photon rapidity
// par[6] - iProcQGP, process number
//
// output:
// d^{2}N/(dp_t dy) (1/GeV)
//---------------------------------------------------
  Double_t pT=x[0];
  Double_t tau0=par[0],t0=par[1],tauCQGP=par[2],yNucl=par[3],tc=par[4],y=par[5];
  Int_t iProcQGP=Int_t(par[6]);

  TF2 frateQGP("frateQGP",&rateQGP,tau0,tauCQGP,-yNucl,yNucl,6);
  frateQGP.SetParameters(pT,y,tau0,t0,tc,iProcQGP);
  frateQGP.SetParNames("transverse momentum","rapidity","initial time","initial temperature","critical temperature","process number");
  return TMath::TwoPi()*pT*frateQGP.Integral(tau0,tauCQGP,-yNucl,yNucl,1e-6);
}   
         
// -----------------------------------------------------------------------------------------------------
static Double_t rateMixQ(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - yprime, space rapidity
// par[0] - p_T (GeV), photon transverse momentum 
// par[1] - y, photon rapidity in the c.m.s. A+A
// par[2] - T_c (GeV), critical temperature
// par[3] - iProcQGP, process number, 0, 1, 2
//
// output:
// EdR/d^3p = EdN/d^4xd^3p (fm fm^{-4}GeV^{-2}) 
//---------------------------------------------------
  Double_t yprime=x[0];
  Double_t pT=par[0],y=par[1],tc=par[2];
  Int_t iProcQGP=Int_t(par[3]),nFl=3;

  Double_t e=pT*TMath::CosH(yprime-y),t=tc;

  const Double_t alpha=1./137.;
// factor to convert from GeV^2 to fm^{-4}GeV^{-2}: (1/hc)**4=(1/0.197)**4=659.921
  const Double_t factor=659.921;
  Double_t alphaS=3.*TMath::TwoPi()/((33.-2.*nFl)*TMath::Log(8.*t/tc));
  const Double_t abc[3]={0.0338, 0.0281, 0.0135}; // a, b, c for nF=3
  Double_t rate=1.;

  switch (iProcQGP) {

    case 0:
// 1-loop
      rate=factor*abc[iProcQGP]*alpha*alphaS*TMath::Exp(-e/t)*t*t*TMath::Log(0.2317*e/alphaS/t);
    break ;

    case 1:
// 2-loop: bremsstrahlung
      rate=factor*abc[iProcQGP]*alpha*alphaS*TMath::Exp(-e/t)*t*t;
    break ;

    case 2:
// 2-loop: annihilation with scattering
      rate=factor*abc[iProcQGP]*alpha*alphaS*TMath::Exp(-e/t)*t*e;
    break ;
  
    default:
      printf("NO iProcQGP=%i \n",iProcQGP);
  }

  return rate;
}   

// -----------------------------------------------------------------------------------------------------
static Double_t fromMixQ(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - p_T (GeV), photon p_T
// par[0] - lamQGP
// par[1] - yNucl, rapidity of projectile nucleus
// par[2] - T_c (GeV), critical temperature
// par[3] - y, photon rapidity
// par[4] - iProcQGP, process number
//
// output:
// d^{2}N/(dp_t dy) (1/GeV)
//---------------------------------------------------
  Double_t pT=x[0];
  Double_t lamQGP=par[0],yNucl=par[1],tc=par[2],y=par[3];
  Int_t iProcQGP=Int_t(par[4]);

  TF1 frateMixQ("frateMixQ",&rateMixQ,-yNucl,yNucl,4);
  frateMixQ.SetParameters(pT,y,tc,iProcQGP);
  frateMixQ.SetParNames("transverse momentum","rapidity","critical temperature","process number");
  return TMath::TwoPi()*pT*lamQGP*frateMixQ.Integral(-yNucl,yNucl);
}   
         
// -----------------------------------------------------------------------------------------------------
static Double_t rateMixH(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - yprime, space rapidity
// par[0] - p_T (GeV), photon transverse momentum 
// par[1] - y, photon rapidity in the c.m.s. A+A
// par[2] - T_c (GeV), critical temperature
// par[3] - iProcHHG, process number
//
// output:
// EdR/d^3p = EdN/d^4xd^3p (fm^{-4}GeV^{-2}) 
//---------------------------------------------------
  Double_t yprime=x[0];
  Double_t pT=par[0],y=par[1],tc=par[2];
  Int_t iProcHHG=Int_t(par[3]);

  Double_t e=pT*TMath::CosH(yprime-y),t=tc;
  const Double_t mPi=0.135;
  Double_t xx=t/mPi,yy=e/mPi;
  Double_t f,rate=1.,emin,factor;
  const Double_t mOm=0.783,width=0.00844,br=0.085,e0=0.379,pi=TMath::Pi();
  const Double_t hc=0.197,hc4=hc*hc*hc*hc; // GeV*fm

  switch (iProcHHG) {

    case 0:
// pi rho --> pi gamma
      f=-2.447+0.796*xx+(0.0338+0.0528*xx)*yy+(-21.447+8.2179*xx)/yy+(1.52436-0.38562*xx)/(yy*yy);
      rate=t*t*TMath::Exp(-e/t+f);
    break ;

    case 1:
// pi pi --> rho gamma
      f=-12.055+4.387*xx+(0.3755+0.00826*xx)*yy+(-0.00777+0.000279*xx)*yy*yy+(5.7869-1.0258*xx)/yy+(-1.979+0.58*xx)/(yy*yy);
      rate=t*t*TMath::Exp(-e/t+f);
    break ;

    case 2:
// rho --> pi pi gamma
      f=-6.295+1.6459*xx+(-0.4015+0.089*xx)*yy+(-0.954+2.05777*xx)/yy;
      rate=t*t*TMath::Exp(-e/t+f);
    break ;
  
    case 3:
// omega --> pi gamma
      emin=mOm*(e*e+e0*e0)/(2.*e*e0);
      factor=(3.*mOm*width*br)/(16.*pi*pi*pi*e0);
      rate=factor*t*(emin+t)*TMath::Exp(-emin/t)/e/hc4;
    break ;
  
    default:
      printf("NO iProcHHG=%i \n",iProcHHG);
  }

  return rate;
}   

// -----------------------------------------------------------------------------------------------------
static Double_t fromMixH(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - p_T (GeV), photon p_T
// par[0] - lamHHG
// par[1] - yNucl, rapidity of projectile nucleus
// par[2] - T_c (GeV), critical temperature
// par[3] - y, photon rapidity
// par[4] - iProcHHG, process number
//
// output:
// d^{2}N/(dp_t dy) (1/GeV)
//---------------------------------------------------
  Double_t pT=x[0];
  Double_t lamHHG=par[0],yNucl=par[1],tc=par[2],y=par[3];
  Int_t iProcHHG=Int_t(par[4]);

  TF1 frateMixH("frateMixH",&rateMixH,-yNucl,yNucl,4);
  frateMixH.SetParameters(pT,y,tc,iProcHHG);
  frateMixH.SetParNames("transverse momentum","rapidity","critical temperature","process number");
  return TMath::TwoPi()*pT*lamHHG*frateMixH.Integral(-yNucl,yNucl);
}   
         
// -----------------------------------------------------------------------------------------------------
static Double_t rateHHG(const Double_t *x, const Double_t *par) {
//---------------------------------------------------
// input:
// x[0] - tau (fm), proper time
// x[1] - yprime, space rapidity
// par[0] - p_T (GeV), photon transverse momentum 
// par[1] - y, photon rapidity in the c.m.s. A+A
// par[2] - tauCHHG (fm), start of HHG 
// par[3] - T_c (GeV), critical temperature
// par[4] - iProcHHG, process number
//
// output:
// EdR/d^3p = EdN/d^4xd^3p (fm^{-4}GeV^{-2}) 
//---------------------------------------------------
  Double_t tau=x[0],yprime=x[1];
  Double_t pT=par[0],y=par[1],tauCHHG=par[2],tc=par[3];
  Int_t iProcHHG=Int_t(par[4]);

  Double_t e=pT*TMath::CosH(yprime-y),t=tc*TMath::Power(tauCHHG/tau,1./3.);
  const Double_t mPi=0.135;
  Double_t xx=t/mPi,yy=e/mPi;
  Double_t f,rate=1.,emin,factor;
  const Double_t mOm=0.783,width=0.00844,br=0.085,e0=0.379,pi=TMath::Pi();
  const Double_t hc=0.197,hc4=hc*hc*hc*hc; // GeV*fm

  switch (iProcHHG) {

    case 0:
// pi rho --> pi gamma
      f=-2.447+0.796*xx+(0.0338+0.0528*xx)*yy+(-21.447+8.2179*xx)/yy+(1.52436-0.38562*xx)/(yy*yy);
      rate=t*t*TMath::Exp(-e/t+f);
    break ;

    case 1:
// pi pi --> rho gamma
      f=-12.055+4.387*xx+(0.3755+0.00826*xx)*yy+(-0.00777+0.000279*xx)*yy*yy+(5.7869-1.0258*xx)/yy+(-1.979+0.58*xx)/(yy*yy);
      rate=t*t*TMath::Exp(-e/t+f);
    break ;

    case 2:
// rho --> pi pi gamma
      f=-6.295+1.6459*xx+(-0.4015+0.089*xx)*yy+(-0.954+2.05777*xx)/yy;
      rate=t*t*TMath::Exp(-e/t+f);
    break ;
  
    case 3:
// omega --> pi gamma
      emin=mOm*(e*e+e0*e0)/(2.*e*e0);
      factor=(3.*mOm*width*br)/(16.*pi*pi*pi*e0);
      rate=factor*t*(emin+t)*TMath::Exp(-emin/t)/e/hc4;
    break ;

    default:
      printf("NO iProcHHG=%i \n",iProcHHG);
  }
  return tau*rate;
}   

// -----------------------------------------------------------------------------------------------------
static Double_t fromHHG(const Double_t *x, const Double_t *par) {
// Thermal photon spectrum from Hot Hadron Gas (HHG)
//  F.D.Steffen, nucl-th/9909035
//  T.Peitzmann and M.H.Thoma, Phys.Rep., 364, 175 (2002), section 2.2.2 
//---------------------------------------------------
// input:
// x[0] - p_T (GeV), photon p_T
// par[0] - tauCHHG (fm), start of HHG 
// par[1] - tauF (fm), freeze-out proper time 
// par[2] - yNucl, rapidity of projectile nucleus
// par[3] - T_c (GeV), critical temperature
// par[4] - y, photon rapidity
// par[5] - iProcHHG, process number
//
// output:
// d^{2}N/(dp_t dy) (1/GeV)
//---------------------------------------------------
  Double_t pT=x[0];
  Double_t tauCHHG=par[0],tauF=par[1],yNucl=par[2],tc=par[3],y=par[4],iProcHHG=par[5];

  TF2 frateHHG("frateHHG",&rateHHG,tauCHHG,tauF,-yNucl,yNucl,5);
  frateHHG.SetParameters(pT,y,tauCHHG,tc,iProcHHG);
  frateHHG.SetParNames("transverse momentum","rapidity","start of HHG","criti temperature","process number");
  return TMath::TwoPi()*pT*frateHHG.Integral(tauCHHG,tauF,-yNucl,yNucl,1e-6);
}   

// -----------------------------------------------------------------------------------------------------
static Double_t fOverlapAB(const Double_t *x, const Double_t *par)
{
//-------------------------------------------------------------------------
// overlap area at the impact parameter b
// input:
// x[0] - impact parameter b < RA+RB
// par[0] - radius of A
// par[1] - radius of B
//-------------------------------------------------------------------------

  Double_t b=x[0], ra=par[0], rb=par[1];
  if(rb>ra) {ra=par[1]; rb=par[0];} // ra > rb

  if(b>(ra+rb)) {
    return 0.;
  }

  if(b<=(ra-rb)) {
    return TMath::Pi()*rb*rb;
  }

  Double_t p=0.5*(b+ra+rb), S=TMath::Sqrt(p*(p-b)*(p-ra)*(p-rb)), h=2.*S/b;
  Double_t sA=ra*ra*TMath::ASin(h/ra)-h*TMath::Sqrt(ra*ra-h*h);
  Double_t sB=rb*rb*TMath::ASin(h/rb)-h*TMath::Sqrt(rb*rb-h*h);
  if(ra>rb && b*b<ra*ra-rb*rb) sB=TMath::Pi()*rb*rb-sB;

  return sA+sB;

}

//_____________________________________________________________________________
AliGenThermalPhotons::AliGenThermalPhotons()
    :AliGenerator(-1),
        fMinImpactParam(0.),
        fMaxImpactParam(0.),
        fTau0(0.),
        fT0(0.),
        fTc(0.),
        fTf(0.),
        fGhhg(0),
        fSumPt()
{
    //
    // Default constructor
    //
    SetCutVertexZ();
    SetPtRange();
    SetYRange();
    fAProjectile = 208;
    fATarget     = 208;
    fEnergyCMS  = 5500.;
}

//_____________________________________________________________________________
AliGenThermalPhotons::AliGenThermalPhotons(Int_t npart)
    :AliGenerator(npart),
        fMinImpactParam(0.),
        fMaxImpactParam(0.),
        fTau0(0.1),
        fT0(0.650),
        fTc(0.170),
        fTf(0.100),
        fGhhg(8),
        fSumPt()
{
  // 
  // Standard constructor
  //

    fName="ThermalPhotons";
    fTitle="Direct thermal photons in 1+1 Bjorken hydrodynamics";

    SetCutVertexZ();
    SetPtRange();
    SetYRange();
    fAProjectile = 208;
    fATarget     = 208;
    fEnergyCMS  = 5500.;
}

//_____________________________________________________________________________
AliGenThermalPhotons::~AliGenThermalPhotons()
{
  //
  // Standard destructor
  //
    delete fSumPt;
}

//_____________________________________________________________________________
void AliGenThermalPhotons::Init()
{
    // Initialisation
  const Double_t step=0.1; 
  Int_t nPt=Int_t((fPtMax-fPtMin)/step);

  fSumPt = new TH1F("fSumPt","thermal #gamma from QGP",nPt,fPtMin,fPtMax);

  Double_t yRap=0.;
  const Int_t nCo=3,nFl=3; //  number of colors for QGP
  Double_t gQGP=2.*(nCo*nCo-1.)+(7./8.)*4.*nCo*nFl; //  number of degrees of freedom in QGP
  Double_t yNucl=TMath::ACosH(fEnergyCMS/2.);
  Double_t tauCQGP=TMath::Power(fT0/fTc,3.)*fTau0,tauCHHG=gQGP*tauCQGP/fGhhg,tauF=tauCHHG*TMath::Power(fTc/fTf,3.);
  Double_t lambda1=tauCQGP*(gQGP/(gQGP-fGhhg)),lambda2=-fGhhg/(gQGP-fGhhg);
  Double_t lamQGP=(tauCHHG-tauCQGP)*(lambda1+0.5*lambda2*(tauCHHG+tauCQGP)),lamHHG=0.5*(tauCHHG-tauCQGP)*(tauCHHG+tauCQGP)-lamQGP;

  Double_t pt,w;
// photons from pure QGP phase
  for(Int_t j=0; j<3; j++) {
    TF1 func("func",&fromQGP,fPtMin,fPtMax,7);
    func.SetParameters(fTau0,fT0,tauCQGP,yNucl,fTc,yRap,j);
    func.SetParNames("nuclear radius","initial time","initial temperature","end of pure QGP","rapidity of projectile nucleus","critical temperature","photon rapidity","process number");
    for(Int_t i=0; i<nPt; i++) {
      pt=fPtMin+(i+0.5)*step;
      w=func.Eval(pt);
      fSumPt->AddBinContent(i+1,w);
    }
  }

// photons from mixed QGP phase
  for(Int_t j=0; j<3; j++) {
    TF1 func("func",&fromMixQ,fPtMin,fPtMax,5);
    func.SetParameters(lamQGP,yNucl,fTc,yRap,j);
    func.SetParNames("lamQGP","rapidity of projectile nucleus","critical temperature","photon rapidity","process number");
    for(Int_t i=0; i<nPt; i++) {
      pt=fPtMin+(i+0.5)*step;
      w=func.Eval(pt);
      fSumPt->AddBinContent(i+1,w);
    }
  }

// photons from mixed HHG phase
  for(Int_t j=0; j<4; j++) {
    TF1 func("func",&fromMixH,fPtMin,fPtMax,5);
    func.SetParameters(lamHHG,yNucl,fTc,yRap,j);
    func.SetParNames("lamQGP","rapidity of projectile nucleus","critical temperature","photon rapidity","process number");
    for(Int_t i=0; i<nPt; i++) {
      pt=fPtMin+(i+0.5)*step;
      w=func.Eval(pt);
      fSumPt->AddBinContent(i+1,w);
    }
  }
  
// photons from pure HHG phase
  for(Int_t j=0; j<4; j++) {
    TF1 func("func",&fromHHG,fPtMin,fPtMax,6);
    func.SetParameters(tauCHHG,tauF,yNucl,fTc,yRap,j);
    func.SetParNames("nuclear radius","start HHG","freeze-out proper time","rapidity of projectile nucleus","critical temperature","photon rapidity","process number");
    for(Int_t i=0; i<nPt; i++) {
      pt=fPtMin+(i+0.5)*step;
      w=func.Eval(pt);
      fSumPt->AddBinContent(i+1,w);
    }
  }
  
}

//_____________________________________________________________________________
void AliGenThermalPhotons::Generate()
{
  //
  // Generate thermal photons of a event 
  //

    Float_t polar[3]= {0,0,0};
    Float_t origin[3];
    Float_t p[3];
    Float_t random[6];
    Int_t nt;

    for (Int_t j=0;j<3;j++) origin[j]=fOrigin[j];
/*
    if(fVertexSmear==kPerEvent) {
      Vertex();
      for (j=0;j<3;j++) origin[j]=fVertex[j];
    }
*/
    TArrayF eventVertex;
    eventVertex.Set(3);
    eventVertex[0] = origin[0];
    eventVertex[1] = origin[1];
    eventVertex[2] = origin[2];

    Int_t nGam;
    Float_t impPar,area,pt,rapidity,phi,ww;
    Float_t r0=1.3,ra=r0*TMath::Power(fAProjectile,1./3.),rb=r0*TMath::Power(fATarget,1./3.);

    TF1 *funcOL = new TF1("funcOL",fOverlapAB,fMinImpactParam,fMaxImpactParam,2); 
    funcOL->SetParameters(ra,rb);
    funcOL->SetParNames("radiusA","radiusB");

    impPar=TMath::Sqrt(fMinImpactParam*fMinImpactParam+(fMaxImpactParam*fMaxImpactParam-fMinImpactParam*fMinImpactParam)*gRandom->Rndm());
    area=funcOL->Eval(impPar);

      ww=area*(fYMax-fYMin)*fSumPt->GetBinWidth(1)*fSumPt->GetSumOfWeights();
      nGam=Int_t(ww);
      if(gRandom->Rndm() < (ww-(Float_t)nGam)) nGam++;

      if(nGam) {
        for(Int_t i=0; i<nGam; i++) {
          pt=fSumPt->GetRandom();
          Rndm(random,2);
          rapidity=(fYMax-fYMin)*random[0]+fYMin;
          phi=2.*TMath::Pi()*random[1];
          p[0]=pt*TMath::Cos(phi);
          p[1]=pt*TMath::Sin(phi);
          p[2]=pt*TMath::SinH(rapidity);

	  if(fVertexSmear==kPerTrack) {
            Rndm(random,6);
	    for (Int_t j=0;j<3;j++) {
	      origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    }
	  }

	  PushTrack(fTrackIt,-1,22,p,origin,polar,0,kPPrimary,nt,1.);
        }
      }

    delete funcOL;
// Header
    AliGenEventHeader* header = new AliGenEventHeader("ThermalPhotons");
// Event Vertex
    header->SetPrimaryVertex(eventVertex);
    header->SetNProduced(fNpart);
    gAlice->SetGenEventHeader(header);

}

void AliGenThermalPhotons::SetPtRange(Float_t ptmin, Float_t ptmax) {
    AliGenerator::SetPtRange(ptmin, ptmax);
}

void AliGenThermalPhotons::SetYRange(Float_t ymin, Float_t ymax) {
    AliGenerator::SetYRange(ymin, ymax);
}
