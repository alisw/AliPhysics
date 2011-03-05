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
// Generator of prompt photons for the reaction A+B, sqrt(S)
//
// main assumptions:
// 1. flat rapidity distribution
// 2. all existing p+p(pbar) data at y_{c.m.} can be described by the function
//           F(x_T) = (sqrt(s))^5 Ed^3sigma/d^3p, x_T = 2p_t/sqrt(s)
//           all data points cover the region x_T: 0.01 - 0.6
//    see Nucl.Phys.A783:577-582,2007, hep-ex/0609037
// 3. binary scaling: for A+B at the impact parameter b
//    Ed^3N^{AB}(b)/d^3p = Ed^3sigma^{pp}/d^3p A B T_{AB}(b),
//    T_{AB}(b) - nuclear overlapping fuction, calculated in the Glauber approach,
//                nuclear density is parametrized by a Woods-Saxon with nuclear radius
//                R_A = 1.19 A^{1/3} - 1.61 A^{-1/3} fm and surface thickness a=0.54 fm
// 4. nuclear effects (Cronin, shadowing, ...) are ignored
//
// input parameters:
//       fAProjectile, fATarget - number of nucleons in a nucleus A and B
//       fMinImpactParam - minimal impct parameter, fm
//       fMaxImpactParam - maximal impct parameter, fm
//       fEnergyCMS - sqrt(S) per nucleon pair, AGeV
//
//       fYMin - minimal rapidity of photons 
//       fYMax - maximal rapidity of photons
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

#include "AliConst.h"
#include "AliGenEventHeader.h"
#include "AliGenPromptPhotons.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliGenPromptPhotons)

TF1*    AliGenPromptPhotons::fgDataPt        = NULL;
TF1*    AliGenPromptPhotons::fgWSzA          = NULL;
TF1*    AliGenPromptPhotons::fgWSzB          = NULL;
TF1*    AliGenPromptPhotons::fgTA            = NULL;
TF1*    AliGenPromptPhotons::fgTB            = NULL;
TF1*    AliGenPromptPhotons::fgTAxTB         = NULL;
TF1*    AliGenPromptPhotons::fgTAB           = NULL;

//_____________________________________________________________________________
AliGenPromptPhotons::AliGenPromptPhotons()
    :AliGenerator(-1),
        fAProjectile(0.),
        fATarget(0.),
        fEnergyCMS(0.),
        fMinImpactParam(0.),
        fMaxImpactParam(0.)
{
    //
    // Default constructor
    //
    SetCutVertexZ();
    SetPtRange();
    SetYRange();
}

//_____________________________________________________________________________
AliGenPromptPhotons::AliGenPromptPhotons(Int_t npart)
    :AliGenerator(npart),
        fAProjectile(208),
        fATarget(208),
        fEnergyCMS(5500.),
        fMinImpactParam(0.),
        fMaxImpactParam(0.)
{
  // 
  // Standard constructor
  //

    fName="PromptPhotons";
    fTitle="Prompt photons from pp data fit + T_AB";

    SetCutVertexZ();
    SetPtRange();
    SetYRange();
}

//_____________________________________________________________________________
AliGenPromptPhotons::~AliGenPromptPhotons()
{
  //
  // Standard destructor
  //
    delete fgDataPt;
    delete fgWSzA;
    delete fgWSzB;
    delete fgTA;
    delete fgTB;
    delete fgTAxTB;
    delete fgTAB;
}

//_____________________________________________________________________________
void AliGenPromptPhotons::Init()
{
  // Initialisation 
  fgDataPt = new TF1("fgDataPt",FitData,fPtMin,fPtMax,1);
  fgDataPt->SetParameter(0,fEnergyCMS);

  const Double_t d=0.54;  // surface thickness (fm)
  const Double_t ra = 1.19*TMath::Power(fAProjectile,1./3.)-1.61/TMath::Power(fAProjectile,1./3.);
  const Double_t eps=0.01; // cut WS at ramax: WS(ramax)/WS(0)=eps
  const Double_t ramax=ra+d*TMath::Log((1.-eps+TMath::Exp(-ra/d))/eps);

  TF1 fWSforNormA("fWSforNormA",&WSforNorm,0,ramax,2);
  fWSforNormA.SetParameters(ra,d);
  fWSforNormA.SetParNames("RA","d");
  const Double_t ca=1./fWSforNormA.Integral(0.,ramax);

  const Double_t rb=1.19*TMath::Power(fATarget,1./3.)-1.61/TMath::Power(fATarget,1./3.);
  const Double_t rbmax=rb+d*TMath::Log((1.-eps+TMath::Exp(-rb/d))/eps);

  TF1 fWSforNormB("fWSforNormB",&WSforNorm,0,rbmax,2);
  fWSforNormB.SetParameters(rb,d);
  fWSforNormB.SetParNames("RB","d");
  const Double_t cb=1./fWSforNormB.Integral(0.,rbmax);

  fgWSzA = new TF1("fgWSzA",WSz,0.,ramax,4);
  fgWSzA->SetParameter(0,ra);
  fgWSzA->SetParameter(1,d);
  fgWSzA->SetParameter(2,ca);

  fgTA = new TF1("fgTA",TA,0.,ramax,1);
  fgTA->SetParameter(0,ramax);

  fgWSzB = new TF1("fgWSzB",WSz,0.,rbmax,4);
  fgWSzB->SetParameter(0,rb);
  fgWSzB->SetParameter(1,d);
  fgWSzB->SetParameter(2,cb);

  fgTB = new TF1("fgTB",TB,0.,TMath::Pi(),3);
  fgTB->SetParameter(0,rbmax);

  fgTAxTB = new TF1("fgTAxTB",TAxTB,0,ramax,2);
  fgTAxTB->SetParameter(0,rbmax);

  fgTAB = new TF1("fgTAB",TAB,0.,ramax+rbmax,2);
  fgTAB->SetParameter(0,ramax);
  fgTAB->SetParameter(1,rbmax);

}

//_____________________________________________________________________________
void AliGenPromptPhotons::Generate()
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
    Float_t b,pt,rapidity,phi,ww;

    b=TMath::Sqrt(fMinImpactParam*fMinImpactParam+(fMaxImpactParam*fMaxImpactParam-fMinImpactParam*fMinImpactParam)*gRandom->Rndm());

    Double_t dsdyPP=fgDataPt->Integral(fPtMin,fPtMax); // pb, fm^2 = 1e10 pb
    ww=dsdyPP*1.0e-10*(fYMax-fYMin)*fAProjectile*fATarget*fgTAB->Eval(b);
    nGam=Int_t(ww);
    if(gRandom->Rndm() < (ww-(Float_t)nGam)) nGam++;

      if(nGam) {
        for(Int_t i=0; i<nGam; i++) {
          pt=fgDataPt->GetRandom();
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

// Header
    AliGenEventHeader* header = new AliGenEventHeader("PromptPhotons");
// Event Vertex
    header->SetPrimaryVertex(eventVertex);
    header->SetNProduced(fNpart);
    gAlice->SetGenEventHeader(header);

}

void AliGenPromptPhotons::SetPtRange(Float_t ptmin, Float_t ptmax) {
    AliGenerator::SetPtRange(ptmin, ptmax);
}

void AliGenPromptPhotons::SetYRange(Float_t ymin, Float_t ymax) {
    AliGenerator::SetYRange(ymin, ymax);
}

//**********************************************************************************
Double_t AliGenPromptPhotons::FitData(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - p_t (GeV).
// par[0]=sqrt(s_NN) (GeV),
//
// output:
// d^{2}#sigma/(dp_t dy) (pb/GeV)
//---------------------------------------------------
//
// d^{2}#sigma/(dp_t dy) = (2 pi p_t) Ed^{3}#sigma/d^{3}p 
//
// data presentation: Nucl.Phys.A783:577-582,2007, hep-ex/0609037, fig.3
// F(x_t)=(#sqrt{s})^{5} Ed^{3}#sigma/d^{3}p
//---------------------------------------------------
// approximate tabulation of F(x_t)
  const Int_t nMax=4;
  const Double_t log10x[nMax]={ -2., -1., -0.6, -0.3};
  const Double_t log10F[nMax]={ 19., 13.,  9.8,   7.};

  const Double_t xT=2.*x[0]/par[0];
  const Double_t log10xT=TMath::Log10(xT);
  Int_t i=0;
  while(log10xT>log10x[i] && i<nMax) i++;
  if(i==0) i=1;
  if(i==nMax) i=nMax-1;
  const Double_t alpha=(log10F[i]-log10F[i-1])/(log10x[i]-log10x[i-1]);
  const Double_t beta=log10F[i-1]-alpha*log10x[i-1];

  return (TMath::TwoPi()*x[0])*TMath::Power(10.,alpha*log10xT+beta)/TMath::Power(par[0],5.);
}

//**********************************************************************************
Double_t AliGenPromptPhotons::WSforNorm(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - r (fm)
// par[0] - R (fm), radius
// par[1] - d (fm), surface thickness
//
// output: 
// 4 pi r**2 /(1+exp((r-R)/d))
//
// Wood Saxon (WS) C/(1+exp((r-RA)/d)) (nuclons/fm^3) 
// To get the normalization A = (Integral 4 pi r**2 dr WS):
// C = A / (Integral 4 pi r**2 dr 1/(1+exp((r-RA)/d)) )
// Thus me need 4 pi r**2 /(1+exp((r-RA)/d)) (1/fm)
//---------------------------------------------------
  const Double_t r=x[0];
  const Double_t bigR=par[0],d=par[1];

  return 4.*TMath::Pi()*r*r/(1.+TMath::Exp((r-bigR)/d));
}

//**********************************************************************************
Double_t AliGenPromptPhotons::WSz(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - z (fm)
// par[0] - R (fm), radius
// par[1] - d (fm), surface thickness
// par[2] - C (nucleons/fm**2), normalization factor 
// par[3] - b (fm), impact parameter
//
// output:
//  Wood Saxon Parameterisation
//  as a function of z for fixed b (1/fm^3)
//---------------------------------------------------
  const Double_t z=x[0];
  const Double_t bigR=par[0],d=par[1],C=par[2],b=par[3];

  return C/(1.+TMath::Exp((TMath::Sqrt(b*b+z*z)-bigR)/d));
}

//**********************************************************************************
Double_t AliGenPromptPhotons::TA(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - b (fm), impact parameter
// par[0] - RAMAX (fm), max. value of projectile radius
//
// output:
// nuclear thickness function T_A(b) (1/fm^2)
//---------------------------------------------------
  const Double_t b=x[0];
  const Double_t ramax=par[0];

  fgWSzA->SetParameter(3,b);

  return 2.*fgWSzA->Integral(0.,TMath::Sqrt(ramax*ramax-b*b));
}

//**********************************************************************************
Double_t AliGenPromptPhotons::TB(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - phi (rad)
// par[0] - RBMAX (fm), max. value of target radius
// par[1] - b (fm), impact parameter
// par[2] - s (fm)
//
// output:
//  nuclear thickness function T_B(phi)=T_B(sqtr(s**2+b**2-2*s*b*cos(phi)))
//---------------------------------------------------
  const Double_t phi=x[0];
  const Double_t rbmax=par[0],b=par[1],s=par[2];

  const Double_t w=TMath::Sqrt(s*s+b*b-2.*s*b*TMath::Cos(phi));

  fgWSzB->SetParameter(3,w);

  return 2.*fgWSzB->Integral(0.,TMath::Sqrt(rbmax*rbmax-w*w));;
}

//**********************************************************************************
Double_t AliGenPromptPhotons::TAxTB(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - s (fm)
// par[0] - RBMAX (fm), max. value of target radius
// par[1] - b (fm), impact parameter
//
// output:
//  s * TA(s) * 2 * Integral(0,phiMax) TB(phi(s,b))
//---------------------------------------------------
  const Double_t s  = x[0];
  const Double_t rbmax=par[0],b=par[1];

  if(s==0.) return 0.;

  fgTB->SetParameter(1,b);
  fgTB->SetParameter(2,s);

  Double_t phiMax;
  if(b<rbmax && s<(rbmax-b)) {
    phiMax=TMath::Pi();
  } else {
    phiMax=TMath::ACos((s*s+b*b-rbmax*rbmax)/(2.*s*b));
  }

  return s*fgTA->Eval(s)*2.*fgTB->Integral(0.,phiMax);
}

// ---------------------------------------------------------------------------------
Double_t AliGenPromptPhotons::TAB(const Double_t* x, const Double_t* par) {
//---------------------------------------------------
// input:
// x[0] - b (fm), impact parameter
// par[0] - RAMAX (fm), max. value of projectile radius
// par[1] - RAMAX (fm), max. value of target radius
//
// output:
// overlap function TAB(b) (1/fm**2)
//---------------------------------------------------
  const Double_t b=x[0];
  const Double_t ramax=par[0],rbmax=par[1];

  Double_t sMin=0.,sMax=ramax;
  if(b>rbmax) sMin=b-rbmax;
  if(b<(ramax-rbmax)) sMax=b+rbmax;

  fgTAxTB->SetParameter(1,b);

  return fgTAxTB->Integral(sMin,sMax);
}
