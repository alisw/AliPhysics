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

///////////////////////////////////////////////////////////////////
// Parameterisation of pi, K, n and p eta and pt distributions   //
// eta: according to HIJING (shadowing + quenching)              //
// pT : according to CDF measurement at 1.8 TeV                  //
// Author: andreas.morsch@cern.ch                                //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TArrayF.h>
#include <TF1.h>
#include <TPDGCode.h>

#include "AliConst.h"
#include "AliGenEventHeader.h"
#include "AliGenHIJINGparaBa.h"
#include "AliRun.h"

ClassImp(AliGenHIJINGparaBa)


static Double_t ptpi(Double_t *px, Double_t *)
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

//_____________________________________________________________________________
static Double_t ptscal(Double_t pt, Int_t np)
{
    //    SCALING EN MASSE PAR RAPPORT A PTPI
    //     MASS PI,K,ETA,RHO,OMEGA,ETA',PHI
    const Double_t khm[10] = {.13957,.493,.5488,.769,.7826,.958,1.02,0,0,0};
    //     VALUE MESON/PI AT 5 GEV
    const Double_t kfmax[10]={1.,0.3,0.55,1.0,1.0,1.0,1.0,0,0,0};
    np--;
    Double_t f5=TMath::Power(((
	sqrt(100.018215)+2.)/(sqrt(100.+khm[np]*khm[np])+2.0)),12.3);
    Double_t fmax2=f5/kfmax[np];
    // PIONS
    Double_t ptpion=100.*ptpi(&pt, (Double_t*) 0);
    Double_t fmtscal=TMath::Power(((
	sqrt(pt*pt+0.018215)+2.)/ (sqrt(pt*pt+khm[np]*khm[np])+2.0)),12.3)/ 
	fmax2;
    return fmtscal*ptpion;
}

//_____________________________________________________________________________
static Double_t ptka( Double_t *px, Double_t *)
{
    //
    // pt parametrisation for k
    //
    return ptscal(*px,2);
}


//_____________________________________________________________________________
static Double_t etapic( Double_t *py, Double_t *)
{
  //
  // eta parametrisation for pi
  //
    const Double_t ka1    = 4913.;
    const Double_t ka2    = 1819.;
    const Double_t keta1  = 0.22;
    const Double_t keta2  = 3.66;
    const Double_t kdeta1 = 1.47;
    const Double_t kdeta2 = 1.51;
    Double_t y=TMath::Abs(*py);
    //
    Double_t ex1 = (y-keta1)*(y-keta1)/(2*kdeta1*kdeta1);
    Double_t ex2 = (y-keta2)*(y-keta2)/(2*kdeta2*kdeta2);
    return ka1*TMath::Exp(-ex1)+ka2*TMath::Exp(-ex2);
}

//_____________________________________________________________________________
static Double_t etakac( Double_t *py, Double_t *)
{
    //
    // eta parametrisation for ka
    //
    const Double_t ka1    = 497.6;
    const Double_t ka2    = 215.6;
    const Double_t keta1  = 0.79;
    const Double_t keta2  = 4.09;
    const Double_t kdeta1 = 1.54;
    const Double_t kdeta2 = 1.40;
    Double_t y=TMath::Abs(*py);
    //
    Double_t ex1 = (y-keta1)*(y-keta1)/(2*kdeta1*kdeta1);
    Double_t ex2 = (y-keta2)*(y-keta2)/(2*kdeta2*kdeta2);
    return ka1*TMath::Exp(-ex1)+ka2*TMath::Exp(-ex2);
}

 static Double_t ptbaryon( Double_t *px, Double_t *)
{
// baryons
//                pt-distribution
//____________________________________________________________

  return ptscal(*px,7);  //  7==> Baryon in the PtScal function
}

 static Double_t etabaryon( Double_t *py, Double_t *)
{
// eta-distribution
//____________________________________________________________
    const Float_t  kp0 =  1.10343e+02;
    const Float_t  kp1 =  1.73247e+01;
    const Float_t  kp2 = -7.23808e+00;
    const Float_t  kp3 =  4.48334e-01;
    const Double_t ky = TMath::Abs(*py);
//
    return (kp0+kp1*ky+kp2*ky*ky+kp3*ky*ky*ky)/20.;
}

AliGenHIJINGparaBa::AliGenHIJINGparaBa()
    :AliGenHIJINGpara(),
     fPtba(0),
     fETAba(0)
{
    //
    // Default constructor
    //
    fName="HIGINGparaBa";
    fTitle="HIJING Parametrisation Particle Generator with Baryons";
}

//_____________________________________________________________________________
AliGenHIJINGparaBa::AliGenHIJINGparaBa(Int_t npart)
  :AliGenHIJINGpara(npart),
     fPtba(0),
     fETAba(0)
{
  // 
  // Standard constructor
  //
    fName="HIGINGparaBa";
    fTitle="HIJING Parametrisation Particle Generator with Baryons";
}

//_____________________________________________________________________________
AliGenHIJINGparaBa::~AliGenHIJINGparaBa()
{
  //
  // Standard destructor
  //
    delete fPtba;
    delete fETAba;
}

//_____________________________________________________________________________
void AliGenHIJINGparaBa::Init()
{
  //
  // Initialise the HIJING parametrisation
  //
    Float_t etaMin =-TMath::Log(TMath::Tan(
	TMath::Min((Double_t)fThetaMax/2,TMath::Pi()/2-1.e-10)));
    Float_t etaMax = -TMath::Log(TMath::Tan(
	TMath::Max((Double_t)fThetaMin/2,1.e-10)));
    fPtpi   = new TF1("ptpi",&ptpi,0,20,0);
    fPtka   = new TF1("ptka",&ptka,0,20,0);
    fPtba   = new TF1("ptbaryon",&ptbaryon,0,20,0);
    fETApic = new TF1("etapic",&etapic,etaMin,etaMax,0);
    fETAkac = new TF1("etakac",&etakac,etaMin,etaMax,0);
    fETAba  = new TF1("etabaryon",&etabaryon,etaMin,etaMax,0);

    TF1 etaPic0("etapic(-7,7)",&etapic,    -7, 7, 0);
    TF1 etaKac0("etakac(-7,7)",&etakac,    -7, 7, 0);
    TF1 etaBar0("etabar(-7,7)",&etabaryon, -7, 7, 0);
    
    TF1 ptPic0("ptpi(0,15)",  &ptpi,     0., 15., 0);
    TF1 ptKac0("ptka(0,15)",  &ptka,     0., 15., 0);
    TF1 ptBar0("ptbar(0,15)", &ptbaryon, 0., 15., 0);

    Float_t intETApi  = etaPic0.Integral(-0.5, 0.5);
    Float_t intETAka  = etaKac0.Integral(-0.5, 0.5);
    Float_t intETAba  = etaBar0.Integral(-0.5, 0.5);

    Float_t scalePi   = 6979./(intETApi/1.5);
    Float_t scaleKa   =  657./(intETAka/2.0);
    Float_t scaleBa   =  364./(intETAba/2.0);

//  Fraction of events corresponding to the selected pt-range    
    Float_t intPt    = (0.837*ptPic0.Integral(0, 15)+
			0.105*ptKac0.Integral(0, 15)+
                        0.058*ptBar0.Integral(0, 15));
    Float_t intPtSel = (0.837*ptPic0.Integral(fPtMin, fPtMax)+
			0.105*ptKac0.Integral(fPtMin, fPtMax)+
	                0.058*ptBar0.Integral(fPtMin, fPtMax));
    Float_t ptFrac   = intPtSel/intPt;

//  Fraction of events corresponding to the selected eta-range    
    Float_t intETASel  = (scalePi*etaPic0.Integral(etaMin, etaMax)+
			  scaleKa*etaKac0.Integral(etaMin, etaMax)+
	                  scaleBa*etaBar0.Integral(etaMin, etaMax));
//  Fraction of events corresponding to the selected phi-range    
    Float_t phiFrac    = (fPhiMax-fPhiMin)/2/TMath::Pi();

    fParentWeight = Float_t(fNpart)/(intETASel*ptFrac*phiFrac);
    
    printf("%s: The number of particles in the selected kinematic region corresponds to %f percent of a full event \n", 
	   ClassName(),100.*fParentWeight);

// Issue warning message if etaMin or etaMax are outside the alowed range 
// of the parametrization
    if (etaMin < -8.001 || etaMax > 8.001) {
	printf("\n \n WARNING FROM AliGenHIJINGParaBa !");
	printf("\n YOU ARE USING THE PARAMETERISATION OUTSIDE ");	
	printf("\n THE ALLOWED PSEUDORAPIDITY RANGE (-8. - 8.)");	    
	printf("\n YOUR LIMITS: %f %f \n \n ", etaMin, etaMax);
    }
}

//_____________________________________________________________________________
void AliGenHIJINGparaBa::Generate()
{
  //
  // Generate one trigger
  //

  
    const Float_t kBorne1 = 0.837;
    const Float_t kBorne2 = kBorne1+0.105;
    
    Float_t polar[3]= {0,0,0};
    //
    const Int_t kPions[3]   = {kPi0, kPiPlus, kPiMinus};
    const Int_t kKaons[4]   = {kK0Long, kK0Short, kKPlus, kKMinus};
    const Int_t kBaryons[4] = {kProton, kProtonBar, kNeutron, kNeutronBar};
    //
    Float_t origin[3];
    Float_t pt, pl, ptot;
    Float_t phi, theta;
    Float_t p[3];
    Int_t i, part, nt, j;
    //
    TF1 *ptf;
    TF1 *etaf;
    //
    Float_t random[6];
    //
    for (j=0;j<3;j++) origin[j]=fOrigin[j];

    if(fVertexSmear == kPerEvent) {
	Float_t dv[3];
	dv[2] = 1.e10;
	while(TMath::Abs(dv[2]) > fCutVertexZ*fOsigma[2]) {
	    Rndm(random,6);
	    for (j=0; j < 3; j++) {
		dv[j] = fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    }
	}
	for (j=0; j < 3; j++) origin[j] += dv[j];
    } // if kPerEvent
    TArrayF eventVertex;
    eventVertex.Set(3);
    eventVertex[0] = origin[0];
    eventVertex[1] = origin[1];
    eventVertex[2] = origin[2];

    for(i=0;i<fNpart;i++) {
	while(1) {
	    Rndm(random,3);
	    if(random[0] < kBorne1) {
		part  = kPions[Int_t (random[1]*3)];
		ptf   = fPtpi;
		etaf  = fETApic;
	    } else if (random[0] < kBorne2) {
		part  = kKaons[Int_t (random[1]*4)];
		ptf   = fPtka;
		etaf  = fETAkac;
	    } else {
		part  = kBaryons[Int_t (random[1]*4)];
		ptf   = fPtba;
		etaf  = fETAba;
	    }
	    
	    phi=fPhiMin+random[2]*(fPhiMax-fPhiMin);
	    theta=2*TMath::ATan(TMath::Exp(-etaf->GetRandom()));
	    if(theta<fThetaMin || theta>fThetaMax) continue;
	    pt=ptf->GetRandom();
	    pl=pt/TMath::Tan(theta);
	    ptot=TMath::Sqrt(pt*pt+pl*pl);
	    if(ptot<fPMin || ptot>fPMax) continue;
	    p[0]=pt*TMath::Cos(phi);
	    p[1]=pt*TMath::Sin(phi);
	    p[2]=pl;
	    if(fVertexSmear==kPerTrack) {
		Rndm(random,6);
		for (j=0;j<3;j++) {
		    origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
			TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
		}
	    }
	    PushTrack(fTrackIt,-1,part,p,origin,polar,0,kPPrimary,nt,fParentWeight);
	    break;
	} // while(1)
    } // Particle loop
// Header
    AliGenEventHeader* header = new AliGenEventHeader("HIJINGparam");
// Event Vertex
    header->SetPrimaryVertex(eventVertex);
    gAlice->SetGenEventHeader(header); 
}



