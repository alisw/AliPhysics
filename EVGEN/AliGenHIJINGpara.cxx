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

/*
$Log$
Revision 1.4  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.3  2000/10/02 21:28:06  fca
Removal of useless dependecies via forward declarations

Revision 1.2  2000/07/11 18:24:55  fca
Coding convention corrections + few minor bug fixes

Revision 1.1  2000/06/09 20:20:30  morsch
Same class as previously in AliSimpleGen.cxx
All coding rule violations except RS3 corrected (AM)

*/

// Parameterisation of pi and K, eta and pt distributions
// used for the ALICE TDRs.
// eta: according to HIJING (shadowing + quenching)
// pT : according to CDF measurement at 1.8 TeV
// Author: andreas.morsch@cern.ch


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

#include "AliGenHIJINGpara.h"
#include "TF1.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliPDG.h"

ClassImp(AliGenHIJINGpara)

AliGenHIJINGpara::AliGenHIJINGpara(const AliGenHIJINGpara & para)
{
// copy constructor
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
AliGenHIJINGpara::AliGenHIJINGpara()
  :AliGenerator()
{
    //
    // Default constructor
    //
    fPtpi = 0;
    fPtka = 0;
    fETApic = 0;
    fETAkac = 0;
}

//_____________________________________________________________________________
AliGenHIJINGpara::AliGenHIJINGpara(Int_t npart)
  :AliGenerator(npart)
{
  // 
  // Standard constructor
  //
    fName="HIGINGpara";
    fTitle="HIJING Parametrisation Particle Generator";
    fPtpi = 0;
    fPtka = 0;
    fETApic = 0;
    fETAkac = 0;
}

//_____________________________________________________________________________
AliGenHIJINGpara::~AliGenHIJINGpara()
{
  //
  // Standard destructor
  //
    delete fPtpi;
    delete fPtka;
    delete fETApic;
    delete fETAkac;
}

//_____________________________________________________________________________
void AliGenHIJINGpara::Init()
{
  //
  // Initialise the HIJING parametrisation
  //
    Float_t etaMin =-TMath::Log(TMath::Tan(
	TMath::Min((Double_t)fThetaMax/2,TMath::Pi()/2-1.e-10)));
    Float_t etaMax = -TMath::Log(TMath::Tan(
	TMath::Max((Double_t)fThetaMin/2,1.e-10)));
    fPtpi = new TF1("ptpi",&ptpi,0,20,0);
    fPtka = new TF1("ptka",&ptka,0,20,0);
    fETApic = new TF1("etapic",&etapic,etaMin,etaMax,0);
    fETAkac = new TF1("etakac",&etakac,etaMin,etaMax,0);
    TF1 *etaPic0 = new TF1("etapic",&etapic,-7,7,0);
    TF1 *etaKac0 = new TF1("etakac",&etakac,-7,7,0);
    Float_t intETApi  = etaPic0->Integral(-0.5, 0.5);
    Float_t intETAka  = etaKac0->Integral(-0.5, 0.5);
    Float_t scalePi=7316/(intETApi/1.5);
    Float_t scaleKa= 684/(intETAka/2.0);
    
    Float_t intPt  = (0.877*etaPic0->Integral(0, 15)+
		      0.123*etaKac0->Integral(0, 15));
    Float_t intPtSel = (0.877*etaPic0->Integral(fPtMin, fPtMax)+
			0.123*etaKac0->Integral(fPtMin, fPtMax));
    Float_t ptFrac = intPtSel/intPt;
    
    
    Float_t intETASel  = (scalePi*etaPic0->Integral(etaMin, etaMax)+
			  scaleKa*etaKac0->Integral(etaMin, etaMax));
    Float_t phiFrac = (fPhiMax-fPhiMin)/2/TMath::Pi();
    fParentWeight = Float_t(fNpart)/intETASel*ptFrac*phiFrac;
    
    printf("\n The number of particles in the selected kinematic region corresponds to %f percent of a full event\n ", 100.*fParentWeight);
    
}

//_____________________________________________________________________________
void AliGenHIJINGpara::Generate()
{
  //
  // Generate one trigger
  //

  
    const Float_t kRaKpic=0.14;
    const Float_t kBorne=1/(1+kRaKpic);
    Float_t polar[3]= {0,0,0};
    //
    const Int_t kPions[3] = {kPi0, kPiPlus, kPiMinus};
    const Int_t kKaons[4] = {kK0Long, kK0Short, kKPlus, kKMinus};
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
    if(fVertexSmear==kPerEvent) {
	Rndm(random,6);
	for (j=0;j<3;j++) {
	    origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	}
    }
    for(i=0;i<fNpart;i++) {
	while(1) {
	    Rndm(random,3);
	    if(random[0]<kBorne) {
		part=kPions[Int_t (random[1]*3)];
		ptf=fPtpi;
	      etaf=fETApic;
	    } else {
		part=kKaons[Int_t (random[1]*4)];
		ptf=fPtka;
		etaf=fETAkac;
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
	    gAlice->SetTrack(fTrackIt,-1,part,p,origin,polar,0,kPPrimary,nt,fParentWeight);
	    break;
	}
    }
}

AliGenHIJINGpara& AliGenHIJINGpara::operator=(const  AliGenHIJINGpara& rhs)
{
// Assignment operator
    return *this;
}


