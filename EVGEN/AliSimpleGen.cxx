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
Revision 1.5  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////
//                                                               //
//    Generate the final state of the interaction as the input   //
//    to the MonteCarlo                                          //
//
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

#include "AliSimpleGen.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliGenHIJINGpara)

//_____________________________________________________________________________
static Double_t ptpi(Double_t *px, Double_t *)
{
  //
  //     PT-PARAMETERIZATION CDF, PRL 61(88) 1819
  //     POWER LAW FOR PT > 500 MEV
  //     MT SCALING BELOW (T=160 MEV)
  //
  const Double_t p0 = 1.3;
  const Double_t xn = 8.28;
  const Double_t xlim=0.5;
  const Double_t t=0.160;
  const Double_t xmpi=0.139;
  const Double_t b=1.;
  Double_t y, y1, xmpi2, ynorm, a;
  Double_t x=*px;
  //
  y1=TMath::Power(p0/(p0+xlim),xn);
  xmpi2=xmpi*xmpi;
  ynorm=b*(TMath::Exp(-sqrt(xlim*xlim+xmpi2)/t));
  a=ynorm/y1;
  if (x > xlim)
    y=a*TMath::Power(p0/(p0+x),xn);
  else
    y=b*TMath::Exp(-sqrt(x*x+xmpi2)/t);
  return y*x;
}

//_____________________________________________________________________________
static Double_t ptscal(Double_t pt, Int_t np)
{
  //    SCALING EN MASSE PAR RAPPORT A PTPI
  //     MASS PI,K,ETA,RHO,OMEGA,ETA',PHI
  const Double_t hm[10] = {.13957,.493,.5488,.769,.7826,.958,1.02,0,0,0};
  //     VALUE MESON/PI AT 5 GEV
  const Double_t fmax[10]={1.,0.3,0.55,1.0,1.0,1.0,1.0,0,0,0};
  np--;
  Double_t f5=TMath::Power(((sqrt(100.018215)+2.)/(sqrt(100.+hm[np]*hm[np])+2.0)),12.3);
  Double_t fmax2=f5/fmax[np];
  // PIONS
  Double_t ptpion=100.*ptpi(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
				 (sqrt(pt*pt+hm[np]*hm[np])+2.0)),12.3)/ fmax2;
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
  const Double_t a1    = 4913.;
  const Double_t a2    = 1819.;
  const Double_t eta1  = 0.22;
  const Double_t eta2  = 3.66;
  const Double_t deta1 = 1.47;
  const Double_t deta2 = 1.51;
  Double_t y=TMath::Abs(*py);
  //
  Double_t ex1 = (y-eta1)*(y-eta1)/(2*deta1*deta1);
  Double_t ex2 = (y-eta2)*(y-eta2)/(2*deta2*deta2);
  return a1*TMath::Exp(-ex1)+a2*TMath::Exp(-ex2);
}

//_____________________________________________________________________________
static Double_t etakac( Double_t *py, Double_t *)
{
  //
  // eta parametrisation for ka
  //
  const Double_t a1    = 497.6;
  const Double_t a2    = 215.6;
  const Double_t eta1  = 0.79;
  const Double_t eta2  = 4.09;
  const Double_t deta1 = 1.54;
  const Double_t deta2 = 1.40;
  Double_t y=TMath::Abs(*py);
  //
  Double_t ex1 = (y-eta1)*(y-eta1)/(2*deta1*deta1);
  Double_t ex2 = (y-eta2)*(y-eta2)/(2*deta2*deta2);
  return a1*TMath::Exp(-ex1)+a2*TMath::Exp(-ex2);
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
  Float_t etaMin = -TMath::Log(TMath::Tan(TMath::Min((Double_t)fThetaMax/2,TMath::Pi()/2-1.e-10)));
  Float_t etaMax = -TMath::Log(TMath::Tan(TMath::Max((Double_t)fThetaMin/2,              1.e-10)));
  fPtpi = new TF1("ptpi",&ptpi,0,20,0);
  fPtka = new TF1("ptka",&ptka,0,20,0);
  fETApic = new TF1("etapic",&etapic,etaMin,etaMax,0);
  fETAkac = new TF1("etakac",&etakac,etaMin,etaMax,0);
  TF1 *ETApic0 = new TF1("etapic",&etapic,-7,7,0);
  TF1 *ETAkac0 = new TF1("etakac",&etakac,-7,7,0);
  Float_t IntETApi  = ETApic0->Integral(-0.5, 0.5);
  Float_t IntETAka  = ETAkac0->Integral(-0.5, 0.5);
  Float_t scalePi=7316/(IntETApi/1.5);
  Float_t scaleKa= 684/(IntETAka/2.0);

  Float_t IntPt  = (0.877*ETApic0->Integral(0, 15)+
		    0.123*ETAkac0->Integral(0, 15));
  Float_t IntPtSel = (0.877*ETApic0->Integral(fPtMin, fPtMax)+
		      0.123*ETAkac0->Integral(fPtMin, fPtMax));
  Float_t PtFrac = IntPtSel/IntPt;
  

  Float_t IntETASel  = (scalePi*ETApic0->Integral(etaMin, etaMax)+
			scaleKa*ETAkac0->Integral(etaMin, etaMax));
  Float_t PhiFrac = (fPhiMax-fPhiMin)/2/TMath::Pi();
  fParentWeight = Float_t(fNpart)/IntETASel*PtFrac*PhiFrac;
  
  printf("\n The number of particles in the selected kinematic region corresponds to %f percent of a full event\n ", 100.*fParentWeight);
  
}

//_____________________________________________________________________________
void AliGenHIJINGpara::Generate()
{
  //
  // Generate one trigger
  //

  
  const Float_t raKpic=0.14;
  const Float_t borne=1/(1+raKpic);
  Float_t polar[3]= {0,0,0};
  //
  const Int_t pions[3] = {kPi0, kPiPlus, kPiMinus};
  const Int_t kaons[4] = {kK0Long, kK0Short, kKPlus, kKMinus};
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
  if(fVertexSmear==perEvent) {
    gMC->Rndm(random,6);
    for (j=0;j<3;j++) {
      origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }
  for(i=0;i<fNpart;i++) {
    while(1) {
      gMC->Rndm(random,3);
      if(random[0]<borne) {
	part=pions[Int_t (random[1]*3)];
	ptf=fPtpi;
	etaf=fETApic;
      } else {
	part=kaons[Int_t (random[1]*4)];
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
      if(fVertexSmear==perTrack) {
	gMC->Rndm(random,6);
	for (j=0;j<3;j++) {
	  origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	}
      }
      gAlice->SetTrack(fTrackIt,-1,part,p,origin,polar,0,"Primary",nt,fParentWeight);
      break;
    }
  }
}
  
ClassImp(AliGenFixed)

//_____________________________________________________________________________
AliGenFixed::AliGenFixed()
  :AliGenerator()
{
  //
  // Default constructor
  //
  fIpart = 0;
}

//_____________________________________________________________________________
AliGenFixed::AliGenFixed(Int_t npart)
  :AliGenerator(npart)
{
  //
  // Standard constructor
  //
  fName="Fixed";
  fTitle="Fixed Particle Generator";
  // Generate Proton by default
  fIpart=kProton;
}

//_____________________________________________________________________________
void AliGenFixed::Generate()
{
  //
  // Generate one trigger
  //
  Float_t polar[3]= {0,0,0};
  Float_t p[3] = {fPMin*TMath::Cos(fPhiMin)*TMath::Sin(fThetaMin),
		  fPMin*TMath::Sin(fPhiMin)*TMath::Sin(fThetaMin),
		  fPMin*TMath::Cos(fThetaMin)};
  Int_t i, nt;
  //
  for(i=0;i<fNpart;i++) {
    gAlice->SetTrack(fTrackIt,-1,fIpart,p,fOrigin.GetArray(),polar,0,"Primary",nt);
  }
}
  
//_____________________________________________________________________________
void AliGenFixed::SetSigma(Float_t, Float_t, Float_t)
{
  //
  // Set the interaction point sigma
  //
  printf("Vertex smearing not implemented for fixed generator\n");
}


ClassImp(AliGenBox)

//_____________________________________________________________________________
AliGenBox::AliGenBox()
    :AliGenerator()
{
  //
  // Default constructor
  //
  fIpart=0;
}

//_____________________________________________________________________________
AliGenBox::AliGenBox(Int_t npart)
  :AliGenerator(npart)
{
  //
  // Standard constructor
  //
  fName="Box";
  fTitle="Box particle generator";
  // Generate Proton by default
  fIpart=kProton;
}

//_____________________________________________________________________________
void AliGenBox::Generate()
{
  //
  // Generate one trigger
  //
  
  Float_t polar[3]= {0,0,0};
  //
  Float_t origin[3];
  Float_t p[3];
  Int_t i, j, nt;
  Float_t pmom, theta, phi;
  //
  Float_t random[6];
  //
  for (j=0;j<3;j++) origin[j]=fOrigin[j];
  if(fVertexSmear==perEvent) {
    gMC->Rndm(random,6);
    for (j=0;j<3;j++) {
      origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }
  for(i=0;i<fNpart;i++) {
    gMC->Rndm(random,3);
    pmom=fPMin+random[0]*(fPMax-fPMin);
    theta=fThetaMin+random[1]*(fThetaMax-fThetaMin);
    phi=fPhiMin+random[2]*(fPhiMax-fPhiMin);
    p[0] = pmom*TMath::Cos(phi)*TMath::Sin(theta);
    p[1] = pmom*TMath::Sin(phi)*TMath::Sin(theta);
    p[2] = pmom*TMath::Cos(theta);
    if(fVertexSmear==perTrack) {
      gMC->Rndm(random,6);
      for (j=0;j<3;j++) {
	origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	  TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
      }
    }
    gAlice->SetTrack(fTrackIt,-1,fIpart,p,origin,polar,0,"Primary",nt);
  }
}


