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

#include <TArrayF.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliDecayer.h"
#include "AliGenEventHeader.h"
#include "AliGenHIJINGpara.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliGenHIJINGpara)



//_____________________________________________________________________________
static Double_t ptpi(Double_t *px, Double_t *)
{
  //
  //     PT-PARAMETERIZATION CDF, PRL 61(88) 1819
  //     POWER LAW FOR PT > 500 MEV
  //     MT SCALING BELOW (T=160 MEV)
  //
  const Double_t kp0    = 1.3;
  const Double_t kxn    = 8.28;
  const Double_t kxlim  = 0.5;
  const Double_t kt     = 0.160;
  const Double_t kxmpi  = 0.139;
  const Double_t kb     = 1.;
  Double_t y, y1, xmpi2, ynorm, a;
  Double_t x = *px;
  //
  y1 = TMath::Power(kp0 / (kp0 + kxlim), kxn);
  xmpi2 = kxmpi * kxmpi;
  ynorm = kb * (TMath::Exp(-sqrt(kxlim * kxlim + xmpi2) / kt ));
  a = ynorm / y1;
  if (x > kxlim)
    y = a * TMath::Power(kp0 / (kp0 + x), kxn);
  else
    y = kb* TMath::Exp(-sqrt(x * x + xmpi2) / kt);
  
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
    :AliGenerator(),
	fNt(-1),
	fNpartProd(0),
	fPi0Decays(kFALSE),
	fPtWgtPi(0.),
	fPtWgtKa(0.),
	fPtpi(0),
	fPtka(0),
	fETApic(0),
	fETAkac(0),
	fDecayer(0)
{
    //
    // Default constructor
    //
    SetCutVertexZ();
    SetPtRange();
}

//_____________________________________________________________________________
AliGenHIJINGpara::AliGenHIJINGpara(Int_t npart)
    :AliGenerator(npart),
	fNt(-1),
	fNpartProd(npart),
	fPi0Decays(kFALSE),
	fPtWgtPi(0.),
	fPtWgtKa(0.),
	fPtpi(0),
	fPtka(0),
	fETApic(0),
	fETAkac(0),
	fDecayer(0)
{
  // 
  // Standard constructor
  //
    fName="HIJINGpara";
    fTitle="HIJING Parametrisation Particle Generator";
    SetCutVertexZ();
    SetPtRange();
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
    fPtpi   = new TF1("ptpi",&ptpi,0,20,0);
    gROOT->GetListOfFunctions()->Remove(fPtpi);
    fPtka   = new TF1("ptka",&ptka,0,20,0);
    gROOT->GetListOfFunctions()->Remove(fPtka);
    fPtpi->SetNpx(1000);
    fPtka->SetNpx(1000);
    fETApic = new TF1("etapic",&etapic,etaMin,etaMax,0);
    gROOT->GetListOfFunctions()->Remove(fETApic);
    fETAkac = new TF1("etakac",&etakac,etaMin,etaMax,0);
    gROOT->GetListOfFunctions()->Remove(fETAkac);

    TF1 etaPic0("etaPic0",&etapic,-7,7,0);
    TF1 etaKac0("etaKac0",&etakac,-7,7,0);

    TF1 ptPic0("ptPic0",&ptpi,0.,15.,0);
    TF1 ptKac0("ptKac0",&ptka,0.,15.,0);

    Float_t intETApi  = etaPic0.Integral(-0.5, 0.5);
    Float_t intETAka  = etaKac0.Integral(-0.5, 0.5);
    Float_t scalePi   = 7316/(intETApi/1.5);
    Float_t scaleKa   =  684/(intETAka/2.0);

//  Fraction of events corresponding to the selected pt-range    
    Float_t intPt    = (0.877*ptPic0.Integral(0, 15)+
			0.123*ptKac0.Integral(0, 15));
    Float_t intPtSel = (0.877*ptPic0.Integral(fPtMin, fPtMax)+
			0.123*ptKac0.Integral(fPtMin, fPtMax));
    Float_t ptFrac   = intPtSel/intPt;

//  Fraction of events corresponding to the selected eta-range    
    Float_t intETASel  = (scalePi*etaPic0.Integral(etaMin, etaMax)+
			  scaleKa*etaKac0.Integral(etaMin, etaMax));
//  Fraction of events corresponding to the selected phi-range    
    Float_t phiFrac    = (fPhiMax-fPhiMin)/2/TMath::Pi();

    
    fParentWeight = (intETASel*ptFrac*phiFrac) / Float_t(fNpart);
    
    if (fAnalog != 0) {
	fPtWgtPi = (fPtMax - fPtMin) / fPtpi->Integral(0., 20.);
	fPtWgtKa = (fPtMax - fPtMin) / fPtka->Integral(0., 20.);
	fParentWeight = (intETASel*phiFrac) / Float_t(fNpart);
    }
    
    
    AliInfo(Form("The number of particles in the selected kinematic region corresponds to %f percent of a full event", 
		 100./ fParentWeight));

// Issue warning message if etaMin or etaMax are outside the alowed range 
// of the parametrization
    if (etaMin < -8.001 || etaMax > 8.001) {
	AliWarning("\nYOU ARE USING THE PARAMETERISATION OUTSIDE ");	
	AliWarning("THE ALLOWED PSEUDORAPIDITY RANGE (-8. - 8.)");	    
	AliWarning(Form("YOUR LIMITS: %f %f \n ", etaMin, etaMax));
    }
//
//
    if (fPi0Decays && gMC)
	fDecayer = gMC->GetDecayer();
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
    Float_t pt, pl, ptot, wgt;
    Float_t phi, theta;
    Float_t p[3];
    Int_t i, part, j;
    //
    TF1 *ptf;
    TF1 *etaf;
    //
    Float_t random[6];
    //
    for (j=0;j<3;j++) origin[j]=fOrigin[j];

    if(fVertexSmear == kPerEvent) {
	Vertex();
	for (j=0; j < 3; j++) origin[j] = fVertex[j];
    } // if kPerEvent
    TArrayF eventVertex;
    eventVertex.Set(3);
    eventVertex[0] = origin[0];
    eventVertex[1] = origin[1];
    eventVertex[2] = origin[2];
     
    for(i=0;i<fNpart;i++) {
	while(1) {
	    Rndm(random,4);
	    if(random[0]<kBorne) {
		part=kPions[Int_t (random[1]*3)];
		ptf=fPtpi;
		etaf=fETApic;
		wgt = fPtWgtPi;
	    } else {
		part=kKaons[Int_t (random[1]*4)];
		ptf=fPtka;
		etaf=fETAkac;
		wgt = fPtWgtKa;
	    }
	    phi=fPhiMin+random[2]*(fPhiMax-fPhiMin);
	    theta=2*TMath::ATan(TMath::Exp(-etaf->GetRandom()));
	    if(theta<fThetaMin || theta>fThetaMax) continue;
	 
	    if (fAnalog == 0) { 
		pt = ptf->GetRandom();
	    } else {
		pt = fPtMin + random[3] * (fPtMax - fPtMin);
	    }
	    
	    
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

	    if (fAnalog == 0) { 
		wgt = fParentWeight;
	    } else {
		wgt *= (fParentWeight * ptf->Eval(pt));
	    }
	    
	    
	    if (part == kPi0 && fPi0Decays){
//
//          Decay pi0 if requested
		PushTrack(0,-1,part,p,origin,polar,0,kPPrimary,fNt,wgt);
		KeepTrack(fNt);
		DecayPi0(origin, p);
	    } else {
      // printf("fNt %d", fNt);
		PushTrack(fTrackIt,-1,part,p,origin,polar,0,kPPrimary,fNt,wgt);

		KeepTrack(fNt);
	    }

	    break;
	}
	SetHighWaterMark(fNt);
    }
//

// Header
    AliGenEventHeader* header = new AliGenEventHeader("HIJINGparam");
// Event Vertex
    header->SetPrimaryVertex(eventVertex);
    header->SetNProduced(fNpartProd);
    gAlice->SetGenEventHeader(header); 
}

void AliGenHIJINGpara::SetPtRange(Float_t ptmin, Float_t ptmax) {
    AliGenerator::SetPtRange(ptmin, ptmax);
}

void AliGenHIJINGpara::DecayPi0(Float_t* orig, Float_t * p) 
{
//
//    Decay the pi0
//    and put decay products on the stack
//
    static TClonesArray *particles;
    if(!particles) particles = new TClonesArray("TParticle",1000);
//    
    const Float_t kMass = TDatabasePDG::Instance()->GetParticle(kPi0)->Mass();
    Float_t       e     = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]+ kMass * kMass);
//
//  Decay the pi0    
    TLorentzVector pmom(p[0], p[1], p[2], e);
    fDecayer->Decay(kPi0, &pmom);
    
//
// Put decay particles on the stack
//
    Float_t polar[3] = {0., 0., 0.};
    Int_t np = fDecayer->ImportParticles(particles);
    fNpartProd += (np-1);
    Int_t nt;    
    for (Int_t i = 1; i < np; i++)
    {
	TParticle* iParticle =  (TParticle *) particles->At(i);
	p[0] = iParticle->Px();
	p[1] = iParticle->Py();
	p[2] = iParticle->Pz();
	Int_t part = iParticle->GetPdgCode();

	PushTrack(fTrackIt, fNt, part, p, orig, polar, 0, kPDecay, nt, fParentWeight);
	KeepTrack(nt);
    }
    fNt = nt;
}

void AliGenHIJINGpara::Draw( const char * /*opt*/)
{
    //
    // Draw the pT and y Distributions
    //
     TCanvas *c0 = new TCanvas("c0","Canvas 0",400,10,600,700);
     c0->Divide(2,1);
     c0->cd(1);
     fPtpi->Draw();
     fPtpi->GetHistogram()->SetXTitle("p_{T} (GeV)");     
     c0->cd(2);
     fPtka->Draw();
     fPtka->GetHistogram()->SetXTitle("p_{T} (GeV)");     

}
