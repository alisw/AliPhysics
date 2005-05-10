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

//
//  Generator for slow nucleons in pA interactions. 
//  Source is modelled by a relativistic Maxwell distributions.
//  This class cooparates with AliCollisionGeometry if used inside AliGenCocktail.
//  In this case the number of slow nucleons is determined from the number of wounded nuclei
//  using a realisation of AliSlowNucleonModel.
//  Original code by  Ferenc Sikler  <sikler@rmki.kfki.hu>
// 

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "AliCollisionGeometry.h"
#include "AliGenSlowNucleons.h"
#include "AliSlowNucleonModel.h"

 ClassImp(AliGenSlowNucleons)
    
 AliGenSlowNucleons::AliGenSlowNucleons():AliGenerator(-1)
{
// Default constructor
    fSlowNucleonModel = 0;
    fCollisionGeometry = 0;
}

AliGenSlowNucleons::AliGenSlowNucleons(Int_t npart)
    :AliGenerator(npart)
{
// Constructor
    fName  = "SlowNucleons";
    fTitle = "Generator for gray particles in pA collisions";
    SetPmax();
    SetTarget();
    SetNominalCmsEnergy();
    SetCharge();
    SetTemperature();
    SetBetaSource();
    fSlowNucleonModel = new AliSlowNucleonModel();
    fCollisionGeometry = 0;
    fDebug = 0;
}

AliGenSlowNucleons::AliGenSlowNucleons(const AliGenSlowNucleons & sn):
    AliGenerator(sn)
{
// Copy constructor
    sn.Copy(*this);
}

//____________________________________________________________
AliGenSlowNucleons::~AliGenSlowNucleons()
{
// Destructor
    delete  fSlowNucleonModel;
}

void AliGenSlowNucleons::SetProtonDirection(Float_t dir) {
// Set direction of the proton to change between pA (1) and Ap (-1)
  fProtonDirection = dir / TMath::Abs(dir);
}

void AliGenSlowNucleons::Init()
{
  //
  // Initialization
  //
    Float_t kMass  = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
    fMomentum = fCMS/2. * fZTarget / fATarget;
    fBeta     = fMomentum / TMath::Sqrt(kMass * kMass + fMomentum * fMomentum);
    if (fDebug) {
	fDebugHist1 = new TH2F("DebugHist1", "nu vs N_slow", 100, 0., 100., 20, 0., 20.);
	fDebugHist2 = new TH2F("DebugHist2", "b  vs N_slow", 100, 0., 100., 15, 0., 15.);
    }
}

void AliGenSlowNucleons::FinishRun()
{
// End of run action
// Show histogram for debugging if requested.
    if (fDebug) {
	TCanvas *c = new TCanvas("c","Canvas 1",400,10,600,700);
	c->Divide(2,1);
	c->cd(1);
	fDebugHist1->Draw();
	c->cd(2);
	fDebugHist2->Draw();
    }
}


void AliGenSlowNucleons::Generate()
{
  //
  // Generate one event
  //
  //
  // Communication with Gray Particle Model 
  // 
    if (fCollisionGeometry) {
	Float_t b   = fCollisionGeometry->ImpactParameter();
	Int_t  nn   = fCollisionGeometry->NN();
	Int_t  nwn  = fCollisionGeometry->NwN();
	Int_t  nnw  = fCollisionGeometry->NNw();
	Int_t  nwnw = fCollisionGeometry->NwNw();

	fSlowNucleonModel->GetNumberOfSlowNucleons(fCollisionGeometry, fNgp, fNgn, fNbp, fNbn);
	if (fDebug) {
	    printf("Nucleons %d %d %d %d \n", fNgp, fNgn, fNbp, fNbn);
	    fDebugHist1->Fill(Float_t(fNgp + fNgn + fNbp + fNbn), fCollisionGeometry->NwN(), 1.);
	    fDebugHist2->Fill(Float_t(fNgp + fNgn + fNbp + fNbn), b, 1.);
	    printf("AliGenSlowNucleons: Impact parameter from Collision Geometry %f %d %d %d %d\n", 
		   b, nn, nwn, nnw, nwnw);
	}
    }     

   //
    Float_t p[3];
    Float_t origin[3] = {0., 0., 0.};
    Float_t polar [3] = {0., 0., 0.};    
    Int_t nt, i, j;
    Int_t kf;
    
    if(fVertexSmear == kPerEvent) {
	Vertex();
	for (j=0; j < 3; j++) origin[j] = fVertex[j];
    } // if kPerEvent
//
//  Gray protons
//
    fCharge = 1;
    kf = kProton;    
    for(i = 0; i < fNgp; i++) {
	GenerateSlow(fCharge, fTemperatureG, fBetaSourceG, p);
	PushTrack(fTrackIt, -1, kf, p, origin, polar,
		 0., kPNoProcess, nt, 1.);
	KeepTrack(nt);
    }
//
//  Gray neutrons
//
    fCharge = 0;
    kf = kNeutron;    
    for(i = 0; i < fNgn; i++) {
	GenerateSlow(fCharge, fTemperatureG, fBetaSourceG, p);
	PushTrack(fTrackIt, -1, kf, p, origin, polar,
		 0., kPNoProcess, nt, 1.);
	KeepTrack(nt);
    }
//
//  Black protons
//
    fCharge = 1;
    kf = kProton;    
    for(i = 0; i < fNbp; i++) {
	GenerateSlow(fCharge, fTemperatureB, fBetaSourceB, p);
	PushTrack(fTrackIt, -1, kf, p, origin, polar,
		 0., kPNoProcess, nt, 1.);
	KeepTrack(nt);
    }
//
//  Black neutrons
//
    fCharge = 0;
    kf = kNeutron;    
    for(i = 0; i < fNbn; i++) {
	GenerateSlow(fCharge, fTemperatureB, fBetaSourceB, p);
	PushTrack(fTrackIt, -1, kf, p, origin, polar,
		 0., kPNoProcess, nt, 1.);
	KeepTrack(nt);
    }
}




void AliGenSlowNucleons::GenerateSlow(Int_t charge, Double_t T, Double_t beta, Float_t* q)

{
/* 
   Emit a slow nucleon with "temperature" T [GeV], 
   from a source moving with velocity beta         
   Three-momentum [GeV/c] is given back in q[3]    
*/

 Double_t m, pmax, p, f, theta, phi;
 TDatabasePDG * pdg = TDatabasePDG::Instance();
 const Double_t kMassProton  = pdg->GetParticle(kProton) ->Mass();
 const Double_t kMassNeutron = pdg->GetParticle(kNeutron)->Mass();
 
 /* Select nucleon type */
 if(charge == 0) m = kMassNeutron;
 else m = kMassProton;

 /* Momentum at maximum of Maxwell-distribution */

 pmax = TMath::Sqrt(2*T*(T+sqrt(T*T+m*m)));

 /* Try until proper momentum                                  */
 /* for lack of primitive function of the Maxwell-distribution */
 /* a brute force trial-accept loop, normalized at pmax        */

 do
 {
     p = Rndm() * fPmax;
     f = Maxwell(m, p, T) / Maxwell(m , pmax, T);
 }
 while(f < Rndm());

 /* Spherical symmetric emission */
 theta = TMath::ACos(2. * Rndm() - 1.);
 phi   = 2. * TMath::Pi() * Rndm();

 /* Determine momentum components in system of the moving source */
 q[0] = p * TMath::Sin(theta) * TMath::Cos(phi);
 q[1] = p * TMath::Sin(theta) * TMath::Sin(phi);
 q[2] = p * TMath::Cos(theta);

 /* Transform to system of the target nucleus                             */
 /* beta is passed as negative, because the gray nucleons are slowed down */
 Lorentz(m, -beta, q);

 /* Transform to laboratory system */
 Lorentz(m, fBeta, q);
 q[2] *= fProtonDirection; 
}

Double_t AliGenSlowNucleons::Maxwell(Double_t m, Double_t p, Double_t T)
{
/* Relativistic Maxwell-distribution */
    Double_t ekin;
    ekin = TMath::Sqrt(p*p+m*m)-m;
    return (p*p * exp(-ekin/T));
}


void AliGenSlowNucleons::Lorentz(Double_t m, Double_t beta, Float_t* q)
{
/* Lorentz transform in the direction of q[2] */
 
    Double_t gamma  = 1/sqrt(1-beta*beta); 
    Double_t energy = sqrt(m*m + q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    q[2] = gamma * (q[2] + beta*energy);
}

	  
AliGenSlowNucleons& AliGenSlowNucleons::operator=(const  AliGenSlowNucleons& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliGenSlowNucleons::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}







