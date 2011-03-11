/**************************************************************************
 * Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
// Afterburner to generate (anti)deuterons simulating coalescence of
// (anti)nucleons as in A. J. Baltz et al., Phys. lett B 325(1994)7.
// If the nucleon generator does not provide a spatial description of 
// the source the afterburner can provide one.
//
// There two options for the source: a thermal picture where all nucleons
// are placed randomly and homogeneously in an spherical volume, or
// an expansion one where they are projected onto its surface.
//
// An (anti)deuteron will form if there is a pair of (anti)proton-(anti)neutron
// with momentum difference less than ~ 300MeV and relative distance less than
// ~ 2.1fm. Only 3/4 of these clusters are accepted due to the deuteron spin.
//
// When there are more than one pair fulfilling the coalescence conditions,
// the selected pair can be the one with the first partner, with
// the lowest relative momentum, the lowest relative distance or both.
//////////////////////////////////////////////////////////////////////////

// Author: Eulogio Serradilla <eulogio.serradilla@ciemat.es>,
//         Arturo Menchaca <menchaca@fisica.unam.mx>
//

#include "Riostream.h"
#include "TParticle.h"
#include "AliRun.h"
#include "AliStack.h"
#include "TMath.h"
#include "TMCProcess.h"
#include "TList.h"
#include "TVector3.h"
#include "AliMC.h"
#include "TArrayF.h"
#include "AliCollisionGeometry.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenCocktailEntry.h"
#include "AliGenCocktailAfterBurner.h"
#include "AliGenDeuteron.h"

ClassImp(AliGenDeuteron)

AliGenDeuteron::AliGenDeuteron(Int_t sign, Double_t pmax, Double_t rmax, Int_t cluster)
 :fSign(1)
 ,fPmax(pmax)
 ,fRmax(rmax)
 ,fSpinProb(0.75)
 ,fMaxRadius(1000.)
 ,fClusterType(cluster)
 ,fModel(0)
 ,fTimeLength(2.5)
 ,fB(0)
 ,fR(0)
 ,fPsiR(0)
 ,fCurStack(0)
 ,fNtrk(0)
{
//
// constructor
//
	fSign = sign > 0 ? 1:-1;
	
}

AliGenDeuteron::~AliGenDeuteron()
{
//
// Default destructor
//
}

void AliGenDeuteron::Init()
{
//
// Standard AliGenerator initializer
//
}

void AliGenDeuteron::Generate()
{
//
// Look for coalescence of (anti)nucleons in the nucleon source
// provided by the generator cocktail
//
	AliInfo(Form("sign: %d, Pmax: %g GeV/c, Rmax: %g fm, cluster: %d",fSign, fPmax, fRmax, fClusterType));
	if(fModel!=kNone) AliInfo(Form("model: %d, time: %g fm/c", fModel, fTimeLength));
	
	AliGenCocktailAfterBurner* gener = (AliGenCocktailAfterBurner*)gAlice->GetMCApp()->Generator();
	
	// find nuclear radius, only for first generator and projectile=target
	Bool_t collisionGeometry=0;
	if(fModel != kNone && gener->FirstGenerator()->Generator()->ProvidesCollisionGeometry())
	{
		TString name;
		Int_t ap, zp, at, zt;
		gener->FirstGenerator()->Generator()->GetProjectile(name,ap,zp);
		gener->FirstGenerator()->Generator()->GetTarget(name,at,zt);
		if(ap != at) AliWarning("projectile and target have different size");
		fR = 1.29*TMath::Power(at, 1./3.);
		collisionGeometry = 1;
	}
	
	for(Int_t ns = 0; ns < gener->GetNumberOfEvents(); ++ns)
	{
		gener->SetActiveEventNumber(ns);
		
		if(fModel != kNone && collisionGeometry)
		{
			fPsiR = (gener->GetCollisionGeometry(ns))->ReactionPlaneAngle();
			fB = (gener->GetCollisionGeometry(ns))->ImpactParameter();
			
			if(fB >= 2.*fR) continue; // no collision
		}
		
		// primary vertex
		TArrayF primVtx;
		gener->GetActiveEventHeader()->PrimaryVertex(primVtx);
		TVector3 r0(primVtx[0],primVtx[1],primVtx[2]);
		
		// emerging nucleons from the collision
		fCurStack = gener->GetStack(ns);
		if(!fCurStack)
		{
			AliWarning("no event stack");
			return;
		}
		
		TList* protons = new TList();
		protons->SetOwner(kFALSE);
		TList* neutrons = new TList();
		neutrons->SetOwner(kFALSE);
		
		for (Int_t i=0; i < fCurStack->GetNtrack(); ++i)
		{
			TParticle* iParticle = fCurStack->Particle(i);
			
			if(iParticle->TestBit(kDoneBit)) continue;
			
			// select only nucleons within the freeze-out volume
			TVector3 r(iParticle->Vx(),iParticle->Vy(),iParticle->Vz());
			if((r-r0).Mag() > fMaxRadius*1.e-13) continue;
			
			Int_t pdgCode = iParticle->GetPdgCode();
			if(pdgCode == fSign*2212)// (anti)proton
			{
				FixProductionVertex(iParticle);
				protons->Add(iParticle);
			}
			else if(pdgCode == fSign*2112) // (anti)neutron
			{
				FixProductionVertex(iParticle);
				neutrons->Add(iParticle);
			}
		}
		
		// look for clusters
		if(fClusterType==kFirstPartner)
		{
			FirstPartner(protons, neutrons);
		}
		else
		{
			WeightMatrix(protons, neutrons);
		}
		
		protons->Clear("nodelete");
		neutrons->Clear("nodelete");
		
		delete protons;
		delete neutrons;
	}
	
	AliInfo("DONE");
}

Double_t AliGenDeuteron::GetCoalescenceProbability(const TParticle* nucleon1, const TParticle* nucleon2) const
{
//
// Coalescence conditions as in
// A. J. Baltz et al., Phys. lett B 325(1994)7
//
// returns < 0 if coalescence is not possible
// otherwise returns a coalescence probability
//
	TVector3 v1(nucleon1->Vx(), nucleon1->Vy(), nucleon1->Vz());
	TVector3 p1(nucleon1->Px(), nucleon1->Py(), nucleon1->Pz());
	
	TVector3 v2(nucleon2->Vx(), nucleon2->Vy(), nucleon2->Vz());
	TVector3 p2(nucleon2->Px(), nucleon2->Py(), nucleon2->Pz());
	
	Double_t deltaP = (p2-p1).Mag();       // relative momentum
	if( deltaP >= fPmax)       return -1.;
	
	Double_t deltaR = (v2-v1).Mag();       // relative distance (cm)
	if(deltaR >= fRmax*1.e-13) return -1.;
	
	if(Rndm() > fSpinProb)    return -1.;  // spin
	
	if(fClusterType == kLowestMomentum) return 1. - deltaP/fPmax;
	if(fClusterType == kLowestDistance) return 1. - 1.e+13*deltaR/fRmax;
	
	return 1. - 1.e+13*(deltaP*deltaR)/(fRmax*fPmax);
}

void AliGenDeuteron::FirstPartner(const TList* protons, TList* neutrons)
{
//
// Clusters are made with the first nucleon pair that fulfill
// the coalescence conditions, starting with the protons
//
	TIter p_next(protons);
	while(TParticle* n0 = (TParticle*) p_next())
	{
		TParticle* partner = 0;
		TIter n_next(neutrons);
		while(TParticle* n1 = (TParticle*) n_next() )
		{
			if(GetCoalescenceProbability(n0, n1) < 0 ) continue; // with next neutron
			partner = n1;
			break;
		}
		
		if(partner == 0) continue; // with next proton
		
		PushDeuteron(n0, partner);
		
		// Remove from the list for the next iteration
		neutrons->Remove(partner);
	}
}

void AliGenDeuteron::WeightMatrix(const TList* protons, const TList* neutrons)
{
//
// Build all possible nucleon pairs with their own probability
// and select only those with the highest coalescence probability
//
	Int_t nMaxPairs = protons->GetSize()*neutrons->GetSize();
	
	TParticle** cProton = new TParticle*[nMaxPairs];
	TParticle** cNeutron = new TParticle*[nMaxPairs];
	Double_t*  cWeight = new Double_t[nMaxPairs];
	
	// build all pairs with probability > 0
	Int_t cIdx = -1;
	TIter p_next(protons);
	while(TParticle* n1 = (TParticle*) p_next())
	{
		TIter n_next(neutrons);
		while(TParticle* n2 = (TParticle*) n_next() )
		{
			Double_t weight = this->GetCoalescenceProbability(n1,n2);
			if(weight < 0) continue;
			++cIdx;
			cProton[cIdx]  = n1;
			cNeutron[cIdx] = n2;
			cWeight[cIdx]  = weight;
		}
		n_next.Reset();
	}
	p_next.Reset();
	
	Int_t nPairs = cIdx + 1;
	
	// find the interacting pairs:
	// remove repeated pairs and select only
	// the pair with the highest coalescence probability
	
	Int_t nMaxIntPair = TMath::Min(protons->GetSize(), neutrons->GetSize());
	
	TParticle** iProton = new TParticle*[nMaxIntPair];
	TParticle** iNeutron = new TParticle*[nMaxIntPair];
	Double_t* iWeight = new Double_t[nMaxIntPair];
	
	Int_t iIdx = -1;
	while(true)
	{
		Int_t j = -1;
		Double_t wMax = 0;
		for(Int_t i=0; i < nPairs; ++i)
		{
			if(cWeight[i] > wMax)
			{
				wMax=cWeight[i];
				j = i;
			}
		}
		
		if(j == -1 ) break; // end
		
		// Save the interacting pair
		++iIdx;
		iProton[iIdx]  = cProton[j];
		iNeutron[iIdx] = cNeutron[j];
		iWeight[iIdx]  = cWeight[j];
		
		// invalidate all combinations with these pairs for the next iteration
		for(Int_t i=0; i < nPairs; ++i)
		{
			if(cProton[i] == iProton[iIdx]) cWeight[i] = -1.;
			if(cNeutron[i] == iNeutron[iIdx]) cWeight[i] = -1.;
		}
	}
	
	Int_t nIntPairs = iIdx + 1;
	
	delete[] cProton;
	delete[] cNeutron;
	delete[] cWeight;
	
	// Add the (anti)deuterons to the current event stack
	for(Int_t i=0; i<nIntPairs; ++i)
	{
		TParticle* n1 = iProton[i];
		TParticle* n2 = iNeutron[i];
		PushDeuteron(n1,n2);
	}
	
	delete[] iProton;
	delete[] iNeutron;
	delete[] iWeight;
}

void AliGenDeuteron::PushDeuteron(TParticle* parent1, TParticle* parent2)
{
//
// Create an (anti)deuteron from parent1 and parent2,
// add to the current stack and set kDoneBit for the parents
//
	const Double_t kDeuteronMass = 1.87561282;
	const Int_t kDeuteronPdg = 1000010020;
	
	// momentum
	TVector3 p1(parent1->Px(), parent1->Py(), parent1->Pz());
	TVector3 p2(parent2->Px(), parent2->Py(), parent2->Pz());
	TVector3 pN = p1+p2;
	
	// production vertex same as the parent1's
	TVector3 vN(parent1->Vx(), parent1->Vy(), parent1->Vz());
	
	// E^2 = p^2 + m^2
	Double_t energy = TMath::Sqrt(pN.Mag2() + kDeuteronMass*kDeuteronMass);
	
	// Add a new (anti)deuteron to current event stack
	fCurStack->PushTrack(1, fCurStack->Particles()->IndexOf(parent1), fSign*kDeuteronPdg,
	                 pN.X(), pN.Y(), pN.Z(), energy,
	                 vN.X(), vN.Y(), vN.Z(), parent1->T(),
	                 0., 0., 0., kPNCapture, fNtrk, 1., 0);
	
	// Set kDoneBit for the parents
	parent1->SetBit(kDoneBit);
	parent2->SetBit(kDoneBit);
}

void AliGenDeuteron::FixProductionVertex(TParticle* i)
{
//
// Replace current generator nucleon spatial distribution
// with a custom distribution according to the selected model
//
	if(fModel == kNone || fModel > kExpansion) return;
	
	// semi-axis from collision geometry (fm)
	Double_t a = fTimeLength + TMath::Sqrt(fR*fR - fB*fB/4.);
	Double_t b = fTimeLength + fR - fB/2.;
	Double_t c = fTimeLength;
	
	Double_t xx = 0;
	Double_t yy = 0;
	Double_t zz = 0;
	
	if(fModel == kThermal)
	{
		// uniformly ditributed in the volume on an ellipsoid
		// random (r,theta,phi) unit sphere
		Double_t r = TMath::Power(Rndm(),1./3.);
		Double_t theta = TMath::ACos(2.*Rndm()-1.);
		Double_t phi = 2.*TMath::Pi()*Rndm();
		
		// transform coordenates
		xx = a*r*TMath::Sin(theta)*TMath::Cos(phi);
		yy = b*r*TMath::Sin(theta)*TMath::Sin(phi);
		zz = c*r*TMath::Cos(theta);
	}
	else if(fModel == kExpansion)
	{
		// project into the surface of an ellipsoid
		xx = a*TMath::Sin(i->Theta())*TMath::Cos(i->Phi());
		yy = b*TMath::Sin(i->Theta())*TMath::Sin(i->Phi());
		zz = c*TMath::Cos(i->Theta());
	}
	
	// rotate by the reaction plane angle
	Double_t x = xx*TMath::Cos(fPsiR)+yy*TMath::Sin(fPsiR);
	Double_t y = -xx*TMath::Sin(fPsiR)+yy*TMath::Cos(fPsiR);
	Double_t z = zz;
	
	// translate by the production vertex (cm)
	i->SetProductionVertex(i->Vx() + 1.e-13*x, i->Vy() + 1.e-13*y, i->Vz() + 1.e-13*z, i->T());
}

