/**************************************************************************
 * Copyright(c) 2009-2014, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////
//
// Afterburner to generate light nuclei for event generators
// such as PYTHIA and PHOJET
//
// Light nuclei are generated whenever a cluster of nucleons is found
// within a sphere of radius p0 (coalescence momentum), i.e. have the
// same momentum.
//
// By default it starts with He4 nuclei which are the most stable,
// then He3 nuclei, tritons and finally deuterons. It can also generate
// a single nucleus species by disabling the others.
//
// Sample code for PYTHIA:
//
//    AliGenLightNuclei* gener = new AliGenLightNuclei();
//
//    AliGenPythia* pythia = new AliGenPythia(-1);
//    pythia->SetCollisionSystem("p+", "p+");
//    pythia->SetEnergyCMS(7000);
//
//    gener->UsePerEventRates();
//    gener->AddGenerator(pythia, "PYTHIA", 1);
//    gener->SetCoalescenceMomentum(0.100); // default (GeV/c)
//
//////////////////////////////////////////////////////////////////////

//
// Author: Eulogio Serradilla <eulogio.serradilla@cern.ch>
//

#include "TMath.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TList.h"
#include "TVector3.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliRun.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliGenLightNuclei.h"

ClassImp(AliGenLightNuclei)

AliGenLightNuclei::AliGenLightNuclei()
:AliGenCocktail()
,fP0(0.100)
,fGenDeuterons(kTRUE)
,fGenTritons(kTRUE)
,fGenHe3Nuclei(kTRUE)
,fGenHe4Nuclei(kTRUE)
{
//
// default constructor
//
}

AliGenLightNuclei::~AliGenLightNuclei()
{
//
// default destructor
//
}

void AliGenLightNuclei::Generate()
{
//
// delegate the particle generation to the cocktail
// and modify the stack adding the light nuclei
//
	AliGenCocktail::Generate();
	
	// find the nucleons and anti-nucleons
	
	TList* protons      = new TList();
	TList* neutrons     = new TList();
	TList* antiprotons  = new TList();
	TList* antineutrons = new TList();
	
	for (Int_t i=0; i < fStack->GetNprimary(); ++i)
	{
		TParticle* iParticle = fStack->Particle(i);
		
		if(iParticle->GetStatusCode() != 1) continue;
		
		switch(iParticle->GetPdgCode())
		{
			case kProton:
				protons->Add(iParticle);
				break;
			case kNeutron:
				neutrons->Add(iParticle);
				break;
			case kProtonBar:
				antiprotons->Add(iParticle);
				break;
			case kNeutronBar:
				antineutrons->Add(iParticle);
				break;
			default:
				break;
		}
	}
	
	// do not delete content
	protons->SetOwner(kFALSE);
	neutrons->SetOwner(kFALSE);
	antiprotons->SetOwner(kFALSE);
	antineutrons->SetOwner(kFALSE);
	
	// first try with He4 nuclei which are the most stable
	
	if(fGenHe4Nuclei)
	{
		this->GenerateHe4Nuclei(protons, neutrons);
		this->GenerateHe4Nuclei(antiprotons, antineutrons);
	}
	
	// then He3 nuclei
	
	if(fGenHe3Nuclei)
	{
		this->GenerateHe3Nuclei(protons, neutrons);
		this->GenerateHe3Nuclei(antiprotons, antineutrons);
	}
	
	// then tritons
	
	if(fGenTritons)
	{
		this->GenerateTritons(protons, neutrons);
		this->GenerateTritons(antiprotons, antineutrons);
	}
	
	// and finally deuterons
	
	if(fGenDeuterons)
	{
		this->GenerateDeuterons(protons, neutrons);
		this->GenerateDeuterons(antiprotons, antineutrons);
	}
	
	protons->Clear("nodelete");
	neutrons->Clear("nodelete");
	antiprotons->Clear("nodelete");
	antineutrons->Clear("nodelete");
	
	delete protons;
	delete neutrons;
	delete antiprotons;
	delete antineutrons;
}

Bool_t AliGenLightNuclei::Coalescence(const TParticle* n1, const TParticle* n2) const
{
//
// returns true if the nucleons are inside of an sphere of radius p0
// (assume the nucleons are in the same place e.g. PYTHIA, PHOJET,...)
//
	Double_t deltaP = this->GetPcm(n1->Px(), n1->Py(), n1->Pz(), n1->GetMass(), 
	                               n2->Px(), n2->Py(), n2->Pz(), n2->GetMass());
	
	return ( deltaP < fP0);
}

TParticle* AliGenLightNuclei::FindPartner(const TParticle* n0, const TList* nucleons, const TParticle* nx) const
{
//
// find the first nucleon partner within a sphere of radius p0
// centered at n0 and exclude nucleon nx
//
	TIter n_next(nucleons);
	while(TParticle* n1 = dynamic_cast<TParticle*>( n_next()) )
	{
		if(n1 == 0) continue;
		if(n1 == nx) continue;
		if(n1->GetStatusCode() == kCluster) continue;
		if( !this->Coalescence(n0, n1) ) continue;
		return n1;
	}
	
	return 0;
}

Int_t AliGenLightNuclei::GenerateDeuterons(const TList* protons, const TList* neutrons)
{
//
// a deuteron is generated from a pair of p-n nucleons
// (the center of the sphere is one of the nucleons)
//
	Int_t npart = 0;
	
	TIter p_next(protons);
	while(TParticle* n0 = dynamic_cast<TParticle*>(p_next()) )
	{
		if(n0 == 0) continue;
		if(n0->GetStatusCode() == kCluster) continue;
		
		TParticle* partner = this->FindPartner(n0, neutrons, 0);
		
		if(partner == 0) continue;
		
		this->PushDeuteron(n0, partner);
		
		++npart;
	}
	
	return npart;
}

Int_t AliGenLightNuclei::GenerateTritons(const TList* protons, const TList* neutrons)
{
//
// a triton is generated from a cluster of p-n-n nucleons with same momentum
// (triangular configuration)
//
	Int_t npart = 0;
	
	TIter p_next(protons);
	while(TParticle* n0 = dynamic_cast<TParticle*>(p_next()) )
	{
		if(n0 == 0) continue;
		if(n0->GetStatusCode() == kCluster) continue;
		
		TParticle* partner1 = this->FindPartner(n0, neutrons, 0);
		
		if(partner1 == 0) continue;
		
		TParticle* partner2 = this->FindPartner(n0, neutrons, partner1);
		
		if(partner2 == 0) continue;
		
		// check that the partners coalesce between themselves
		
		if(!this->Coalescence(partner1, partner2)) continue;
		
		this->PushTriton(n0, partner1, partner2);
		
		++npart;
	}
	
	return npart;
}

Int_t AliGenLightNuclei::GenerateHe3Nuclei(const TList* protons, const TList* neutrons)
{
//
// a He3 nucleus is generated from a cluster of p-n-p nucleons with same momentum
// (triangular configuration)
//
	Int_t npart = 0;
	
	TIter p_next(protons); // center of the sphere
	while(TParticle* n0 = dynamic_cast<TParticle*>(p_next()) )
	{
		if(n0 == 0) continue;
		if(n0->GetStatusCode() == kCluster) continue;
		
		TParticle* partner1 = this->FindPartner(n0, neutrons, 0);
		
		if(partner1 == 0) continue;
		
		TParticle* partner2 = this->FindPartner(n0, protons, n0);
		
		if(partner2 == 0) continue;
		
		// check that the partners coalesce between themselves
		
		if(!this->Coalescence(partner1, partner2)) continue;
		
		this->PushHe3Nucleus(n0, partner1, partner2);
		
		++npart;
	}
	
	return npart;
}

Int_t AliGenLightNuclei::GenerateHe4Nuclei(const TList* protons, const TList* neutrons)
{
//
// a He4 nucleus is generated from a cluster of p-n-p-n nucleons with same momentum
// (tetrahedron configuration)
//
	Int_t npart = 0;
	
	TIter p_next(protons); // center of the sphere
	while(TParticle* n0 = dynamic_cast<TParticle*>(p_next()) )
	{
		if(n0 == 0) continue;
		if(n0->GetStatusCode() == kCluster) continue;
		
		TParticle* partner1 = this->FindPartner(n0, neutrons, 0);
		
		if(partner1 == 0) continue;
		
		TParticle* partner2 = this->FindPartner(n0, protons, n0);
		
		if(partner2 == 0) continue;
		
		TParticle* partner3 = this->FindPartner(n0, neutrons, partner1);
		
		if(partner3 == 0) continue;
		
		// check that the partners coalesce between themselves
		
		if(!this->Coalescence(partner1, partner2)) continue;
		if(!this->Coalescence(partner1, partner3)) continue;
		if(!this->Coalescence(partner2, partner3)) continue;
		
		this->PushHe4Nucleus(n0, partner1, partner2, partner3 );
		
		++npart;
	}
	
	return npart;
}

void AliGenLightNuclei::PushDeuteron(TParticle* parent1, TParticle* parent2)
{
//
// push a deuteron to the particle stack
//
	Int_t pdg = ( parent1->GetPdgCode() > 0 ) ? kDeuteron : kAntiDeuteron;
	this->PushNucleus(pdg, 1.87561282, parent1, parent2);
}

void AliGenLightNuclei::PushTriton(TParticle* parent1, TParticle* parent2, TParticle* parent3)
{
//
// push a triton to the particle stack
//
	Int_t pdg = ( parent1->GetPdgCode() > 0 ) ? kTriton : kAntiTriton;
	this->PushNucleus(pdg, 2.80925, parent1, parent2, parent3);
}

void AliGenLightNuclei::PushHe3Nucleus(TParticle* parent1, TParticle* parent2, TParticle* parent3)
{
//
// push a He3 nucleus to the particle stack
//
	Int_t pdg = ( parent1->GetPdgCode() > 0 ) ? kHe3Nucleus : kAntiHe3Nucleus;
	this->PushNucleus(pdg, 2.80923, parent1, parent2, parent3);
}

void AliGenLightNuclei::PushHe4Nucleus(TParticle* parent1, TParticle* parent2, TParticle* parent3, TParticle* parent4)
{
//
// push a He4 nucleus to the particle stack
//
	Int_t pdg = ( parent1->GetPdgCode() > 0 ) ? kAlpha : kAntiAlpha;
	this->PushNucleus(pdg, 3.727417, parent1, parent2, parent3, parent4);
}

void AliGenLightNuclei::PushNucleus(Int_t pdg, Double_t mass, TParticle* parent1, TParticle* parent2, TParticle* parent3, TParticle* parent4)
{
//
// push a nucleus to the stack and tag the parents with the kCluster status code
//
	Int_t ntrk;
	
	// momentum
	TVector3 p1(parent1->Px(), parent1->Py(), parent1->Pz());
	TVector3 p2(parent2->Px(), parent2->Py(), parent2->Pz());
	TVector3 p3(0, 0, 0);
	TVector3 p4(0, 0, 0);
	if(parent3 != 0) p3.SetXYZ(parent3->Px(), parent3->Py(), parent3->Pz());
	if(parent4 != 0) p4.SetXYZ(parent4->Px(), parent4->Py(), parent4->Pz());
	
	// momentum
	TVector3 pN = p1+p2+p3+p4;
	
	// E^2 = p^2 + m^2
	Double_t energy = TMath::Sqrt(pN.Mag2() + mass*mass);
	
	// production vertex same as the parent1'
	TVector3 vN(parent1->Vx(), parent1->Vy(), parent1->Vz());
	
	Double_t weight = 1;
	Int_t is = 1; // final state particle
	
	// add a new nucleus to current event stack
	fStack->PushTrack(1, -1, pdg,
	                 pN.X(), pN.Y(), pN.Z(), energy,
	                 vN.X(), vN.Y(), vN.Z(), parent1->T(),
	                 0., 0., 0., kPNCapture, ntrk, weight, is);
	
	// change the status code of the parents
	parent1->SetStatusCode(kCluster);
	parent2->SetStatusCode(kCluster);
	if(parent3 != 0) parent3->SetStatusCode(kCluster);
	if(parent4 != 0) parent4->SetStatusCode(kCluster);
	
	// set no transport for the parents
	parent1->SetBit(kDoneBit);
	parent2->SetBit(kDoneBit);
	if(parent3 != 0) parent3->SetBit(kDoneBit);
	if(parent4 != 0) parent4->SetBit(kDoneBit);
}

Double_t AliGenLightNuclei::GetS(Double_t p1x, Double_t p1y, Double_t p1z, Double_t m1, Double_t p2x, Double_t p2y, Double_t p2z, Double_t m2) const
{
//
// square of the center of mass energy
//
	Double_t E1 = TMath::Sqrt( p1x*p1x + p1y*p1y + p1z*p1z + m1*m1);
	Double_t E2 = TMath::Sqrt( p2x*p2x + p2y*p2y + p2z*p2z + m2*m2);
	
	return (E1+E2)*(E1+E2) - ((p1x+p2x)*(p1x+p2x) + (p1y+p2y)*(p1y+p2y) + (p1z+p2z)*(p1z+p2z));
}

Double_t AliGenLightNuclei::GetPcm(Double_t p1x, Double_t p1y, Double_t p1z, Double_t m1, Double_t p2x, Double_t p2y, Double_t p2z, Double_t m2) const
{
//
// momentum in the CM frame for 2 particles
//
	Double_t s = this->GetS(p1x, p1y, p1z, m1, p2x, p2y, p2z, m2);
	
	return TMath::Sqrt( (s-(m1-m2)*(m1-m2))*(s-(m1+m2)*(m1+m2)) )/(2.*TMath::Sqrt(s));
}
