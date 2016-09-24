/**************************************************************************
 * Copyright(c) 2009-2016, ALICE Experiment at CERN, All rights reserved. *
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
// whose momenta in their CM frame is less than p0 (coalescence momentum).
//
// Sample code for PYTHIA:
//
//    AliGenCocktail* gener = new AliGenCocktail();
//
//    AliGenPythia* pythia = new AliGenPythia(-1);
//    pythia->SetCollisionSystem("p+", "p+");
//    pythia->SetEnergyCMS(7000);
//
//    AliGenLightNuclei* aft = new AliGenLightNuclei();
//    aft->SetNucleusPdgCode(AliGenLightNuclei::kDeuteron); // default
//    aft->SetCoalescenceMomentum(0.100); // default (GeV/c)
//    
//    gener->AddGenerator(pythia, "PYTHIA", 1);
//    gener->AddGenerator(aft, "deuteron", 1);
//
// Notice that the order in which the afterburner is added is
// important and more than one afterburner can be added
// to generate different nucleus species.
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
#include "TMath.h"
#include "TLorentzVector.h"
#include "AliGenLightNuclei.h"

ClassImp(AliGenLightNuclei)

AliGenLightNuclei::AliGenLightNuclei()
:AliGenerator()
,fPdg(kDeuteron)
,fP0(0.100)
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
// modify current stack adding light nuclei
//
	if(fStack == 0)
	{
		  AliRunLoader* rl = AliRunLoader::Instance();
		  if(rl != 0)  fStack = rl->Stack();
	}
	
	if(fStack == 0) AliFatal("no stack");
	
	// find nucleons and anti-nucleons
	TList* protons      = new TList();
	TList* neutrons     = new TList();
	TList* lambdas      = new TList();
	TList* antiprotons  = new TList();
	TList* antineutrons = new TList();
	TList* antilambdas  = new TList();
	
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
			case kLambda0:
				lambdas->Add(iParticle);
				break;
			case kProtonBar:
				antiprotons->Add(iParticle);
				break;
			case kNeutronBar:
				antineutrons->Add(iParticle);
				break;
			case kLambda0Bar:
				antilambdas->Add(iParticle);
				break;
			default:
				break;
		}
	}
	
	// do not delete content
	protons->SetOwner(kFALSE);
	neutrons->SetOwner(kFALSE);
	lambdas->SetOwner(kFALSE);
	antiprotons->SetOwner(kFALSE);
	antineutrons->SetOwner(kFALSE);
	antilambdas->SetOwner(kFALSE);
	
	if(TMath::Abs(fPdg) == kDeuteron)
	{
		this->GenerateNuclei(kDeuteron, 1.87561282, protons, neutrons);
		this->GenerateNuclei(-kDeuteron, 1.87561282, antiprotons, antineutrons);
	}
	else if(TMath::Abs(fPdg) == kLambdaN)
	{
		this->GenerateNuclei(kLambdaN, 2.054, lambdas, neutrons);
		this->GenerateNuclei(-kLambdaN, 2.054, antilambdas, antineutrons);
	}
	else if(TMath::Abs(fPdg) == kHDibarion)
	{
		this->GenerateNuclei(kHDibarion, 2.231, lambdas, lambdas);
		this->GenerateNuclei(-kHDibarion, 2.231, antilambdas, antilambdas);
	}
	else if(TMath::Abs(fPdg) == kTriton)
	{
		this->GenerateNuclei(kTriton, 2.80925, protons, neutrons, neutrons);
		this->GenerateNuclei(-kTriton, 2.80925, antiprotons, antineutrons, antineutrons);
	}
	else if(TMath::Abs(fPdg) == kHyperTriton )
	{
		this->GenerateNuclei(kHyperTriton, 2.99131, protons, neutrons, lambdas);
		this->GenerateNuclei(-kHyperTriton, 2.99131, antiprotons, antineutrons, antilambdas);
	}
	else if(TMath::Abs(fPdg) == kLambdaNN )
	{
		this->GenerateNuclei(kLambdaNN, 2.982, lambdas, neutrons, neutrons);
		this->GenerateNuclei(-kLambdaNN, 2.982, antilambdas, antineutrons, antineutrons);
	}
	else if(TMath::Abs(fPdg) == kLambdaLN )
	{
		this->GenerateNuclei(kLambdaLN, 3.171, lambdas, lambdas, neutrons);
		this->GenerateNuclei(-kLambdaLN, 3.171, antilambdas, antilambdas, antineutrons);
	}
	else if(TMath::Abs(fPdg) == kLambdaLP )
	{
		this->GenerateNuclei(kLambdaLP, 3.17, lambdas, lambdas, protons);
		this->GenerateNuclei(-kLambdaLP, 3.17, antilambdas, antilambdas, antiprotons);
	}
	else if(TMath::Abs(fPdg) == kHe3Nucleus)
	{
		this->GenerateNuclei(kHe3Nucleus, 2.80923, protons, neutrons, protons);
		this->GenerateNuclei(-kHe3Nucleus, 2.80923, antiprotons, antineutrons, antiprotons);
	}
	else if(TMath::Abs(fPdg) == kAlpha )
	{
		this->GenerateNuclei(kAlpha, 3.727417, protons, neutrons, protons, neutrons);
		this->GenerateNuclei(-kAlpha, 3.727417, antiprotons, antineutrons, antiprotons, antineutrons);
	}
	else
	{
		AliFatal(Form("Unknown nucleus PDG: %d", fPdg));
	}
	
	protons->Clear("nodelete");
	neutrons->Clear("nodelete");
	lambdas->Clear("nodelete");
	antiprotons->Clear("nodelete");
	antineutrons->Clear("nodelete");
	antilambdas->Clear("nodelete");
	
	delete protons;
	delete neutrons;
	delete lambdas;
	delete antiprotons;
	delete antineutrons;
	delete antilambdas;
}

Bool_t AliGenLightNuclei::Coalescence(const TLorentzVector& p1, const TLorentzVector& p2) const
{
//
// returns true if the 2 nucleons have momentum < p0 in their CM frame
//
	TLorentzVector p1cm(p1);
	TLorentzVector p2cm(p2);
	
	TLorentzVector p = p1 + p2;
	TVector3 b = -p.BoostVector();
	
	p1cm.Boost(b);
	p2cm.Boost(b);
	
	return (p1cm.Vect().Mag() < fP0) && (p2cm.Vect().Mag() < fP0);
}

Bool_t AliGenLightNuclei::Coalescence(const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3) const
{
//
// returns true if the 3 nucleons have momentum < p0 in their CM frame
//
	TLorentzVector p1cm(p1);
	TLorentzVector p2cm(p2);
	TLorentzVector p3cm(p3);
	
	TLorentzVector p = p1 + p2 + p3;
	TVector3 b = -p.BoostVector();
	
	p1cm.Boost(b);
	p2cm.Boost(b);
	p3cm.Boost(b);
	
	return (p1cm.Vect().Mag() < fP0) && (p2cm.Vect().Mag() < fP0) && (p3cm.Vect().Mag() < fP0);
}

Bool_t AliGenLightNuclei::Coalescence(const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& p4) const
{
//
// returns true if the 4 nucleons have momentum < p0 in their CM frame
//
	TLorentzVector p1cm(p1);
	TLorentzVector p2cm(p2);
	TLorentzVector p3cm(p3);
	TLorentzVector p4cm(p4);
	
	TLorentzVector p = p1 + p2 + p3 + p4;
	TVector3 b = -p.BoostVector();
	
	p1cm.Boost(b);
	p2cm.Boost(b);
	p3cm.Boost(b);
	p4cm.Boost(b);
	
	return (p1cm.Vect().Mag() < fP0) && (p2cm.Vect().Mag() < fP0) && (p3cm.Vect().Mag() < fP0) && (p4cm.Vect().Mag() < fP0);
}

Int_t AliGenLightNuclei::GenerateNuclei(Int_t pdg, Double_t mass, const TList* l1, const TList* l2)
{
//
// a nucleus with A=2 is generated from the first n1-n2 nucleon cluster
// that fulfill the coalescence condition
//
	Int_t npart = 0;
	
	TIter n1iter(l1);
	while(TParticle* n1 = dynamic_cast<TParticle*>(n1iter()))
	{
		if(n1 == 0) continue;
		if(n1->GetStatusCode() == kCluster) continue;
		
		TLorentzVector p1(n1->Px(), n1->Py(), n1->Pz(), n1->Energy());
		
		TIter n2iter(l2);
		if(l2 == l1) n2iter = n1iter;
		
		while(TParticle* n2 = dynamic_cast<TParticle*>( n2iter()))
		{
			if(n2 == 0) continue;
			if(n2 == n1) continue;
			if(n2->GetStatusCode() == kCluster) continue;
			
			TLorentzVector p2(n2->Px(), n2->Py(), n2->Pz(), n2->Energy());
			
			if(!this->Coalescence(p1, p2)) continue;
			
			this->PushNucleus(pdg, mass, n1, n2);
			
			++npart;
			
			break;
		}
	}
	
	return npart;
}

Int_t AliGenLightNuclei::GenerateNuclei(Int_t pdg, Double_t mass, const TList* l1, const TList* l2, const TList* l3)
{
//
// a nucleus with A=3 is generated from the first n1-n2-n3 nucleon cluster
// that fulfill the coalescence condition
//
	Int_t npart = 0;
	
	TIter n1iter(l1);
	while(TParticle* n1 = dynamic_cast<TParticle*>(n1iter()))
	{
		if(n1 == 0) continue;
		if(n1->GetStatusCode() == kCluster) continue;
		
		TLorentzVector p1(n1->Px(), n1->Py(), n1->Pz(), n1->Energy());
		
		TIter n2iter(l2);
		if(l2 == l1) n2iter = n1iter;
		
		while(TParticle* n2 = dynamic_cast<TParticle*>( n2iter()) )
		{
			if(n2 == 0) continue;
			if(n2 == n1) continue;
			if(n2->GetStatusCode() == kCluster) continue;
			
			TLorentzVector p2(n2->Px(), n2->Py(), n2->Pz(), n2->Energy());
			
			TIter n3iter(l3);
			if(l3 == l1) n3iter = n1iter;
			if(l3 == l2) n3iter = n2iter;
			
			while(TParticle* n3 = dynamic_cast<TParticle*>( n3iter()) )
			{
				if(n3 == 0) continue;
				if(n3 == n1) continue;
				if(n3 == n2) continue;
				if(n3->GetStatusCode() == kCluster) continue;
				
				TLorentzVector p3(n3->Px(), n3->Py(), n3->Pz(), n3->Energy());
				
				if(!this->Coalescence(p1, p2, p3)) continue;
				
				this->PushNucleus(pdg, mass, n1, n2, n3);
				
				++npart;
				
				break;
			}
			
			if(n2->GetStatusCode() == kCluster) break;
		}
	}
	
	return npart;
}

Int_t AliGenLightNuclei::GenerateNuclei(Int_t pdg, Double_t mass, const TList* l1, const TList* l2, const TList* l3, const TList* l4)
{
//
// a nucleus with A=4 is generated from the first n1-n2-n3-n4 nucleon cluster
// that fulfill the coalescence condition
//
	Int_t npart = 0;
	
	TIter n1iter(l1);
	while(TParticle* n1 = dynamic_cast<TParticle*>(n1iter()))
	{
		if(n1 == 0) continue;
		if(n1->GetStatusCode() == kCluster) continue;
		
		TLorentzVector p1(n1->Px(), n1->Py(), n1->Pz(), n1->Energy());
		
		TIter n2iter(l2);
		if(l2 == l1) n2iter = n1iter;
		
		while(TParticle* n2 = dynamic_cast<TParticle*>( n2iter()))
		{
			if(n2 == 0) continue;
			if(n2->GetStatusCode() == kCluster) continue;
			
			TLorentzVector p2(n2->Px(), n2->Py(), n2->Pz(), n2->Energy());
			
			TIter n3iter(l3);
			if(l3 == l1) n3iter = n1iter;
			if(l3 == l2) n3iter = n2iter;
			
			while(TParticle* n3 = dynamic_cast<TParticle*>( n3iter()))
			{
				if(n3 == 0) continue;
				if(n3 == n1) continue;
				if(n3 == n2) continue;
				if(n3->GetStatusCode() == kCluster) continue;
				
				TLorentzVector p3(n3->Px(), n3->Py(), n3->Pz(), n3->Energy());
				
				TIter n4iter(l4);
				if(l4 == l1) n4iter = n1iter;
				if(l4 == l2) n4iter = n2iter;
				if(l4 == l3) n4iter = n3iter;
			
				while(TParticle* n4 = dynamic_cast<TParticle*>( n4iter()))
				{
					if(n4 == 0) continue;
					if(n4 == n1) continue;
					if(n4 == n2) continue;
					if(n4 == n3) continue;
					if(n4->GetStatusCode() == kCluster) continue;
					
					TLorentzVector p4(n4->Px(), n4->Py(), n4->Pz(), n4->Energy());
					
					if(!this->Coalescence(p1, p2, p3, p4)) continue;
					
					this->PushNucleus(pdg, mass, n1, n2, n3, n4);
					
					++npart;
					
					break;
				}
				
				if(n3->GetStatusCode() == kCluster) break;
			}
			
			if(n2->GetStatusCode() == kCluster) break;
		}
	}
	
	return npart;
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
