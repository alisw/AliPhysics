#include "Riostream.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliGenDeuteron.h"
#include "AliGenCocktailAfterBurner.h"
#include "TMath.h"
#include "TMCProcess.h"
#include "TList.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "AliMC.h"
/*
A very simple model based on the article "Physics Letters B325 (1994) 7-12". 
The only addition is to model the freeze-out at which the
light nuclei can form for the energies of the LHC, since PYTHIA does not give
any detailed spatial description of the generated particles at the level of a
few fermis.

There are two options to place the nucleons: a thermal picture where all
nucleons are placed randomly and homogeneously in an spherical volume of a
few fermis, and expansion one where they are projected on its surface.

A (anti)deuteron will form if there is a pair of (anti)proton-(anti)neutron
with momentum difference less than ~ 300MeV and relative distance less than a
~ 2.1fm. Only 3/4 of this clusters are accepted due to the deuteron spin.
*/

ClassImp(AliGenDeuteron)


AliGenDeuteron::AliGenDeuteron(Int_t sign, Double_t pmax, Double_t rmax, Double_t rsrc, Int_t model):
fDeuteronMass(1.87561282), //pdg
fPmax(pmax), // GeV/c
fRmax(rmax*1.e-13), //cm
fRsrc(rsrc*1.e-13), // cm
fSpinProb(0.75),
fModel(model),
fSign(1)
{
	fSign = sign > 0 ? 1:-1;
}

AliGenDeuteron::~AliGenDeuteron()
{
}

void AliGenDeuteron::Init()
{
	// Standard AliGenerator Initializer
}

void AliGenDeuteron::Generate()
{
// Modify the stack of each event, adding the new particles at the end and
// not transporting the parents.

	Info("Generate","freeze-out model : %d (0 expansion, 1 thermal)",fModel);
	Info("Generate","relative momentum: %g GeV/c",fPmax);
	Info("Generate","relative distance: %g cm",fRmax);
	Info("Generate","source radius    : %g cm",fRsrc);
	Info("Generate","spin probability : %g ",fSpinProb);
	Info("Generate","sign             : %d ",fSign);
	
	TRandom3 rnd(0);
	TRandom3 rnd2(0);
	
	Int_t ntr;
	
	// Get the cocktail generator
	AliGenCocktailAfterBurner* gener = (AliGenCocktailAfterBurner*)gAlice->GetMCApp()->Generator();
	
	// Loop over events
	for(Int_t ns = 0; ns < gener->GetNumberOfEvents(); ++ns)
	{
		gener->SetActiveEventNumber(ns);
		
		AliStack* stack = gener->GetStack(ns); // Stack of event ns
		if(!stack)
		{
			Info("Generate", "no stack, exiting");
			return;
		}
		
		// find emerging nucleons from the collision of current event
		
		TList* protons = new TList();
		protons->SetOwner(kFALSE);
		TList* neutrons = new TList();
		neutrons->SetOwner(kFALSE);
		
		// FIXME: primary vertex
		//TVector3 r0(AliGenerator::fVertex[0],AliGenerator::fVertex[1],AliGenerator::fVertex[2]);
		
		// workaround for primary vertex
		TParticle* prim = stack->Particle(0);
		TVector3 r0(prim->Vx(),prim->Vy(),prim->Vz()); // primary vertex
		
		for (Int_t i=0; i<stack->GetNtrack(); ++i)
		{
			TParticle* iParticle = stack->Particle(i);
			
			if(iParticle->TestBit(kDoneBit)) continue;
			// select only nucleons within the freeze-out volume
			TVector3 r(iParticle->Vx(),iParticle->Vy(),iParticle->Vz());
			if((r-r0).Mag()>fRsrc) continue;
			
			Int_t pdgCode = iParticle->GetPdgCode();
			if(pdgCode == fSign*2212)// (anti)proton
			{
				FixProductionVertex(iParticle,rnd2);
				protons->Add(iParticle);
			}
			else if(pdgCode == fSign*2112) // (anti)neutron
			{
				FixProductionVertex(iParticle,rnd2);
				neutrons->Add(iParticle);
			}
		}
		
		// coalescence conditions
		// A.J. Baltz et al., Phys. lett B 325(1994)7
		
		TIter p_next(protons);
		while(TParticle* n0 = (TParticle*) p_next())
		{
			TParticle* iProton = 0;
			TParticle* jNeutron = 0;
			
			// Select the first nucleon
			iProton = n0;
			
			TVector3 v0(n0->Vx(), n0->Vy(), n0->Vz());
			TVector3 p0(n0->Px(), n0->Py(), n0->Pz()), p1(0,0,0);
			
			// See if there is a neibourgh...
			TIter n_next(neutrons);
			while(TParticle* n1 = (TParticle*) n_next() )
			{
				TVector3 v(n1->Vx(), n1->Vy(), n1->Vz());
				TVector3 p(n1->Px(), n1->Py(), n1->Pz());
				
				if((p-p0).Mag() > fPmax) continue;
				if((v-v0).Mag() > fRmax) continue;
				
				jNeutron = n1;
				p1 = p;
				break;
			}
			
			if(jNeutron == 0) continue; // with next proton
			
			if(rnd.Rndm() > fSpinProb) continue;
			
			// neutron captured!
			
			TVector3 pN = p0+p1;
			TVector3 vN(iProton->Vx(), iProton->Vy(), iProton->Vz()); // production vertex = p = n = collision vertex
			
			// E^2 = p^2 + m^2
			Double_t EN = TMath::Sqrt(pN.Mag2()+fDeuteronMass*fDeuteronMass);
			
			// Add a new (anti)deuteron to the current event stack
			stack->PushTrack(1, stack->Particles()->IndexOf(iProton), fSign*1000010020,
				         pN.X(), pN.Y(), pN.Z(), EN,
				         vN.X(), vN.Y(), vN.Z(), iProton->T(),
				         0., 0., 0., kPNCapture, ntr,1.,0);
			
			//Info("Generate","neutron capture (NEW DEUTERON)");
			
			// Set bit not to transport for the nucleons
			iProton->SetBit(kDoneBit);
			jNeutron->SetBit(kDoneBit);
			
			// Remove from the list
			neutrons->Remove(jNeutron);
		}
		
		protons->Clear("nodelete");
		neutrons->Clear("nodelete");
		
		delete protons;
		delete neutrons;
	}
	
	Info("Generate","Coalescence afterburner: DONE");
}

// create the freeze-out nucleon distribution around the collision vertex
void AliGenDeuteron::FixProductionVertex(TParticle* i, TRandom3& rnd)
{
	Double_t x,y,z;
	
	if(fModel == kThermal) // homogeneous volume
	{
		// random (r,theta,phi)
		Double_t r = fRsrc*TMath::Power(rnd.Rndm(),1./3.);
		Double_t theta = TMath::ACos(2.*rnd.Rndm()-1.);
		Double_t phi = 2.*TMath::Pi()*rnd.Rndm();
	
		// transform coordenates
		x = r*TMath::Sin(theta)*TMath::Cos(phi);
		y = r*TMath::Sin(theta)*TMath::Sin(phi);
		z = r*TMath::Cos(theta);
	}
	else // projection into the surface of an sphere
	{
		x = fRsrc*TMath::Sin(i->Theta())*TMath::Cos(i->Phi());
		y = fRsrc*TMath::Sin(i->Theta())*TMath::Sin(i->Phi());
		z = fRsrc*TMath::Cos(i->Theta());
	}
	
	// assume nucleons in the freeze-out have the same production vertex
	i->SetProductionVertex(i->Vx()+x, i->Vy()+y, i->Vz()+z, i->T());
}
