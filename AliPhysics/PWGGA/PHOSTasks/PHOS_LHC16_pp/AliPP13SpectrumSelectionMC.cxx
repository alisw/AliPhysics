// --- Custom header files ---
#include "AliPP13SpectrumSelectionMC.h"

// --- AliRoot header files ---
#include <AliPHOSAodCluster.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13SpectrumSelectionMC);


//________________________________________________________________
TLorentzVector AliPP13SpectrumSelectionMC::ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const
{
    Float_t energy = c1->E();

    TLorentzVector p;
    c1->GetMomentum(p, eflags.vtxBest);
    p *= fWeights->Nonlinearity(energy);
	return p;
}

//________________________________________________________________
Bool_t AliPP13SpectrumSelectionMC::IsPrimary(const AliAODMCParticle * particle, Double_t rcut) const
{
	// Look what particle left vertex (e.g. with vertex with radius <1 cm)
	Double_t r2 = particle->Xv() * particle->Xv() + particle->Yv() * particle->Yv()	;
	return r2 < rcut * rcut;
}
