// --- Custom header files ---
#include "PhysPhotonSelectionMC.h"

// --- AliRoot header files ---
#include <AliPHOSAodCluster.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(PhysPhotonSelectionMC);


//________________________________________________________________
TLorentzVector PhysPhotonSelectionMC::ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const
{
    Float_t energy = c1->E();

    TLorentzVector p;
    c1->GetMomentum(p, eflags.vtxBest);
    p *= Nonlinearity(energy);
	return p;
}

//________________________________________________________________
Float_t PhysPhotonSelectionMC::Nonlinearity(Float_t x) const
{
	return fGlobalEnergyScale * (1. + fNonA * TMath::Exp(-x / 2. * x / fNonSigma / fNonSigma));
}


