// --- Custom header files ---
#include "AliPP13PhysPhotonSelectionMC.h"

// --- AliRoot header files ---
#include <AliPHOSAodCluster.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13PhysPhotonSelectionMC);


//________________________________________________________________
TLorentzVector AliPP13PhysPhotonSelectionMC::ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const
{
    Float_t energy = c1->E();

    TLorentzVector p;
    c1->GetMomentum(p, eflags.vtxBest);
    p *= Nonlinearity(energy);
	return p;
}

//________________________________________________________________
Float_t AliPP13PhysPhotonSelectionMC::Nonlinearity(Float_t x) const
{
	return fGlobalEnergyScale * (1. + fNonA * TMath::Exp(-x / 2. * x / fNonSigma / fNonSigma));
}


