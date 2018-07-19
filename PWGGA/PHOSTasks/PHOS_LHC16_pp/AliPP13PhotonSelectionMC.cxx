// --- Custom header files ---
#include "AliPP13PhotonSelectionMC.h"

// --- ROOT system ---
#include <TH2F.h>
#include <TH3F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13PhotonSelectionMC);


//________________________________________________________________
TLorentzVector AliPP13PhotonSelectionMC::ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const
{
	TLorentzVector p;
	c1->GetMomentum(p, eflags.vtxBest);

	// NB: Apply nonlinearity Correction Here
	Float_t energy = c1->E();
	p *= fWeights->Nonlinearity(energy);
	return p;
}