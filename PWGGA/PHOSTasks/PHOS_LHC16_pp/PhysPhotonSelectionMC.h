#ifndef PHYSPHOTONSELECTIONMC_H
#define PHYSPHOTONSELECTIONMC_H

// --- Custom header files ---
#include "PhysPhotonSelection.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

class PhysPhotonSelectionMC : public PhysPhotonSelection
{
public:
	PhysPhotonSelectionMC(): PhysPhotonSelection() {}
	PhysPhotonSelectionMC(const char * name, const char * title, ClusterCuts cuts, 
		Float_t nona = 0., Float_t nonsigma = 1., Float_t genergy = 1.):
		PhysPhotonSelection(name, title, cuts),
		fNonA(nona),
		fNonSigma(nonsigma),
		fGlobalEnergyScale(genergy)
	{
		fCuts.fTimingCut = 99999; // No timing cut in MC
	}

protected:
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;
	virtual Float_t Nonlinearity(Float_t x) const;

	PhysPhotonSelectionMC(const PhysPhotonSelectionMC &);
	PhysPhotonSelectionMC & operator = (const PhysPhotonSelectionMC &);

	// Parameters of nonlinearity parametrization
	Float_t fNonA;
	Float_t fNonSigma;
	Float_t fGlobalEnergyScale;

private:
	ClassDef(PhysPhotonSelectionMC, 2)
};
#endif