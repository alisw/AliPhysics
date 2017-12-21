#ifndef ALIPP13PHYSPHOTONSELECTIONMC_H
#define ALIPP13PHYSPHOTONSELECTIONMC_H

// --- Custom header files ---
#include "AliPP13PhysPhotonSelection.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13PhysPhotonSelectionMC : public AliPP13PhysPhotonSelection
{
public:
	AliPP13PhysPhotonSelectionMC(): AliPP13PhysPhotonSelection() {}
	AliPP13PhysPhotonSelectionMC(const char * name, const char * title, AliPP13ClusterCuts cuts, 
		Float_t nona = 0., Float_t nonsigma = 1., Float_t genergy = 1.):
		AliPP13PhysPhotonSelection(name, title, cuts),
		fNonA(nona),
		fNonSigma(nonsigma),
		fGlobalEnergyScale(genergy)
	{
		fCuts.fTimingCut = 99999; // No timing cut in MC
	}

protected:
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;
	virtual Float_t Nonlinearity(Float_t x) const;

	AliPP13PhysPhotonSelectionMC(const AliPP13PhysPhotonSelectionMC &);
	AliPP13PhysPhotonSelectionMC & operator = (const AliPP13PhysPhotonSelectionMC &);

	// Parameters of nonlinearity parametrization
	Float_t fNonA;
	Float_t fNonSigma;
	Float_t fGlobalEnergyScale;

private:
	ClassDef(AliPP13PhysPhotonSelectionMC, 2)
};
#endif