#ifndef ALIPP13PHYSPHOTONSELECTIONMC_H
#define ALIPP13PHYSPHOTONSELECTIONMC_H

// --- Custom header files ---
#include "AliPP13PhysPhotonSelection.h"
#include "AliPP13SelectionWeights.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13PhysPhotonSelectionMC : public AliPP13PhysPhotonSelection
{
public:
	AliPP13PhysPhotonSelectionMC(): AliPP13PhysPhotonSelection() {}
	AliPP13PhysPhotonSelectionMC(const char * name, const char * title,
			AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13PhysPhotonSelection(name, title, cuts, w)
	{
		fCuts.fTimingCut = 99999; // No timing cut in MC
	}

protected:
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;
	
	AliPP13PhysPhotonSelectionMC(const AliPP13PhysPhotonSelectionMC &);
	AliPP13PhysPhotonSelectionMC & operator = (const AliPP13PhysPhotonSelectionMC &);

private:
	ClassDef(AliPP13PhysPhotonSelectionMC, 2)
};
#endif