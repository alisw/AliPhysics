#ifndef ALIPP13SPECTRUMSELECTIONMC_H
#define ALIPP13SPECTRUMSELECTIONMC_H

// --- Custom header files ---
#include "AliPP13SpectrumSelection.h"
#include "AliPP13SelectionWeights.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13SpectrumSelectionMC : public AliPP13SpectrumSelection
{
public:
	AliPP13SpectrumSelectionMC(): AliPP13SpectrumSelection() {}
	AliPP13SpectrumSelectionMC(const char * name, const char * title,
			AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13SpectrumSelection(name, title, cuts, w)
	{
		fCuts.fTimingCut = 99999; // No timing cut in MC
	}

protected:
	virtual Bool_t IsPrimary(const AliAODMCParticle * particle, Double_t rcut=1.) const;
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;
	
	AliPP13SpectrumSelectionMC(const AliPP13SpectrumSelectionMC &);
	AliPP13SpectrumSelectionMC & operator = (const AliPP13SpectrumSelectionMC &);

private:
	ClassDef(AliPP13SpectrumSelectionMC, 2)
};
#endif
