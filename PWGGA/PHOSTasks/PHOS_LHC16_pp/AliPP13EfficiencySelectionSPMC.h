#ifndef ALIPP13EFFICIENCYSELECTIONSPMC_H
#define ALIPP13EFFICIENCYSELECTIONSPMC_H


#include <map>

// --- Custom header files ---
#include "AliPP13SelectionWeights.h"
#include "AliPP13EfficiencySelectionMC.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>
#include <TF1.h>

// --- AliRoot header files ---
#include <AliAODMCParticle.h>
#include <AliVCluster.h>
#include <AliStack.h>
#include <AliLog.h>


class AliPP13EfficiencySelectionSPMC: public AliPP13EfficiencySelectionMC
{
public:
	AliPP13EfficiencySelectionSPMC(): AliPP13EfficiencySelectionMC() {}
	AliPP13EfficiencySelectionSPMC(const char * name, const char * title, AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13EfficiencySelectionMC(name, title, cuts, w)
	{
	}

	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);

protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	AliPP13EfficiencySelectionSPMC(const AliPP13EfficiencySelectionSPMC &);
	AliPP13EfficiencySelectionSPMC & operator = (const AliPP13EfficiencySelectionSPMC &);
	ClassDef(AliPP13EfficiencySelectionSPMC, 2)
};
#endif