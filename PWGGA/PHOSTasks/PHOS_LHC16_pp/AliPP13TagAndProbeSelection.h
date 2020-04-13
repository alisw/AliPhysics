#ifndef ALIPP13TAGANDPROBESELECTION_H
#define ALIPP13TAGANDPROBESELECTION_H

// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliLog.h>

class AliPP13TagAndProbeSelection: public AliPP13PhysicsSelection
{
public:
	AliPP13TagAndProbeSelection():
		AliPP13PhysicsSelection(),
		fTimingCut(999999),
		fMassEnergyAll(),
		fMassEnergyTOF()
	{
	}

	AliPP13TagAndProbeSelection(const char * name, const char * title, AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w),
		fTimingCut(cuts.fTimingCut),  // Copy timing cut information
		fMassEnergyAll(),
		fMassEnergyTOF()
	{
		fCuts.fTimingCut = 99999; // Don't use timecut for cluster selection
	}

	virtual void InitSelectionHistograms();
	
protected:
	virtual void SelectTwoParticleCombinations(const TObjArray & photonCandidates, const EventFlags & eflags);
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	AliPP13TagAndProbeSelection(const AliPP13TagAndProbeSelection &);
	AliPP13TagAndProbeSelection & operator = (const AliPP13TagAndProbeSelection &);
private:
	Float_t fTimingCut;
	AliPP13DetectorHistogram * fMassEnergyAll[2]; //!
	AliPP13DetectorHistogram * fMassEnergyTOF[2]; //!

	ClassDef(AliPP13TagAndProbeSelection, 2)
};
#endif
