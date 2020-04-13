#ifndef ALIPP13TRIGGEREFFICIENCY_H
#define ALIPP13TRIGGEREFFICIENCY_H

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

class AliPP13TriggerEfficiency: public AliPP13PhysicsSelection
{
	enum {kTRUs = 7};
public:
	AliPP13TriggerEfficiency():
		AliPP13PhysicsSelection(),
		fNevents(0),
		fTotalMassEnergyAll(),
		fTotalMassEnergyTrigger()
	{
	}

	AliPP13TriggerEfficiency(const char * name, const char * title, AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w),
		fNevents(0),
		fTotalMassEnergyAll(),
		fTotalMassEnergyTrigger()
	{
	}

	~AliPP13TriggerEfficiency()
	{
		// for (Int_t tru = 0; tru < kTRUs; ++tru)
		// {
		// 	for (Int_t i = 0; i < 2; ++i)
		// 	{
		// 		delete fMassEnergyAll[tru][i];
		// 		delete fMassEnergyTrigger[tru][i];
		// 	}
		// }

		for (Int_t i = 0; i < 2; ++i)
		{
			delete fTotalMassEnergyAll[i];
			delete fTotalMassEnergyTrigger[i];
		}

	}

	virtual void InitSelectionHistograms();

protected:
	virtual void SelectTwoParticleCombinations(const TObjArray & photonCandidates, const EventFlags & eflags);
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	virtual void SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags);

	AliPP13TriggerEfficiency(const AliPP13TriggerEfficiency &);
	AliPP13TriggerEfficiency & operator = (const AliPP13TriggerEfficiency &);
private:
	TH1 * fNevents; //!
	AliPP13DetectorHistogram *  fTotalMassEnergyAll[2]; //!
	AliPP13DetectorHistogram *  fTotalMassEnergyTrigger[2]; //!

	// AliPP13DetectorHistogram * fMassEnergyAll[kTRUs][2]; //!
	// AliPP13DetectorHistogram * fMassEnergyTrigger[kTRUs][2]; //!
	ClassDef(AliPP13TriggerEfficiency, 2)
};
#endif
