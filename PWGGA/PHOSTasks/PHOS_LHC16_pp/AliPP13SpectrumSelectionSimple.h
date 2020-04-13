#ifndef ALIPP13SPECTRUMSELECTIONSIMPLE_H
#define ALIPP13SPECTRUMSELECTIONSIMPLE_H

// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"
	

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13SpectrumSelectionSimple : public AliPP13PhysicsSelection
{
public:
	AliPP13SpectrumSelectionSimple():
		AliPP13PhysicsSelection(),
		fMassPt(),
		fMixMassPt()
	{}

	AliPP13SpectrumSelectionSimple(const char * name, const char * title, AliPP13ClusterCuts cuts,
			AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w),
		fMassPt(),
		fMixMassPt()
	{}


	virtual void InitSelectionHistograms();

protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	AliPP13SpectrumSelectionSimple(const AliPP13SpectrumSelectionSimple &);
	AliPP13SpectrumSelectionSimple & operator = (const AliPP13SpectrumSelectionSimple &);

private:
	TH2 * fMassPt; //!
	TH2 * fMixMassPt; //!
	ClassDef(AliPP13SpectrumSelectionSimple, 2)
};
#endif
