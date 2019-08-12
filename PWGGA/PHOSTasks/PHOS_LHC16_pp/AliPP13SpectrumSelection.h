#ifndef ALIPP13SPECTRUMSELECTION_H
#define ALIPP13SPECTRUMSELECTION_H

// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"
	

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13SpectrumSelection : public AliPP13PhysicsSelection
{
public:
	AliPP13SpectrumSelection():
		AliPP13PhysicsSelection(),
		fInvariantMass(),
		fClusters(0)
	{}

	AliPP13SpectrumSelection(const char * name, const char * title, AliPP13ClusterCuts cuts,
			AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w),
		fInvariantMass(),
		fClusters(0)
	{}

	virtual ~AliPP13SpectrumSelection()
	{
		// NB: Don't use this 
		// delete [] fInvariantMass;

		for(Int_t i = 0; i < 2; ++i)
			if (fInvariantMass[i])
				delete fInvariantMass[i];

		// Don't delete fClusters, as ROOT will take 
		// care of it.
	}
	
	virtual void InitSelectionHistograms();

protected:
	virtual void FillClusterHistograms(const AliVCluster * clus, const EventFlags & eflags);
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	AliPP13SpectrumSelection(const AliPP13SpectrumSelection &);
	AliPP13SpectrumSelection & operator = (const AliPP13SpectrumSelection &);

private:
	AliPP13DetectorHistogram * fInvariantMass[2]; //!
	TH1 * fClusters; //!
	ClassDef(AliPP13SpectrumSelection, 2)
};
#endif
