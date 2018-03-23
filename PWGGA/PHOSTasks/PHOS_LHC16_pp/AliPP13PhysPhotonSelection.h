#ifndef ALIPP13PHYSPHOTONSELECTION_H
#define ALIPP13PHYSPHOTONSELECTION_H

// --- Custom header files ---
#include "AliPP13PhotonSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"
	

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13PhysPhotonSelection : public AliPP13PhotonSelection
{
public:
	AliPP13PhysPhotonSelection():
		AliPP13PhotonSelection(),
		fInvariantMass(),
		fClusters(0)
	{}

	AliPP13PhysPhotonSelection(const char * name, const char * title, AliPP13ClusterCuts cuts,
			AliPP13SelectionWeights * w):
		AliPP13PhotonSelection(name, title, cuts, w),
		fInvariantMass(),
		fClusters(0)
	{}

	virtual ~AliPP13PhysPhotonSelection()
	{
		// NB: Don't use this 
		// delete [] fInvariantMass;

		for(Int_t i = 0; i < 2; ++i)
			if (fInvariantMass[i]) delete fInvariantMass[i];

		// Don't delete fClusters, as ROOT will take 
		// care of it.
	}
	
	virtual void InitSelectionHistograms();

protected:
	virtual void FillClusterHistograms(const AliVCluster * clus, const EventFlags & eflags);
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	AliPP13PhysPhotonSelection(const AliPP13PhysPhotonSelection &);
	AliPP13PhysPhotonSelection & operator = (const AliPP13PhysPhotonSelection &);

private:
	AliPP13DetectorHistogram * fInvariantMass[2]; //!
	TH1 * fClusters; //!
	ClassDef(AliPP13PhysPhotonSelection, 2)
};
#endif