#ifndef PHYSPHOTONSELECTION_H
#define PHYSPHOTONSELECTION_H

// --- Custom header files ---
#include "PhotonSelection.h"
#include "DetectorHistogram.h"

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

class PhysPhotonSelection : public PhotonSelection
{
public:
	PhysPhotonSelection():
		PhotonSelection(),
		fInvariantMass(),
		fClusters(0)
	{}

	PhysPhotonSelection(const char * name, const char * title, ClusterCuts cuts):
		PhotonSelection(name, title, cuts),
		fInvariantMass(),
		fClusters(0)
	{}

	virtual ~PhysPhotonSelection()
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

	PhysPhotonSelection(const PhysPhotonSelection &);
	PhysPhotonSelection & operator = (const PhysPhotonSelection &);

private:
	DetectorHistogram * fInvariantMass[2]; //!
	TH1 * fClusters; //!
	ClassDef(PhysPhotonSelection, 2)
};
#endif