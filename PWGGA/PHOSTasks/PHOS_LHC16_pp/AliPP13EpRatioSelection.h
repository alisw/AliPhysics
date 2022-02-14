#ifndef ALIPP13EPRATIOSELECTION_H
#define ALIPP13EPRATIOSELECTION_H

// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"
	

// --- ROOT system ---
#include <TObjArray.h>
#include <TVector3.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13EpRatioSelection : public AliPP13PhysicsSelection
{
public:
	AliPP13EpRatioSelection():
		AliPP13PhysicsSelection(),
		fEpE(),
		fEpPt(),
		fPosition(),
		fPIDCriteria(),
		fTPCSignal()
	{}

	AliPP13EpRatioSelection(const char * name, const char * title, AliPP13ClusterCuts cuts,
		AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w),
		fEpE(),
		fEpPt(),
		fPosition(),
		fPIDCriteria(),
		fTPCSignal()
	{}

	virtual ~AliPP13EpRatioSelection()
	{
		// NB: Don't use this 
		// delete [] fInvariantMass;
		for(Int_t i = 0; i < 2; ++i)
		{
			if (fEpE[i]) 
				delete fEpE[i];

			if (fEpPt[i]) 
				delete fEpPt[i];
		}

		for(Int_t i = 0; i < 4; ++i)
		{
			if(fPosition[i])
				delete fPosition[i];
		}	
		
		// Don't delete fClusters, as ROOT will take 
		// care of it.
	}
	
	virtual void InitSelectionHistograms();

	// NB: It actually doesn't fill Pi0Mass
	virtual void FillHistograms(TObjArray * clusArray, TList * pool, const EventFlags & eflags); 

protected:
	TVector3 LocalPosition(const AliVCluster * clus) const;
	virtual void FillClusterHistograms(const AliVCluster * clus, const EventFlags & eflags);

	AliPP13EpRatioSelection(const AliPP13EpRatioSelection &);
	AliPP13EpRatioSelection & operator = (const AliPP13EpRatioSelection &);
private:
	AliPP13DetectorHistogram * fEpE[2]; //!
	AliPP13DetectorHistogram * fEpPt[2]; //!
	AliPP13DetectorHistogram * fPosition[4]; //!
	TH2 * fPIDCriteria[9]; //!
	TH2 * fTPCSignal[4]; //!
	ClassDef(AliPP13EpRatioSelection, 2)
};
#endif
