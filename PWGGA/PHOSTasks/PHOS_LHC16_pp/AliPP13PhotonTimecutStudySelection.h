#ifndef ALIPP13PHOTONTIMECUTSELECTION_H
#define ALIPP13PHOTONTIMECUTSELECTION_H

// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"

// --- ROOT header files ---
#include "TH2F.h"


class AliPP13PhotonTimecutStudySelection : public AliPP13PhysicsSelection
{
public:

	AliPP13PhotonTimecutStudySelection():
		AliPP13PhysicsSelection(),
		fTimingCutPair(999999),
		fMassPt(),
		fMassPtMainMain(),
		fMassPtMainPileup(),
		fMassPtPileupPileup()
	{}

	AliPP13PhotonTimecutStudySelection(const char * name, const char * title, 
			AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w), 
		fTimingCutPair(cuts.fTimingCut), // Use timing cut for pair of clusters
		fMassPt(),
		fMassPtMainMain(),
		fMassPtMainPileup(),
		fMassPtPileupPileup()

	{
		fCuts.fTimingCut = 99999; // Pass infinite (99999) cut to select all clusters 
	}
	virtual void InitSelectionHistograms();

protected:

	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	AliPP13PhotonTimecutStudySelection(const AliPP13PhotonTimecutStudySelection &);
	AliPP13PhotonTimecutStudySelection & operator = (const AliPP13PhotonTimecutStudySelection &);
	virtual Bool_t IsMainBC(const AliVCluster * clus) const;

	// This one should't be used for selection,
	// it should be applied for combinations.
	//
	Float_t fTimingCutPair;
private:

	AliPP13DetectorHistogram * fMassPt[2];             //!
	AliPP13DetectorHistogram * fMassPtMainMain[2];     //!
	AliPP13DetectorHistogram * fMassPtMainPileup[2];   //!
	AliPP13DetectorHistogram * fMassPtPileupPileup[2]; //!

	ClassDef(AliPP13PhotonTimecutStudySelection, 1)
};
#endif
