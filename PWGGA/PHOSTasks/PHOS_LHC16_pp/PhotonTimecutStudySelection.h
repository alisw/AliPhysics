#ifndef PHOTONTIMECUTSELECTION_H
#define PHOTONTIMECUTSELECTION_H

// --- Custom header files ---
#include "PhotonSelection.h"
#include "DetectorHistogram.h"

// --- ROOT header files ---
#include "TH2F.h"


class PhotonTimecutStudySelection : public PhotonSelection
{
public:

	PhotonTimecutStudySelection():
		PhotonSelection(),
		fTimingCutPair(999999),
		fMassPt(),
		fMassPtMainMain(),
		fMassPtMainPileup(),
		fMassPtPileupPileup()
	{}

	PhotonTimecutStudySelection(const char * name, const char * title, ClusterCuts cuts):
		PhotonSelection(name, title, cuts), 
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
	PhotonTimecutStudySelection(const PhotonTimecutStudySelection &);
	PhotonTimecutStudySelection & operator = (const PhotonTimecutStudySelection &);
	virtual Bool_t IsMainBC(const AliVCluster * clus) const;

	// This one should't be used for selection,
	// it should be applied for combinations.
	//
	Float_t fTimingCutPair;
private:

	DetectorHistogram * fMassPt[2];             //!
	DetectorHistogram * fMassPtMainMain[2];     //!
	DetectorHistogram * fMassPtMainPileup[2];   //!
	DetectorHistogram * fMassPtPileupPileup[2]; //!

	ClassDef(PhotonTimecutStudySelection, 1)
};
#endif