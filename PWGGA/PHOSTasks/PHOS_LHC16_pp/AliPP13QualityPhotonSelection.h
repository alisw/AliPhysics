#ifndef ALIPP13QUALITYPHOTONSELECTION_H
#define ALIPP13QUALITYPHOTONSELECTION_H

// --- Custom header files ---
#include "AliPP13DetectorHistogram.h"
#include "AliPP13PhysicsSelection.h"
#include "AliPP13SelectionWeights.h"

// --- ROOT system ---
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <TObjArray.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliLog.h>

class AliPP13QualityPhotonSelection : public AliPP13PhysicsSelection
{
public:
	AliPP13QualityPhotonSelection():
		AliPP13PhysicsSelection(),
		fClusterNXZ(),
		fClusterEXZ(),
		fClusterTime(0),
		fClusterTimeWide(0),
		fClusterEvsT(0),
		fClusterEvsTWide(0),
		fClusterTimeMap(0),
		fClusterIdN(),
		fClusterIdE(),
		fMassPtA(),
		fAsymmetry(0),
		fZvertex(0),
		fNcellsE(0),
		fShapeE(0)
	{
	}

	AliPP13QualityPhotonSelection(const char * name, const char * title, AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13PhysicsSelection(name, title, cuts, w),
		fClusterNXZ(),
		fClusterEXZ(),
		fClusterTime(0),
		fClusterTimeWide(0),
		fClusterEvsT(0),
		fClusterEvsTWide(0),
		fClusterTimeMap(0),
		fClusterIdN(),
		fClusterIdE(),
		fMassPtA(),
		fAsymmetry(0),
		fZvertex(0),
		fNcellsE(0),
		fShapeE(0)
	{
	}

	~AliPP13QualityPhotonSelection()
	{
		// NB: Don't use "delete []", to supress warnings.
		//

		for(Int_t i = 0; i < 2; ++i)
		{
			if (fClusterNXZ[i]) delete fClusterNXZ[i];
			if (fClusterEXZ[i]) delete fClusterEXZ[i];
		}

		if (fClusterTime)    delete fClusterTime;
		if (fClusterTimeWide)    delete fClusterTime;
		if (fClusterEvsT)    delete fClusterEvsT;
		if (fClusterEvsT)    delete fClusterEvsTWide;
		if (fClusterTimeMap) delete fClusterTimeMap;
	}

	virtual void InitSelectionHistograms();
	virtual Bool_t SelectEvent(const EventFlags & flgs);

protected:
	virtual void SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags);
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	AliPP13QualityPhotonSelection(const AliPP13QualityPhotonSelection &);
	AliPP13QualityPhotonSelection & operator = (const AliPP13QualityPhotonSelection &);

private:
	virtual Int_t AbsId(Int_t x, Int_t z, Int_t sm) const;

	AliPP13DetectorHistogram * fClusterNXZ[2]; //!
	AliPP13DetectorHistogram * fClusterEXZ[2]; //!
	AliPP13DetectorHistogram * fClusterTime; //!
	AliPP13DetectorHistogram * fClusterTimeWide; //!
	AliPP13DetectorHistogram * fClusterEvsT; //!
	AliPP13DetectorHistogram * fClusterEvsTWide; //!
	AliPP13DetectorHistogram * fClusterTimeMap; //!
	TH1F * fClusterIdN[2]; //!
	TH1F * fClusterIdE[2]; //!

	TH3F * fMassPtA[2]; //!
	TH1F * fAsymmetry; //!
	TH1F * fZvertex; //!
	TH2F * fNcellsE; //!
	TH2F * fShapeE;  //!

	ClassDef(AliPP13QualityPhotonSelection, 2)
};
#endif