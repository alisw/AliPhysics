#ifndef ALIPP13NONLINEARITYSELECTION_H
#define ALIPP13NONLINEARITYSELECTION_H

// --- Custom header files ---
#include "AliPP13PhotonSelection.h"
#include "AliPP13SelectionWeights.h"
#include "AliPP13DetectorHistogram.h"

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliLog.h>

// TODO: Merge nonlinearity selection wiht nonlinearity Study
//

class AliPP13NonlinearitySelection: public AliPP13PhotonSelection
{
public:
	AliPP13NonlinearitySelection():
		AliPP13PhotonSelection(),
		fMassPt()
	{
	}

	AliPP13NonlinearitySelection(const char * name, const char * title, 
		AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w, Bool_t isMC = kTRUE):
		AliPP13PhotonSelection(name, title, cuts, w),
		fMassPt()
	{
		// NB: Disable timing cuts for MC
		if(isMC)
			fCuts.fTimingCut = 999999;
	}

	virtual void InitSelectionHistograms();

	~AliPP13NonlinearitySelection()
	{
		for(Int_t i = 0; i < 2; ++i)
		{
			if(fMassPt[i])
				delete fMassPt[i];
		}
	}
	
protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	virtual void FillPi0Mass(TObjArray * clusArray, TList * pool, const EventFlags & eflags);

	AliPP13NonlinearitySelection(const AliPP13NonlinearitySelection &);
	AliPP13NonlinearitySelection & operator = (const AliPP13NonlinearitySelection &);

private:
	AliPP13DetectorHistogram * fMassPt[2]; //!
	ClassDef(AliPP13NonlinearitySelection, 2)
};
#endif