#ifndef TAGANDPROBESELECTION_H
#define TAGANDPROBESELECTION_H

// --- Custom header files ---
#include "PhotonSelection.h"
#include "DetectorHistogram.h"

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliLog.h>

class TagAndProbeSelection: public PhotonSelection
{
public:
	TagAndProbeSelection():
		PhotonSelection(),
		fTimingCut(999999),
		fMassEnergyAll(),
		fMassEnergyTOF()
	{
	}

	TagAndProbeSelection(const char * name, const char * title, ClusterCuts cuts):
		PhotonSelection(name, title, cuts),
		fTimingCut(cuts.fTimingCut),  // Copy timing cut information
		fMassEnergyAll(),
		fMassEnergyTOF()
	{
		fCuts.fTimingCut = 99999; // Don't use timecut for cluster selection
	}

	virtual void InitSelectionHistograms();
	
protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	virtual void FillPi0Mass(TObjArray * clusArray, TList * pool, const EventFlags & eflags);

	TagAndProbeSelection(const TagAndProbeSelection &);
	TagAndProbeSelection & operator = (const TagAndProbeSelection &);
private:
	Float_t fTimingCut;
	DetectorHistogram * fMassEnergyAll[2]; //!
	DetectorHistogram * fMassEnergyTOF[2]; //!

	ClassDef(TagAndProbeSelection, 2)
};
#endif