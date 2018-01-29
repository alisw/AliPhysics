#ifndef ALIPP13PHOTONSELECTIONMC_H
#define ALIPP13PHOTONSELECTIONMC_H

// --- Custom header files ---
#include "AliPP13PhotonSelection.h"
#include "AliPP13DetectorHistogram.h"
#include "AliPP13SelectionWeights.h"

// --- ROOT system ---
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13PhotonSelectionMC : public AliPP13PhotonSelection
{
public:
	AliPP13PhotonSelectionMC():
		AliPP13PhotonSelection()
	{}

	AliPP13PhotonSelectionMC(const char * name, const char * title, 
		AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13PhotonSelection(name, title, cuts, w)
	{
		fCuts.fTimingCut = 99999; // No timing cut in MC
	}

	virtual ~AliPP13PhotonSelectionMC()
	{
	}
	
protected:
	TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;
	AliPP13PhotonSelectionMC(const AliPP13PhotonSelectionMC &);
	AliPP13PhotonSelectionMC & operator = (const AliPP13PhotonSelectionMC &);

private:
	ClassDef(AliPP13PhotonSelectionMC, 2)
};
#endif