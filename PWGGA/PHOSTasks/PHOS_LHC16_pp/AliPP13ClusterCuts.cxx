// --- Custom libraries ---
#include "AliPP13ClusterCuts.h"


Bool_t AliPP13ClusterCuts::AcceptCluster(AliVCluster * clus) const
{
	if (!clus->IsPHOS())
		return kFALSE;

	if (clus->GetType() != AliVCluster::kPHOSNeutral)
		return kFALSE; // don't use CPV

	// if (clus->E() < 0.1)
	// 	return kFALSE;

	// if (clus->E() > 1. && clus->GetNCells() < fNCellsCut)
	// 	return kFALSE;

	// if (clus->E() > 1. && clus->GetM02() < 0.1)
	// 	return kFALSE;

	if (clus->GetNCells() < fNCellsCut)
		return kFALSE;

	if (clus->E() < fClusterMinE)
		return kFALSE;

	if (TMath::Abs(clus->GetTOF()) > fTimingCut)
		return kFALSE;

	if (clus->GetDistanceToBadChannel() < fMinimalDistance)
		return kFALSE;

	return kTRUE;
}

Bool_t AliPP13ClusterCuts::AcceptPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags) const
{
	(void) eflags;

	Double_t asym = TMath::Abs( (c1->E() - c2->E()) / (c1->E() + c2->E()) );
	if (asym > fAsymmetryCut)
		return kFALSE;

	// Accept the pair
	return kTRUE;
}

AliPP13ClusterCuts AliPP13ClusterCuts::GetClusterCuts(Int_t ctype)
{
	// This parameter is redundant,
	// but it will be easier to extend this for other periods/datasets
	(void) ctype;
	AliPP13ClusterCuts cuts;
	cuts.fNCellsCut = 3;
	cuts.fClusterMinE = 0.3;
	cuts.fTimingCut = 12.5e-9;
	cuts.fAsymmetryCut = 1.0;
	cuts.fNContributors = 1;
	cuts.fMinimalDistance = 0;
	return cuts;
}
