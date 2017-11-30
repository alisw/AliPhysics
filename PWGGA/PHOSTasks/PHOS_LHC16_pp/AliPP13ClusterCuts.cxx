// --- Custom libraries ---
#include "AliPP13ClusterCuts.h"


Bool_t AliPP13ClusterCuts::AcceptCluster(AliVCluster * clus) const
{
	if (clus->GetNCells() < fNCellsCut)
		return kFALSE;

	if (clus->E() < fClusterMinE)
		return kFALSE;

	if (TMath::Abs(clus->GetTOF()) > fTimingCut)
		return kFALSE;

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
	return cuts;
}
