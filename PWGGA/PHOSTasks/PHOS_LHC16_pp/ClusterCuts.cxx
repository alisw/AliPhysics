// --- Custom libraries ---
#include "ClusterCuts.h"


Bool_t ClusterCuts::AcceptCluster(AliVCluster * clus) const
{
	if (clus->GetNCells() < fNCellsCut)
		return kFALSE;

	if (clus->E() < fClusterMinE)
		return kFALSE;

	if (TMath::Abs(clus->GetTOF()) > fTimingCut)
		return kFALSE;

	return kTRUE;
}

ClusterCuts ClusterCuts::GetClusterCuts(Int_t ctype)
{
	// This parameter is redundant,
	// but it will be easier to extend this for other periods/datasets
	(void) ctype;
	ClusterCuts cuts;
	cuts.fNCellsCut = 3;
	cuts.fClusterMinE = 0.3;
	cuts.fTimingCut = 12.5e-9;
	cuts.fAsymmetryCut = 1.0;
	return cuts;
}
