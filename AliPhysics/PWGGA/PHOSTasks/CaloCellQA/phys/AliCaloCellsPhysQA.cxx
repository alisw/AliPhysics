#include <AliCaloCellsPhysQA.h>

ClassImp(AliCaloCellsPhysQA);

Int_t AliCaloCellsPhysQA::CheckClusterGetSM(AliVCluster * clus)
{
	// Reject CPV clusters (new in luster
	if (clus->GetType() != AliVCluster::kPHOSNeutral)
		return -1;

	if (clus->GetNCells() < fNCells) 
		return -1;

	if (clus->E() < fEmin)
		return -1;

	if (TMath::Abs(clus->GetTOF()) > fTimingCut)
		return -1;

	return AliCaloCellsQA::CheckClusterGetSM(clus);
}