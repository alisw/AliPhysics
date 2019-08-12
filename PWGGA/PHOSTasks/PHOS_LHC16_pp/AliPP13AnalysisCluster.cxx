// --- Custom header files ---
#include <AliPP13AnalysisCluster.h>

// --- ROOT system ---

// --- AliRoot header files ---


ClassImp(AliPP13AnalysisCluster);

//________________________________________________________________
AliPP13AnalysisCluster::AliPP13AnalysisCluster(): 
	AliAODCaloCluster(),
	// fCluster(0),
	fTRU(-1),
	fTRUCh(-1),
	fTRUChX(-1),
	fTRUChZ(-1),
	fModule(-1),
	fTrigger(kFALSE) 
{}

//________________________________________________________________
AliPP13AnalysisCluster::AliPP13AnalysisCluster(const AliAODCaloCluster & c):
	AliAODCaloCluster(c),
	// fCluster(c),
	fTRU(-1),
	fTRUCh(-1),
	fTRUChX(-1),
	fTRUChZ(-1),
	fModule(-1),
	fTrigger(kFALSE) 	
{}
