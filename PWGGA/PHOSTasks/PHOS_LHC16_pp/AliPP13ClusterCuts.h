#ifndef ALIPP13CLUSTERCUTS_H
#define ALIPP13CLUSTERCUTS_H

// --- Custom libraries ---
#include <AliPP13SelectionWeights.h>

// --- AliRoot header files ---
#include <AliVCluster.h>

// NB: There is no need to derive it from TObject 
// as it should be a lightweight class


struct AliPP13ClusterCuts
{

	enum PredefinedSet{kStandardPHOS};

	Bool_t AcceptCluster(AliVCluster * clus) const;
	Bool_t AcceptPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags) const;
	static AliPP13ClusterCuts GetClusterCuts(Int_t ctype = kStandardPHOS);

	Float_t fClusterMinE;
	Float_t fAsymmetryCut;
	Float_t fTimingCut;
	Float_t fMinimalDistance;
	Int_t fNCellsCut;
    Int_t fNContributors;
};

#endif
