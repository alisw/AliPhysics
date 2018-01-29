#ifndef ALIPP13CLUSTERCUTS_H
#define ALIPP13CLUSTERCUTS_H

// --- AliRoot header files ---
#include <AliVCluster.h>

// NB: There is no need to derive it from TObject 
// as it should be a lightweight class


// TODO: Derive a new version for photon analysis
//       
struct AliPP13ClusterCuts
{

	enum PredefinedSet{kStandardPHOS};

	Bool_t AcceptCluster(AliVCluster * clus) const;
	static AliPP13ClusterCuts GetClusterCuts(Int_t ctype = kStandardPHOS);

	Float_t fClusterMinE;
	Float_t fAsymmetryCut;
	Float_t fTimingCut;
	Int_t fNCellsCut;
    Int_t fNContributors;
};

#endif