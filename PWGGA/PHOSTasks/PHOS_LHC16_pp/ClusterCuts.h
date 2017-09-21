#ifndef CLUSTERCUTS_H
#define CLUSTERCUTS_H

// --- AliRoot header files ---
#include <AliVCluster.h>

// NB: There is no need to derive it from TObject 
// as it should be a lightweight class
struct ClusterCuts
{

	enum PredefinedSet{kStandardPHOS};

	Bool_t AcceptCluster(AliVCluster * clus) const;
	static ClusterCuts GetClusterCuts(Int_t ctype = kStandardPHOS);

	Float_t fClusterMinE;
	Float_t fAsymmetryCut;
	Float_t fTimingCut;
	Int_t fNCellsCut;
};

#endif