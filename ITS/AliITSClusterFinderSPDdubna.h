#ifndef ALIITSCLUSTERFINDERSPDDUBNA_H
#define ALIITSCLUSTERFINDERSPDDUBNA_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////

#include "AliITSClusterFinder.h"

class AliITSMapA1;
class AliITSClusterFinderSPDdubna : public AliITSClusterFinder{
 public:
    AliITSClusterFinderSPDdubna(AliITSsegmentation *seg,TClonesArray *dig,
				TClonesArray *recp);
    AliITSClusterFinderSPDdubna();
    virtual ~AliITSClusterFinderSPDdubna();
    // copy ctor
    AliITSClusterFinderSPDdubna(const AliITSClusterFinderSPDdubna &source);
    // = opperator
    AliITSClusterFinderSPDdubna& operator=(const AliITSClusterFinderSPDdubna
					   &source);

    void SetDx(Float_t dx=1.) {// set dx
	fDx=dx;}
    void SetDz(Float_t dz=0.) {// set dz
	fDz=dz;}
    void SetNCells(Int_t minc=0) {// set ncells
	fMinNCells=minc;}

    // Search for clusters
    void FindRawClusters(Int_t mod=0);
    void Find1DClusters(Int_t mod);
    void GroupClusters();
    void TracksInCluster();
    void SelectClusters() {// selects clusters
    }
    void  GetRecPoints();
 private:
    TClonesArray       *fClusters;      // clusters
    Int_t               fNclusters;     // num of clusters
    Float_t             fDz;            // dz
    Float_t             fDx;            // dx
    Int_t               fMinNCells;     // min num of cells in the cluster
  
    ClassDef(AliITSClusterFinderSPDdubna,1)  // SPD clustering - Boris B.
	                                     // algo based on Piergiorgio's
	                                     // algo
};
#endif
