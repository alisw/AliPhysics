#ifndef ALIITSCLUSTERFINDERSPD_H
#define ALIITSCLUSTERFINDERSPD_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////

#include "AliITSClusterFinder.h"

class AliITSMapA1;

class AliITSClusterFinderSPD : public AliITSClusterFinder{
 public:
    AliITSClusterFinderSPD(AliITSsegmentation *segmentation,
			   TClonesArray *digits, TClonesArray *recpoints);
    AliITSClusterFinderSPD();
    virtual ~AliITSClusterFinderSPD(){// destructor
    }
    // copy constructor
    AliITSClusterFinderSPD(const AliITSClusterFinderSPD &source);
    // assignment operator
    AliITSClusterFinderSPD& operator=(const AliITSClusterFinderSPD &source);
  
    virtual void SetDx(Float_t dx=1.) {// set dx
	fDx=dx;}
    virtual void SetDz(Float_t dz=0.) {// set dz
	fDz=dz;}
    // Search for clusters
    virtual void FindRawClusters(Int_t module); 
    void  DigitToPoint(Int_t nclus, Float_t *xcenter, Float_t *zcenter,
		       Float_t *errxcenter,Float_t *errzcenter,
		       Int_t *tr1clus, Int_t *tr2clus, Int_t *tr3clus);
    void  ClusterFinder(Int_t ndigits,Int_t digx[],Int_t digz[],
			Int_t digtr1[],Int_t digtr2[],Int_t digtr3[],
			Int_t digtr4[],
			Int_t &nclus,
			Float_t xcenter[],Float_t zcenter[],
			Float_t errxcenter[],Float_t errzcenter[],  
			Int_t tr1clus[],Int_t tr2clus[], Int_t tr3clus[],
			Int_t module);
 private:
    TClonesArray       *fClusters;      // clusters
    Int_t               fNclusters;     // num of clusters
    Float_t             fDz;            // dz
    Float_t             fDx;            // dx
    Int_t               fMinNCells;     // min num of cells in the cluster
  
    ClassDef(AliITSClusterFinderSPD,1)  // SPD clustering
};
#endif
