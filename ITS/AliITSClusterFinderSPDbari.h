#ifndef ALIITSCLUSTERFINDERSPDBARI_H
#define ALIITSCLUSTERFINDERSPDBARI_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////

#include "AliITSClusterFinder.h"

class AliITSMapA1;

class AliITSClusterFinderSPDbari :
  public AliITSClusterFinder

{
public:
  AliITSClusterFinderSPDbari
  (AliITSsegmentation *segmentation,
   TClonesArray *digits, TClonesArray *recpoints);
  AliITSClusterFinderSPDbari();
  virtual ~AliITSClusterFinderSPDbari(){
    // destructor
  }
  AliITSClusterFinderSPDbari(const AliITSClusterFinderSPDbari &source); // copy constructor
  AliITSClusterFinderSPDbari& operator=(const AliITSClusterFinderSPDbari &source); // assignment operator
  
  
  virtual void SetDx(Float_t dx=1.) {
    // set dx
    fDx=dx;
  }
  virtual void SetDz(Float_t dz=0.) {
    // set dz
    fDz=dz;
  }

  // Search for clusters
  virtual void FindRawClusters(); 
  void  DigitToPoint(Int_t nclus, Float_t *xcenter, Float_t *zcenter,
		   Float_t *errxcenter,Float_t *errzcenter,
		  Int_t *tr1clus, Int_t *tr2clus, Int_t *tr3clus);
  void  ClusterFinder(Int_t ndigits,Int_t digx[],Int_t digz[],
         	      Int_t digtr1[],Int_t digtr2[],Int_t digtr3[],
		      Int_t &nclus,
		      Float_t xcenter[],Float_t zcenter[],
		      Float_t errxcenter[],Float_t errzcenter[],  
		      Int_t tr1clus[],Int_t tr2clus[], Int_t tr3clus[]);  
  
  
  
private:
  
  TClonesArray       *fClusters;      // clusters
  Int_t               fNclusters;     // num of clusters
  Float_t             fDz;            // dz
  Float_t             fDx;            // dx
  
  Int_t               fMinNCells;     // min num of cells in the cluster
  
  ClassDef(AliITSClusterFinderSPDbari,1)  // SPD clustering based
                                          // on Nico Di Bari algorithm
    };
#endif

