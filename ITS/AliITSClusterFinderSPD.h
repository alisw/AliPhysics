#ifndef ALIITSCLUSTERFINDERSPD_H
#define ALIITSCLUSTERFINDERSPD_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////

#include "AliITSClusterFinder.h"

class AliITSMapA1;

class AliITSClusterFinderSPD :
  public AliITSClusterFinder

{
public:
  AliITSClusterFinderSPD
       (AliITSsegmentation *seg, TClonesArray *dig, TClonesArray *recp);
  AliITSClusterFinderSPD();
  virtual ~AliITSClusterFinderSPD();
  
  void SetDx(Float_t dx=1.) {
    // set dx
    fDx=dx;
  }
  void SetDz(Float_t dz=0.) {
    // set dz
    fDz=dz;
  }
  void SetNCells(Int_t minc=0) {
    // set ncells
    fMinNCells=minc;
  }
  
  // Search for clusters
  void FindRawClusters(Int_t mod=0);
  void  Find1DClusters();
  void  GroupClusters();
  void  TracksInCluster();
  void  SelectClusters() {
    // selects clusters
  }
  void  GetRecPoints();
  
  private:

  AliITSClusterFinderSPD(const AliITSClusterFinderSPD &source); // copy ctor
  AliITSClusterFinderSPD& operator=(const AliITSClusterFinderSPD &source); 

private:
  
  TClonesArray       *fClusters;      // clusters
  Int_t               fNclusters;     // num of clusters
  Float_t             fDz;            // dz
  Float_t             fDx;            // dx
  
  Int_t               fMinNCells;     // min num of cells in the cluster
  
  ClassDef(AliITSClusterFinderSPD,1)  // SPD clustering - Boris B. algo based
                                      // on Piergiorgio's algo
    };
#endif







