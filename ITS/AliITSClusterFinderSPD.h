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
  (AliITSsegmentation *segmentation,
   TClonesArray *digits, TClonesArray *recpoints);
  AliITSClusterFinderSPD();
  virtual ~AliITSClusterFinderSPD();
  AliITSClusterFinderSPD(const AliITSClusterFinderSPD &source); // copy constructor
  AliITSClusterFinderSPD& operator=(const AliITSClusterFinderSPD &source); // assignment operator
  
  virtual void SetMap();
  virtual void SetDx(Float_t dx=1.) {
    // set dx
    fDx=dx;
  }
  virtual void SetDz(Float_t dz=0.) {
    // set dz
    fDz=dz;
  }
  virtual void SetNCells(Int_t minc=0) {
    // set ncells
    fMinNCells=minc;
  }
  
  // Search for clusters
  virtual void FindRawClusters();
  void  Find1DClusters();
  void  GroupClusters();
  void  SelectClusters() {
    // selects clusters
  }
  void  GetRecPoints();
  
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







