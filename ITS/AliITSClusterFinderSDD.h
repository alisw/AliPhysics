#ifndef ALIITSCLUSTERFINDERSDD_H
#define ALIITSCLUSTERFINDERSDD_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////

#include "AliITSClusterFinder.h"

class AliITSClusterFinderSDD :
 public AliITSClusterFinder

{
public:
  AliITSClusterFinderSDD
  (AliITSsegmentation *seg, 
   AliITSresponse *response, TClonesArray *digits,TClonesArray *recpoints);
  AliITSClusterFinderSDD();
  virtual ~AliITSClusterFinderSDD(){
    // destructor
  }
  AliITSClusterFinderSDD(const AliITSClusterFinderSDD &source); // copy constructor
  AliITSClusterFinderSDD& operator=(const AliITSClusterFinderSDD &source); // assignment operator
  
  virtual void SetMap();
  virtual void SetCutAmplitude(Float_t thres=1.2) {
    // set cut amplitude
    fCutAmplitude=thres;
  }
  virtual void SetDAnode(Float_t danode=4.2) {
    // setDAnode
    fDAnode=danode;
  }
  virtual void SetDTime(Float_t dtime=75) {
    // SetDTime
    fDTime=dtime;
  }
  virtual void SetMinPeak(Int_t minpeak=7) {
    // SetMinPeak
    fMinPeak=minpeak;
  }
  virtual void SetNCells(Int_t minc=4) {
    // setNCells
    fMinNCells=minc;
  }
  
  void FillMap();
  
  // Search for clusters
  virtual void FindRawClusters();
  void  Find1DClusters();
  void  GroupClusters();
  void  SelectClusters();
  void  GetRecPoints();
  
private:
  
  TClonesArray       *fClusters;      // clusters
  Int_t               fNclusters;     // num of clusters
  AliITSMapA2        *fMap;           // map
  Float_t             fCutAmplitude;  // cut amplitude
  Float_t             fDAnode;        // fDanode
  Float_t             fDTime;         // fDtime
  
  Int_t               fMinPeak;       // min peak
  Int_t               fMinNCells;     // min num of cells
  
  
  ClassDef(AliITSClusterFinderSDD,1) // SDD clustering - Piergiorgio C. algo
    };
#endif







