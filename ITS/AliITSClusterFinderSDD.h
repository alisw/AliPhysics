#ifndef ALIITSCLUSTERFINDERSDD_H
#define ALIITSCLUSTERFINDERSDD_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////


#include "AliITSClusterFinder.h"

class AliITSMapA2;
class TFile;

class AliITSClusterFinderSDD :
 public AliITSClusterFinder

{
public:
  AliITSClusterFinderSDD
  (AliITSsegmentation *seg, 
   AliITSresponse *response, TClonesArray *digits,TClonesArray *recpoints);
  AliITSClusterFinderSDD();
  virtual ~AliITSClusterFinderSDD();
  AliITSClusterFinderSDD(const AliITSClusterFinderSDD &source); // copy constructor
  AliITSClusterFinderSDD& operator=(const AliITSClusterFinderSDD &source); // assignment operator
  
  virtual void SetCutAmplitude(Int_t thres=0) {
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
  virtual void SetMinNCells(Int_t minc=6) {
    // setNCells
    fMinNCells=minc;
  }
   virtual void SetMaxNCells(Int_t maxc=10) {
    // setNCells
    fMaxNCells=maxc;
  }
   virtual void SetTimeCorr(Float_t timec=70.) {
    // setNCells
    fTimeCorr=timec;
  }
  
  // Search for clusters
  virtual void FindRawClusters();
  void  Find1DClusters();
  void  GroupClusters();
  void  SelectClusters();
  void  GetRecPoints();
  
private:
  
  TClonesArray       *fClusters;      // clusters
  Int_t               fNclusters;     // num of clusters
  Float_t             fDAnode;        // fDanode
  Float_t             fDTime;         // fDtime
  Float_t             fTimeCorr;      // Correction factor along time coord
  
  Int_t               fCutAmplitude;  // cut amplitude
  Int_t               fMinPeak;       // min peak
  Int_t               fMinNCells;     // min num of cells
  Int_t               fMaxNCells;     // max num of cells
  ClassDef(AliITSClusterFinderSDD,1) // SDD clustering - Piergiorgio C. algo
    };
#endif







