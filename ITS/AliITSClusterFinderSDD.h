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
  virtual Int_t CutAmplitude() {
    // get cut amplitude
    return fCutAmplitude;
  }
  virtual void SetDAnode(Float_t danode=4.2) {
    // setDAnode
    fDAnode=danode;
  }
  virtual Float_t DAnode() {
    // get DAnode
    return fDAnode;
  }
  virtual void SetDTime(Float_t dtime=75) {
    // SetDTime
    fDTime=dtime;
  }
  virtual Float_t DTime() {
    // get DTime
    return fDTime;
  }
  virtual void SetMinPeak(Int_t minpeak=10) {
    // SetMinPeak
    fMinPeak=minpeak;
  }
  virtual Int_t MinPeak() {
    // get MinPeak
    return fMinPeak;
  }
  virtual void SetMinCharge(Int_t mincharge=30) {
    // SetMinCharge
    fMinCharge=mincharge;
  }
  virtual Int_t MinCharge() {
    // get MinCharge
    return fMinCharge;
  }
  virtual void SetMinNCells(Int_t minc=3) {
    // setNCells
    fMinNCells=minc;
  }
  virtual Int_t MinNCells() {
    // get MinNCells
    return fMinNCells;
  }
   virtual void SetMaxNCells(Int_t maxc=10) {
    // setNCells
    fMaxNCells=maxc;
  }
  virtual Int_t MaxNCells() {
    // get MaxNCells
    return fMaxNCells;
  }
   virtual void SetTimeCorr(Float_t timec=23.) {
    // setNCells
    fTimeCorr=timec;
  }
  virtual Float_t TimeCorr() {
    // get Time Correction (ns)
    return fTimeCorr;
  }
  
  // Search for clusters
  virtual void FindRawClusters();
  void  Find1DClusters();
  void  Find1DClustersE();
  void  GroupClusters();
  void  SelectClusters();
  void  GetRecPoints();
  void ResolveClusters(); // Boris........ 
  void ResolveClustersE(); // Ernesto 
  Int_t SearchPeak( Float_t *spect, Int_t xdim, Int_t zdim, Int_t *peakX, Int_t
                    *peakZ, Float_t *peakAmp, Float_t minpeak ); // Ernesto
  Int_t noLinearFit( Int_t xdim, Int_t zdim, Float_t *param, Float_t *spe, Int_t
                     *niter, Float_t *chir );
  void minim( Int_t xdim, Int_t zdim, Float_t *param, Float_t *prm0, Float_t *steprm, Float_t *chisqr, 
		      Float_t *spe, Float_t *speFit );
  Float_t chisq( Int_t xdim, Int_t zdim, Float_t *spe, Float_t *speFit );
  void PeakFunc( Int_t xdim, Int_t zdim, Float_t *par, Float_t *spe, Float_t
                 *Integral=0 ); 

  virtual void Print();

private:
  
  TClonesArray       *fClusters;      // clusters
  Int_t               fNclusters;     // num of clusters
  Float_t             fDAnode;        // fDanode
  Float_t             fDTime;         // fDtime
  Float_t             fTimeCorr;      // Correction factor along time coord
  
  Int_t               fCutAmplitude;  // cut amplitude
  Int_t               fMinPeak;       // min peak
  Int_t               fMinCharge;     // min charge
  Int_t               fMinNCells;     // min num of cells
  Int_t               fMaxNCells;     // max num of cells
  ClassDef(AliITSClusterFinderSDD,1) // SDD clustering - Piergiorgio C. algo
    };
#endif







