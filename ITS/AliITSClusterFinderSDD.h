#ifndef ALIITSCLUSTERFINDERSDD_H
#define ALIITSCLUSTERFINDERSDD_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////
/*
  $Id$
*/

#include "AliITSClusterFinder.h"

class AliITSMapA2;
class TFile;

class AliITSClusterFinderSDD : public AliITSClusterFinder{
 public:
    AliITSClusterFinderSDD
	(AliITSsegmentation *seg,AliITSresponse *response,
	 TClonesArray *digits,TClonesArray *recpoints);
    AliITSClusterFinderSDD();
    virtual ~AliITSClusterFinderSDD();

    virtual void  SetCutAmplitude(Float_t nsigma=4);
    virtual Int_t CutAmplitude() const {// get cut amplitude
	return fCutAmplitude;}
    virtual void SetDAnode(Float_t danode=4.2) {// setDAnode
	fDAnode=danode;}
    virtual Float_t DAnode() const {// get DAnode
	return fDAnode;}
    virtual void SetDTime(Float_t dtime=75) {// SetDTime
	fDTime=dtime;}
    virtual Float_t DTime() const {// get DTime
	return fDTime;}
    virtual void SetMinPeak(Int_t minpeak=10) {// SetMinPeak
	fMinPeak=minpeak;}
    virtual Int_t MinPeak() const {// get MinPeak
	return fMinPeak;}
    virtual void SetMinCharge(Int_t mincharge=30) {// SetMinCharge
	fMinCharge=mincharge;}
    virtual Int_t MinCharge() const {// get MinCharge
	return fMinCharge;}
    virtual void SetMinNCells(Int_t minc=3) {// setNCells
	fMinNCells=minc;}
    virtual Int_t MinNCells() const {// get MinNCells
	return fMinNCells;}
    virtual void SetMaxNCells(Int_t maxc=10) {// setNCells
	fMaxNCells=maxc;}
    virtual Int_t MaxNCells() const {// get MaxNCells
	return fMaxNCells;}
    virtual void SetTimeCorr(Float_t timec=19.3) {// setNCells
	fTimeCorr=timec;}
    virtual Float_t TimeCorr() const{// get Time Correction (ns)
	return fTimeCorr;}

    // Search for clusters
    virtual void FindRawClusters(Int_t mod=0);
    void  Find1DClusters();
    void  Find1DClustersE();
    void  GroupClusters();
    void  SelectClusters();
    void  GetRecPoints();
    void  ResolveClusters(); // Boris........ 
    void  ResolveClustersE(); // Ernesto 
    Int_t SearchPeak(Float_t *spect,Int_t xdim,Int_t zdim,Int_t *peakX,
		     Int_t *peakZ,Float_t *peakAmp,Float_t minpeak); // Ernesto
    Int_t NoLinearFit( Int_t xdim, Int_t zdim, Float_t *param, Float_t *spe,
		       Int_t *niter, Float_t *chir );
    void  Minim( Int_t xdim, Int_t zdim, Float_t *param, Float_t *prm0,
		Float_t *steprm, Float_t *chisqr,Float_t *spe,Float_t *speFit);
    Float_t ChiSqr ( Int_t xdim, Int_t zdim, Float_t *spe, Float_t *speFit ) const;
    void  PeakFunc( Int_t xdim, Int_t zdim, Float_t *par, Float_t *spe,
		   Float_t *Integral=0 );
    void  Print() const;

 private:
    AliITSClusterFinderSDD(const AliITSClusterFinderSDD &source); // copy ctor
    AliITSClusterFinderSDD& operator=(const AliITSClusterFinderSDD &source);
 private:
    Int_t               fModule;        //! ITS current module 
    TClonesArray       *fClusters;      //! clusters
    Int_t               fNclusters;     //! num of clusters
    Float_t             fDAnode;        //! fDanode
    Float_t             fDTime;         //! fDtime
    Float_t             fTimeCorr;      //! Correction factor along time coord
    Int_t               fCutAmplitude;  //! cut amplitude
    Int_t               fMinPeak;       //! min peak
    Int_t               fMinCharge;     //! min charge
    Int_t               fMinNCells;     //! min num of cells
    Int_t               fMaxNCells;     //! max num of cells

    ClassDef(AliITSClusterFinderSDD,1) // SDD clustering - Piergiorgio C. algo
    };

#endif
