#ifndef ALITPCSEED_H
#define ALITPCSEED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------
//   TPC seed
//   Class  needed for TPC parallel tracking 
//
//   Origin: 
//-------------------------------------------------------

#include <TError.h>
#include <TBits.h>

#include "AliTPCtrack.h"
#include "AliComplexCluster.h"
#include "AliPID.h"

class TFile;
class AliTPCParam;
class AliTPCseed;
class AliTPCclusterMI;
class AliTPCTrackerPoint;
class AliESD;
class AliTPCCalPad;
class TClonesArray;

class AliTPCseed : public AliTPCtrack {
  public:  
     AliTPCseed();
     virtual ~AliTPCseed();
     AliTPCseed(const AliTPCtrack &t);
     AliTPCseed(const AliTPCseed &s, Bool_t clusterOwner = kFALSE);
     //AliTPCseed(const AliTPCseed &t, Double_t a);
     AliTPCseed(Double_t xr, Double_t alpha, const Double_t xx[5], 
                const Double_t cc[15], Int_t i);   
     AliTPCseed &operator = (const AliTPCseed & param);  
     static Int_t  RefitTrack(AliTPCseed* seed, AliExternalTrackParam * in, AliExternalTrackParam * out);
     Bool_t RefitTrack(AliTPCseed* seed, Bool_t out);
     Int_t Compare(const TObject *o) const;
     void Reset(Bool_t all = kTRUE);
     Int_t GetProlongation(Double_t xr, Double_t &y, Double_t & z) const;
     virtual Double_t GetPredictedChi2(const AliCluster *cluster2) const;
     virtual Bool_t Update(const AliCluster* c2, Double_t chi2, Int_t i);
     AliTPCTrackerPoint * GetTrackPoint(Int_t i);
     AliTPCclusterMI * GetClusterFast(Int_t irow){ return fClusterPointer[irow];}
     void SetClusterPointer(Int_t irow, AliTPCclusterMI* cl) {fClusterPointer[irow]=cl;}
     Double_t GetDensityFirst(Int_t n);
     Double_t GetSigma2C() const {
       Double_t cnv=GetBz()*kB2C;
       return GetSigma1Pt2()*cnv*cnv;
     }
     void GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2);
     
     void Modify(Double_t factor);
     void SetClusterIndex2(Int_t row, Int_t index) {
       fIndex[row] = index;
     }
     Int_t  GetClusterIndex2(Int_t row) const {
       return fIndex[row];
     }
     Int_t GetClusterSector(Int_t row) const {
       Int_t pica = -1;
       if (fIndex[row]>=0) pica =  ((fIndex[row]&0xff000000)>>24);
       return pica;
     }
    
     Double_t GetYat(Double_t x) const;

     void SetErrorY2(Float_t sy2){fErrorY2=sy2;}
     void SetErrorZ2(Float_t sz2){fErrorZ2=sz2;}
     Float_t  CookdEdx(Double_t low=0.05, Double_t up=0.70, Int_t i1=0, Int_t i2=159, Bool_t onlyused = kFALSE);
     Float_t  CookShape(Int_t type);
     //  Float_t CookShape2(Int_t type,Bool_t withQ);
     void CookPID();
     Bool_t IsActive() const { return !(fRemoval);}
     void Desactivate(Int_t reason){ fRemoval = reason;} 
     AliTPCclusterMI* GetClusterPointer(Int_t i) const {return fClusterPointer[i];}
     Int_t GetSector() const {return fSector;}
     Float_t GetCurrentSigmaY2() const {return fCurrentSigmaY2;}
     Float_t GetCurrentSigmaZ2() const {return fCurrentSigmaZ2;}
     Int_t GetRelativeSector() const {return fRelativeSector;}
     Char_t GetCircular() const {return fCircular;}

     void SetCurrentSigmaY2(Float_t s) {fCurrentSigmaY2=s;}
     void SetCurrentSigmaZ2(Float_t s) {fCurrentSigmaZ2=s;}
     void SetRelativeSector(Int_t r) {fRelativeSector=r;}
     void SetCircular(Char_t c) {fCircular=c;}
     void SetIsSeeding(Bool_t s) {fIsSeeding=s;}
     void SetSeedType(Int_t s) {fSeedType=s;}
     void SetSeed1(Int_t s) {fSeed1=s;}
     void SetSeed2(Int_t s) {fSeed2=s;}
     void SetESD(AliESDtrack* esd) {fEsd=esd;}
     void SetBSigned(Bool_t s) {fBSigned=s;}
     void SetSort(Int_t s) {fSort=s;}
     void SetOverlapLabel(Int_t i, Int_t l) {fOverlapLabels[i]=l;}
     void SetCurrentCluster(AliTPCclusterMI* cl) {fCurrentCluster=cl;}
     void SetNoCluster(Int_t n) {fNoCluster=n;}
     void SetRow(Int_t n) {fRow=n;}
     void SetSector(Int_t n) {fSector=n;}
     void SetCurrentClusterIndex1(Int_t n) {fCurrentClusterIndex1=n;}
     void SetInDead(Bool_t s) {fInDead=s;}

     Double_t TPCrPID(Int_t i) const {return fTPCr[i];}
     Double_t* TPCrPIDs() {return fTPCr;}
     Bool_t GetIsSeeding() const {return fIsSeeding;}
     Int_t GetSeedType() const {return fSeedType;}
     Int_t GetSeed1() const {return fSeed1;}
     Int_t GetSeed2() const {return fSeed2;}
     AliESDtrack* GetESD() {return fEsd;}
     Float_t GetSDEDX(Int_t i) const {return fSDEDX[i];}
     Float_t GetDEDXregion(Int_t i) const {return fDEDX[i];}
     Int_t GetNCDEDX(Int_t i) const {return fNCDEDX[i];}
     Int_t GetNCDEDXInclThres(Int_t i) const {return fNCDEDXInclThres[i];}
     Bool_t GetBSigned() const {return fBSigned;}
     Int_t GetSort() const {return fSort;}
     Int_t GetOverlapLabel(Int_t i) const {return fOverlapLabels[i];}
     AliTPCclusterMI* GetCurrentCluster() const {return fCurrentCluster;}
     Int_t GetNoCluster() const {return fNoCluster;}
     Int_t GetRow() const {return fRow;}
     Int_t GetCurrentClusterIndex1() const {return fCurrentClusterIndex1;}
     Bool_t GetInDead() const {return fInDead;}
     Float_t GetErrorY2() const {return fErrorY2;}
     Float_t GetErrorZ2() const {return fErrorZ2;}
  Float_t GetCMeanSigmaY2p30() const {return fCMeanSigmaY2p30;}
  Float_t GetCMeanSigmaZ2p30() const {return fCMeanSigmaZ2p30;}
  Float_t GetCMeanSigmaY2p30R() const {return fCMeanSigmaY2p30R;}
  Float_t GetCMeanSigmaZ2p30R() const {return fCMeanSigmaZ2p30R;}
     //
     //

  Float_t  CookdEdxNorm(Double_t low=0.05, Double_t up=0.70, Int_t type=0, Int_t i1=0, Int_t i2=159, Bool_t shapeNorm=kTRUE, Int_t posNorm=0, Int_t padNorm=0,Int_t returnVal=0);

  Float_t  CookdEdxAnalytical(Double_t low=0.05, Double_t up=0.70, Int_t type=0, Int_t i1=0, Int_t i2=159, Int_t returnVal=0, Int_t rowThres = 2, Int_t mode=0);

 static   void GetError(AliTPCclusterMI* cluster, AliExternalTrackParam * param, 
			 Double_t& erry, Double_t &errz);
 static   void GetShape(AliTPCclusterMI* cluster, AliExternalTrackParam * param, 
			 Double_t& rmsy, Double_t &rmsz);
  static   Double_t GetQCorrGeom(Float_t ty, Float_t tz);
  static   Double_t GetQCorrShape(Int_t ipad, Int_t type,Float_t z, Float_t ty, Float_t tz, Float_t q, Float_t thr);
  //
  Float_t GetTPCClustInfo(Int_t nNeighbours, Int_t type, Int_t row0, Int_t row1);
  
 private:
     //     AliTPCseed & operator = (const AliTPCseed &)
     //  {::Fatal("= operator","Not Implemented\n");return *this;}
     AliESDtrack * fEsd; //!
     AliTPCclusterMI*   fClusterPointer[160];  // array of cluster pointers  - 
     Bool_t             fClusterOwner;         // indicates the track is owner of cluster
     //---CURRENT VALUES
     Int_t fRow;                 // current row number  
     Int_t fSector;              // current sector number
     Int_t fRelativeSector;      // index of current relative sector
     Float_t fCurrentSigmaY2;    //!expected current cluster sigma Y
     Float_t fCurrentSigmaZ2;    //!expected current cluster sigma Z
     Float_t fCMeanSigmaY2p30;   //! current mean sigma Y2 - mean30%
     Float_t fCMeanSigmaZ2p30;   //! current mean sigma Z2 - mean30%
     Float_t fCMeanSigmaY2p30R;   //! current relative mean sigma Y2 - mean30%
     Float_t fCMeanSigmaZ2p30R;   //! current relative mean sigma Z2 - mean30%
     Float_t fErrorY2;           //!sigma of current cluster 
     Float_t fErrorZ2;           //!sigma of current cluster    
     AliTPCclusterMI * fCurrentCluster; //!pointer to the current cluster for prolongation
     Int_t   fCurrentClusterIndex1; //! index of the current cluster
     Bool_t  fInDead;            //! indicate if the track is in dead zone
     Bool_t  fIsSeeding;         //!indicates if it is proces of seeading
     Int_t   fNoCluster;         //!indicates number of rows without clusters
     Int_t   fSort;              //!indicate criteria for sorting
     Bool_t  fBSigned;        //indicates that clusters of this trackes are signed to be used
     //
     //
     Float_t fDEDX[5];            // dedx according padrows
     Float_t fSDEDX[4];           // sdedx according padrows
     Int_t   fNCDEDX[4];          // number of clusters for dedx measurment
     Int_t   fNCDEDXInclThres[4]; // number of clusters for dedx measurment including sub-threshold clusters
     Double_t fTPCr[AliPID::kSPECIES];   // rough PID according TPC   
     //
     Int_t   fSeedType;         //seeding type
     Int_t   fSeed1;            //first row for seeding
     Int_t   fSeed2;            //last row for seeding
     Int_t   fOverlapLabels[12];  //track labels and the length of the  overlap     
     Float_t fMAngular;           // mean angular factor
     Char_t   fCircular;           // indicates curlin track
     AliTPCTrackerPoint  fTrackPoints[160];  //track points - array track points
     ClassDef(AliTPCseed,5)  
};





#endif


