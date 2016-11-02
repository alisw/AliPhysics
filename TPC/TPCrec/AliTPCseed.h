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

#include <TVectorFfwd.h>
#include "AliTPCtrack.h"
#include "AliComplexCluster.h"
#include "AliPID.h"
#include "AliVTPCseed.h"
#include "AliTPCreco.h"

class TFile;
class AliTPCParam;
class AliTPCseed;
class AliTPCclusterMI;
class AliESD;
class AliTPCCalPad;
class TClonesArray;

class AliTPCseed : public AliTPCtrack, public AliVTPCseed {
  public:  
  enum {
    kInDead=BIT(15)            //! indicate if the track is in dead zone
    ,kIsSeeding=BIT(16)         //! indicates if it is proces of seeading
    ,kBSigned=BIT(17)           //indicates that clusters of this trackes are signed to be used
  };
     AliTPCseed();
     virtual ~AliTPCseed();
     virtual TObject* Clone(const char* newname = "") const;
     AliTPCseed(const AliTPCtrack &t);
     AliTPCseed(const AliTPCseed &s, Bool_t clusterOwner = kFALSE);
     //AliTPCseed(const AliTPCseed &t, Double_t a);
     AliTPCseed(Double_t xr, Double_t alpha, const Double_t xx[5], 
                const Double_t cc[15], Int_t i);   
     AliTPCseed &operator = (const AliTPCseed & param);  
     void Clear(Option_t* = "");
     static Int_t  RefitTrack(AliTPCseed* seed, AliExternalTrackParam * in, AliExternalTrackParam * out);
     Bool_t RefitTrack(AliTPCseed* seed, Bool_t out);
     Int_t Compare(const TObject *o) const;
     void Reset(Bool_t all = kTRUE);
     Int_t GetProlongation(Double_t xr, Double_t &y, Double_t & z) const;
     virtual Double_t GetPredictedChi2(const AliCluster *cluster2) const;
     virtual Bool_t Update(const AliCluster* c2, Double_t chi2, Int_t i);
     //
     const AliTPCTrackerPoints::Point* GetTrackPoint(Int_t i) const { return fTrackPointsArr.GetPoint(i); }
     //
     AliTPCclusterMI * GetClusterFast(Int_t irow){return fClusterPointer ? ((AliTPCclusterMI*)fClusterPointer[irow]):0;}
     AliTPCclusterMI * GetClusterFast(Int_t irow) const { return fClusterPointer ? fClusterPointer[irow]:0;}
     const AliTPCclusterMI** GetClusters() const {return (const AliTPCclusterMI**)fClusterPointer;}
     void  SetClustersArrayTMP(AliTPCclusterMI** arr) {fClusterPointer = arr; fNClStore = arr ? kMaxRow : 0;}
     void SetClusterPointer(Int_t irow, AliTPCclusterMI* cl);
     Double_t GetDensityFirst(Int_t n);
     Double_t GetSigma2C() const {
       Double_t cnv=GetBz()*kB2C;
       return GetSigma1Pt2()*cnv*cnv;
     }
     void GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2);
     void GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable);
     
     void Modify(Double_t factor);
     void  SetClusterIndex2(Int_t row, Int_t index) {fIndex[row] = index;} //RS: no check on index range
     Int_t GetClusterIndex2(Int_t row) const { return fIndex[row];} // RS: no check on index range
     Int_t GetClusterSector(Int_t row) const {
       Int_t pica = -1;
       if (fIndex[row]>=0) pica =  ((fIndex[row]&0xff000000)>>24);
       return pica;
     }
    
     Double_t GetYat(Double_t x) const;

     void SetErrorY2(Float_t sy2){fErrorY2=sy2;}
     void SetErrorZ2(Float_t sz2){fErrorZ2=sz2;}
     void SetErrorY2Syst(Float_t sy2){fErrorY2Syst=sy2;}
     void SetErrorZ2Syst(Float_t sz2){fErrorZ2Syst=sz2;}
     Float_t  CookdEdx(Double_t low=0.05, Double_t up=0.70, Int_t i1=0, Int_t i2=kMaxRow, Bool_t onlyused = kFALSE);
     Float_t  CookShape(Int_t type);
     //  Float_t CookShape2(Int_t type,Bool_t withQ);
     void CookPID();
     Bool_t IsActive() const { return !(fRemoval);}
     void Desactivate(Int_t reason){ fRemoval = reason;} 
     AliTPCclusterMI* GetClusterPointer(Int_t i) const {return GetClusterFast(i);}
     Int_t GetSector() const {return fSector;}
     Float_t GetCurrentSigmaY2() const {return fCurrentSigmaY2;}
     Float_t GetCurrentSigmaZ2() const {return fCurrentSigmaZ2;}
     Int_t GetRelativeSector() const {return fRelativeSector;}
     Char_t GetCircular() const {return fCircular;}

     void SetCurrentSigmaY2(Float_t s) {fCurrentSigmaY2=s;}
     void SetCurrentSigmaZ2(Float_t s) {fCurrentSigmaZ2=s;}
     void SetRelativeSector(Int_t r) {fRelativeSector=r;}
     void SetCircular(Char_t c) {fCircular=c;}
     void SetIsSeeding(Bool_t s) {SetBit(kIsSeeding,s);}
     void SetSeedType(Int_t s) {fSeedType=s;}
     void SetSeed1(Int_t s) {fSeed1=s;}
     void SetSeed2(Int_t s) {fSeed2=s;}
     void SetESD(AliESDtrack* esd) {fEsd=esd;}
     void SetBSigned(Bool_t s) {SetBit(kBSigned,s);}
     void SetSort(Int_t s) {fSort=s;}
     void SetOverlapLabel(Int_t i, Int_t l) {fOverlapLabels[i]=l;}
     void SetCurrentCluster(AliTPCclusterMI* cl) {fCurrentCluster=cl;}
     void SetNoCluster(Int_t n) {fNoCluster=n;}
     void SetRow(Int_t n) {fRow=n;}
     void SetSector(Int_t n) {fSector=n;}
     void SetCurrentClusterIndex1(Int_t n) {fCurrentClusterIndex1=n;}
     void SetInDead(Bool_t s) {SetBit(kInDead,s);}

     Double_t TPCrPID(Int_t i) const {return fTPCr[i];}
     Double_t* TPCrPIDs() {return fTPCr;}
     Bool_t GetIsSeeding() const {return TestBit(kIsSeeding);}
     Int_t GetSeedType() const {return fSeedType;}
     Int_t GetSeed1() const {return fSeed1;}
     Int_t GetSeed2() const {return fSeed2;}
     AliESDtrack* GetESD() {return fEsd;}
     Float_t GetSDEDX(Int_t i) const {return fSDEDX[i];}
     Float_t GetDEDXregion(Int_t i) const {return fDEDX[i];}
     Int_t GetNCDEDX(Int_t i) const {return fNCDEDX[i];}
     Int_t GetNCDEDXInclThres(Int_t i) const {return fNCDEDXInclThres[i];}
     Bool_t GetBSigned() const {return TestBit(kBSigned);}
     Int_t GetSort() const {return fSort;}
     Int_t GetOverlapLabel(Int_t i) const {return fOverlapLabels[i];}
     AliTPCclusterMI* GetCurrentCluster() const {return fCurrentCluster;}
     Int_t GetNoCluster() const {return fNoCluster;}
     Int_t GetRow() const {return fRow;}
     Int_t GetCurrentClusterIndex1() const {return fCurrentClusterIndex1;}
     Bool_t GetInDead() const {return TestBit(kInDead);}
     Float_t GetErrorY2() const {return fErrorY2;}
     Float_t GetErrorZ2() const {return fErrorZ2;}
     Float_t GetErrorY2Syst() const {return fErrorY2Syst;}
     Float_t GetErrorZ2Syst() const {return fErrorZ2Syst;}
  Float_t GetCMeanSigmaY2p30() const {return fCMeanSigmaY2p30;}
  Float_t GetCMeanSigmaZ2p30() const {return fCMeanSigmaZ2p30;}
  Float_t GetCMeanSigmaY2p30R() const {return fCMeanSigmaY2p30R;}
  Float_t GetCMeanSigmaZ2p30R() const {return fCMeanSigmaZ2p30R;}
     //
  const Double_t *GetLSCovY() const {return fLSCovY;}
  const Double_t *GetLSCovZ() const {return fLSCovZ;}

     //

  Float_t  CookdEdxNorm(Double_t low=0.05, Double_t up=0.70, Int_t type=0, Int_t i1=0, Int_t i2=kMaxRow, Bool_t shapeNorm=kTRUE, Int_t posNorm=0, Int_t padNorm=0,Int_t returnVal=0);

  Float_t  CookdEdxAnalytical(Double_t low=0.05, Double_t up=0.70, Int_t type=0, Int_t i1=0, Int_t i2=kMaxRow, Int_t returnVal=0, Int_t rowThres = 2, Int_t mode=0, TVectorT<float> *returnVec = NULL);

 static   void GetError(AliTPCclusterMI* cluster, AliExternalTrackParam * param, 
			 Double_t& erry, Double_t &errz);
 static   void GetShape(AliTPCclusterMI* cluster, AliExternalTrackParam * param, 
			 Double_t& rmsy, Double_t &rmsz);
  static   Double_t GetQCorrGeom(Float_t ty, Float_t tz);
  static   Double_t GetQCorrShape(Int_t ipad, Int_t type,Float_t z, Float_t ty, Float_t tz, Float_t q, Float_t thr);
  //
  // Float_t GetTPCClustInfo(Int_t nNeighbours, Int_t type, Int_t row0, Int_t row1, TVectorT<float> *returnVec);
  void    SetPoolID(Int_t id) {fPoolID = id;}
  Int_t   GetPoolID()  const {return fPoolID;}
  Int_t   GetNumberOfClustersIndices();  // Should be in AliTPCtrack

  // AliVVTPCseed interface

  void CopyToTPCseed( AliTPCseed &s) const { s = *this; }
  void SetFromTPCseed( const AliTPCseed* seed ) { *this=*seed; }

  void    SetShared(int i)      {fTrackPointsArr.SetShared(i);}
  Bool_t  IsShared(int i) const {return fTrackPointsArr.IsShared(i);}
  //
  Bool_t  GetClusterOwner() const {return fClusterOwner;}
  void    SetClusterOwner(Bool_t v) {fClusterOwner = v;}
  //
  void    TagSuppressSharedClusters();

 private:
     //     AliTPCseed & operator = (const AliTPCseed &)
     //  {::Fatal("= operator","Not Implemented\n");return *this;}
     AliESDtrack * fEsd; //!
     Int_t fNClStore;                       // size of stored cluster pointers array
     AliTPCclusterMI**  fClusterPointer;    //[fNClStore] array of cluster pointers  - 
     Bool_t             fClusterOwner;         // indicates the track is owner of cluster
     //---CURRENT VALUES
     Short_t fRow;               // current row number  
     Char_t  fSector;            // current sector number
     Char_t  fRelativeSector;    // index of current relative sector
     Float_t fCurrentSigmaY2;    //!expected current cluster sigma Y
     Float_t fCurrentSigmaZ2;    //!expected current cluster sigma Z
     Float_t fCMeanSigmaY2p30;   //! current mean sigma Y2 - mean30%
     Float_t fCMeanSigmaZ2p30;   //! current mean sigma Z2 - mean30%
     Float_t fCMeanSigmaY2p30R;   //! current relative mean sigma Y2 - mean30%
     Float_t fCMeanSigmaZ2p30R;   //! current relative mean sigma Z2 - mean30%
     Float_t fErrorY2;           //!sigma of current cluster 
     Float_t fErrorZ2;           //!sigma of current cluster    
     Float_t fErrorY2Syst;        //!syst sigma of current cluster 
     Float_t fErrorZ2Syst;        //!syst sigma of current cluster    
     AliTPCclusterMI * fCurrentCluster; //!pointer to the current cluster for prolongation
     Int_t   fCurrentClusterIndex1; //! index of the current cluster
     UChar_t  fNoCluster;         //!indicates number of rows without clusters
     Char_t   fSort;              //!indicate criteria for sorting
     //
     //
     Float_t fDEDX[9];            // dedx according padrows
     Float_t fSDEDX[4];           // sdedx according padrows
     UChar_t fNCDEDX[4];          // number of clusters for dedx measurment
     UChar_t  fNCDEDXInclThres[4]; // number of clusters for dedx measurment including sub-threshold clusters
     Double_t fTPCr[AliPID::kSPECIES];   // rough PID according TPC   
     Double_t fLSCovY[5];      // sum x^i/sigmaYi^2, i=0,4 for LSM matrix
     Double_t fLSCovZ[3];      // sum x^i/sigmaZi^2, i=0,2 for LSM matrix
     //
     UChar_t   fSeedType;         //seeding type
     UChar_t   fSeed1;            //first row for seeding
     UChar_t   fSeed2;            //last row for seeding
     Char_t   fCircular;           // indicates curlin track
     Int_t   fOverlapLabels[12];  //track labels and the length of the  overlap     
     Float_t fMAngular;           // mean angular factor
     Int_t   fPoolID;              //! id in the pool
     AliTPCTrackerPoints fTrackPointsArr;  // track points - array track points
     ClassDef(AliTPCseed,11)
};





#endif


