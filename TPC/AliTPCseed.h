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

#include "AliTPCtrack.h"
#include "AliComplexCluster.h"
#include "AliPID.h"

class TFile;
class AliTPCParam;
class AliTPCseed;
class AliTPCclusterMI;
class AliTPCTrackerPoint;
class AliESD;   
class TClonesArray;

class AliTPCseed : public AliTPCtrack {
  friend class AliTPCtrackerMI;
  public:  
     AliTPCseed();
     virtual ~AliTPCseed();
     AliTPCseed(const AliTPCtrack &t);
     AliTPCseed(const AliTPCseed &s, Bool_t clusterOwner = kFALSE);
     //AliTPCseed(const AliTPCseed &t, Double_t a);
     AliTPCseed(Double_t xr, Double_t alpha, const Double_t xx[5], 
                const Double_t cc[15], Int_t i);     
     Int_t Compare(const TObject *o) const;
     void Reset(Bool_t all = kTRUE);
     Int_t GetProlongation(Double_t xr, Double_t &y, Double_t & z) const;
     virtual Double_t GetPredictedChi2(const AliCluster *cluster2) const;
     virtual Bool_t Update(const AliCluster* c2, Double_t chi2, Int_t i);
     AliTPCTrackerPoint * GetTrackPoint(Int_t i);
     AliTPCclusterMI * GetClusterFast(Int_t irow){ return fClusterPointer[irow];}
     void RebuildSeed(); // rebuild seed to be ready for storing
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
     void CookPID();
     Double_t Bethe(Double_t bg);     // return bethe-bloch
     //     void CookdEdx2(Double_t low=0.05, Double_t up=0.70);
     Bool_t IsActive() const { return !(fRemoval);}
     void Desactivate(Int_t reason){ fRemoval = reason;} 
     //
     //
 private:
     //     AliTPCseed & operator = (const AliTPCseed &)
     //  {::Fatal("= operator","Not Implemented\n");return *this;}
     AliESDtrack * fEsd; //!
     AliTPCclusterMI*   fClusterPointer[160];  // array of cluster pointers  - 
     Bool_t             fClusterOwner;         // indicates the track is owner of cluster
     TClonesArray * fPoints;              //!array with points along the track
     TClonesArray * fEPoints;             //! array with exact points - calculated in special macro not used in tracking
     //---CURRENT VALUES
     Int_t fRow;                 //!current row number  
     Int_t fSector;              //!current sector number
     Int_t fRelativeSector;      //! index of current relative sector
     Float_t fCurrentSigmaY2;    //!expected current cluster sigma Y
     Float_t fCurrentSigmaZ2;    //!expected current cluster sigma Z
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
     Float_t fDEDX[4];         // dedx according padrows
     Float_t fSDEDX[4];        // sdedx according padrows
     Int_t   fNCDEDX[4];       // number of clusters for dedx measurment
     Double_t fTPCr[AliPID::kSPECIES];   // rough PID according TPC   
     //
     Int_t   fSeedType;         //seeding type
     Int_t   fSeed1;            //first row for seeding
     Int_t   fSeed2;            //last row for seeding
     Int_t   fOverlapLabels[12];  //track labels and the length of the  overlap     
     Float_t fMAngular;           // mean angular factor
     Char_t   fCircular;           // indicates curlin track
     AliTPCTrackerPoint  fTrackPoints[160];  //track points - array track points
     ClassDef(AliTPCseed,1)  
};





#endif


