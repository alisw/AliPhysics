#ifndef ALITPCTRACKERMI_H
#define ALITPCTRACKERMI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                       TPC trackerMI
//
//   Origin: 
//-------------------------------------------------------
#include "AliTracker.h"
#include "AliTPCtrack.h"
#include "AliTPCClustersArray.h"

#include "AliTPCreco.h" 
#include "Rtypes.h"

class TFile;
class AliTPCParam;
class AliTPCseed;
class AliTPCclusterMI;
class AliTPCTrackPoint;



class AliTPCclusterTracks {
 public: 
  AliTPCclusterTracks();
  Float_t fDistance[3];   // distance to the 3 nerest track if there overlap with cluster
  Short_t fTrackIndex[3]; // indexes of the  tracks overlapped with clusters 
};

class AliTPCseed;

class AliTPCKalmanSegment: public TObject {
  //
  // class to store tracklet parameters
  // needed to probabilistically define track beginning and track end  
 public:
  AliTPCKalmanSegment();
  void Init(AliTPCseed* seed);             // in initialization initial entrance integral chi2, fNCFoundable and fNC stored 
  void Finish(AliTPCseed* seed);           // in finish state vector stored and chi2 and fNC... calculated
  void GetState(Double_t &x, Double_t & alpha, Double_t state[5]);        
  void GetCovariance(Double_t covariance[15]);
  void GetStatistic(Int_t & nclusters, Int_t & nfoundable, Float_t & chi2); 
 private:
  Float_t fX;               //  x - state
  Float_t fAlpha;           //  Rotation angle the local (TPC sector) 
  Float_t fState[5];        // state vector
  Float_t fChi2;            // chi2 - for given tracklet
  Float_t fCovariance[15];  // 15 elements of covariance matrix
  Int_t fNCFoundable;       // number of foundable clusters on tracklet (effect of dead zone)
  Int_t fNC;                // number of accepted clusters for tracklet  
  //  Int_t fN;                 // total number of padrows for given tracklet
  ClassDef(AliTPCKalmanSegment,1) 
}; 


class AliTPCseed : public AliTPCtrack {
   public:
     AliTPCseed();
     virtual ~AliTPCseed();
     AliTPCseed(const AliTPCtrack &t);
     AliTPCseed(const AliKalmanTrack &t, Double_t a);
     Int_t Compare(const TObject *o) const;
     void Reset();
     Int_t GetProlongation(Double_t xr, Double_t &y, Double_t & z) const;
     virtual Double_t GetPredictedChi2(const AliTPCclusterMI *cluster) const;
     virtual Int_t Update(const AliTPCclusterMI* c, Double_t chi2, UInt_t i);
     AliTPCTrackPoint * GetTrackPoint(Int_t i);
     void RebuildSeed(); // rebuild seed to be ready for storing
     AliTPCseed(UInt_t index, const Double_t xx[5], 
                const Double_t cc[15], Double_t xr, Double_t alpha);
     void SetClusterIndex(Int_t index){
       fClusterIndex[fRow] = index;
     }
     void SetErrorY2(Float_t sy2){fErrorY2=sy2;}
     void SetErrorZ2(Float_t sz2){fErrorZ2=sz2;}
     void CookdEdx(Double_t low=0.05, Double_t up=0.70);
     Bool_t IsActive(){ return !(fRemoval);}
     void Desactivate(Int_t reason){ fRemoval = reason;} 
     
     //     Float_t GetRadius(){ return (1-fP2)/fP4;}  
     Int_t fRelativeSector;  // ! index of current relative sector
     Int_t   fClusterIndex[200];  //array of cluster indexes
     Float_t fClusterDensity[16]; //array with cluster densities 
    
     Int_t fRemoval;               //reason - why was track removed - 0 - means still active     
     TClonesArray * fPoints;              // array with points along the track   
     TClonesArray * fEPoints;             // array with exact points - calculated in special macro not used in tracking
     Int_t fRow;                 //!current row number  
     Int_t fSector;              //!current sector number
     Float_t fCurrentSigmaY;     //!expected current cluster sigma Y
     Float_t fCurrentSigmaZ;     //!expected current cluster sigma Z
     AliTPCclusterMI * fCurrentCluster; //!pointer to the current cluster for prolongation
     Int_t   fCurrentClusterIndex1; //! index of the current cluster
     Int_t   fCurrentClusterIndex2; //! index of the current cluster
    
     Float_t fErrorY2;   //!sigma of current cluster 
     Float_t fErrorZ2;   //!sigma of current cluster    
     Int_t   fNFoundable;      //number of foundable clusters - dead zone taken to the account
     Bool_t  fInDead;         // indicate if the track is in dead zone
     Int_t   fFirstPoint;    // first cluster position
     Int_t   fLastPoint;     // last  cluster position     
     Int_t   fNShared;       // number of shared points
     Bool_t  fIsSeeding;     //indicates if it is proces of seeading
     Bool_t  fStopped;      // indicate that track cann't be prolongate anymore (for secondaries)
   private:
     Float_t fSdEdx;           // sigma of dedx
     Float_t fMAngular;        // mean angular factor
     AliTPCTrackPoint   ** fTrackPoints;  //!track points - array track points
     Float_t fDEDX[4];         // dedx according padrows
     Float_t fSDEDX[4];        // sdedx according padrows
     Int_t   fNCDEDX[4];       // number of clusters for dedx measurment
   
     ClassDef(AliTPCseed,1)  
};




class AliTPCtrackerMI : public AliTracker {
public:
   AliTPCtrackerMI():AliTracker(),fkNIS(0),fkNOS(0) {
      fInnerSec=fOuterSec=0; fSeeds=0; 
   }
   AliTPCtrackerMI(const AliTPCParam *par, Int_t eventn=0);
  ~AliTPCtrackerMI();

   Int_t ReadSeeds(const TFile *in);
   Int_t LoadClusters();
   void UnloadClusters();

   void LoadInnerSectors();
   void LoadOuterSectors();
   AliCluster * GetCluster (int) const {return 0;}
   AliTPCclusterMI *GetClusterMI(Int_t index) const;
   Int_t Clusters2Tracks(const TFile *in, TFile *out);
   //   Int_t PropagateBack(const TFile *in, TFile *out);

   virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const; 
   void RotateToLocal(AliTPCseed *seed);
   virtual Double_t ErrY2(AliTPCseed* seed, AliTPCclusterMI * cl = 0);
   virtual Double_t ErrZ2(AliTPCseed* seed, AliTPCclusterMI * cl = 0);   

   Double_t f1(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t f2(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t f3(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2); 
public:
//**************** Internal tracker class ********************** 
   class AliTPCRow {
   public:
     AliTPCRow() {fN=0; fClusterTracks=0;}
     ~AliTPCRow();
     void InsertCluster(const AliTPCclusterMI *c, UInt_t index);
     operator int() const {return fN;}
     const AliTPCclusterMI* operator[](Int_t i) const {return fClusters[i];}
     UInt_t GetIndex(Int_t i) const {return fIndex[i];}
     Int_t Find(Double_t z) const; 
     AliTPCclusterMI *  FindNearest(Double_t y, Double_t z, Double_t roady, Double_t roadz) const;
     void SetX(Double_t x) {fX=x;}
     Double_t GetX() const {return fX;}
     AliTPCclusterTracks *  GetClusterTracks(Int_t index){ return ( (index<fN) && fClusterTracks!=0)? &(fClusterTracks[index]):0;}
     void UpdateClusterTrack(Int_t clindex, Int_t trindex,AliTPCseed * seed); 
     void MakeClusterTracks();
     void ClearClusterTracks();
     Float_t fDeadZone;  // the width of the dead zone
   private:
     Int_t fN;                                          //number of clusters 
     const AliTPCclusterMI *fClusters[kMaxClusterPerRow]; //pointers to clusters
     UInt_t fIndex[kMaxClusterPerRow];                  //indeces of clusters
     Double_t fX;                                 //X-coordinate of this row
     AliTPCclusterTracks * fClusterTracks;        // array of cluster tracks - for overlap calculation
   private:
     AliTPCRow(const AliTPCRow& r);            //dummy copy constructor
     AliTPCRow &operator=(const AliTPCRow& r); //dummy assignment operator
   };

//**************** Internal tracker class ********************** 
   class AliTPCSector {
   public:
     AliTPCSector() { fN=0; fRow = 0; }
    ~AliTPCSector() { delete[] fRow; }
     AliTPCRow& operator[](Int_t i) const { return *(fRow+i); }
     Int_t GetNRows() const { return fN; }
     void Setup(const AliTPCParam *par, Int_t flag);
     Double_t GetX(Int_t l) const {return fRow[l].GetX();}
     Double_t GetMaxY(Int_t l) const {
         return GetX(l)*TMath::Tan(0.5*GetAlpha());
     } 
     Double_t GetAlpha() const {return fAlpha;}
     Double_t GetAlphaShift() const {return fAlphaShift;}     
     Int_t GetRowNumber(Double_t x) const {
        //return pad row number for this x
       Double_t r;
       if (fN < 64){
         r=fRow[fN-1].GetX();
         if (x > r) return fN;
         r=fRow[0].GetX();
         if (x < r) return -1;
         return Int_t((x-r)/fPadPitchLength + 0.5);}
       else{    
           r=fRow[fN-1].GetX();
           if (x > r) return fN;
           r=fRow[0].GetX();
           if (x < r) return -1;
          Double_t r1=fRow[64].GetX();
          if(x<r1){       
            return Int_t((x-r)/f1PadPitchLength + 0.5);}
          else{
            return (Int_t((x-r1)/f2PadPitchLength + 0.5)+64);} 
       }
     }
     Double_t GetPadPitchWidth()  const {return fPadPitchWidth;}
     Double_t GetPadPitchLength() const {return fPadPitchLength;}
     Double_t GetPadPitchLength(Float_t x) const {return (x<200) ? fPadPitchLength:f2PadPitchLength ;}

   private:
     Int_t fN;                        //number of pad rows 
     AliTPCRow *fRow;                    //array of pad rows
     Double_t fAlpha;                    //opening angle
     Double_t fAlphaShift;               //shift angle;
     Double_t fPadPitchWidth;            //pad pitch width
     Double_t fPadPitchLength;           //pad pitch length
     Double_t f1PadPitchLength;           //pad pitch length
     Double_t f2PadPitchLength;           //pad pitch length
    
   private:
     AliTPCSector(const AliTPCSector &s);           //dummy copy contructor
     AliTPCSector& operator=(const AliTPCSector &s);//dummy assignment operator
   };

   Float_t OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t &sum2);
   void  SignShared(AliTPCseed * s1, AliTPCseed * s2);
   void  RemoveOverlap(TObjArray * arr, Float_t factor, Int_t removalindex, Bool_t shared=kFALSE);
   void  RemoveUsed(TObjArray * arr, Float_t factor, Int_t removalindex);

private:
   Float_t  GetSigmaY(AliTPCseed * seed);
   Float_t  GetSigmaZ(AliTPCseed * seed);

   void MakeSeeds(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2);
   void MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2);

   TObjArray * MakeSeedsSectors(Int_t sec1, Int_t sec2);   // make seeds from all sectors
   void MakeSeedsAll();
   Int_t FollowProlongation(AliTPCseed& t, Int_t rf=0);
   //Int_t FollowProlongationFast(AliTPCseed& t, Int_t step);   
   Int_t FollowBackProlongation(AliTPCseed& t, Int_t rf=0);

   Int_t FollowToNext(AliTPCseed& t, Int_t nr);
   Int_t UpdateClusters(AliTPCseed& t, Int_t trindex,  Int_t nr);
   Int_t FollowToNextCluster( Int_t trindex, Int_t nr);

   virtual Int_t PropagateBack (const TFile *, TFile *){return 0;}
   void ParallelTracking(Int_t rfirst, Int_t rlast);
   void SetSampledEdx(AliTPCseed *t, Float_t q, Int_t i) {;}
   Int_t UpdateTrack(AliTPCseed *t,AliTPCclusterMI* c, Double_t chi2, UInt_t i); //update trackinfo

   //   Int_t FollowBackProlongation(AliTPCseed &s, const AliTPCtrack &t);

   AliTPCtrackerMI(const AliTPCtrackerMI& r);           //dummy copy constructor
   AliTPCtrackerMI &operator=(const AliTPCtrackerMI& r);//dummy assignment operator

   const Int_t fkNIS;        //number of inner sectors
   AliTPCSector *fInnerSec;  //array of inner sectors;
   const Int_t fkNOS;        //number of outer sectors
   AliTPCSector *fOuterSec;  //array of outer sectors;

   Int_t fN;               //number of loaded sectors
   AliTPCSector *fSectors; //pointer to loaded sectors;

   Int_t fEventN;                      //event number
   AliTPCClustersArray fClustersArray; //array of TPC clusters
   Int_t fNtracks;                     //current number of tracks
   TObjArray *fSeeds;                  //array of track seeds
   //   TObjArray * fTrackPointPool;        // ! pool with track points
   const AliTPCParam *fParam;          //pointer to the parameters
};

#endif


