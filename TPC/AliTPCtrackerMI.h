#ifndef ALITPCTRACKERMI_H
#define ALITPCTRACKERMI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------
//                       TPC tracker
//   Parallel tracker 
//
//   Origin: 
//-------------------------------------------------------

#include <TError.h>
#include "AliTracker.h"
#include "AliTPCreco.h"
#include "AliPID.h"
#include "AliTPCclusterMI.h"


class TFile;
class AliTPCParam;
class AliTPCseed;
class AliTPCTrackerPoint;
class AliESD;   
class TTree;
class AliESDkink;
class TTreeSRedirector;
class AliTrackPoint;

class AliTPCtrackerMI : public AliTracker {
public:
  AliTPCtrackerMI():AliTracker(),fkNIS(0),fkNOS(0) {
    fInnerSec=fOuterSec=0; fSeeds=0; 
  }
  AliTPCtrackerMI(const AliTPCParam *par); 
  virtual ~AliTPCtrackerMI();
  //
  void SetIteration(Int_t iteration){fIteration = iteration;}
  virtual Int_t Clusters2Tracks (AliESD *esd);
  virtual Int_t RefitInward (AliESD *esd);
  virtual Int_t LoadClusters (TTree * tree);
  Int_t  LoadClusters();
  void   UnloadClusters();
  void   Transform(AliCluster * cluster);
  //
  void SetIO();  //set default IO from folders
  void SetIO(TTree * input, TTree * output, AliESD * event);
  void FillESD(TObjArray* arr);
  void WriteTracks();
  void WriteTracks(TTree * tree);  
  void DeleteSeeds();
  void SetDebug(Int_t debug){ fDebug = debug;}
  void FindKinks(TObjArray * array, AliESD * esd);
  void FindV0s(TObjArray * array, AliESD * esd);
  void UpdateKinkQualityM(AliTPCseed * seed);
  void UpdateKinkQualityD(AliTPCseed * seed);
  Int_t CheckKinkPoint(AliTPCseed*seed, AliTPCseed &mother, AliTPCseed &daughter, AliESDkink &kink);
  Int_t RefitKink(AliTPCseed &mother, AliTPCseed &daughter, AliESDkink &kink);
   Int_t ReadSeeds(const TFile *in);
   TObjArray * GetSeeds(){return fSeeds;}
   //   
   AliCluster * GetCluster(Int_t index) const {return (AliCluster*)GetClusterMI(index);}
   AliTPCclusterMI *GetClusterMI(Int_t index) const;
   Int_t Clusters2Tracks();
   virtual void  CookLabel(AliKalmanTrack *tk,Float_t wrong) const; 
   virtual Int_t   CookLabel(AliTPCseed *t,Float_t wrong, Int_t first,Int_t last ) const; 
   
   void RotateToLocal(AliTPCseed *seed);
  
   
   Int_t FollowProlongation(AliTPCseed& t, Int_t rf=0, Int_t step=1);
   Int_t FollowProlongationFast(AliTPCseed& t, Int_t rf=0, Int_t step=1);
   Bool_t GetTrackPoint(Int_t index, AliTrackPoint &p ) const; 

   Int_t FollowBackProlongation(AliTPCseed& t, Int_t rf);
   Int_t FollowToNext(AliTPCseed& t, Int_t nr);
   Int_t FollowToNextFast(AliTPCseed& t, Int_t nr);
   Int_t UpdateClusters(AliTPCseed& t,  Int_t nr);
   Int_t FollowToNextCluster( AliTPCseed& t, Int_t nr);

   Int_t PropagateBack(TObjArray * arr);
   Int_t PropagateBack(AliESD * event);
   Int_t PropagateBack(AliTPCseed *pt, Int_t row0, Int_t row1);   
   Int_t PropagateForward();
   Int_t PropagateForward2(TObjArray * arr);

   void SortTracks(TObjArray * arr, Int_t mode) const;
  

   virtual Double_t ErrY2(AliTPCseed* seed, AliTPCclusterMI * cl = 0);
   virtual Double_t ErrZ2(AliTPCseed* seed, AliTPCclusterMI * cl = 0);   

   Double_t F1(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t F1old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t F2(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t F2old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 

   Double_t F3(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2); 
   Double_t F3n(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2, 
                Double_t c); 
   Bool_t GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z);

 public:
//**************** Internal tracker class ********************** 
   class AliTPCSector;
   class AliTPCRow {
   public:
     AliTPCRow();
     ~AliTPCRow();
     void InsertCluster(const AliTPCclusterMI *c, UInt_t index);
     void ResetClusters();
     operator int() const {return fN;}
     Int_t GetN() const {return fN;}
     const AliTPCclusterMI* operator[](Int_t i) const {return fClusters[i];}
     UInt_t GetIndex(Int_t i) const {return fIndex[i];}
     inline Int_t Find(Double_t z) const; 
     AliTPCclusterMI *  FindNearest(Double_t y, Double_t z, Double_t roady, Double_t roadz) const;
     AliTPCclusterMI *  FindNearest2(Double_t y, Double_t z, Double_t roady, Double_t roadz, UInt_t & index) const;
     AliTPCclusterMI *  FindNearest3(Double_t y, Double_t z, Double_t roady, Double_t roadz, UInt_t & index) const;
     
     void SetX(Double_t x) {fX=x;}
     Double_t GetX() const {return fX;}
     Float_t GetDeadZone() const {return fDeadZone;}
     void SetDeadZone(Float_t d) {fDeadZone=d;}
     Int_t GetN1() const {return fN1;}
     void SetN1(Int_t n) {fN1=n;}
     Int_t GetN2() const {return fN2;}
     void SetN2(Int_t n) {fN2=n;}
     AliTPCclusterMI* GetClusters1() const {return fClusters1;}
     AliTPCclusterMI* GetClusters2() const {return fClusters2;}
     void SetClusters1(AliTPCclusterMI* cl) {fClusters1=cl;}
     void SetClusters2(AliTPCclusterMI* cl) {fClusters2=cl;}
     void SetCluster1(Int_t i, AliTPCclusterMI cl) {fClusters1[i]=cl;}
     void SetCluster2(Int_t i, AliTPCclusterMI cl) {fClusters2[i]=cl;}
     AliTPCclusterMI* GetCluster1(Int_t i) const {return &fClusters1[i];}
     AliTPCclusterMI* GetCluster2(Int_t i) const {return &fClusters2[i];}
     Short_t GetFastCluster(Int_t i) const {return fFastCluster[i];}
     void SetFastCluster(Int_t i, Short_t cl) {fFastCluster[i]=cl;}

private:  
     AliTPCRow & operator=(const AliTPCRow & );
     AliTPCRow(const AliTPCRow& /*r*/);           //dummy copy constructor
     Float_t fDeadZone;  // the width of the dead zone
     AliTPCclusterMI *fClusters1; //array with clusters 1
     Int_t fN1;  //number of clusters on left side
     AliTPCclusterMI *fClusters2; //array with clusters 2
     Int_t fN2; // number of clusters on right side of the TPC
     Short_t fFastCluster[510];   //index of the nearest cluster at given position
     Int_t fN;                                          //number of clusters 
     const AliTPCclusterMI *fClusters[kMaxClusterPerRow]; //pointers to clusters
                               // indexes for cluster at given position z  
     // AliTPCclusterMI *fClustersArray;                     // 
     UInt_t fIndex[kMaxClusterPerRow];                  //indeces of clusters
     Double_t fX;                                 //X-coordinate of this row

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
    //Int_t GetFirst(){return fFirstRow;}
    Int_t GetRowNumber(Double_t  x) const;
    Double_t GetPadPitchWidth()  const {return fPadPitchWidth;}
    Double_t GetPadPitchLength() const {return fPadPitchLength;}
    Double_t GetPadPitchLength(Float_t x) const {return (x<200) ? fPadPitchLength:f2PadPitchLength ;}
    
   private:
    AliTPCSector & operator=(const AliTPCSector & );
    AliTPCSector(const AliTPCSector &/*s*/);           //dummy copy contructor 
    Int_t fN;                        //number of pad rows 
    //Int_t fFirstRow;                 //offset
    AliTPCRow *fRow;                    //array of pad rows
    Double_t fAlpha;                    //opening angle
    Double_t fAlphaShift;               //shift angle;
    Double_t fPadPitchWidth;            //pad pitch width
    Double_t fPadPitchLength;           //pad pitch length
    Double_t f1PadPitchLength;           //pad pitch length
    Double_t f2PadPitchLength;           //pad pitch length
    
   };

   Float_t OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t &sum2);
   void  SignShared(AliTPCseed * s1, AliTPCseed * s2);
   void  SignShared(TObjArray * arr);

   void  RemoveUsed(TObjArray * arr, Float_t factor1, Float_t factor2,  Int_t removalindex);
   void  RemoveUsed2(TObjArray * arr, Float_t factor1, Float_t factor2, Int_t minimal);
   void  RemoveDouble(TObjArray * arr, Float_t factor1, Float_t factor2,  Int_t removalindex);

   void  StopNotActive(TObjArray * arr, Int_t row0, Float_t th0, Float_t th1, Float_t th2) const;
   void  StopNotActive(AliTPCseed * seed, Int_t row0, Float_t th0, Float_t th1, Float_t th2) const;
   Int_t AcceptCluster(AliTPCseed * seed, AliTPCclusterMI * cluster, Float_t factor, Float_t cory=1., Float_t corz=1.);

private:
  AliTPCtrackerMI(const AliTPCtrackerMI& r);           //dummy copy constructor
  AliTPCtrackerMI &operator=(const AliTPCtrackerMI& r);//dummy assignment operator
   inline AliTPCRow &GetRow(Int_t sec, Int_t row);
   inline Bool_t     IsActive(Int_t sec, Int_t row);
   inline Double_t  GetXrow(Int_t row) const;
   inline Double_t  GetMaxY(Int_t row) const;
   inline Int_t GetRowNumber(Double_t x) const;
   Int_t GetRowNumber(Double_t x[3]) const;
   inline Double_t GetPadPitchLength(Double_t x) const;
   inline Double_t GetPadPitchLength(Int_t row) const;

   Float_t  GetSigmaY(AliTPCseed * seed);
   Float_t  GetSigmaZ(AliTPCseed * seed);
   void GetShape(AliTPCseed * seed, Int_t row);
 
   void ReadSeeds(AliESD *event, Int_t direction);  //read seeds from the event

   void MakeSeeds3(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Int_t ddsec=0); 
   void MakeSeeds5(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1);

   void MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Bool_t bconstrain=kTRUE);
  

   AliTPCseed *MakeSeed(AliTPCseed *t, Float_t r0, Float_t r1, Float_t r2); //reseed
   AliTPCseed *ReSeed(AliTPCseed *t, Float_t r0, Float_t r1, Float_t r2); //reseed
   AliTPCseed *ReSeed(AliTPCseed *t, Int_t r0, Bool_t forward); //reseed


  
   AliTPCseed * ReSeed(AliTPCseed *t);
   Int_t LoadInnerSectors();
   Int_t LoadOuterSectors();
   void UnsignClusters();
   void SignClusters(TObjArray * arr, Float_t fnumber=3., Float_t fdensity=2.);  

   void ParallelTracking(TObjArray * arr, Int_t rfirst, Int_t rlast);
   void Tracking(TObjArray * arr);
   TObjArray * Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4], Float_t dy=-1, Int_t dsec=0);
   TObjArray * Tracking();
   TObjArray * TrackingSpecial();
   void SumTracks(TObjArray *arr1,TObjArray *arr2) const;
   void PrepareForBackProlongation(TObjArray * arr, Float_t fac) const;
   void PrepareForProlongation(TObjArray * arr, Float_t fac) const;

   void SetSampledEdx(AliTPCseed */*t*/, Float_t /*q*/, Int_t /*i*/) {;}
   Int_t UpdateTrack(AliTPCseed *t, Int_t accept); //update trackinfo


   const Int_t fkNIS;        //number of inner sectors
   AliTPCSector *fInnerSec;  //array of inner sectors;
   const Int_t fkNOS;        //number of outer sectors
   AliTPCSector *fOuterSec;  //array of outer sectors;

   Int_t fN;               //number of loaded sectors
   AliTPCSector *fSectors; //pointer to loaded sectors;
   //
   TTree * fInput;       // input tree with clusters
   TTree * fOutput;      // output tree with tracks
   TTree * fSeedTree;    // output tree with seeds - filled in debug mode 1
   TTree * fTreeDebug;   // output with a debug information about track
   AliESD * fEvent;      // output with esd tracks
   Int_t    fDebug;      // debug option        
   Bool_t   fNewIO;      // indicated if we have data using New IO 
   Int_t fNtracks;                     //current number of tracks
   TObjArray *fSeeds;                  //array of track seeds
   Int_t fIteration;                   // indicate iteration - 0 - froward -1 back - 2forward - back->forward
   //   TObjArray * fTrackPointPool;        // ! pool with track points
   //   TObjArray * fSeedPool;              //! pool with seeds
   Double_t fXRow[200];                // radius of the pad row
   Double_t fYMax[200];                // max y for given pad row
   Double_t fPadLength[200];                // max y for given pad row
   const AliTPCParam *fParam;          //pointer to the parameters
   TTreeSRedirector *fDebugStreamer;     //!debug streamer
   ClassDef(AliTPCtrackerMI,1) 
};


AliTPCtrackerMI::AliTPCRow & AliTPCtrackerMI::GetRow(Int_t sec, Int_t row)
{
  //
  return (row>=fInnerSec->GetNRows()) ? fOuterSec[sec][row-fInnerSec->GetNRows()]:fInnerSec[sec][row];
}

Bool_t   AliTPCtrackerMI::IsActive(Int_t sec, Int_t row)
{
  //
  // check if the given sector row is active 
  //
  return (row>=fInnerSec->GetNRows()) ? fOuterSec[sec][row-fInnerSec->GetNRows()].GetN()>0:fInnerSec[sec][row].GetN()>0;
}


Double_t  AliTPCtrackerMI::GetXrow(Int_t row) const {
  //  return (row>=fInnerSec->GetNRows()) ? fOuterSec->GetX(row-fInnerSec->GetNRows()):fInnerSec->GetX(row);
  return fXRow[row];
}

Double_t  AliTPCtrackerMI::GetMaxY(Int_t row) const {
  //return (row>=fInnerSec->GetNRows()) ? fOuterSec->GetMaxY(row-fInnerSec->GetNRows()):fInnerSec->GetMaxY(row);
  return fYMax[row];
}

Int_t AliTPCtrackerMI::GetRowNumber(Double_t x) const
{
  //
  return (x>133.) ? fOuterSec->GetRowNumber(x)+fInnerSec->GetNRows():fInnerSec->GetRowNumber(x);
}

Double_t  AliTPCtrackerMI::GetPadPitchLength(Double_t x) const
{
  //
  return (x>133.) ? fOuterSec->GetPadPitchLength(x):fInnerSec->GetPadPitchLength(x);
  //return fPadLength[row];
}

Double_t  AliTPCtrackerMI::GetPadPitchLength(Int_t row) const
{
  //
  return fPadLength[row];
}



#endif


