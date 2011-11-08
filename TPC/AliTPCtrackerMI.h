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

#include "AliTracker.h"
#include "AliTPCreco.h"
#include "AliTPCclusterMI.h"
#include "AliTPCtrackerSector.h"


class TFile;
class AliTPCParam;
class AliTPCseed;
class AliTPCTrackerPoint;
class AliESDEvent;
class AliESDtrack;
class TTree;
class AliESDkink;
class TTreeSRedirector;
class AliTrackPoint;



class AliTPCtrackerMI : public AliTracker {
public:
  AliTPCtrackerMI();
  AliTPCtrackerMI(const AliTPCParam *par); 
  virtual ~AliTPCtrackerMI();
  //
  void SetIteration(Int_t iteration){fIteration = iteration;}
  virtual Int_t Clusters2Tracks (AliESDEvent *const esd);
  virtual Int_t RefitInward (AliESDEvent *esd);
  virtual Int_t LoadClusters (TTree * const tree);
  virtual Int_t LoadClusters (const TObjArray * arr); // another input
  virtual Int_t LoadClusters (const TClonesArray * arr); // another input
  Int_t  LoadClusters();
  void   UnloadClusters();
  Int_t LoadInnerSectors();
  Int_t LoadOuterSectors();
  virtual void FillClusterArray(TObjArray* array) const;
  void   Transform(AliTPCclusterMI * cluster);
  //
  void FillESD(const TObjArray* arr);
  void DeleteSeeds();
  void SetDebug(Int_t debug){ fDebug = debug;}
  void FindKinks(TObjArray * array, AliESDEvent * esd);
  //
  void FindCurling(const TObjArray * array, AliESDEvent * esd, Int_t iter);     
  void FindSplitted(TObjArray * array, AliESDEvent * esd, Int_t iter);       
  void FindMultiMC(const TObjArray * array, AliESDEvent * esd, Int_t iter);     
  //
  void UpdateKinkQualityM(AliTPCseed * seed);
  void UpdateKinkQualityD(AliTPCseed * seed);
  Int_t CheckKinkPoint(AliTPCseed*seed, AliTPCseed &mother, AliTPCseed &daughter, const AliESDkink &kink);
  Int_t RefitKink(AliTPCseed &mother, AliTPCseed &daughter, const AliESDkink &kink);
   Int_t ReadSeeds(const TFile *in);
   TObjArray * GetSeeds() const {return fSeeds;}
   //   
   AliCluster * GetCluster(Int_t index) const {return (AliCluster*)GetClusterMI(index);}
   AliTPCclusterMI *GetClusterMI(Int_t index) const;
   Int_t Clusters2Tracks();
   virtual void  CookLabel(AliKalmanTrack *tk,Float_t wrong) const; 
   virtual Int_t   CookLabel(AliTPCseed *const t,Float_t wrong, Int_t first,Int_t last ) const; 
   
   void RotateToLocal(AliTPCseed *seed);
   
   Int_t FollowProlongation(AliTPCseed& t, Int_t rf=0, Int_t step=1, Bool_t fromSeeds=0);
   Bool_t GetTrackPoint(Int_t index, AliTrackPoint &p ) const; 

   Int_t FollowBackProlongation(AliTPCseed& t, Int_t rf, Bool_t fromSeeds=0);
   Int_t FollowToNext(AliTPCseed& t, Int_t nr);
   Int_t UpdateClusters(AliTPCseed& t,  Int_t nr);
   Int_t FollowToNextCluster( AliTPCseed& t, Int_t nr);

   Int_t PropagateBack(const TObjArray *const arr);
   Int_t PropagateBack(AliESDEvent * event);
   Int_t PropagateBack(AliTPCseed *const pt, Int_t row0, Int_t row1);   
   Int_t PropagateForward();
   Int_t PropagateForward2(const TObjArray *const arr);

   void SortTracks(TObjArray * arr, Int_t mode) const;
  
   virtual Double_t ErrY2(AliTPCseed* seed, const AliTPCclusterMI * cl = 0);
   virtual Double_t ErrZ2(AliTPCseed* seed, const AliTPCclusterMI * cl = 0);   

   Double_t F1(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 
   Double_t F1old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 
   Double_t F2(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 
   Double_t F2old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 

   Double_t F3(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2) const; 
   Double_t F3n(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2, 
                Double_t c) const; 
   Bool_t GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z) const;

 public:
   void SetUseHLTClusters(Int_t useHLTClusters) {fUseHLTClusters = useHLTClusters;} // set usage from HLT clusters from rec.C options

   Float_t OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t &sum2);
   void  SignShared(AliTPCseed * s1, AliTPCseed * s2);
   void  SignShared(TObjArray * arr);

   void  RemoveUsed2(TObjArray * arr, Float_t factor1, Float_t factor2, Int_t minimal);

   void  StopNotActive(const TObjArray * arr, Int_t row0, Float_t th0, Float_t th1, Float_t th2) const;
   void  StopNotActive(AliTPCseed * seed, Int_t row0, Float_t th0, Float_t th1, Float_t th2) const;
   Int_t AcceptCluster(AliTPCseed * seed, AliTPCclusterMI * cluster);

private:
  Bool_t IsFindable(AliTPCseed & t);
  AliTPCtrackerMI(const AliTPCtrackerMI& r);           //dummy copy constructor
  AliTPCtrackerMI &operator=(const AliTPCtrackerMI& r);//dummy assignment operator
  void AddCovariance(AliTPCseed * seed);               // add covariance
  void AddCovarianceAdd(AliTPCseed * seed);               // add covariance

   inline AliTPCtrackerRow &GetRow(Int_t sec, Int_t row);
   inline Bool_t     IsActive(Int_t sec, Int_t row);
   inline Double_t  GetXrow(Int_t row) const;
   inline Double_t  GetMaxY(Int_t row) const;
   inline Int_t GetRowNumber(Double_t x) const;
   Int_t GetRowNumber(Double_t x[3]) const;
   inline Double_t GetPadPitchLength(Double_t x) const;
   inline Double_t GetPadPitchLength(Int_t row) const;

    void GetShape(AliTPCseed * seed, Int_t row);
 
   void ReadSeeds(const AliESDEvent *const event, Int_t direction);  //read seeds from the event

   void MakeSeeds3(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Int_t ddsec=0); 
   void MakeSeeds5(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1);

   void MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Bool_t bconstrain=kTRUE);
  

   AliTPCseed *MakeSeed(AliTPCseed *const track, Float_t r0, Float_t r1, Float_t r2); //reseed
   AliTPCseed *ReSeed(const AliTPCseed *t, Float_t r0, Float_t r1, Float_t r2); //reseed
   AliTPCseed *ReSeed(AliTPCseed *t, Int_t r0, Bool_t forward); //reseed


  
   AliTPCseed * ReSeed(AliTPCseed *t);
   //Int_t LoadInnerSectors();
   //Int_t LoadOuterSectors();
   void DumpClusters(Int_t iter, TObjArray *trackArray);
   void UnsignClusters();
   void SignClusters(const TObjArray * arr, Float_t fnumber=3., Float_t fdensity=2.);  

   void ParallelTracking(TObjArray *const arr, Int_t rfirst, Int_t rlast);
   void Tracking(TObjArray * arr);
   TObjArray * Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4], Float_t dy=-1, Int_t dsec=0);
   TObjArray * Tracking();
   TObjArray * TrackingSpecial();
   void SumTracks(TObjArray *arr1,TObjArray *&arr2) const;
   void PrepareForBackProlongation(const TObjArray *const arr, Float_t fac) const;
   void PrepareForProlongation(TObjArray *const arr, Float_t fac) const;

   Int_t UpdateTrack(AliTPCseed *t, Int_t accept); //update trackinfo

   void MakeESDBitmaps(AliTPCseed *t, AliESDtrack *esd);

   const Int_t fkNIS;        //number of inner sectors
   AliTPCtrackerSector *fInnerSec;  //array of inner sectors;
   const Int_t fkNOS;        //number of outer sectors
   AliTPCtrackerSector *fOuterSec;  //array of outer sectors;

   Int_t fN;               //number of loaded sectors
   AliTPCtrackerSector *fSectors; //pointer to loaded sectors;
   //
   TTree * fInput;       // input tree with clusters
   TTree * fOutput;      // output tree with tracks
   TTree * fSeedTree;    // output tree with seeds - filled in debug mode 1
   TTree * fTreeDebug;   // output with a debug information about track
   AliESDEvent * fEvent;      // output with esd tracks
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
   const AliTPCParam *fkParam;          //pointer to the parameters
   TTreeSRedirector *fDebugStreamer;     //!debug streamer
   Int_t  fUseHLTClusters;              // use HLT clusters instead of offline clusters
   ClassDef(AliTPCtrackerMI,2) 
};


AliTPCtrackerRow & AliTPCtrackerMI::GetRow(Int_t sec, Int_t row)
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


