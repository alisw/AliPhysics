
#ifndef ALITPCTRACKER_H
#define ALITPCTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------
//                       TPC tracker
//   Parallel tracker 
//
//   Origin: 
//-------------------------------------------------------

#include <TArrayI.h>
#include <TMatrixD.h>
#include "AliTracker.h"
#include "AliTPCreco.h"
#include "AliTPCclusterMI.h"
#include "AliTPCtrackerSector.h"
#include "AliESDfriend.h"
#include "AliTPCseed.h"


class TFile;
class AliTPCParam;
class AliTPCseed;
class AliESDEvent;
class AliESDtrack;
class TTree;
class AliESDkink;
class TTreeSRedirector;
class AliTrackPoint;
class AliDCSSensorArray;
class AliDCSSensor;
class TGraphErrors;


class AliTPCtracker : public AliTracker {
public:
  enum EStreamFlags{ // flags to store addition data/code debugging infomation which is not stored in ESD but in specaial TPCdebug.root file
    kStreamTransform          =0x00001,    // flag:stream cluster transformation 
    kStreamErrParam           =0x00002,    // flag:stream in debug mode cluster and track extrapolation at given row together with error nad shape estimate
    kStreamFilterClusterInfo  =0x00004,    // flag:stream TPC data ouliers filtering information
    kStreamClDump             =0x00008,    // flag:stream clusters at the end of process (signed with useage flags)
    kStreamRemoveUsed         =0x00010,    // flag:stream information about TPC tracks which were descarded (double track removal)
    kStreamRemoveShort        =0x00020,    // flag:stream information about TPC tracks which were discarded (short track removal) 
    kStreamSplitted2          =0x00040,    // flag:stream information about discarded TPC tracks pair algorithm 
    kStreamFillESD            =0x00080,    // flag: stream track information in FillESD function (after track Iteration 0)
    kStreamCPropagateBack     =0x00100,    // flag: stream track information in PropagateBack function (after tracking Iteration 1)
    kStreamRecoverBack        =0x00200,    // flag: stream track information for track  failing PropagateBack function and recovered back
    kStreamRefitInward        =0x00400,    // flag: stream track information in RefitInward function (after tracking Iteration 2)
    kStreamRecoverIn          =0x00800,    // flag: stream track information for track  failing in RefitInward function and recovered back
    kStreamUpdateTrack        =0x01000,    // flag: stream track/cluster infroamtion in track update method
    //
    kStreamCrosstalkMatrix    =0x02000,    // flag: stream crosstalk matrix as used in the reconstruction at given region of TPC
    kStreamXtalk              =0x04000,    // flag: stream crosstalk correction as applied to cluster
    kStreamIonTail            =0x08000,    // flag: stream ion tail correction  as applied to cluster
    kStreamFindMultiMC        =0x10000,    // flag: stream MC infomation about the multiple find track (ONLY for MC data)
    kStreamFindCurling        =0x20000,    // flag: stream track infroamtion in the FindCurling tracks method
    kStreamFindKinks          =0x40000,    // flag: stream track infroamtion in the FindKinks method
    kStreamSeeddEdx           =0x80000,    // flag: stream TPC dEdx intermediate information  AliTPCseed::CookdEdxNorm (to check and validate methods used in calibration
    kStreamOuterDet           =0x100000    // flag: stream matching with outer detectors 
  };
  enum {kMaxFriendTracks=2000};

  AliTPCtracker();
  AliTPCtracker(const AliTPCParam *par); 
  virtual ~AliTPCtracker();
  //
  void SetIteration(Int_t iteration){fIteration = iteration;}
  virtual Int_t Clusters2TracksHLT(AliESDEvent *const esd, const AliESDEvent *hltEvent);
  virtual Int_t Clusters2Tracks (AliESDEvent *const esd);
  virtual Int_t RefitInward (AliESDEvent *esd);
  virtual Int_t LoadClusters (TTree * const tree);
  virtual Int_t LoadClusters (const TObjArray * arr); // another input
  virtual Int_t LoadClusters (const TClonesArray * arr); // another input
  void    FilterOutlierClusters();   // filter outlier clusters  
  virtual Int_t PostProcess(AliESDEvent *esd); 
  Int_t  LoadClusters();
  void   UnloadClusters();
  Int_t LoadInnerSectors();
  Int_t LoadOuterSectors();
  virtual void FillClusterArray(TObjArray* array) const;
  void Transform(AliTPCclusterMI * cluster);
  void ApplyTailCancellation();
  void ApplyXtalkCorrection();
  void CalculateXtalkCorrection();
  void GetTailValue(Float_t ampfactor,Double_t &ionTailMax,Double_t &ionTailTotal,TGraphErrors **graphRes,Float_t *indexAmpGraphs,AliTPCclusterMI *cl0,AliTPCclusterMI *cl1);
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
   void SetSeeds(TObjArray * seeds) { fSeeds = seeds;}
   //   
   AliCluster * GetCluster(Int_t index) const {return (AliCluster*)GetClusterMI(index);}
   AliTPCclusterMI *GetClusterMI(Int_t index) const;
   Int_t Clusters2Tracks();
   virtual void  CookLabel(AliKalmanTrack *tk,Float_t wrong) const; 
   virtual Int_t   CookLabel(AliTPCseed *const t,Float_t wrong, Int_t first,Int_t last ) const; 
   
   void RotateToLocal(AliTPCseed *seed);
   
   Int_t FollowProlongation(AliTPCseed& t, Int_t rf=0, Int_t step=1, Bool_t fromSeeds=0);
   Bool_t GetTrackPoint(Int_t index, AliTrackPoint &p ) const; 
   void   AddSystCovariance(AliTPCseed* t);

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
  
   void    ErrY2Z2(AliTPCseed* seed, const AliTPCclusterMI *cl, double &erry2, double &errz2);
   virtual Double_t ErrY2(AliTPCseed* seed, const AliTPCclusterMI * cl);
   virtual Double_t ErrZ2(AliTPCseed* seed, const AliTPCclusterMI * cl);   

   Double_t F1(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 
   Double_t F1old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 
   Double_t F2(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 
   Double_t F2old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3) const; 

   Double_t F3(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2) const; 
   Double_t F3n(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2, 
                Double_t c) const; 
   Bool_t GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z) const;
   Bool_t GetProlongationLine(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z) const;
   //
   void ResetSeedsPool();
   void MarkSeedFree( TObject* seed );
   TObject *&NextFreeSeed();
   //
   void FillSeedClusterStatCache(const AliTPCseed* seed);
   void GetCachedSeedClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2) const;
   void GetSeedClusterStatistic(const AliTPCseed* seed, Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2) const;
   //
 public:
   void SetUseHLTClusters(Int_t useHLTClusters) {fUseHLTClusters = useHLTClusters;} // set usage from HLT clusters from rec.C options

   inline void SetTPCtrackerSectors(AliTPCtrackerSector *innerSec, AliTPCtrackerSector *outerSec); // set the AliTPCtrackerSector arrays from outside (toy MC)

   Float_t OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t &sum2);
   void  SignShared(AliTPCseed * s1, AliTPCseed * s2);
   void  SignShared(TObjArray * arr);

   void  RemoveUsed2(TObjArray * arr, Float_t factor1, Float_t factor2, Int_t minimal);

   Int_t AcceptCluster(AliTPCseed * seed, AliTPCclusterMI * cluster);

   Bool_t IsTPCHVDipEvent(AliESDEvent const *esdEvent);

   // public for ToyMC usage
   void MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Bool_t bconstrain=kTRUE); 
   void MakeSeeds2Dist(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Bool_t bconstrain=kTRUE); 
   void MakeSeeds3(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Int_t ddsec=0); 
   void MakeSeeds3Dist(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1, Int_t ddsec=0); 
   void SumTracks(TObjArray *arr1,TObjArray *&arr2);
   void SignClusters(const TObjArray * arr, Float_t fnumber=3., Float_t fdensity=2.);  
   //
   Bool_t DistortX(const AliTPCseed* seed, double& x, int row);
   //
   virtual Bool_t OwnsESDObjects() const {return kTRUE;} //RS TPC owns the seeds stored in the friends
   virtual void   CleanESDFriendsObjects(AliESDEvent* esd);
   virtual void   CleanESDTracksObjects(TObjArray* trcList);
   //
   Double_t GetDistortionX(double x, double y, double z, int sec, int row);
   Double_t GetYSectEdgeDist(int sec, int row, double ymax, double z);
   static Int_t GetTrackSector(double alpha);

private:
  Bool_t IsFindable(AliTPCseed & t);
  AliTPCtracker(const AliTPCtracker& r);           //dummy copy constructor
  AliTPCtracker &operator=(const AliTPCtracker& r);//dummy assignment operator
  void AddCovariance(AliTPCseed * seed);               // add covariance
  void AddCovarianceAdd(AliTPCseed * seed);               // add covariance

   inline AliTPCtrackerRow &GetRow(Int_t sec, Int_t row);
   inline Bool_t     IsActive(Int_t sec, Int_t row);
   inline Double_t  GetXrow(Int_t row) const;
   inline Double_t  GetMaxY(Int_t row) const;
   inline Int_t GetRowNumber(Double_t x) const;
   inline Int_t GetRowNumber(const AliTPCseed* seed) const;
   Int_t GetRowNumber(Double_t x[3]) const;
   inline Double_t GetPadPitchLength(Double_t x) const;
   inline Double_t GetPadPitchLength(Int_t row) const;

    void GetShape(AliTPCseed * seed, Int_t row);
 
   void ReadSeeds(const AliESDEvent *const event, Int_t direction);  //read seeds from the event

   void MakeSeeds5(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1);
   void MakeSeeds5Dist(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t cuts[4], Float_t deltay = -1);
  

   AliTPCseed *MakeSeed(AliTPCseed *const track, Float_t r0, Float_t r1, Float_t r2); //reseed
   AliTPCseed *ReSeed(const AliTPCseed *t, Float_t r0, Float_t r1, Float_t r2); //reseed
   AliTPCseed *ReSeed(AliTPCseed *t, Int_t r0, Bool_t forward); //reseed


  
   AliTPCseed * ReSeed(AliTPCseed *t);
   //Int_t LoadInnerSectors();
   //Int_t LoadOuterSectors();
   void DumpClusters(Int_t iter, TObjArray *trackArray);
   void UnsignClusters();

   void FillClusterOccupancyInfo();

   void ParallelTracking(TObjArray *const arr, Int_t rfirst, Int_t rlast);
   void Tracking(TObjArray * arr);
   TObjArray * Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4], Float_t dy=-1, Int_t dsec=0);
   TObjArray * Tracking();
   TObjArray * TrackingSpecial();
   void PrepareForBackProlongation(const TObjArray *const arr, Float_t fac) const;
   void PrepareForProlongation(TObjArray *const arr, Float_t fac) const;

   Int_t UpdateTrack(AliTPCseed *t, Int_t accept); //update trackinfo

   void MakeESDBitmaps(AliTPCseed *t, AliESDtrack *esd);

   Int_t PropagateToRowHLT(AliTPCseed *pt, int nrow);
   void TrackFollowingHLT(TObjArray *const arr);
   TObjArray * MakeSeedsHLT(const AliESDEvent *hltEvent);

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
   const AliESDEvent * fEventHLT;      // input with HLT tracks
   Int_t    fDebug;      // debug option        
   Bool_t   fNewIO;      // indicated if we have data using New IO 
   Int_t fNtracks;                     //current number of tracks
   TObjArray *fSeeds;                  //array of track seeds
   Int_t fIteration;                   // indicate iteration - 0 - froward -1 back - 2forward - back->forward
   //   TObjArray * fTrackPointPool;        // ! pool with track points
   Double_t fXRow[200];                // radius of the pad row
   Double_t fYMax[200];                // max y for given pad row
   Double_t fPadLength[200];                // max y for given pad row
   const AliTPCParam *fkParam;          //pointer to the parameters
   TTreeSRedirector *fDebugStreamer;     //!debug streamer
   Int_t  fUseHLTClusters;              // use HLT clusters instead of offline clusters
   //
   Double_t fClExtraRoadY;              //! extra additiom to Y road for FindCluster
   Double_t fClExtraRoadZ;              //! extra addition Z road for FindCluster
   Double_t fExtraClErrYZ2;             //! extra cl.error Y^2+Z^2
   Double_t fExtraClErrY2;              //! extra cl.error Y^2
   Double_t fExtraClErrZ2;              //! extra cl.error Z^2
   Double_t fPrimaryDCAZCut;            //! special cut on DCAz for primaries tracking only, disables secondaries seeding
   Double_t fPrimaryDCAYCut;            //! special cut on DCAy for primaries tracking only, disables secondaries seeding
   Bool_t   fDisableSecondaries;        //! special flag to disable secondaries seeding
   TObjArray * fCrossTalkSignalArray;  // for 36 sectors 
   AliTPCclusterMI** fClPointersPool;  //! pool of cluster pointers for seeds stored in friends
   AliTPCclusterMI** fClPointersPoolPtr; //! pointer on the current free slot in the pool
   Int_t fClPointersPoolSize; // number of seeds holding the cluster arrays
   TClonesArray* fSeedsPool;            //! pool of seeds
   TClonesArray* fHelixPool;            //! pool of helises
   TClonesArray* fETPPool;              //! pool of helises
   TArrayI fFreeSeedsID;                //! array of ID's of freed seeds
   Int_t fNFreeSeeds;                   //! number of seeds freed in the pool
   Int_t fLastSeedID;                   //! id of the pool seed on which is returned by the NextFreeSeed method
   //
   Bool_t fClStatFoundable[kMaxRow];    //! cached info on foundable clusters of the seed
   Bool_t fClStatFound[kMaxRow];        //! cached info on found clusters of the seed   
   Bool_t fClStatShared[kMaxRow];       //! cached info on shared clusters of the seed   
   //
   Int_t fAccountDistortions;           //! flag to account for distortions. RS: to set!
   //
   ClassDef(AliTPCtracker,5) 
};


AliTPCtrackerRow & AliTPCtracker::GetRow(Int_t sec, Int_t row)
{
  //
  return (row>=fInnerSec->GetNRows()) ? fOuterSec[sec][row-fInnerSec->GetNRows()]:fInnerSec[sec][row];
}

Bool_t   AliTPCtracker::IsActive(Int_t sec, Int_t row)
{
  //
  // check if the given sector row is active 
  //
  return (row>=fInnerSec->GetNRows()) ? fOuterSec[sec][row-fInnerSec->GetNRows()].GetN()>0:fInnerSec[sec][row].GetN()>0;
}


Double_t  AliTPCtracker::GetXrow(Int_t row) const {
  //  return (row>=fInnerSec->GetNRows()) ? fOuterSec->GetX(row-fInnerSec->GetNRows()):fInnerSec->GetX(row);
  return fXRow[row];
}

Double_t  AliTPCtracker::GetMaxY(Int_t row) const {
  //return (row>=fInnerSec->GetNRows()) ? fOuterSec->GetMaxY(row-fInnerSec->GetNRows()):fInnerSec->GetMaxY(row);
  return fYMax[row];
}

Int_t AliTPCtracker::GetRowNumber(Double_t x) const
{
  //
  if (x<133.) return TMath::Max(fInnerSec->GetRowNumber(x),0);
  return TMath::Min(fOuterSec->GetRowNumber(x)+fInnerSec->GetNRows(),158);
}

Int_t AliTPCtracker::GetRowNumber(const AliTPCseed* t) const
{
  // last row of the seed or row corresponding to x
  int row = t->GetRow();
  /*
  int rowX = GetRowNumber(t->GetX());
  if (TMath::Abs(row-rowX)>1 && row>-1) {
    printf("MISMATCH Seed %d\n",t->GetPoolID());
    printf("GetROW: %d %d %f %s\n",row,rowX,t->GetX(),row==rowX ? "":"!!!!!!!!!**********!!!!!!!!!!!************!!!!!!!!!!!");
  }
  */
  return (row>-1&&row<159) ? row : GetRowNumber(t->GetX());
}

Double_t  AliTPCtracker::GetPadPitchLength(Double_t x) const
{
  //
  return (x>133.) ? fOuterSec->GetPadPitchLength(x):fInnerSec->GetPadPitchLength(x);
  //return fPadLength[row];
}

Double_t  AliTPCtracker::GetPadPitchLength(Int_t row) const
{
  //
  return fPadLength[row];
}

void  AliTPCtracker::SetTPCtrackerSectors(AliTPCtrackerSector *innerSec, AliTPCtrackerSector *outerSec)
{
  //
  fInnerSec = innerSec;
  fOuterSec = outerSec;
}

#endif


