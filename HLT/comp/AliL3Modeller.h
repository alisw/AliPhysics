// @(#) $Id$

#ifndef AliL3_Modeller
#define AliL3_Modeller


#include "AliL3RootTypes.h"

class AliL3TrackArray;
class AliL3MemHandler;
class AliL3DigitRowData;
class AliL3ModelTrack;

struct Cluster {
  UInt_t fCharge; // Charge
  UInt_t fPad; // Pad
  UInt_t fTime; // Time
  UInt_t fSigmaY2; // SigmaY2
  UInt_t fSigmaZ2; // SigmaZ2
};

struct Digit {
  Short_t fCharge; // Charge
  Bool_t fUsed; // Flag if used
};

struct ClusterRegion {
  Int_t fMintime; // Min time
  Int_t fMaxtime; // Max time
};

class AliL3Modeller {
  
 public:
  
  AliL3Modeller();
  virtual ~AliL3Modeller();
  
  virtual void FindClusters();
  void Init(Int_t slice,Int_t patch,Char_t *trackdata,Char_t *path,Bool_t houghtracks,Bool_t binary=kTRUE);
  void CheckForOverlaps(Float_t dangle=-1,Int_t *rowrange=0);
  void CalculateCrossingPoints();
  void RemoveBadTracks();
  void WriteRemaining();
  
  void SetInputData(AliL3DigitRowData *digits) {fRowData = digits;}
  void SetTrackThreshold(Int_t i=0) {fTrackThreshold=i;}
  void SetOverlap(Int_t p=6,Int_t t=8) {fPadOverlap=p;fTimeOverlap=t;}
  void SetSearchRange(Int_t p=1,Int_t t=2) {fPadSearch=p;fTimeSearch=t;}
  void SetInnerSearchRange(Int_t p,Int_t t) {fInnerPadSearch=p; fInnerTimeSearch=t;}
  void SetOuterSearchRange(Int_t p,Int_t t) {fOuterPadSearch=p; fOuterTimeSearch=t;}
  void SetMaxClusterRange(Int_t p,Int_t t) {fMaxPads=p; fMaxTimebins=t;}
  void Debug() {fDebug=kTRUE;}
  virtual Bool_t SetFitRange(AliL3ModelTrack *track,Int_t */*p*/,Int_t */*t*/) {return kFALSE;}
  virtual void SetNmaxOverlaps(Int_t /*i*/) {return;}
  virtual void SetChiSqMax(Float_t /*f*/,Int_t /*p*/) {return;}
  
  AliL3TrackArray *GetTracks() {return fTracks;}
    
 protected:
  
  AliL3TrackArray *fTracks; //! Array of tracks
  AliL3DigitRowData *fRowData;//! Row data
  Digit *fRow; //! Current row
  Char_t fPath[1024]; // Path to the files

  Bool_t fDebug; // Flag to switch on/off the debugging
  Int_t fNClusters; // Number of clusters
  Int_t fMaxClusters; // Max clusters (?)
  Int_t fCurrentPadRow; // Current pad row
  Int_t fMaxPads; // Max pads (?)
  Int_t fMaxTimebins; // Max time bins (?)
  Int_t fPadSearch; // Pad search (?)
  Int_t fTimeSearch; // Time search (?)
  Int_t fInnerPadSearch; // Inner pad search (?)
  Int_t fInnerTimeSearch; // Inner time search (?)
  Int_t fOuterPadSearch; // Outer Pad search (?)
  Int_t fOuterTimeSearch; // Outer time search (?)
  
  Int_t fSlice; // Slice
  Int_t fPatch; // Patch
  
  void FillCluster(AliL3ModelTrack *track,Cluster *cluster,Int_t row,Int_t npads);
  void FillZeros(AliL3DigitRowData *digPt,Bool_t reversesign=kFALSE);
  void LocateCluster(AliL3ModelTrack *track,ClusterRegion *region,Int_t &padmin,Int_t &padmax);
  void GetTrackID(Int_t pad,Int_t time,Int_t *trackID);
  
 private:
  Bool_t fHoughTracks; // Flag to switch on/off Hough tracks
  Bool_t CheckCluster(Int_t hitpad,Int_t hittime);
  Float_t fPadOverlap; // Pad overlap (?)
  Float_t fTimeOverlap; // Time overlap (?)
  Int_t fTrackThreshold; //minimum weigth track need in order to be included.(=Nhits/weight)
  AliL3MemHandler *fMemHandler; //! Pointer to the memory handler

  void CalcClusterWidth(Cluster *cl,Float_t &sigmaY2,Float_t &sigmaZ2);
  
  ClassDef(AliL3Modeller,1) //Modeller class
    
};

#endif
