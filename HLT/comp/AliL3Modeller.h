// @(#) $Id$

#ifndef AliL3_Modeller
#define AliL3_Modeller


#include "AliL3RootTypes.h"

class AliL3TrackArray;
class AliL3MemHandler;
class AliL3DigitRowData;
class AliL3ModelTrack;

struct Cluster {
  UInt_t fCharge;
  UInt_t fPad;
  UInt_t fTime;
  UInt_t fSigmaY2;
  UInt_t fSigmaZ2;
};

struct Digit {
  Short_t fCharge;
  Bool_t fUsed;
};

struct ClusterRegion {
  Int_t mintime;
  Int_t maxtime;
};

class AliL3Modeller {
  
 private:
  Bool_t fHoughTracks;
  Bool_t CheckCluster(Int_t hitpad,Int_t hittime);
  Float_t fPadOverlap;
  Float_t fTimeOverlap;
  Int_t fTrackThreshold; //minimum weigth track need in order to be included.(=Nhits/weight)
  AliL3MemHandler *fMemHandler; //!

  void CalcClusterWidth(Cluster *cl,Float_t &sigmaY2,Float_t &sigmaZ2);
  
 protected:
  
  AliL3TrackArray *fTracks; //!
  AliL3DigitRowData *fRowData;//!
  Digit *fRow; //!
  Char_t fPath[1024];

  Bool_t fDebug;
  Int_t fNClusters;
  Int_t fMaxClusters;
  Int_t fCurrentPadRow;
  Int_t fMaxPads;
  Int_t fMaxTimebins;
  Int_t fPadSearch;
  Int_t fTimeSearch;
  Int_t fInnerPadSearch;
  Int_t fInnerTimeSearch;
  Int_t fOuterPadSearch;
  Int_t fOuterTimeSearch;
  
  Int_t fSlice;
  Int_t fPatch;
  
  void FillCluster(AliL3ModelTrack *track,Cluster *cluster,Int_t row,Int_t npads);
  void FillZeros(AliL3DigitRowData *digPt,Bool_t reversesign=kFALSE);
  void LocateCluster(AliL3ModelTrack *track,ClusterRegion *region,Int_t &padmin,Int_t &padmax);
  void GetTrackID(Int_t pad,Int_t time,Int_t *trackID);
  
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
  virtual void SetFitRange(Int_t /*p*/,Int_t /*t*/) {return;}
  virtual void SetNmaxOverlaps(Int_t /*i*/) {return;}
  virtual void SetChiSqMax(Float_t /*f*/) {return;}
  
  AliL3TrackArray *GetTracks() {return fTracks;}
    
  ClassDef(AliL3Modeller,1) //Modeller class
    
};

#endif
