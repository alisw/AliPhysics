// @(#) $Id$
// Original: AliHLTMerger.h,v 1.8 2004/06/11 16:06:33 loizides 
#ifndef ALIHLTTPCMERGER_H
#define ALIHLTTPCMERGER_H
//#define PI 3.14159265358979312

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCTrack;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCVertex;
class AliHLTTPCTrackArray;

class AliHLTTPCMerger {
 private:
  Double_t fMaxY;    //maxy
  Double_t fMaxZ;    //maxz
  Double_t fMaxKappa;//maxkappa
  Double_t fMaxPsi;  //maxpsi
  Double_t fMaxTgl;  //maxtgl
  Char_t fTrackType; //track type to merge
  
  AliHLTTPCTrackArray **fInTrack;//!
  Int_t fNIn; //ntracks

  AliHLTTPCTrackArray *fOutTrack;//!

  void SetArray(Int_t nin);
  void DeleteArray();

 protected:
  Int_t fCurrentTracks; //current number
  Int_t fSlice;         //slice
  AliHLTTPCVertex *fVertex; //!
  Bool_t f2Global; //global
  Bool_t Is2Global(Bool_t is){f2Global=is;return f2Global;}
  void InitMerger(Int_t ntrackarrays,Char_t *tracktype="AliHLTTPCTrack");
  
 public:
  AliHLTTPCMerger();
  virtual ~AliHLTTPCMerger();

  Int_t GetNIn() const {return fNIn;}
  AliHLTTPCTrackArray *GetInTracks(Int_t in){return fInTrack[in];}
  AliHLTTPCTrackArray *GetOutTracks(){return fOutTrack;}

  Bool_t Is2Global() const {return f2Global;}
  void SetVertex(AliHLTTPCVertex *vertex){fVertex=vertex;}
  void Reset();
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr); //Fill tracks in fTrackArray[fCurrentTracks] 
  Double_t GetAngle(Double_t a1,Double_t a2);
  void* GetNtuple() const;
  void* GetNtuple(char *varlist) const;
  Bool_t WriteNtuple(char *filename,void* nt) const;
  void FillNtuple(void* nt,Float_t *data) const ;
  void FillNtuple(void* nt,AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  void AddAllTracks();//Copy all Tracks to Output Array
  void SortGlobalTracks(AliHLTTPCTrack **tracks, Int_t ntrack);
  virtual void SortTracks(AliHLTTPCTrack **tracks, Int_t ntrack) const;
  virtual void AddTrack(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrack *track);
  virtual AliHLTTPCTrack * MultiMerge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrack **tracks, Int_t ntrack);
  AliHLTTPCTrack * MergeTracks(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrack *t0,AliHLTTPCTrack *t1);
  virtual Bool_t IsTrack(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  Bool_t IsRTrack(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  Double_t TrackDiff(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  void Print();
  void PrintDiff(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  void PrintTrack(AliHLTTPCTrack *track);
  
  ClassDef(AliHLTTPCMerger,1) //Merging base class
};

#endif
