// @(#) $Id$
#ifndef ALIL3MERGER_H
#define ALIL3MERGER_H
//#define PI 3.14159265358979312

#include "AliHLTRootTypes.h"

class AliHLTTrack;
class AliHLTTrackSegmentData;
class AliHLTVertex;
class AliHLTTrackArray;

class AliHLTMerger {
 private:
  Double_t fMaxY;    //maxy
  Double_t fMaxZ;    //maxz
  Double_t fMaxKappa;//maxkappa
  Double_t fMaxPsi;  //maxpsi
  Double_t fMaxTgl;  //maxtgl
  Char_t fTrackType; //track type to merge
  
  AliHLTTrackArray **fInTrack;//!
  Int_t fNIn; //ntracks

  AliHLTTrackArray *fOutTrack;//!

  void SetArray(Int_t nin);
  void DeleteArray();

 protected:
  Int_t fCurrentTracks; //current number
  Int_t fSlice;         //slice
  AliHLTVertex *fVertex; //!
  Bool_t f2Global; //global
  Bool_t Is2Global(Bool_t is){f2Global=is;return f2Global;}
  void InitMerger(Int_t ntrackarrays,Char_t *tracktype="AliHLTTrack");
  
 public:
  AliHLTMerger();
  virtual ~AliHLTMerger();

  Int_t GetNIn() const {return fNIn;}
  AliHLTTrackArray *GetInTracks(Int_t in){return fInTrack[in];}
  AliHLTTrackArray *GetOutTracks(){return fOutTrack;}

  Bool_t Is2Global() const {return f2Global;}
  void SetVertex(AliHLTVertex *vertex){fVertex=vertex;}
  void Reset();
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void FillTracks(Int_t ntracks, AliHLTTrackSegmentData* tr); //Fill tracks in fTrackArray[fCurrentTracks] 
  Double_t GetAngle(Double_t a1,Double_t a2);
  void* GetNtuple() const;
  void* GetNtuple(char *varlist) const;
  Bool_t WriteNtuple(char *filename,void* nt) const;
  void FillNtuple(void* nt,Float_t *data) const ;
  void FillNtuple(void* nt,AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  void AddAllTracks();//Copy all Tracks to Output Array
  void SortGlobalTracks(AliHLTTrack **tracks, Int_t ntrack);
  virtual void SortTracks(AliHLTTrack **tracks, Int_t ntrack) const;
  virtual void AddTrack(AliHLTTrackArray *mergedtrack,AliHLTTrack *track);
  virtual AliHLTTrack * MultiMerge(AliHLTTrackArray *mergedtrack,AliHLTTrack **tracks, Int_t ntrack);
  AliHLTTrack * MergeTracks(AliHLTTrackArray *mergedtrack,AliHLTTrack *t0,AliHLTTrack *t1);
  virtual Bool_t IsTrack(AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  Bool_t IsRTrack(AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  Double_t TrackDiff(AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  void Print();
  void PrintDiff(AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  void PrintTrack(AliHLTTrack *track);
  
  ClassDef(AliHLTMerger,1) //Merging base class
};

typedef AliHLTMerger AliL3Merger; // for backward compatibility

#endif
