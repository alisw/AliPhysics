// @(#) $Id$
#ifndef ALIL3MERGER_H
#define ALIL3MERGER_H
//#define PI 3.14159265358979312

#include "AliL3RootTypes.h"

class AliL3Track;
class AliL3TrackSegmentData;
class AliL3Vertex;
class AliL3TrackArray;

class AliL3Merger {
 private:
  Double_t fMaxY;    //maxy
  Double_t fMaxZ;    //maxz
  Double_t fMaxKappa;//maxkappa
  Double_t fMaxPsi;  //maxpsi
  Double_t fMaxTgl;  //maxtgl
  Char_t fTrackType; //track type to merge
  
  AliL3TrackArray **fInTrack;//!
  Int_t fNIn; //ntracks

  AliL3TrackArray *fOutTrack;//!

  void SetArray(Int_t nin);
  void DeleteArray();

 protected:
  Int_t fCurrentTracks; //current number
  Int_t fSlice;         //slice
  AliL3Vertex *fVertex; //!
  Bool_t f2Global; //global
  Bool_t Is2Global(Bool_t is){f2Global=is;return f2Global;}
  void InitMerger(Int_t ntrackarrays,Char_t *tracktype="AliL3Track");
  
 public:
  AliL3Merger();
  virtual ~AliL3Merger();

  Int_t GetNIn() const {return fNIn;}
  AliL3TrackArray *GetInTracks(Int_t in){return fInTrack[in];}
  AliL3TrackArray *GetOutTracks(){return fOutTrack;}

  Bool_t Is2Global() const {return f2Global;}
  void SetVertex(AliL3Vertex *vertex){fVertex=vertex;}
  void Reset();
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void FillTracks(Int_t ntracks, AliL3TrackSegmentData* tr); //Fill tracks in fTrackArray[fCurrentTracks] 
  Double_t GetAngle(Double_t a1,Double_t a2);
  void* GetNtuple() const;
  void* GetNtuple(char *varlist) const;
  Bool_t WriteNtuple(char *filename,void* nt) const;
  void FillNtuple(void* nt,Float_t *data) const ;
  void FillNtuple(void* nt,AliL3Track *innertrack,AliL3Track *outertrack);
  void AddAllTracks();//Copy all Tracks to Output Array
  void SortGlobalTracks(AliL3Track **tracks, Int_t ntrack);
  virtual void SortTracks(AliL3Track **tracks, Int_t ntrack) const;
  virtual void AddTrack(AliL3TrackArray *mergedtrack,AliL3Track *track);
  virtual AliL3Track * MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack);
  AliL3Track * MergeTracks(AliL3TrackArray *mergedtrack,AliL3Track *t0,AliL3Track *t1);
  virtual Bool_t IsTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  Bool_t IsRTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  Double_t TrackDiff(AliL3Track *innertrack,AliL3Track *outertrack);
  void Print();
  void PrintDiff(AliL3Track *innertrack,AliL3Track *outertrack);
  void PrintTrack(AliL3Track *track);
  
  ClassDef(AliL3Merger,1) //Merging base class
};

#endif
