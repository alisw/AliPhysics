#ifndef ALIL3MERGER_H
#define ALIL3MERGER_H
#define PI 3.14159265358979312

#include "AliL3RootTypes.h"

class AliL3Track;
class AliL3TrackSegmentData;
class AliL3Vertex;
class AliL3Transform;
class AliL3TrackArray;

class AliL3Merger {
 private:
  Double_t fMaxY;
  Double_t fMaxZ;
  Double_t fMaxKappa;
  Double_t fMaxPsi;
  Double_t fMaxTgl;
  void SetArray(Int_t nin);
  void DeleteArray();

  AliL3TrackArray **fInTrack;//!
  Int_t fNIn;

  AliL3TrackArray *fOutTrack;//!

 protected:
  Int_t fCurrentTracks;
  Int_t fSlice;
  AliL3Vertex *fVertex;//!
  AliL3Transform *fTransformer;//!  
  Bool_t f2Global;
  Bool_t Is2Global(Bool_t is){f2Global=is;return f2Global;}

 public:
  AliL3Merger();
  AliL3Merger(Int_t ntrackarrays);
  virtual ~AliL3Merger();

  Int_t GetNIn(){return fNIn;}
  AliL3TrackArray *GetInTracks(Int_t in){return fInTrack[in];}
  AliL3TrackArray *GetOutTracks(){return fOutTrack;}

  Bool_t Is2Global(){return f2Global;}
  void SetTransformer(AliL3Transform *trans){fTransformer = trans;}
  void SetVertex(AliL3Vertex *vertex){fVertex=vertex;}
  void Reset();
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void FillTracks(Int_t ntracks, AliL3TrackSegmentData* tr); //Fill tracks in fTrackArray[fCurrentTracks] 
  Double_t GetAngle(Double_t a1,Double_t a2);
  void* GetNtuple();
  void* GetNtuple(char *varlist);
  Bool_t WriteNtuple(char *filename,void* nt);
  void FillNtuple(void* nt,Float_t *data);
  void FillNtuple(void* nt,AliL3Track *innertrack,AliL3Track *outertrack);
  void AddAllTracks();//Copy all Tracks to Output Array
  void SortGlobalTracks(AliL3Track **tracks, Int_t ntrack);
  void SortTracks(AliL3Track **tracks, Int_t ntrack);
  void AddTrack(AliL3TrackArray *mergedtrack,AliL3Track *track);
  AliL3Track * MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack);
  AliL3Track * MergeTracks(AliL3TrackArray *mergedtrack,AliL3Track *t0,AliL3Track *t1);
  Bool_t IsTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  Bool_t IsRTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  Double_t TrackDiff(AliL3Track *innertrack,AliL3Track *outertrack);
  void Print();
  void PrintDiff(AliL3Track *innertrack,AliL3Track *outertrack);
  void PrintTrack(AliL3Track *track);
//  Int_t WriteTracks(Char_t *file);
//  Int_t WriteInTracks(Char_t *file);
//  Int_t WriteAllTracks(Char_t *file);
  
  ClassDef(AliL3Merger,1) 
};

#endif
