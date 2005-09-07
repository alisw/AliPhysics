// @(#) $Id$
#ifndef ALIHLTTPCMERGER_H
#define ALIHLTTPCMERGER_H
#define PI 3.14159265358979312

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCTrack;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCVertex;
class AliHLTTPCTrackArray;

class AliHLTTPCMerger {
 private:
  Double_t fMaxY;
  Double_t fMaxZ;
  Double_t fMaxKappa;
  Double_t fMaxPsi;
  Double_t fMaxTgl;
  void SetArray(Int_t nin);
  void DeleteArray();
  Char_t fTrackType;
  
  AliHLTTPCTrackArray **fInTrack;//!
  Int_t fNIn;

  AliHLTTPCTrackArray *fOutTrack;//!

 protected:
  Int_t fCurrentTracks;
  Int_t fSlice;
  AliHLTTPCVertex *fVertex;//!
  Bool_t f2Global;
  Bool_t Is2Global(Bool_t is){f2Global=is;return f2Global;}
  void InitMerger(Int_t ntrackarrays,Char_t *tracktype="AliHLTTPCTrack");
  
 public:
  AliHLTTPCMerger();
  virtual ~AliHLTTPCMerger();

  Int_t GetNIn(){return fNIn;}
  AliHLTTPCTrackArray *GetInTracks(Int_t in){return fInTrack[in];}
  AliHLTTPCTrackArray *GetOutTracks(){return fOutTrack;}

  Bool_t Is2Global(){return f2Global;}
  void SetVertex(AliHLTTPCVertex *vertex){fVertex=vertex;}
  void Reset();
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr); //Fill tracks in fTrackArray[fCurrentTracks] 
  Double_t GetAngle(Double_t a1,Double_t a2);
  void* GetNtuple();
  void* GetNtuple(char *varlist);
  Bool_t WriteNtuple(char *filename,void* nt);
  void FillNtuple(void* nt,Float_t *data);
  void FillNtuple(void* nt,AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  void AddAllTracks();//Copy all Tracks to Output Array
  void SortGlobalTracks(AliHLTTPCTrack **tracks, Int_t ntrack);
  virtual void SortTracks(AliHLTTPCTrack **tracks, Int_t ntrack);
  virtual void AddTrack(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrack *track);
  virtual AliHLTTPCTrack * MultiMerge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrack **tracks, Int_t ntrack);
  AliHLTTPCTrack * MergeTracks(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrack *t0,AliHLTTPCTrack *t1);
  virtual Bool_t IsTrack(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  Bool_t IsRTrack(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  Double_t TrackDiff(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  void Print();
  void PrintDiff(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack);
  void PrintTrack(AliHLTTPCTrack *track);
//  Int_t WriteTracks(Char_t *file);
//  Int_t WriteInTracks(Char_t *file);
//  Int_t WriteAllTracks(Char_t *file);
  
  ClassDef(AliHLTTPCMerger,1) //Merging base class
};

#endif
