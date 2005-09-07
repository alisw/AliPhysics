// @(#) $Id$

#ifndef ALIHLTTPCTRACKARRAY_H
#define ALIHLTTPCTRACKARRAY_H

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCConfMapTrack;
class AliHLTTPCTrack;
class AliHLTTPCTrackSegmentData;

class AliHLTTPCTrackArray{
 private:
  void DeleteArray();

  Char_t fTrackType;
  Int_t fSize;
  Bool_t *fIsPresent;//!
  Int_t fNAbsent;

  AliHLTTPCTrack **fTrack;//!
  Int_t fNTracks;

  UInt_t WriteConfMapTracks(AliHLTTPCTrackSegmentData* tr); 

 public:
  AliHLTTPCTrackArray();
  AliHLTTPCTrackArray(Int_t ntrack);
  AliHLTTPCTrackArray(char* tracktype,Int_t ntrack);
  AliHLTTPCTrackArray(char* tracktype);
  virtual ~AliHLTTPCTrackArray();
  Int_t GetTrackType(){return fTrackType;}
  Int_t GetSize(){return fSize;}
  Bool_t SetSize(Int_t newsize=2000);

  Int_t GetNPresent(){return (fNTracks- fNAbsent);}

  Int_t GetNTracks(){return fNTracks;}
  AliHLTTPCTrack *NextTrack();
  AliHLTTPCTrack *GetCheckedTrack(Int_t t){if(fIsPresent[t]) return fTrack[t]; return 0;}
  AliHLTTPCTrack *GetTrack(Int_t t){return fTrack[t];}

  void Remove(Int_t track); 
  void RemoveLast() {fNTracks--;}
  void Compress();
  void Reset();
  void QSort();
  void QSort( AliHLTTPCTrack **a, Int_t first, Int_t last);
  Int_t TrackCompare(AliHLTTPCTrack *a, AliHLTTPCTrack *b);


  void FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr,Int_t slice); //Fill tracks and transform
  void FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr); //Fill tracks
  UInt_t WriteTracks(AliHLTTPCTrackSegmentData* tr); //Write tracks
  UInt_t WriteTracks(UInt_t & ntracks,AliHLTTPCTrackSegmentData* tr); //Write tracks
  UInt_t GetOutSize();
  UInt_t GetOutCount(){return (UInt_t) GetNPresent();}
  void AddTracks(AliHLTTPCTrackArray *newtrack,Bool_t remove_old=kTRUE,Int_t slice=-1);//add all Tracks to this 
  void AddLast(AliHLTTPCTrack *track);

  ClassDef(AliHLTTPCTrackArray,1) //Track array class
};

#endif
