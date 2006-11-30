// @(#) $Id$

#ifndef ALIL3TRACKARRAY_H
#define ALIL3TRACKARRAY_H

#include "AliHLTRootTypes.h"

class AliHLTConfMapTrack;
class AliHLTTrack;
class AliHLTTrackSegmentData;

class AliHLTTrackArray {

 private:

  Char_t fTrackType; //track type
  Int_t fSize; //size of arra
  Bool_t *fIsPresent;//!
  Int_t fNAbsent; //ntracks absent

  AliHLTTrack **fTrack;//!
  Int_t fNTracks; //ntracks in

  UInt_t WriteConfMapTracks(AliHLTTrackSegmentData* tr); 
  void DeleteArray();

 public:
  AliHLTTrackArray();
  AliHLTTrackArray(Int_t ntrack);
  AliHLTTrackArray(char* tracktype,Int_t ntrack);
  AliHLTTrackArray(char* tracktype);
  virtual ~AliHLTTrackArray();
  Int_t GetTrackType(){return fTrackType;}
  Int_t GetSize() const {return fSize;}
  Bool_t SetSize(Int_t newsize=2000);

  Int_t GetNPresent() const {return (fNTracks- fNAbsent);}
  Int_t GetNTracks() const {return fNTracks;}
  AliHLTTrack *NextTrack();
  AliHLTTrack *GetCheckedTrack(Int_t t){if(fIsPresent[t]) return fTrack[t]; return 0;}
  AliHLTTrack *GetTrack(Int_t t){return fTrack[t];}

  void Remove(Int_t track); 
  void RemoveLast() {fNTracks--;}
  void Compress();
  void Reset();
  void QSort();
  void QSort( AliHLTTrack **a, Int_t first, Int_t last);
  Int_t TrackCompare(AliHLTTrack *a, AliHLTTrack *b) const;

  void FillTracks(Int_t ntracks, AliHLTTrackSegmentData* tr,Int_t slice); //Fill tracks and transform
  void FillTracks(Int_t ntracks, AliHLTTrackSegmentData* tr); //Fill tracks
  UInt_t WriteTracks(AliHLTTrackSegmentData* tr); //Write tracks
  UInt_t WriteTracks(UInt_t & ntracks,AliHLTTrackSegmentData* tr); //Write tracks
  UInt_t GetOutSize();
  UInt_t GetOutCount(){return (UInt_t) GetNPresent();}
  void AddTracks(AliHLTTrackArray *newtrack,Bool_t remove_old=kTRUE,Int_t slice=-1);//add all Tracks to this 
  void AddLast(AliHLTTrack *track);

  ClassDef(AliHLTTrackArray,1) //Track array class
};

typedef AliHLTTrackArray AliL3TrackArray; // for backward compatibility

#endif
