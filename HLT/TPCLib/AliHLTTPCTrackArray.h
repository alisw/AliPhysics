// -*- Mode: C++ -*-
// $Id$
// Original: AliHLTTrackArray.h,v 1.7 2004/06/11 16:06:33 loizides 
#ifndef ALIHLTTPCTRACKARRAY_H
#define ALIHLTTPCTRACKARRAY_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCTrackArray.h
/// @author Uli Frankenfeld, maintained by Matthias Richter
/// @date   
/// @brief  Array of AliHLTTPCTracks
///

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCConfMapTrack;
class AliHLTTPCTrack;
struct AliHLTTPCTrackSegmentData;
struct AliHLTTPCTrackSegmentDataV1;
struct AliHLTExternalTrackParam;
/**
 * @class AliHLTTPCTrackArray
 * Array of AliHLTTrack objects.
 * The class implements a dynamic array and handler methods.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCTrackArray {
 public:
  /** default constructor */
  AliHLTTPCTrackArray();
  /**
   * constructor 
   * @param ntrack     initial size
   */
  AliHLTTPCTrackArray(Int_t ntrack);
  /**
   * constructor 
   * @param tracktype  string describing type, one of
   *   - AliHLTTPCTrack        -> 't'
   *   - AliHLTTPCConfMapTrack -> 'c'
   *   - AliHLTTPCHoughTrack   -> 'h'
   *   - AliHLTTPCModelTrack   -> 'm'
   * @param ntrack     initial size
   */
  AliHLTTPCTrackArray(const char* tracktype,Int_t ntrack);
  /**
   * constructor 
   * @param tracktype  string describing type, one of
   *   - AliHLTTPCTrack        -> 't'
   *   - AliHLTTPCConfMapTrack -> 'c'
   *   - AliHLTTPCHoughTrack   -> 'h'
   *   - AliHLTTPCModelTrack   -> 'm'
   */
  AliHLTTPCTrackArray(const char* tracktype);
  /** destructor */
  virtual ~AliHLTTPCTrackArray();

  /**
   * Get type of track.
   * @return one of
   *   - 't' -> AliHLTTPCTrack        
   *   - 'c' -> AliHLTTPCConfMapTrack 
   *   - 'h' -> AliHLTTPCHoughTrack   
   *   - 'm' -> AliHLTTPCModelTrack   
   */
  Int_t GetTrackType(){return fTrackType;}

  /**
   * Get size of the array.
   * @return size of the array
   */
  Int_t GetSize() const {return fSize;}

  /**
   * Set size.
   * If the current size is smaller, the array is grown to the new size.
   * @return kTRUE if the array was grown, kFALSE otherwise
   */
  Bool_t SetSize(Int_t newsize=0);

  Int_t GetNPresent() const {return (fNTracks- fNAbsent);}
  Int_t GetNTracks() const {return fNTracks;}

  /**
   * Return pointer to next free track object.
   * The array is grown if necessary.
   */
  AliHLTTPCTrack *NextTrack();
  AliHLTTPCTrack *GetCheckedTrack(Int_t t) const {if(fIsPresent[t]) return fTrack[t]; return 0;}
  AliHLTTPCTrack *GetTrack(Int_t t) const {return fTrack[t];}

  void Remove(Int_t track); 
  void RemoveLast() {fNTracks--;}
  void Compress();
  void Reset();
  void QSort();
  void QSort( AliHLTTPCTrack **a, Int_t first, Int_t last);
  Int_t TrackCompare(AliHLTTPCTrack *a, AliHLTTPCTrack *b) const;

  /**
   * Fill track array from track segment array.
   * Old method excluding buffer protection kept for backward compatibility.
   * Forwarded to FillTracksChecked().
   * @param ntracks      size of the input array
   * @param tr           array of AliHLTTrackSegmentData
   * @param slice        slice no to transform the tracks to
   * @param bTransform   transform to global coordinates if 1
   */
  int FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr, Int_t slice=-1, Int_t bTransform=1);

  /**
   * Fill track array from track segment array.
   * Reads the track from an array of AliHLTTrackSegmentData. The coordinates
   * are transformed to global coordinates if the slice parameter is specified.
   * In that case to internal slice variable is set to zero.
   * 
   * The sizeInByte parameter allows an additional buffer protection if
   * non-zero. The size of the AliHLTTPCTrackSegmentData is not fixed due to
   * variable array at the end of the structure. The pointer to the next
   * entry must be set according to the variable array.
   * @param ntracks      size of the input array
   * @param tr           array of AliHLTTrackSegmentData
   * @param sizeInByte   additional size protection
   * @param slice        slice no to transform the tracks to
   * @param bTransform   transform to global coordinates if 1
   */
  int FillTracksChecked(AliHLTTPCTrackSegmentData* tr, Int_t ntracks, unsigned int sizeInByte, Int_t slice=-1, Int_t bTransform=1);
  int FillTracksChecked(AliHLTExternalTrackParam* tr, Int_t ntracks, unsigned int sizeInByte, Int_t slice=-1, Int_t bTransform=1);

  /**
   * Fill array from version1 structure.
   * The version 1 of ALiHLTTPCTrackSegmentData was valid until July 2008
   * revision 27415.
   *
   * Similar behavior like FillTracksChecked.
   */
  int FillTracksVersion1(AliHLTTPCTrackSegmentDataV1* tr, Int_t ntracks, unsigned int sizeInByte, Int_t slice=-1, Int_t bTransform=1);

  UInt_t WriteTracks(AliHLTTPCTrackSegmentData* tr); //Write tracks
  UInt_t WriteTracks(UInt_t & ntracks,AliHLTTPCTrackSegmentData* tr); //Write tracks
  UInt_t GetOutSize();
  UInt_t GetOutCount(){return (UInt_t) GetNPresent();}
  void AddTracks(AliHLTTPCTrackArray *newtrack,Bool_t remove_old=kTRUE,Int_t slice=-1);//add all Tracks to this 
  void AddLast(AliHLTTPCTrack *track);

  AliHLTTPCTrack* operator[](int index);

 private:
  /** copy constructor prohibited */
  AliHLTTPCTrackArray(const AliHLTTPCTrackArray&);
  /** assignment operator prohibited */
  AliHLTTPCTrackArray& operator=(const AliHLTTPCTrackArray&);

  Char_t fTrackType; //track type
  Int_t fSize; //size of array
  Bool_t *fIsPresent;//!
  Int_t fNAbsent; //ntracks absent

  AliHLTTPCTrack **fTrack;//!
  Int_t fNTracks; //ntracks in

  UInt_t WriteConfMapTracks(AliHLTTPCTrackSegmentData* tr); 
  void DeleteArray();

  ClassDef(AliHLTTPCTrackArray,2) //Track array class
};

#endif
