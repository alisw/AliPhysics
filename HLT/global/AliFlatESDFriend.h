#ifndef ALIFLATESDFRIEND_H
#define ALIFLATESDFRIEND_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */

#include "Rtypes.h"
#include "AliFlatESDMisc.h"
#include "AliVVeventFriend.h"
#include "AliFlatESDFriendTrack.h"

class AliVVVZEROfriend;
class AliVVTZEROfriend;

//_____________________________________________________________________________
class AliFlatESDFriend : public AliVVfriend {
public:
  AliFlatESDFriend();
  ~AliFlatESDFriend() {}

  // Implementation of virtual methods of AliVVfriend

  Int_t GetNumberOfTracks() const { return fNTracks; }
  Int_t GetEntriesInTracks() const { return fNTrackEntries; }
  const AliVVfriendTrack* GetTrack(Int_t i) const {return GetFlatTrack(i); }
  
  AliVVVZEROfriend *GetVZEROfriend(){ return NULL; }
  AliVVTZEROfriend *GetTZEROfriend(){ return NULL; }

  void Ls() const;

  void Reset();
  
  Bool_t TestSkipBit() const { return (fBitFlags!=0); }

  Int_t GetNclustersTPC(UInt_t sector) const { return (sector<72)?fNclustersTPC[sector]:0; }
  Int_t GetNclustersTPCused(UInt_t sector) const { return (sector<72)?fNclustersTPCused[sector]:0; }
  
  //virtual void AddTrack(const AliVVfriendTrack *t) {}
  //virtual void AddTrackAt(const AliVVfriendTrack* /*t*/, Int_t /*i*/) {}
  //virtual void SetVZEROfriend(AliESDVZEROfriend* /*obj*/) {}
  //virtual void SetTZEROfriend(AliESDTZEROfriend * obj) {}
  //void SetSkipBit(Bool_t skip){}

  // Own methods   
  
  void SetSkipBit(Bool_t skip){ fBitFlags = skip; }
  void SetNclustersTPC(UInt_t sector, Int_t occupancy ) { if (sector<72) fNclustersTPC[sector]=occupancy; }
  void SetNclustersTPCused(UInt_t sector, Int_t occupancy ) {if (sector<72) fNclustersTPCused[sector]=occupancy; }

  const AliFlatESDFriendTrack  *GetFlatTrack( Int_t i ) const { 
    const Long64_t *table = reinterpret_cast<const Long64_t*> (fContent + fTrackTablePointer);
    if( i<0 || i>fNTracks || table[i]<0 ) return NULL;
    return reinterpret_cast<const AliFlatESDFriendTrack*>( fContent + table[i] );
  }

  // -- Size methods

  ULong64_t  GetSize()  const { return fContent - reinterpret_cast<const Byte_t*>(this) + fContentSize; }

  void Reinitialize()
  {
    new (this) AliFlatESDFriend(AliFlatESDReinitialize);
  }

private: 

  AliFlatESDFriend(const AliFlatESDFriend&);
  AliFlatESDFriend& operator=(const AliFlatESDFriend& );  

  // special constructor, to be called by placement new,
  // when accessing information after reinterpret_cast
  // so that vtable is generated, but values are not overwritten
  AliFlatESDFriend(AliFlatESDSpecialConstructorFlag);
 
  AliFlatESDFriendTrack  *GetFlatTrackNonConst( Int_t i ){ 
    const Long64_t *table = reinterpret_cast<const Long64_t*> (fContent + fTrackTablePointer);
    if( i<0 || i>fNTracks || table[i]<0 ) return NULL;
    return reinterpret_cast<AliFlatESDFriendTrack*>( fContent + table[i] );
  }

  size_t fContentSize;     // Size of fContent
  UInt_t fBitFlags; // bit flags
  Int_t fNTracks;                   // Number of tracks in vector
  Int_t fNTrackEntries;             // Number of non-empty track friends in vector
  Int_t fNclustersTPC[72]; //cluster occupancy per sector per sector
  Int_t fNclustersTPCused[72]; //number of clusters used in tracking per sector
 
  // Pointers to specific data in fContent
  
  size_t fTrackTablePointer;         // position of the first track in fContent
  size_t fTracksPointer;         // position of the first track in fContent

  // -- Variable Size Object

  Byte_t fContent[1];                  // Variale size object, which contains all data
};

#endif
