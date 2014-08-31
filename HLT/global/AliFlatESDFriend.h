#ifndef ALIFLATESDFRIEND_H
#define ALIFLATESDFRIEND_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */

#include "Rtypes.h"
#include "AliVMisc.h"
#include "AliVfriendEvent.h"
#include "AliFlatESDFriendTrack.h"


class AliESDfriend;
//class AliESDVZEROfriend;
//class AliESDTZEROfriend;


//_____________________________________________________________________________
class AliFlatESDFriend : public AliVfriendEvent {
public:
  AliFlatESDFriend();
  ~AliFlatESDFriend() {}
 
   // constructor and method for reinitialisation of virtual table
   AliFlatESDFriend( AliVConstructorReinitialisationFlag );
   void Reinitialize(){ new (this) AliFlatESDFriend(AliVReinitialize); }

  // Implementation of virtual methods of AliVfriend

  Int_t GetNumberOfTracks() const { return fNTracks; }
  const AliVfriendTrack* GetTrack(Int_t i) const {return GetFlatTrack(i); }
  Int_t GetEntriesInTracks() const { return fNTrackEntries; }
 
  //AliESDVZEROfriend *GetVZEROfriend(){ return NULL; }
  //AliESDTZEROfriend *GetTZEROfriend(){ return NULL; }

  void Ls() const;
  void Reset();

  // bit manipulation for filtering
  void SetSkipBit(Bool_t skip){ fBitFlags = skip; }
  Bool_t TestSkipBit() const { return (fBitFlags!=0); }

  //TPC cluster occupancy
  Int_t GetNclustersTPC(UInt_t sector) const { return (sector<72)?fNclustersTPC[sector]:0; }
  Int_t GetNclustersTPCused(UInt_t sector) const { return (sector<72)?fNclustersTPCused[sector]:0; }
  
  // -- Own methods  -- 

  // Set methods

  Int_t SetFromESDfriend( size_t allocatedMemorySize, const AliESDfriend *esdFriend );
    
  void SetNclustersTPC(UInt_t sector, Int_t occupancy ) { if (sector<72) fNclustersTPC[sector]=occupancy; }
  void SetNclustersTPCused(UInt_t sector, Int_t occupancy ) {if (sector<72) fNclustersTPCused[sector]=occupancy; }

  Int_t SetTracksStart( AliFlatESDFriendTrack* &t, Long64_t* &table, Int_t nTracks, size_t freeMem );
  void  SetTracksEnd( Int_t nTracks, Int_t nTrackEntries, size_t tracksSize );

  // other methods

  const AliFlatESDFriendTrack  *GetFlatTrack( Int_t i ) const { 
    const Long64_t *table = reinterpret_cast<const Long64_t*> (fContent + fTrackTablePointer);
    if( i<0 || i>fNTracks || table[i]<0 ) return NULL;
    return reinterpret_cast<const AliFlatESDFriendTrack*>( fContent + table[i] );
  }

  // -- Size methods

  ULong64_t  GetSize()  const { return fContent - reinterpret_cast<const Byte_t*>(this) + fContentSize; }


private: 

  AliFlatESDFriend(const AliFlatESDFriend&);
  AliFlatESDFriend& operator=(const AliFlatESDFriend& );  

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
  
  size_t fTrackTablePointer;     // position of the first track in fContent
  size_t fTracksPointer;         // position of the first track in fContent

  // -- Variable Size Object

  Byte_t fContent[1];                  // Variale size object, which contains all data
};


inline Int_t AliFlatESDFriend::SetTracksStart( AliFlatESDFriendTrack* &t, Long64_t* &table, Int_t nTracks, size_t freeMem)
{
  fNTracks = 0;
  fNTrackEntries = 0;
  if( nTracks*sizeof(Long64_t)  > freeMem ) return -1;
  fTrackTablePointer = fContentSize;
  fContentSize += nTracks*sizeof(Long64_t);
  fTracksPointer = fContentSize;
  table = reinterpret_cast< Long64_t* >( fContent + fTrackTablePointer );
  t = reinterpret_cast< AliFlatESDFriendTrack* >( fContent + fTracksPointer );
  return 0;
}

inline void AliFlatESDFriend::SetTracksEnd( Int_t nTracks, Int_t nTrackEntries, size_t tracksSize )
{
  if( nTracks<0 ) return;
  Long64_t *table = reinterpret_cast< Long64_t*> (fContent + fTrackTablePointer);
  for( int i=0; i<nTracks; i++ ) table[i]+=fTracksPointer;
  fNTracks = nTracks;
  fNTrackEntries = nTrackEntries;
  fContentSize += tracksSize;
}

#endif
