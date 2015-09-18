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
#include "AliFlatESDVZEROFriend.h"
#include "AliESDVZEROfriend.h"

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

  AliVVZEROfriend* GetVVZEROfriend() { return GetFlatVZEROFriendNonConst(); }
  Int_t GetESDVZEROfriend( AliESDVZEROfriend &v ) const;

 
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

  Int_t  SetVZEROFriend( const AliESDVZEROfriend *v, size_t freeMem );

  Int_t SetTracksStart( AliFlatESDFriendTrack* &t, Long64_t* &table, Int_t nTracks, size_t freeMem );
  void  SetTracksEnd( Int_t nTracks, Int_t nTrackEntries, size_t tracksSize );

  // other methods

  const AliFlatESDVZEROFriend *GetFlatVZEROFriend() const { 
    return (fVZEROFriendPointer>0) ?reinterpret_cast<const AliFlatESDVZEROFriend*>( fContent + fVZEROFriendPointer ) :NULL; 
  }

  AliFlatESDVZEROFriend *GetFlatVZEROFriendNonConst()  { 
    return (fVZEROFriendPointer>0) ?reinterpret_cast<AliFlatESDVZEROFriend*>( fContent + fVZEROFriendPointer ) :NULL; 
  }

  const AliFlatESDFriendTrack  *GetFlatTrack( Int_t i ) const { 
    const Long64_t *table = reinterpret_cast<const Long64_t*> (fContent + fTrackTablePointer);
    if( i<0 || i>=fNTracks || table[i]<0 ) return NULL;
    return reinterpret_cast<const AliFlatESDFriendTrack*>( fContent + table[i] );
  }
  
  AliFlatESDFriendTrack  *GetFlatTrackNonConst( Int_t i ){ 
    const Long64_t *table = reinterpret_cast<const Long64_t*> (fContent + fTrackTablePointer);
    if( i<0 || i>=fNTracks || table[i]<0 ) return NULL;
    return reinterpret_cast<AliFlatESDFriendTrack*>( fContent + table[i] );
  }

  const AliFlatESDFriendTrack  *GetFirstTrackEntry() const { 
    return reinterpret_cast<const AliFlatESDFriendTrack*>( fContent + fTracksPointer );
  }
 
  AliFlatESDFriendTrack  *GetFirstTrackEntryNonConst(){ 
    return reinterpret_cast<AliFlatESDFriendTrack*>( fContent + fTracksPointer );
  }
 

  // -- Size methods

  ULong64_t  GetSize()  const { return fContent - reinterpret_cast<const Byte_t*>(this) + fContentSize; }

  static ULong64_t EstimateSize(AliESDfriend* esdFriend );

private: 

  AliFlatESDFriend(const AliFlatESDFriend&);
  AliFlatESDFriend& operator=(const AliFlatESDFriend& );  

  size_t fContentSize;     // Size of fContent
  UInt_t fBitFlags; // bit flags
  Int_t fNTracks;                   // Number of tracks in vector
  Int_t fNTrackEntries;             // Number of non-empty track friends in vector
  Int_t fNclustersTPC[72]; //cluster occupancy per sector per sector
  Int_t fNclustersTPCused[72]; //number of clusters used in tracking per sector
 
  // Pointers to specific data in fContent
  
  Long_t fVZEROFriendPointer;     // position of flat VZERO friend in fContent
  size_t fTrackTablePointer;     // position of the first track in fContent
  size_t fTracksPointer;         // position of the first track in fContent

  // -- Variable Size Object

  Byte_t fContent[1];                  // Variale size object, which contains all data

  ClassDef(AliFlatESDFriend,0)
};


inline AliFlatESDFriend::AliFlatESDFriend() 
:
  fContentSize(0),
  fBitFlags(0),
  fNTracks(0),
  fNTrackEntries(0),
  fVZEROFriendPointer(-1),
  fTrackTablePointer(0),
  fTracksPointer(0)
{
  // Default constructor
  Reset();
}

inline void AliFlatESDFriend::Reset() 
{
  fBitFlags = 0;
  fNTracks = 0;
  fNTrackEntries = 0; 
  fVZEROFriendPointer = -1;
  fTrackTablePointer = 0;
  fTracksPointer = 0;
  for( int i=0; i<72; i++ ){
    fNclustersTPC[i]=0;
    fNclustersTPCused[i]=0;
  }
  // We set size of the fContent array such, that it reaches the end of the AliFlatESDFriend structure. 
  // To be sure that actual GetSize() is always >= size of the structure. 
  // First, it makes the memory alignment better. Second, just for a case..
  fContentSize = sizeof(AliFlatESDFriend) - (fContent - reinterpret_cast<const Byte_t*>(this));
  for( UInt_t i=0; i<fContentSize; i++ ) fContent[i]=0;
}
 
#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatESDFriend::AliFlatESDFriend(AliVConstructorReinitialisationFlag f)
:
  AliVfriendEvent(f)
{
  //special constructor, used to restore the vtable pointer

  // Reinitialise VZERO information  
  {    
    AliFlatESDVZEROFriend * vzero =  GetFlatVZEROFriendNonConst();
    if( vzero ) vzero->Reinitialize();    
  }

  // track info
  AliFlatESDFriendTrack  *tr = GetFirstTrackEntryNonConst();
  for( int i=0; i<fNTrackEntries; i++ ){
    tr->Reinitialize();
    tr = tr->GetNextTrackNonConst();
  }
}
#pragma GCC diagnostic warning "-Weffc++" 


inline Int_t AliFlatESDFriend::SetTracksStart( AliFlatESDFriendTrack* &t, Long64_t* &table, Int_t nTracks, size_t freeMem)
{
  fNTracks = 0;
  fNTrackEntries = 0;
  size_t memoryAlignment = reinterpret_cast<size_t>(fContent+fContentSize) % 64;
  if (memoryAlignment) memoryAlignment = 64 - memoryAlignment;
  if( nTracks*sizeof(Long64_t)  + memoryAlignment> freeMem ) return -1;
  fContentSize += memoryAlignment;
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


inline Int_t AliFlatESDFriend::GetESDVZEROfriend( AliESDVZEROfriend &v ) const
{
  const AliFlatESDVZEROFriend* flatVZERO = GetFlatVZEROFriend();
  if( !flatVZERO ) return -1;
  flatVZERO->GetESDVZEROfriend( v );
  return 0;
}

#endif
