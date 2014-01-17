#ifndef ALIFLATESDEVENT_H
#define ALIFLATESDEVENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */

#include "Rtypes.h"
#include "AliFlatESDTrack.h"
#include "AliFlatESDV0.h"

class AliESDEvent;
struct AliFlatESDVertex;
class AliESDVertex;
class AliESDV0;

class AliFlatESDEvent {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatESDEvent();   
  AliFlatESDEvent(AliESDEvent *esd);   
  AliFlatESDEvent(AliESDEvent *esd, Bool_t useESDFriends);   
  ~AliFlatESDEvent();  

  // --------------------------------------------------------------------------------
  // -- Fill / Set methods
  void Reset();

  Int_t Fill( const AliESDEvent *esd, const Bool_t useESDFriends = kTRUE, const Bool_t fillV0s=kTRUE );

  void FillPrimaryVertices( const AliESDVertex *vertexSPD,
			    const AliESDVertex *vertexTracks );

 
  AliFlatESDTrack *GetNextTrackPointer(){ return reinterpret_cast<AliFlatESDTrack*>(fContent + fSize); }

  void StoreLastTrack(){ 
    fNTracks++;
    fSize+= GetNextTrackPointer()->GetSize();
    fV0Pointer = fSize;
  }

  AliFlatESDV0 *GetNextV0Pointer(){ return reinterpret_cast<AliFlatESDV0*>(fContent + fSize); }

  void StoreLastV0(){ 
    fNV0s++;
    fSize+= sizeof(AliFlatESDV0);
  }

  // --------------------------------------------------------------------------------
  // -- Getter methods

  AliFlatESDVertex* GetPrimaryVertexSPD(){
    return (fPrimaryVertexMask & 0x1) ? reinterpret_cast<AliFlatESDVertex*>(fContent) : NULL;
  } 

  AliFlatESDVertex* GetPrimaryVertexTracks() { 
    return (fPrimaryVertexMask & 0x2) ? reinterpret_cast<AliFlatESDVertex*>(fContent + CountBits(fPrimaryVertexMask, 0x1)) : NULL;
  } 

  Int_t GetNumberOfV0s() {return fNV0s;}

  Int_t GetNumberOfTracks() {return fNTracks;}
  
  AliFlatESDTrack *GetTracks() {return reinterpret_cast<AliFlatESDTrack*>(fContent + fTracksPointer);}
  
  // --------------------------------------------------------------------------------
  // -- Size methods
  static ULong64_t EstimateSize(AliESDEvent*esd, Bool_t useESDFriends = kTRUE, Bool_t fillV0s=kTRUE);
         ULong64_t GetSize()    {return fContent - reinterpret_cast<Byte_t*>(this) + fSize;}
  
 private:
  AliFlatESDEvent(const AliFlatESDEvent&);
  AliFlatESDEvent& operator=(const AliFlatESDEvent&);

  void FillPrimaryVertex(const AliESDVertex *v, Byte_t flag);
  Int_t FillNextTrack( const AliESDtrack* esdTrack,  AliESDfriendTrack* friendTrack);
  Int_t FillNextV0( const AliESDV0 *v0);
  UInt_t CountBits(Byte_t field, UInt_t mask);

  // --------------------------------------------------------------------------------
  // -- Fixed size member objects
  //    -> Try to align in memory

  Byte_t   fPrimaryVertexMask;            // Bit mask specfifying which primary vertices are present
  Int_t    fNTracks;                   // Number of tracks in vector
  ULong64_t fTracksPointer;            // position of the first track in fContent
  Int_t fNV0s; // Number of v0's
  ULong64_t fV0Pointer;            // position of the first V0 in fContent


  ULong64_t fSize;                      // Size of this object
  
  // --------------------------------------------------------------------------------
  // -- Variable Size Object
  Byte_t fContent[1];                  // Variale size object, which contains all data
  
};
#endif
