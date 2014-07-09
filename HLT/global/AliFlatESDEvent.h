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
#include "AliVVevent.h"
#include "AliFlatESDVertex.h"
#include "AliFlatESDMisc.h"

class AliESDEvent;
class AliESDVertex;
class AliESDV0;
class TString;
class AliVVv0;

class AliFlatESDEvent: public AliVVevent {
//class AliFlatESDEvent {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatESDEvent();   


  AliFlatESDEvent(AliESDEvent *esd);   
  AliFlatESDEvent(AliESDEvent *esd, Bool_t useESDFriends);   
  virtual ~AliFlatESDEvent();  

  // --------------------------------------------------------------------------------
  // -- Fill / Set methods
  void Reset();

  Int_t Fill( const AliESDEvent *esd, const Bool_t useESDFriends = kTRUE, const Bool_t fillV0s=kTRUE );

  void FillPrimaryVertices( const AliESDVertex *vertexSPD,
			    const AliESDVertex *vertexTracks );

 
  AliFlatESDTrack *GetNextTrackPointer(){ 
    return reinterpret_cast<AliFlatESDTrack*> (fContent + fSize);
  }

  void StoreLastTrack(){ 
    fNTracks++;
    fSize+= GetNextTrackPointer()->GetSize();
    fV0Pointer = fSize;
  }

  AliFlatESDV0 *GetNextV0Pointer(){
		return reinterpret_cast<AliFlatESDV0*> (fContent + fSize);
}

  void StoreLastV0(){ 
    fNV0s++;
    fSize+= sizeof(AliFlatESDV0);
  }

  // --------------------------------------------------------------------------------
  // -- Getter methods

   const AliFlatESDVertex* GetPrimaryVertexSPD() const {
    if (fPrimaryVertexMask & 0x1){
		 return reinterpret_cast< AliFlatESDVertex*> (const_cast<Byte_t*>(fContent));
	}
	else return NULL;
  } 

  const  AliFlatESDVertex* GetPrimaryVertexTracks() const { 
    if (fPrimaryVertexMask & 0x2){
		 return reinterpret_cast< AliFlatESDVertex*> (const_cast<Byte_t*>(fContent)) + CountBits(fPrimaryVertexMask, 0x1);
	}
	else return NULL;
   } 

  Int_t GetNumberOfV0s() const {return fNV0s;}

  Int_t GetNumberOfTracks() const {return fNTracks;}
  
  AliFlatESDV0* GetV0s() {
  	return reinterpret_cast<AliFlatESDV0*> (fContent + fV0Pointer);
  }

  AliFlatESDTrack *GetTracks() {
		return reinterpret_cast<AliFlatESDTrack*> (fContent + fTracksPointer);
  }

  const AliVVvertex* GetPrimaryVertex() const {return NULL;}
  const AliVVvertex* GetPrimaryVertexTPC() const {return NULL;}
  AliFlatESDTrack* GetTrack(Int_t /*i*/) const {return NULL;}
  AliVVkink* GetKink(Int_t /*i*/) const {return NULL;}
  AliFlatESDV0* GetV0(Int_t /*i*/) const {return NULL;}
  Int_t GetNumberOfKinks() const {return 0;}
  Int_t GetEventNumberInFile() const {return -1;}
  const AliMultiplicity* GetMultiplicity() const {return NULL;} //by default SPDmult
  Int_t GetRunNumber() const {return -1;}
  TString GetFiredTriggerClasses() const {TString string; return string;}
  TObject* FindListObject(const char* /*name*/) const {return NULL;}
  ULong64_t GetTriggerMask() const {return 0;}
  Double_t GetMagneticField() const {return 0;}
  UInt_t GetTimeStamp() const { return 0;}
  UInt_t GetEventSpecie() const { return 0;}

  
  // --------------------------------------------------------------------------------
  // -- Size methods
  static ULong64_t EstimateSize(AliESDEvent*esd, Bool_t useESDFriends = kTRUE, Bool_t fillV0s=kTRUE);
         ULong64_t GetSize()    {return fContent - reinterpret_cast<Byte_t*>(this) + fSize;}


//
// Initializing function
//
// to be called after event is received via reinterpret_cast from memory

void Reinitialize();

 private:
  AliFlatESDEvent(const AliFlatESDEvent&);
  AliFlatESDEvent& operator=(const AliFlatESDEvent&);
  
// special constructor, to be called by placement new,
// when accessing information after reinterpret_cast
// so that vtable is generated, but values are not overwritten
	AliFlatESDEvent(AliFlatESDSpecialConstructorFlag){}

  void FillPrimaryVertex(const AliESDVertex *v, Byte_t flag);
  Int_t FillNextTrack( const AliESDtrack* esdTrack,  AliESDfriendTrack* friendTrack);
  Int_t FillNextV0( const AliESDV0 *v0);
  UInt_t CountBits(Byte_t field, UInt_t mask) const;

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
