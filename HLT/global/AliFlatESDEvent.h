#ifndef ALIFLATESDEVENT_H
#define ALIFLATESDEVENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */

#include "Rtypes.h"
#include "AliVVevent.h"
#include "AliFlatESDVertex.h"

class AliFlatESDTrack;
class AliFlatESDV0;
class AliFlatESDTrigger;

class AliESDEvent;
class AliESDVertex;
class AliESDtrack;
class AliESDfriendTrack;
class AliESDv0;

class AliFlatESDEvent :public AliVVevent {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatESDEvent();
  ~AliFlatESDEvent();

  // Interface to AliVVEvent

  void Reset(){}

  Double_t  GetMagneticField() const { return fMagneticField; }
  UInt_t    GetPeriodNumber()  const { return fPeriodNumber; }
  Int_t     GetRunNumber()     const { return fRunNumber; }
  UInt_t    GetOrbitNumber()   const { return fOrbitNumber; }
  UShort_t  GetBunchCrossNumber() const { return fBunchCrossNumber; }
  UInt_t    GetTimeStamp()     const { return fTimeStamp; }
  ULong64_t GetTriggerMask()   const { return fTriggerMask; }
  
  TString GetFiredTriggerClasses() const { TString s; return s; } //!!
  UInt_t GetEventSpecie() const { return 0; }  //!!

  Int_t GetNumberOfTracks() const { return fNTracks; }
  Int_t GetNumberOfV0s() const { return fNV0s; }
  Int_t GetNumberOfKinks() const { return 0; }  

  AliVVtrack* GetVVTrack(Int_t /*i*/) const { return NULL; }
  //AliVVtrack* GetTrack(Int_t /*i*/) const { return NULL; }
  AliESDkink* GetKink(Int_t /*i*/) const { return NULL;}

  void ConnectTracks(){} // initialisation, needed by the ESD event, part of VVevent interface

  // --------------------------------------------------------------------------------
  // -- Fill / Set methods
  void Init();

  Int_t FillFromESD( const size_t allocatedMemorySize, const AliESDEvent *esd, const Bool_t useESDFriends = kTRUE, const Bool_t fillV0s=kTRUE );

  void  SetMagneticField( Double_t mf ){ fMagneticField = mf; }
  void  SetPeriodNumber( Int_t n ) { fPeriodNumber = n; }
  void  SetRunNumber( Int_t n ) { fRunNumber = n; }
  void  SetOrbitNumber( UInt_t n ) { fOrbitNumber = n; }
  void  SetBunchCrossNumber( UShort_t n ) { fBunchCrossNumber = n; }
  void  SetTimeStamp( UInt_t timeStamp ){ fTimeStamp = timeStamp; }
  void  SetTriggerMask(ULong64_t n) { fTriggerMask = n; }
  
  AliFlatESDTrigger *FillTriggersStart();
  void FillTriggersEnd( Int_t nTriggerClasses, size_t triggersSize );

  void FillPrimaryVertices( const AliESDVertex *vertexSPD, const AliESDVertex *vertexTracks );

  AliFlatESDTrack *FillTracksStart();
  void FillTracksEnd( Int_t nTracks, size_t tracksSize );

  AliFlatESDV0 *FillV0sStart();
  void FillV0sEnd( Int_t nV0s, size_t v0sSize );  


  // --------------------------------------------------------------------------------
  // -- Getter methods


  const AliFlatESDVertex* GetPrimaryVertexSPD() const ;
  const AliFlatESDVertex* GetPrimaryVertexTracks() const ;

  Int_t GetNumberOfTriggerClasses() const { return fNTriggerClasses; }
   
  const AliFlatESDTrigger *GetTriggerClasses() const { return reinterpret_cast<const AliFlatESDTrigger*>( fContent + fTriggerPointer ); }
  const AliFlatESDTrack   *GetTracks() const { return reinterpret_cast<const AliFlatESDTrack*>( fContent + fTracksPointer ); }
  const AliFlatESDV0      *GetV0s() const { return reinterpret_cast<const AliFlatESDV0*>( fContent + fV0Pointer ); }

  // --------------------------------------------------------------------------------
  // -- Size methods

  ULong64_t  GetSize()  const { return fContent + fContentSize - reinterpret_cast<const Byte_t*>(this); }

  static ULong64_t EstimateSize(AliESDEvent*esd, Bool_t useESDFriends = kTRUE, Bool_t fillV0s=kTRUE );

 private:

  AliFlatESDEvent( const AliFlatESDEvent& );
  AliFlatESDEvent& operator=( const AliFlatESDEvent& );

  Int_t AddTriggerClass(  const char *TriggerClassName, Int_t TriggerIndex, ULong64_t MaxSize );
  void  FillPrimaryVertex(const AliESDVertex *v, Byte_t flag);
  Int_t FillNextTrack( const AliESDtrack* esdTrack,  AliESDfriendTrack* friendTrack);
  Int_t FillNextV0( const AliESDv0 *v0);

  const AliFlatESDVertex* GetFirstPrimaryVertex() const ;
  Byte_t    *GetEndOfContent(){ return fContent + fContentSize; }

  UInt_t CountBits(Byte_t field, UInt_t mask) const;

  // --------------------------------------------------------------------------------

  // -- Fixed size member objects ( try to align in memory )

  size_t     fContentSize;     // Size of fContent
  Double32_t fMagneticField;   // Solenoid Magnetic Field in kG : for compatibility with AliMagF
  UInt_t     fPeriodNumber;    // PeriodNumber
  Int_t      fRunNumber;       // Run Number
  UInt_t     fOrbitNumber;     // Orbit Number
  UInt_t     fTimeStamp;         // Time stamp
  UShort_t   fBunchCrossNumber;  // Bunch Crossing Number
  Byte_t     fPrimaryVertexMask;            // Bit mask specfifying which primary vertices are present
  ULong64_t  fTriggerMask;      // Trigger mask

  UInt_t  fNTriggerClasses;  // N trigger classes
  UInt_t  fNPrimaryVertices; // Number of primary vertices in array
  UInt_t  fNTracks;          // Number of tracks in array
  UInt_t  fNV0s;             // Number of v0's
  
  // Pointers to specific data in fContent
  
  size_t fTriggerPointer;        // position of the first trigger description in fContent
  size_t fPrimaryVertexPointer;  // position of the first primary vertex in fContent
  size_t fTracksPointer;         // position of the first track in fContent
  size_t fV0Pointer;             // position of the first V0 in fContent  

  // --------------------------------------------------------------------------------

  // -- Variable Size Object

  Byte_t fContent[1];                  // Variale size object, which contains all data

  ClassDef(AliFlatESDEvent,0)
};

// Inline implementations 

inline AliFlatESDTrack *AliFlatESDEvent::FillTracksStart()
{
  fNTracks = 0;
  fTracksPointer = fContentSize;
  return reinterpret_cast< AliFlatESDTrack* >( fContent + fContentSize );
}

inline void AliFlatESDEvent::FillTracksEnd( Int_t nTracks, size_t tracksSize )
{
  fNTracks = nTracks;
  fContentSize += tracksSize;
}

inline AliFlatESDV0 *AliFlatESDEvent::FillV0sStart()
{
  fNV0s = 0;
  fV0Pointer = fContentSize;
  return reinterpret_cast< AliFlatESDV0* >( fContent + fContentSize );
}
   
inline void AliFlatESDEvent::FillV0sEnd( Int_t nV0s, size_t v0sSize )
{
  fNV0s = nV0s;
  fContentSize += v0sSize;
}
  
inline AliFlatESDTrigger *AliFlatESDEvent::FillTriggersStart()
{
  fNTriggerClasses = 0;
  fTriggerPointer = fContentSize;
  return reinterpret_cast< AliFlatESDTrigger* >( fContent + fContentSize );
}

inline void AliFlatESDEvent::FillTriggersEnd( Int_t nTriggerClasses, size_t triggersSize )
{
  fNTriggerClasses = nTriggerClasses;
  fContentSize += triggersSize;
}

inline const AliFlatESDVertex* AliFlatESDEvent::GetFirstPrimaryVertex() const 
{
  return reinterpret_cast<const AliFlatESDVertex*>(fContent + fPrimaryVertexPointer);
}
 
inline const AliFlatESDVertex* AliFlatESDEvent::GetPrimaryVertexSPD() const 
{
  return (fPrimaryVertexMask & 0x1) ? GetFirstPrimaryVertex() : NULL;
} 

inline const AliFlatESDVertex* AliFlatESDEvent::GetPrimaryVertexTracks() const 
{ 
  return (fPrimaryVertexMask & 0x2) ? (GetFirstPrimaryVertex() + (fPrimaryVertexMask & 0x1) ) : NULL;
} 

#endif
