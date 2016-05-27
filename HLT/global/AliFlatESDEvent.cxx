/* $Id$ */

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * >> Flat structure representing an ESD <<
 *
 * To be used in the online and offline calibration schema.
 *
 * Class provides interface methods for 
 *   - Filling from AliESDEvent, but also from HLT (to be added)
 *   - Getter methods
 *
 * In the online case, the structure can be directly written into a shared 
 * memory, in the offline case, the size has to be estimated first.
 *
 * 
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 * ************************************************************************
 * Offline usage:
 *
 *  ...
 *  AliESDEvent* esd = ....;
 *  Bool_t fillV0s = kTRUE;
 *
 *  // -- Book memory for AliFlatESDEvent
 *  Int_t memSize = AliFlatESDEvent::EstimateSize(esd, fillV0s);
 *  Byte_t *mem = new Byte_t[memSize];
 *  AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent*>(mem);
 *
 *  // -- Fill AliFlatESDEvent
 *  flatEsd->SetFromESD(memSize, esd, fillV0s);
 *  ...
 *
 **************************************************************************/

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDv0.h"

#include "AliFlatESDEvent.h"
#include "AliFlatESDTrack.h"
#include "AliFlatTPCCluster.h"
#include "AliFlatExternalTrackParam.h"
#include "Riostream.h"

#include "AliFlatESDVZERO.h"
#include "AliFlatMultiplicity.h"
#include "AliFlatESDVertex.h"

#include "AliFlatESDV0.h"
#include "AliFlatESDTrigger.h"

#include "AliESDEvent.h"
#include "AliESDVertex.h"

#include <bitset>

// _______________________________________________________________________________________________________
AliFlatESDEvent::AliFlatESDEvent() 
  :
  fContentSize(0),
  fMagneticField(0),
  fPeriodNumber(0),
  fRunNumber(0),
  fOrbitNumber(0),
  fTimeStamp(0),
  fEventSpecie(0),
  fBunchCrossNumber(0),
  fPrimaryVertexMask(0),
  fTriggerMask(0),
  fTriggerMaskNext50(0),
  fNTriggerClasses(0),
  fNTracks(0),
  fNV0s(0),
  fNtrackletsSPD(0),
  fTriggerPointer(0),
  fVZEROPointer(-1),
  fMultiplicityPointer(-1),
  fPrimaryVertexTracksPointer(0),
  fPrimaryVertexTPCPointer(0),
  fPrimaryVertexSPDPointer(0),
  fTrackTablePointer(0),
  fTracksPointer(0),
  fV0Pointer(0),
  fFriendEvent(NULL)
{
  // Default constructor
  for( int i=0; i<6; i++ ) fNclustersITS[i] = 0;
  fContent[0]=0;
}

#pragma GCC diagnostic ignored "-Weffc++" 
AliFlatESDEvent::AliFlatESDEvent( AliVConstructorReinitialisationFlag /*f*/ ) 
  :
 fFriendEvent(NULL)
{
  // Constructor for reinitialisation of vtable

  // Reinitialise trigger information  
  {
    AliFlatESDTrigger * trigger =  reinterpret_cast< AliFlatESDTrigger*>( fContent + fTriggerPointer ); 
    for( UInt_t i=0; i<fNTriggerClasses; i++ ){
      trigger->Reinitialize();
      trigger = trigger->GetNextTriggerNonConst();
    }
  }

  // Reinitialise VZERO information  
  {
    AliFlatESDVZERO * vzero =  GetFlatVZERONonConst();
    if( vzero ) vzero->Reinitialize();    
  }

  // Reinitialise Multiplicity information
  {
    AliFlatMultiplicity * mult =  GetFlatMultiplicityNonConst();
    if( mult ) mult->Reinitialize();
  }

  // Reinitialise primary vertices

  if( fPrimaryVertexMask & 0x1 ){
    AliFlatESDVertex *vtxSPD = reinterpret_cast<AliFlatESDVertex*>(fContent + fPrimaryVertexSPDPointer);
    vtxSPD->Reinitialize();
  }
  if( fPrimaryVertexMask & 0x2 ){
    AliFlatESDVertex *vtxTracks = reinterpret_cast<AliFlatESDVertex*>(fContent + fPrimaryVertexTracksPointer);
    vtxTracks->Reinitialize();
  }
  if( fPrimaryVertexMask & 0x4 ){
    AliFlatESDVertex *vtxTPC = reinterpret_cast<AliFlatESDVertex*>(fContent + fPrimaryVertexTPCPointer);
    vtxTPC->Reinitialize();
  }

  // Reinitialise tracks 
  {
    AliFlatESDTrack *track = reinterpret_cast<AliFlatESDTrack*>( fContent + fTracksPointer );
    for( UInt_t i=0; i<fNTracks; i++ ){
      track->Reinitialize();
      track = track->GetNextTrackNonConst();
    }
  }

  // Reinitialise V0s
  {
    AliFlatESDV0 *v0 = reinterpret_cast<AliFlatESDV0*>( fContent + fV0Pointer ); 
    for( UInt_t i=0; i < fNV0s; i++){
      v0->Reinitialize();
      v0 = v0->GetNextV0NonConst();
    }
  }
}
#pragma GCC diagnostic warning "-Weffc++" 


Bool_t AliFlatESDEvent::IsTriggerClassFired(const char* name) const
{
  return GetFiredTriggerClasses().Contains(name);
}

TString AliFlatESDEvent::GetFiredTriggerClasses() const 
{ 
  // Fired trigger classes
  TString trclasses; 
  const AliFlatESDTrigger *tr = GetTriggerClasses();
  std::bitset<AliESDRun::kNTriggerClasses>  mask = GetTriggerMask() | (GetTriggerMaskNext50()<<50);
  for(Int_t i = 0; i < GetNumberOfTriggerClasses(); i++) {
    int index = tr->GetTriggerIndex();    
    if( mask.test(index) ){
      trclasses += " ";
      trclasses += tr->GetTriggerClassName();
    }
    tr = tr->GetNextTrigger();
  }
  return trclasses;
}


void AliFlatESDEvent::Reset()
{
  // Init

  fContentSize = 0;
  fMagneticField = 0;
  fPeriodNumber = 0;
  fRunNumber = 0;
  fOrbitNumber = 0;
  fTimeStamp = 0;
  fBunchCrossNumber = 0;
  fPrimaryVertexMask = 0;
  fTriggerMask = 0;
  fTriggerMaskNext50 = 0;
  fNTriggerClasses = 0;
  fNTracks = 0;
  fNV0s = 0;
  fNtrackletsSPD = 0;
  fTriggerPointer = 0;
  fVZEROPointer = -1;
  fMultiplicityPointer = -1;
  fPrimaryVertexTracksPointer = 0;
  fPrimaryVertexTPCPointer = 0;
  fPrimaryVertexSPDPointer = 0;
  fTrackTablePointer = 0;
  fTracksPointer = 0;
  fV0Pointer = 0;
  for( int i=0; i<6; i++ ) fNclustersITS[i] = 0;
}

// _______________________________________________________________________________________________________
  ULong64_t AliFlatESDEvent::EstimateSize(AliESDEvent *esd, Bool_t fillV0s) 
{
  // Estimate upper limit of the object size
  // -> Added objects have to be added here as well
  
  ULong64_t size = sizeof(AliFlatESDEvent);
  size += 2 * sizeof( AliFlatESDVertex );
	// one Long64_t per track for tracks table
  size += esd->GetNumberOfTracks() * ( AliFlatESDTrack::EstimateSize() + sizeof(Long64_t) );
  size += AliESDRun::kNTriggerClasses * sizeof(AliFlatESDTrigger) ;
  if( esd->GetVZEROData() ) size += sizeof(AliFlatESDVZERO) ;
  if( esd->GetMultiplicity() ) size += sizeof(AliFlatMultiplicity) ;
  if( fillV0s ) size += esd->GetNumberOfV0s()*sizeof(AliFlatESDV0);
  return size;
}

Int_t AliFlatESDEvent::SetVZEROData( const AliESDVZERO *vzero, size_t allocatedVZEROMemory )
{
  // fill VZERO info
  fVZEROPointer = -1;
  if( !vzero ) return 0;
  if( allocatedVZEROMemory < sizeof(AliFlatESDVZERO) ) return -1;
  fVZEROPointer = fContentSize;
  AliFlatESDVZERO *flatVZERO = reinterpret_cast<AliFlatESDVZERO*> (fContent + fVZEROPointer );
  flatVZERO->SetFromESDVZERO( *vzero );
  fContentSize += flatVZERO->GetSize();
  return 0;
}

Int_t AliFlatESDEvent::SetMultiplicity( const AliMultiplicity *mult, size_t allocatedMemory )
{
  // fill Multiplicity info
  fMultiplicityPointer = -1;
  if( !mult ) return 0;
  if( allocatedMemory < sizeof(AliFlatMultiplicity) ) return -1;
  fMultiplicityPointer = fContentSize;
  AliFlatMultiplicity *flatMult = reinterpret_cast<AliFlatMultiplicity*> (fContent + fMultiplicityPointer );
  flatMult->SetFromMultiplicity( *mult );
  fContentSize += flatMult->GetSize();
  return 0;
}

Int_t AliFlatESDEvent::SetMultiplicity( const AliFlatMultiplicity *mult, size_t allocatedMemory )
{
  // fill Multiplicity info
  fMultiplicityPointer = -1;
  if( !mult ) return 0;
  if( allocatedMemory < sizeof(AliFlatMultiplicity) ) return -1;
  fMultiplicityPointer = fContentSize;
  AliFlatMultiplicity *flatMult = reinterpret_cast<AliFlatMultiplicity*> (fContent + fMultiplicityPointer );
  *flatMult = *mult;
  fContentSize += flatMult->GetSize();
  return 0;
}

Int_t AliFlatESDEvent::SetPrimaryVertexTracks( const AliESDVertex *vtx, size_t allocatedVtxMemory )
{
  // fill primary vertex tracks
  if( !vtx ) return 0;
  if(!vtx->GetStatus()) return 0;
  if( allocatedVtxMemory < sizeof(AliFlatESDVertex) ) return -1;
  fPrimaryVertexMask |= 0x1;
  fPrimaryVertexTracksPointer = fContentSize;
  AliFlatESDVertex *flatVtx = reinterpret_cast<AliFlatESDVertex*> (fContent + fContentSize);
  flatVtx->SetFromESDVertex( *vtx );
  fContentSize += flatVtx->GetSize();
  return 0;
}

Int_t AliFlatESDEvent::SetPrimaryVertexSPD( const AliESDVertex *vtx, size_t allocatedVtxMemory  )
{
  // fill primary vertex SPD
  if( !vtx ) return 0;
  if(!vtx->GetStatus()) return 0;
  if( allocatedVtxMemory < sizeof(AliFlatESDVertex) ) return -1;
  fPrimaryVertexMask |= 0x2;
  fPrimaryVertexSPDPointer = fContentSize;
  AliFlatESDVertex *flatVtx = reinterpret_cast<AliFlatESDVertex*> (fContent + fContentSize);
  flatVtx->SetFromESDVertex( *vtx );
  fContentSize += flatVtx->GetSize();
  return 0;
}

Int_t AliFlatESDEvent::SetPrimaryVertexTPC( const AliESDVertex *vtx, size_t allocatedVtxMemory )
{
  // fill primary vertex tracks
  if( !vtx ) return 0;
  if(!vtx->GetStatus()) return 0;
  if( allocatedVtxMemory < sizeof(AliFlatESDVertex) ) return -1;
  fPrimaryVertexMask |= 0x4;
  fPrimaryVertexTPCPointer = fContentSize;
  AliFlatESDVertex *flatVtx = reinterpret_cast<AliFlatESDVertex*> (fContent + fContentSize);
  flatVtx->SetFromESDVertex( *vtx );
  fContentSize += flatVtx->GetSize();
  return 0;
}


// _______________________________________________________________________________________________________
Int_t AliFlatESDEvent::SetFromESD( const size_t allocatedMemorySize, const AliESDEvent *esd, const Bool_t fillV0s)
{
  // Fill flat ESD event from  ALiESDEvent
  // - Fill tracks + v0s
  // -> Added objects have to be added here as well
 
  if( allocatedMemorySize < sizeof(AliFlatESDEvent) ) return -1;

  Reset();
  
  if( !esd ) return 0;
  
  Int_t err = 0;
  size_t freeSpace = allocatedMemorySize - GetSize();

  // fill run info
  {
    SetMagneticField( esd->GetMagneticField() );
    SetPeriodNumber( esd->GetPeriodNumber() );
    SetRunNumber( esd->GetRunNumber() );
    SetOrbitNumber( esd->GetOrbitNumber() );
    SetBunchCrossNumber( esd->GetBunchCrossNumber() );
    SetTimeStamp( esd->GetTimeStamp() );
    SetEventSpecie( esd->GetEventSpecie() );
    SetTriggerMask( esd->GetTriggerMask() );
    SetTriggerMaskNext50( esd->GetTriggerMaskNext50() );
  }
 
  // Fill trigger information  
  {
    size_t triggerSize = 0;
    int nTriggers = 0;
    AliFlatESDTrigger *trigger = SetTriggersStart();
    const AliESDRun*  esdRun = esd->GetESDRun();
    if( esdRun ){ 
      for( int index=0; index<AliESDRun::kNTriggerClasses; index++){
        const char* name = esdRun->GetTriggerClass(index);
	if( name && name[0]!='\0' ){
          err = trigger->SetTriggerClass( name, index, freeSpace );
	  if( err!=0 ) return err;
	  nTriggers++;
	  freeSpace -= trigger->GetSize();
	  triggerSize += trigger->GetSize();
	  trigger = trigger->GetNextTriggerNonConst();
        }
      }
    }
    SetTriggersEnd( nTriggers, triggerSize );    
  }

  // fill VZERO info 

  
  err = SetVZEROData( esd->GetVZEROData(), freeSpace );
  if( err!=0 ) return err;
  freeSpace = allocatedMemorySize - GetSize();  
  
  // fill Multiplicity info

  err = SetMultiplicity( esd->GetMultiplicity(), freeSpace );
  if( err!=0 ) return err;
  freeSpace = allocatedMemorySize - GetSize();

  // fill primary vertices

  err = SetPrimaryVertexTracks( esd->GetPrimaryVertexTracks(), freeSpace );
  if( err!=0 ) return err;
  freeSpace = allocatedMemorySize - GetSize();

  err = SetPrimaryVertexTPC( esd->GetPrimaryVertexTPC(), freeSpace );
  if( err!=0 ) return err;
  freeSpace = allocatedMemorySize - GetSize();

  err = SetPrimaryVertexSPD( esd->GetPrimaryVertexSPD(), freeSpace );
  if( err!=0 ) return err;
  freeSpace = allocatedMemorySize - GetSize();

  // fill tracks 
  {
   size_t trackSize = 0;
   int nTracks = 0;
   Long64_t *table = NULL;
   AliFlatESDTrack *flatTrack = NULL;
   err = SetTracksStart( flatTrack, table, esd->GetNumberOfTracks(), freeSpace );
   if( err!=0 ) return err;
   freeSpace = allocatedMemorySize - GetSize();

   for (Int_t idxTrack = 0; idxTrack < esd->GetNumberOfTracks(); ++idxTrack) {
     AliESDtrack *esdTrack = esd->GetTrack(idxTrack);
     table[idxTrack] = -1;
     if (esdTrack) {
       table[idxTrack] = trackSize;
       if( freeSpace<flatTrack->EstimateSize() ) return -1;
       new (flatTrack) AliFlatESDTrack;       
       flatTrack->SetFromESDTrack( esdTrack );
       trackSize += flatTrack->GetSize();
       freeSpace -= flatTrack->GetSize();
       nTracks++;
       flatTrack = flatTrack->GetNextTrackNonConst();
     }
   }
   SetTracksEnd( nTracks, trackSize );
  }

  // fill V0s

  if( fillV0s ){
    size_t v0size = 0;
    int nV0s = 0; 
    AliFlatESDV0 *flatV0 = SetV0sStart();
    for( int i=0; i < esd->GetNumberOfV0s(); i++){
      AliESDv0 *esdV0 = esd->GetV0( i );
      if( !esdV0 || freeSpace < flatV0->GetSize() ) return -1;
      new (flatV0) AliFlatESDV0;
      flatV0->SetNegTrackID( esdV0->GetNindex());
      flatV0->SetPosTrackID( esdV0->GetPindex());
      nV0s++;
      v0size += flatV0->GetSize();
      freeSpace -= flatV0->GetSize(); 
      flatV0 = flatV0->GetNextV0NonConst();
    }
    SetV0sEnd( nV0s, v0size );
  }
  
  return 0;
}

AliVParticle* AliFlatESDEvent::GetTrack(Int_t i) const
{
  return const_cast<AliFlatESDTrack*>(GetFlatTrack(i));
}
  
AliVEvent::EDataLayoutType AliFlatESDEvent::GetDataLayoutType() const 
{
  return AliVEvent::kFlat;
}

void  AliFlatESDEvent::GetESDEvent( AliESDEvent *esd ) const 
{
  // Fill flat ESD event to ALiESDEvent

  if( !esd ) return;  
  
  esd->Reset(); 
  esd->CreateStdContent();

  // fill run info
  {  
    esd->SetMagneticField( GetMagneticField() );
    esd->SetRunNumber( GetRunNumber() );
    esd->SetPeriodNumber( GetPeriodNumber() );
    esd->SetOrbitNumber( GetOrbitNumber() );
    esd->SetBunchCrossNumber( GetBunchCrossNumber() );
    esd->SetTimeStamp( GetTimeStamp() );
    esd->SetEventSpecie( GetEventSpecie() );
  }

  // Fill trigger information  
  {
    for( int i=0; i<GetNumberOfTriggerClasses(); i++ ){
      const AliFlatESDTrigger &trg = GetTriggerClasses()[i];    
      esd->SetTriggerClass( trg.GetTriggerClassName(), trg.GetTriggerIndex() );
    }
    
    esd->SetTriggerMask( GetTriggerMask() );
    esd->SetTriggerMaskNext50( GetTriggerMaskNext50() );
  }

  // fill VZERO info 
  
  {
    AliESDVZERO v;
    if( GetVZEROData( v )>=0 ) esd->SetVZEROData( &v );
  }
  
  // fill Multiplicity info

  {
    AliMultiplicity v;
    if( GetMultiplicity( v )>=0 ) esd->SetMultiplicity( &v );
  }

  // fill primary vertices
  {
    AliESDVertex v;
    if( GetPrimaryVertexSPD( v ) >= 0 ) esd->SetPrimaryVertexSPD( &v );
    if( GetPrimaryVertexTPC( v ) >= 0 ) esd->SetPrimaryVertexTPC( &v );
    if( GetPrimaryVertexTracks( v ) >= 0 ) esd->SetPrimaryVertexTracks( &v );
  }

   // fill tracks 
  {
    const AliFlatESDTrack *flatTrack = GetTracks();
    for( int i=0; i<GetNumberOfTracks(); i++ ){
      AliESDtrack esdtrack;
      flatTrack->GetESDTrack( &esdtrack );
      esd->AddTrack( &esdtrack );
      flatTrack = flatTrack->GetNextTrack();
    }
  }
 
  // fill V0s
  {
    const AliFlatESDV0 *flatV0 = GetV0s();
    for( int i=0; i<GetNumberOfV0s(); i++ ){
      Int_t negID = flatV0->GetNegTrackID();
      Int_t posID = flatV0->GetPosTrackID();
      const AliESDtrack *negTr = esd->GetTrack( negID );
      const AliESDtrack *posTr = esd->GetTrack( posID );
      if( negTr && posTr ){
	AliESDv0 esdV0( *negTr, negID, *posTr, posID );
	esd->AddV0( &esdV0 );
      }
      flatV0 = flatV0->GetNextV0();  
    }
  }  
  
}

//_____________________________________________________________________________
//example static conversion functions, use with care!
//
AliESDEvent* AliFlatESDEvent::MakeESDevent(AliFlatESDEvent* flatEvent)
{
  //produce an AliESDEvent from a AliFlatESDEvent
  //caller then owns the object
  AliESDEvent* esdEvent = new AliESDEvent();
  flatEvent->GetESDEvent(esdEvent);
  return esdEvent;
}

AliFlatESDEvent* AliFlatESDEvent::MakeFlatEvent(AliESDEvent* esdEvent,Bool_t fillV0s)
{
  //produce an AliFlatESDEvent from en AliESDEvent
  //caller then owns the object
  //object has to be destroyed using DestroyFlatEvent
  Int_t flatEventSize = AliFlatESDEvent::EstimateSize(esdEvent,fillV0s);
  Byte_t *flatEventBuffer = new Byte_t[flatEventSize];
  AliFlatESDEvent *flatEvent = reinterpret_cast<AliFlatESDEvent*>(flatEventBuffer);
  Int_t err = flatEvent->SetFromESD(flatEventSize, esdEvent, fillV0s);
  if (err!=0) { delete [] flatEventBuffer; return NULL; }
  return flatEvent;
}

void AliFlatESDEvent::DestroyFlatEvent(AliFlatESDEvent* flatEvent)
{
  //destroy the object created with MakeFlatEvent
  delete [] reinterpret_cast<Byte_t*>(flatEvent);
}

