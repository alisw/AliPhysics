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
 *  Bool_t useESDFriends = kTRUE;
 *
 *  // -- Book memory for AliFlatESDEvent
 *  Byte_t *mem = new Byte_t[AliFlatESDEvent::EstimateSize(esd, useESDFriends)];
 *  AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent*>(mem);
 *
 *  // -- Fill AliFlatESDEvent
 *  flatEsd->Fill(esd, useESDFriends);  
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

#include "AliFlatESDVertex.h"

#include "AliFlatESDV0.h"
#include "AliFlatESDTrigger.h"

#include "AliESDEvent.h"
#include "AliESDVertex.h"

// _______________________________________________________________________________________________________
AliFlatESDEvent::AliFlatESDEvent() :
  AliVVevent(),
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
  fNPrimaryVertices(0),
  fNTracks(0),
  fNV0s(0),
  fTriggerPointer(0),
  fPrimaryVertexTracksPointer(0),
  fPrimaryVertexSPDPointer(0),
  fTrackTablePointer(0),
  fTracksPointer(0),
  fV0Pointer(0)
{
  // Default constructor
  fContent[0]=0;
}


TString AliFlatESDEvent::GetFiredTriggerClasses() const 
{ 
  // Fired trigger classes
  TString trclasses; 
  const AliFlatESDTrigger *tr = GetTriggerClasses();
  ULong64_t mask = GetTriggerMask() | GetTriggerMaskNext50();
  for(Int_t i = 0; i < GetNumberOfTriggerClasses(); i++) {
    int index = tr->GetTriggerIndex();    
    if( mask & (1ull<<index) ){
      trclasses += " ";
      trclasses += tr->GetTriggerClassName();
      trclasses += " ";
    }
  }
  tr = tr->GetNextTrigger();
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
  fNPrimaryVertices = 0;
  fNTracks = 0;
  fNV0s = 0;
  fTriggerPointer = 0;
  fPrimaryVertexTracksPointer = 0;
  fPrimaryVertexSPDPointer = 0;
  fTrackTablePointer = 0;
  fTracksPointer = 0;
  fV0Pointer = 0;
}

// _______________________________________________________________________________________________________
  ULong64_t AliFlatESDEvent::EstimateSize(AliESDEvent *esd, Bool_t fillV0s) 
{
  // Estimate upper limit of the object size
  // -> Added objects have to be added here as well
  
  ULong64_t size = sizeof(AliFlatESDEvent);
  size += 2 * sizeof( AliFlatESDVertex );
  size += esd->GetNumberOfTracks() * ( AliFlatESDTrack::EstimateSize() + sizeof(Long64_t) );
  if( fillV0s ) size += esd->GetNumberOfV0s()*sizeof(AliFlatESDV0);
  return size;
}

Int_t AliFlatESDEvent::SetPrimaryVertexTracks( const AliESDVertex *vtx, size_t allocatedVtxMemory )
{
  // fill primary vertex tracks
  if( !vtx ) return 0;
  if( allocatedVtxMemory < sizeof(AliFlatESDVertex) ) return -1;
  fPrimaryVertexMask |= 0x1;
  fPrimaryVertexTracksPointer = fContentSize;
  AliFlatESDVertex *flatVtx = reinterpret_cast<AliFlatESDVertex*> (fContent + fContentSize);
  flatVtx->Set( *vtx );
  fContentSize += flatVtx->GetSize();
  return 0;
}

Int_t AliFlatESDEvent::SetPrimaryVertexSPD( const AliESDVertex *vtx, size_t allocatedVtxMemory  )
{
  // fill primary vertex SPD
  if( !vtx ) return 0;
  if( allocatedVtxMemory < sizeof(AliFlatESDVertex) ) return -1;
  fPrimaryVertexMask |= 0x2;
  fPrimaryVertexSPDPointer = fContentSize;
  AliFlatESDVertex *flatVtx = reinterpret_cast<AliFlatESDVertex*> (fContent + fContentSize);
  flatVtx->Set( *vtx );
  fContentSize += flatVtx->GetSize();
  return 0;
}


// _______________________________________________________________________________________________________
Int_t AliFlatESDEvent::SetFromESD( const size_t allocatedMemorySize, const AliESDEvent *esd, const Bool_t fillV0s)
{
  // Fill flat ESD event from normal ALiESDEvent
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

  // fill primary vertices

  err = SetPrimaryVertexTracks( esd->GetPrimaryVertexTracks(), freeSpace );
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
       flatTrack->Set( esdTrack );
       trackSize += flatTrack->GetSize();
       freeSpace -= flatTrack->GetSize();
       nTracks++;
       flatTrack = flatTrack->GetNextTrack();
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
    }
    SetV0sEnd( nV0s, v0size );
  }
  
  return 0;
}

