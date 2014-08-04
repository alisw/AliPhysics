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
  fContentSize(0),
  fMagneticField(0),
  fPeriodNumber(0),
  fRunNumber(0),
  fOrbitNumber(0),
  fTimeStamp(0),
  fBunchCrossNumber(0),
  fPrimaryVertexMask(0),
  fTriggerMask(0),
  fNTriggerClasses(0),
  fNPrimaryVertices(0),
  fNTracks(0),
  fNV0s(0),
  fTriggerPointer(0),
  fPrimaryVertexPointer(0),
  fTracksPointer(0),
  fV0Pointer(0)
{
  // Default constructor
  fContent[0]=0;
}

/*
// _______________________________________________________________________________________________________
AliFlatESDEvent::AliFlatESDEvent(AliESDEvent *esd) :
  fSize(0),
  fMagneticField(0),
  fPeriodNumber(0),
  fRunNumber(0),
  fOrbitNumber(0),
  fTimeStamp(0),
  fBunchCrossNumber(0),
  fPrimaryVertexMask(0),
  fTriggerMask(0),
  fNTriggerClasses(0),
  fNPrimaryVertices(0),
  fNTracks(0),
  fNV0s(0),
  fTriggerPointer(0),
  fPrimaryVertexPointer(0),
  fTracksPointer(0),
  fV0Pointer(0)
{ 
  fContent[0]=0;
  Fill(esd);
}

// _______________________________________________________________________________________________________
AliFlatESDEvent::AliFlatESDEvent(AliESDEvent *esd, Bool_t useESDFriends) :
  // Constructor
  fPrimaryVertexMask(0),
  fNTracks(0),
  fTracksPointer(0),
  fNV0s(0),
  fV0Pointer(0),
  fSize(0),
  fMagneticField(0),
  fPeriodNumber(0),
  fRunNumber(0),
  fOrbitNumber(0),
  fBunchCrossNumber(0),
  fTimeStamp(0),
  fContent() 
{ 
  Fill(esd, useESDFriends);
}
*/

// _______________________________________________________________________________________________________
AliFlatESDEvent::~AliFlatESDEvent() 
{
  // Destructor
}

void AliFlatESDEvent::Init()
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
  fNTriggerClasses = 0;
  fNPrimaryVertices = 0;
  fNTracks = 0;
  fNV0s = 0;
  fTriggerPointer = 0;
  fPrimaryVertexPointer = 0;
  fTracksPointer = 0;
  fV0Pointer = 0;
}

// _______________________________________________________________________________________________________
  ULong64_t AliFlatESDEvent::EstimateSize(AliESDEvent *esd, Bool_t useESDFriends, Bool_t fillV0s) 
{
  // Estimate upper limit of the object size
  // -> Added objects have to be added here as well
  
  ULong64_t size = sizeof(AliFlatESDEvent);
  size += 2 * sizeof( AliFlatESDVertex );
  size += esd->GetNumberOfTracks() * AliFlatESDTrack::EstimateSize(useESDFriends);
  if( fillV0s ) size += esd->GetNumberOfV0s()*sizeof(AliFlatESDV0);
  return size;
}

void AliFlatESDEvent::FillPrimaryVertices( const AliESDVertex *vertexSPD,
					    const AliESDVertex *vertexTracks )
{
  // fill primary vertices
  
  fPrimaryVertexMask = 0;
  fContentSize = 0;

  Byte_t flag = 0x1;
  FillPrimaryVertex(vertexSPD, flag);

  flag = 0x2;
  FillPrimaryVertex(vertexTracks, flag);

  fTracksPointer = fContentSize;
  fV0Pointer = fContentSize;
}

void AliFlatESDEvent::FillPrimaryVertex(const AliESDVertex *v, Byte_t flag) 
{
  
  // Fill primary vertex parameters

  if (!v) return;

  AliFlatESDVertex *vtx = reinterpret_cast<AliFlatESDVertex*> (fContent + fContentSize);
  vtx->Set( *v );    
  fPrimaryVertexMask |= flag;
  fContentSize += sizeof(AliFlatESDVertex);
}


Int_t AliFlatESDEvent::FillNextTrack( const AliESDtrack* esdTrack, AliESDfriendTrack* friendTrack)
{
  // fill next track

  AliFlatESDTrack *flatTrack = reinterpret_cast<AliFlatESDTrack*>(fContent+fContentSize);
  //new (flatTrack) AliFlatESDTrack;
  flatTrack->Fill(esdTrack, friendTrack);
  fContentSize += flatTrack->GetSize();
  ++fNTracks;
  return 0;
}

Int_t AliFlatESDEvent::AddTriggerClass(  const char *TriggerClassName, Int_t TriggerIndex, ULong64_t MaxSize )
{
  // add trigger class 
  AliFlatESDTrigger *flatTrigger = reinterpret_cast<AliFlatESDTrigger*>(fContent+fContentSize);
  Int_t err = flatTrigger->SetTriggerClass( TriggerClassName, TriggerIndex, MaxSize );
  if( err==0 ){
    fContentSize+= flatTrigger->GetSize();
    fNTriggerClasses++;
  }
  return err;
}
 

// _______________________________________________________________________________________________________
Int_t AliFlatESDEvent::FillFromESD( const size_t allocatedMemorySize, const AliESDEvent *esd, const Bool_t useESDFriends, const Bool_t fillV0s)
{
  // Fill flat ESD event from normal ALiESDEvent
  // - Fill tracks + friends (if requested)
  // -> Added objects have to be added here as well
 
#ifdef xxx
  if( allocatedMemorySize < sizeof(AliFlatESDEvent) ) return -1;

  Init();
  
  if( !esd ) return 0;
  
  Int_t err = 0;
  size_t freeSpace = allocatedMemorySize - GetSize();

  // fill run info
  {
    SetTriggerMask( esd->GetTriggerMask() );
  }


  // Fill trigger information  
  {
    size_t triggerSize = 0;
    int nTriggers = 0;
    AliFlatESDTrigger *trigger = FillTriggersStart();
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
	  trigger = trigger->GetNextTrigger();
        }
      }
    }
    FillTriggersEnd( nTriggers, triggerSize );    
  }

  // fill primary vertices

  
  err = FillPrimaryVertices( freeSpace, esd->GetPrimaryVertexSPD(), esd->GetPrimaryVertexTracks() );
  
  if( err!=0 ) return err;
  
  freeSpace = allocatedMemorySize - GetSize();

 
 ULong64_t GetTriggerMask() const {return fHeader?fHeader->GetTriggerMask():0;}
  TString   GetFiredTriggerClasses() const {return (fESDRun&&fHeader)?fESDRun->GetFiredTriggerClasses(fHeader->GetTriggerMask()):"";}
  Bool_t    IsTriggerClassFired(const char *name) const {return (fESDRun&&fHeader)?fESDRun->IsTriggerClassFired(fHeader->GetTriggerMask(),name):kFALSE;}

  const AliHLTCTPData* pCTPData=CTPData();
  if (pCTPData) {
    AliHLTUInt64_t mask=pCTPData->ActiveTriggers(trigData);
    for (int index=0; index<gkNCTPTriggerClasses; index++) {
      if ((mask&((AliHLTUInt64_t)0x1<<index)) == 0) continue;
      pESD->SetTriggerClass(pCTPData->Name(index), index);
    }
    pESD->SetTriggerMask(mask);
  }

  
  //allocatedMemorySize - GetSize();

  FillTriggersEnd(nTriggers, triggerSize );
 
  // -- Get ESD friends
  // -------------------------------------------------------
  Bool_t connectESDFriends = useESDFriends;
  AliESDfriend* esdFriend  = NULL;

   if (connectESDFriends) {
    esdFriend = dynamic_cast<AliESDfriend*>(esd->FindListObject("AliESDfriend"));  
    if (!esdFriend) {
      connectESDFriends = kFALSE;
      Printf("WARNING: friends not available, cluster information will not be included");
    }
    else 
      Printf("INFO: friends are available, cluster information will be included");
  }

  // -- Track loop, fill AliFlatESDTrack sub structure
  // -------------------------------------------------------
  for (Int_t idxTrack = 0; idxTrack < esd->GetNumberOfTracks(); ++idxTrack) {
    AliESDtrack       *esdTrack    = esd->GetTrack(idxTrack);
    AliESDfriendTrack *friendTrack = NULL;

    if (esdTrack) {
      if (connectESDFriends){
	friendTrack = esdFriend->GetTrack(idxTrack);
      }
      FillNextTrack( esdTrack, friendTrack);
    }
  }

  // Fill V0s
  
  fV0Pointer = fContentSize;
  fNV0s = 0;
  if( fillV0s ){
    for( int i=0; i < esd->GetNumberOfV0s(); i++){
      AliESDv0 *esdV0 = esd->GetV0( i );
      AliFlatESDV0 *v0 = GetNextV0Pointer();
      if( !v0 ) continue;
      v0->fNegTrackID = esdV0->GetNindex();
      v0->fPosTrackID = esdV0->GetNindex();
      StoreLastV0();      
    }
  }
#endif
  return 0;
}

UInt_t AliFlatESDEvent::CountBits(Byte_t field, UInt_t mask) const {
  // Count bits in field
  UInt_t count = 0; 
  UInt_t reg = 0x0; 
  
  reg |= field;   
  reg &= mask;
  
  for (count = 0; reg; count++)
    reg &= reg - 1; 

  return count;
}
