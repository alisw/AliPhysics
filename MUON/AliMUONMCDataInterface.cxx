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

// $Id$

/// \class AliMUONMCDataInterface
///
/// Easy to use MC data accessor
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONMCDataInterface.h"

#include "AliLog.h"
#include "AliMUONDataManager.h"
#include "AliMUONVStore.h"
#include "AliMUONVHitStore.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include <Riostream.h>
#include <TClonesArray.h>
#include <TParticle.h>

/// \cond CLASSIMP
ClassImp(AliMUONMCDataInterface)
/// \endcond

//_____________________________________________________________________________
AliMUONMCDataInterface::AliMUONMCDataInterface(const char* filename) :
TObject(),
fDataManager(new AliMUONDataManager(filename))
{
  /// ctor
  if (!IsValid())
  {
    AliError(Form("Could not access %s filename. Object is unuseable",filename));
  }
}

//_____________________________________________________________________________
AliMUONMCDataInterface::~AliMUONMCDataInterface()
{
  /// dtor
  delete fDataManager;
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpHits(Int_t event) const
{
  /// Dump all the hits for one event
  
  Int_t ntracks = NumberOfTracks(event);
  
  for ( Int_t i = 0; i < ntracks; ++i ) 
  {
    cout << ">> Track " << i << endl;
    AliMUONVHitStore* hitStore = HitStore(event,i);
    hitStore->Print("","full");
    delete hitStore;
  }
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpKine(Int_t event) const
{
  /// Dump all generated particles for one event
  AliStack* stack = Stack(event);
  
  if ( stack ) 
  {
    Int_t nparticles = (Int_t) stack->GetNtrack();
  
    for (Int_t iparticle=0; iparticle<nparticles; ++iparticle) 
    {
      stack->Particle(iparticle)->Print("");  
    }
  }
  else
  {
    AliError("Could not get stack");
  }
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpTrackRefs(Int_t event) const
{
  /// Dump track references for one event
  Int_t ntrackrefs = NumberOfTrackRefs(event);
  
  for ( Int_t i = 0; i < ntrackrefs; ++i ) 
  {
    TClonesArray* trackRefs = TrackRefs(event,i);
    trackRefs->Print("","*");
    delete trackRefs;
  }
}

//_____________________________________________________________________________
AliMUONVHitStore* 
AliMUONMCDataInterface::HitStore(Int_t event, Int_t track) const
{
  /// Return the hitStore for a given track of one event
  if ( !IsValid() ) return 0x0;
  
  fDataManager->Load("H");

  if ( fDataManager->Load(event) ) return 0x0;
  
  TTree* treeH = fDataManager->Tree("H");
  
  AliMUONVHitStore* hitStore(0x0);
  
  if (treeH)
  {
    hitStore = AliMUONVHitStore::Create(*treeH);
    if ( hitStore )
    {
      hitStore->Connect(*treeH);
      if ( treeH->GetEvent(track) == 0 ) 
      {
        hitStore = 0x0;
      }
    }
  }
  else
  {
    AliError("Could not get TreeH");
  }
  
  fDataManager->Unload("H");
  
  return hitStore;
}

//_____________________________________________________________________________
Bool_t
AliMUONMCDataInterface::IsValid() const
{
  /// Whether we were initialized properly or not
  return fDataManager->IsValid();
}

//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfEvents() const
{
  /// Number of events in the file we're connected to
  if (!IsValid()) return 0;
  return fDataManager->NumberOfEvents();
}

//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfTracks(Int_t event) const
{
  /// Number of tracks in the event
  if (!IsValid()) return 0;
  
  fDataManager->Load("H");

  if ( fDataManager->Load(event) ) return 0;
  
  Int_t rv(0);
  
  TTree* treeH = fDataManager->Tree("H");
  if (treeH)
  {
    rv=static_cast<Int_t>(fDataManager->Tree("H")->GetEntries());
  }
  else
  {
    AliError("Could not get TreeH");
  }

  fDataManager->Unload("H");
  
  return rv;
}

//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfTrackRefs(Int_t event) const
{
  /// Number of track references in the event
  if (!IsValid()) return 0;
  
  fDataManager->Load("TR");
  
  if ( fDataManager->Load(event) ) return 0;
  
  Int_t rv(0);
  
  TTree* treeH = fDataManager->Tree("TR");
  if (treeH)
  {
    rv=static_cast<Int_t>(fDataManager->Tree("TR")->GetEntries());
  }
  else
  {
    AliError("Could not get TreeTR");
  }
  
  fDataManager->Unload("TR");
  
  return rv;
}

//_____________________________________________________________________________
AliStack*
AliMUONMCDataInterface::Stack(Int_t event) const
{
  /// Get the Stack (list of generated particles) for one event
  if (!IsValid()) return 0x0;

  fDataManager->Load("K");

  if ( fDataManager->Load(event) ) return 0x0;
  
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader();
  AliDebug(1,Form("RunLoader=%p",runLoader));
  return runLoader->Stack();
}

//_____________________________________________________________________________
TClonesArray*
AliMUONMCDataInterface::TrackRefs(Int_t event, Int_t track) const
{
  /// Get the track references for a given (generated) track of one event
  if ( !IsValid() ) return 0x0;
  
  fDataManager->Load("TR");
  
  if ( fDataManager->Load(event) ) return 0x0;
  
  TTree* treeTR = fDataManager->Tree("TR");
  
  TClonesArray* rv = 0;
  
  if (treeTR)
  {
    if ( treeTR->GetEvent(track) > 0 ) 
    {
      TBranch* branch = treeTR->GetBranch("MUON");
      branch->SetAddress(&rv);
      branch->GetEvent(track);
    }
  }
  else
  {
    AliError("Could not get TreeTR");
  }
  
  fDataManager->Unload("TR");
  
  return rv;
}
