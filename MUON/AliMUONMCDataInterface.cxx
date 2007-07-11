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

//-----------------------------------------------------------------------------
/// \class AliMUONMCDataInterface
///
/// Easy to use MC data accessor
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONMCDataInterface.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVHitStore.h"
#include "AliMUONVStore.h"
#include "AliMUONVTriggerStore.h"

#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliStack.h"

#include <Riostream.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TParticle.h>

/// \cond CLASSIMP
ClassImp(AliMUONMCDataInterface)
/// \endcond

Int_t AliMUONMCDataInterface::fgInstanceCounter(0);

//_____________________________________________________________________________
AliMUONMCDataInterface::AliMUONMCDataInterface(const char* filename) :
TObject(),
fLoader(0x0),
fHitStore(0x0),
fSDigitStore(0x0),
fDigitStore(0x0),
fTriggerStore(0x0),
fTrackRefs(0x0),
fCurrentEvent(-1),
fIsValid(kFALSE)
{
  /// ctor
  
  ++fgInstanceCounter;
  
  Open(filename);
}

//_____________________________________________________________________________
AliMUONMCDataInterface::~AliMUONMCDataInterface()
{
  /// dtor
  if ( fLoader ) 
  {
    delete fLoader->GetRunLoader();
  }
  --fgInstanceCounter;
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONMCDataInterface::DigitStore(Int_t event)
{
  /// Return a pointer to the digitStore for a given event (or 0 if not found)
  /// Returned pointer should not be deleted
  
  if ( LoadEvent(event) ) return 0x0;
  
  fLoader->LoadDigits();
  
  TTree* treeD = fLoader->TreeD();
  
  if (!treeD)
  {
    AliError("Could not get treeD");
    return 0x0;
  }
  
  if (!fDigitStore)
  {
    fDigitStore = AliMUONVDigitStore::Create(*treeD);
  }
  
  if ( fDigitStore ) 
  {
    fDigitStore->Clear();
    fDigitStore->Connect(*treeD);
    treeD->GetEvent(0);
  }
  
  fLoader->UnloadDigits();
  
  return fDigitStore;
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpDigits(Int_t event, Bool_t sorted)
{
  /// Dump the digits for a given event, sorted if requested.
  DigitStore(event);
  
  if ( fDigitStore ) 
  {
    if ( sorted ) 
    {
      DumpSorted(*fDigitStore);
    }
    else
    {
      fDigitStore->Print();
    }  
  }
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpHits(Int_t event)
{
  /// Dump all the hits for one event
  
  Int_t ntracks = NumberOfTracks(event);
  
  for ( Int_t i = 0; i < ntracks; ++i ) 
  {
    cout << ">> Track " << i << endl;
    HitStore(event,i);
    if ( fHitStore )
    {
      fHitStore->Print("","full");
    }
  }
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpKine(Int_t event)
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
AliMUONMCDataInterface::DumpSDigits(Int_t event, Bool_t sorted)
{
  /// Dump the SDigits for a given event, sorted if requested
  SDigitStore(event);
  
  if ( fSDigitStore ) 
  {
    if ( sorted ) 
    {
      DumpSorted(*fSDigitStore);
    }
    else
    {
      fSDigitStore->Print();
    }
  }
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpSorted(const AliMUONVStore& store) const
{
  /// Dump the given store in sorted order
  
  TIter next(store.CreateIterator());
  TObject* object;
  TList list;
  list.SetOwner(kFALSE);
  
  while ( ( object = next() ) )
  {
    list.Add(object);
  }
  
  list.Sort();
  
  list.Print();
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpTrackRefs(Int_t event)
{
  /// Dump track references for one event
  Int_t ntrackrefs = NumberOfTrackRefs(event);
  
  for ( Int_t i = 0; i < ntrackrefs; ++i ) 
  {
    TrackRefs(event,i);
    if ( fTrackRefs ) 
    {
      fTrackRefs->Print("","*");
    }
  }
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpTrigger(Int_t event)
{
  /// Dump trigger for a given event (trigger is read from TreeD)
  
  TriggerStore(event);

  if ( fTriggerStore ) 
  {
    fTriggerStore->Print();
  }
}

//_____________________________________________________________________________
AliMUONVHitStore* 
AliMUONMCDataInterface::HitStore(Int_t event, Int_t track)
{
  /// Return the hitStore for a given track of one event
  /// Return 0 if event and/or track not found
  /// Returned pointer should not be deleted
  
  if ( !IsValid() ) return 0x0;
  
  if ( LoadEvent(event) ) return 0x0;

  if ( fHitStore) fHitStore->Clear();
  
  fLoader->LoadHits();
  
  TTree* treeH = fLoader->TreeH();
  
  if (!treeH) 
  {
    AliError("Could not get treeH");
    return 0x0;
  }
  
  if ( !fHitStore ) 
  {
    fHitStore = AliMUONVHitStore::Create(*treeH);
    AliDebug(1,"Creating hitStore from treeH");
  }
  
  if ( fHitStore )
  {
    fHitStore->Connect(*treeH);
    if ( treeH->GetEvent(track) == 0 ) 
    {
      AliError(Form("Could not read track %d",track));
      fHitStore->Clear();
    }
  }
  
  fLoader->UnloadHits();

  return fHitStore;
}

//_____________________________________________________________________________
Bool_t
AliMUONMCDataInterface::IsValid() const
{
  /// Whether we were initialized properly or not
  return fIsValid;
}

//_____________________________________________________________________________
Int_t
AliMUONMCDataInterface::LoadEvent(Int_t event)
{
  /// Load event if different from the current one.
  if ( event != fCurrentEvent ) 
  {
    fCurrentEvent = event;
    AliDebug(1,Form("Loading event %d using runLoader %p",event,fLoader->GetRunLoader()));
    if ( event < NumberOfEvents() )
    {
      return fLoader->GetRunLoader()->GetEvent(event);
    }
    else
    {
      return 1;
    }
  }
  return 0;
}


//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfEvents() const
{
  /// Number of events in the file we're connected to
  if (!IsValid()) return 0;
  return fLoader->GetRunLoader()->GetNumberOfEvents();
}

//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfTracks(Int_t event)
{
  /// Number of tracks in the event
  if (!IsValid()) return 0;
  
  if ( LoadEvent(event) ) return 0;
  
  fLoader->LoadHits();
  
  Int_t rv(0);
  
  TTree* treeH = fLoader->TreeH();
  if (treeH)
  {
    rv = static_cast<Int_t>(treeH->GetEntries());
  }
  else
  {
    AliError("Could not get TreeH");
  }

  fLoader->UnloadHits();
  
  return rv;
}

//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfTrackRefs(Int_t event)
{
  /// Number of track references in the event
  if (!IsValid()) return 0;
  
  if ( LoadEvent(event) ) return 0;

  fLoader->GetRunLoader()->LoadTrackRefs();
  
  Int_t rv(0);
  
  TTree* treeTR = fLoader->GetRunLoader()->TreeTR();
  if (treeTR)
  {
    rv = static_cast<Int_t>(treeTR->GetEntries());
  }
  else
  {
    AliError("Could not get TreeTR");
  }
  
  fLoader->GetRunLoader()->UnloadTrackRefs();
  
  return rv;
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::Open(const char* filename)
{
  /// Connect to a given galice.root file
  
  delete fHitStore; 
  fHitStore=0x0;
  delete fSDigitStore;
  fSDigitStore=0x0;
  delete fDigitStore;
  fDigitStore=0x0;
  delete fTrackRefs;
  fTrackRefs=0x0;
  delete fTriggerStore;
  fTriggerStore=0x0;
  
  fCurrentEvent=-1;

  if ( fLoader ) 
  {
    delete fLoader->GetRunLoader();
  }
  
  fLoader = 0x0;
  
  fIsValid = kTRUE;
  
  TString foldername(Form("%s-%d",ClassName(),fgInstanceCounter));
  
  while (AliRunLoader::GetRunLoader(foldername)) 
  {
    delete AliRunLoader::GetRunLoader(foldername);
  }
  
  AliRunLoader* runLoader = AliRunLoader::Open(filename,foldername);
  if (!runLoader) 
  {
    AliError(Form("Cannot open file %s",filename));    
    fIsValid = kFALSE;
  }
  fLoader = runLoader->GetDetectorLoader("MUON");
  if (!fLoader) 
  {
    AliError("Cannot get AliMUONLoader");
    fIsValid = kFALSE;
  }
  
  if (!IsValid())
  {
    AliError(Form("Could not access %s filename. Object is unuseable",filename));
  }
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONMCDataInterface::SDigitStore(Int_t event)
{
  /// Return the SDigit store for a given event.
  /// Return 0 if event not found
  /// Returned pointer should not be deleted
  
  if ( LoadEvent(event) ) return 0x0;
  
  fLoader->LoadSDigits();
  
  TTree* treeS = fLoader->TreeS();
  
  if (!treeS)
  {
    AliError("Could not get treeS");
    return 0x0;
  }
  
  if (!fSDigitStore)
  {
    fSDigitStore = AliMUONVDigitStore::Create(*treeS);
  }
  
  if ( fSDigitStore ) 
  {
    fSDigitStore->Clear();
    fSDigitStore->Connect(*treeS);
    treeS->GetEvent(0);
  }
  
  fLoader->UnloadSDigits();
  
  return fSDigitStore;
}

//_____________________________________________________________________________
AliStack*
AliMUONMCDataInterface::Stack(Int_t event)
{
  /// Get the Stack (list of generated particles) for one event
  /// Returned pointer should not be deleted
  
  if (!IsValid()) return 0x0;

  if ( LoadEvent(event) ) return 0;
  
  fLoader->GetRunLoader()->LoadKinematics();
  
  return fLoader->GetRunLoader()->Stack();
}

//_____________________________________________________________________________
TClonesArray*
AliMUONMCDataInterface::TrackRefs(Int_t event, Int_t track)
{
  /// Get the track references for a given (generated) track of one event
  /// Returned pointer should not be deleted
  
  if ( !IsValid() ) return 0x0;

  if ( LoadEvent(event) ) return 0;
  
  fLoader->GetRunLoader()->LoadTrackRefs();
  
  TTree* treeTR = fLoader->GetRunLoader()->TreeTR();
  
  if ( fTrackRefs ) fTrackRefs->Clear("C");
  
  if (treeTR)
  {
    if ( treeTR->GetEvent(track) > 0 ) 
    {
      TBranch* branch = treeTR->GetBranch("MUON");
      branch->SetAddress(&fTrackRefs);
      branch->GetEvent(track);
    }
  }
  else
  {
    AliError("Could not get TreeTR");
  }
  
  fLoader->GetRunLoader()->UnloadTrackRefs();
  
  return fTrackRefs;
}

//_____________________________________________________________________________
AliMUONVTriggerStore*
AliMUONMCDataInterface::TriggerStore(Int_t event)
{
  /// Return the triggerStore for a given event.
  /// Return 0x0 if event not found.
  /// Returned pointer should not be deleted.
  
  if ( LoadEvent(event) ) return 0x0;
  
  fLoader->LoadDigits();
  
  TTree* treeD = fLoader->TreeD();
  
  if ( !treeD ) 
  {
    AliError("Could not get treeD");
    return 0x0;
  }
  
  if (!fTriggerStore)
  {
    fTriggerStore = AliMUONVTriggerStore::Create(*treeD);
  }
  
  if ( fTriggerStore ) 
  {
    fTriggerStore->Clear();
    fTriggerStore->Connect(*treeD);
    treeD->GetEvent(0);
  }
  
  fLoader->UnloadDigits();
  
  return fTriggerStore;
}
