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
///
// Moved parts of old AliMUONDataInterface interface to AliMUONMCDataInterface
//  Artur Szostak <artursz@iafrica.com> (University of Cape Town)
//-----------------------------------------------------------------------------

#include "AliMUONMCDataInterface.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVHitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONHit.h"
#include "AliMUONVDigit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONGlobalTrigger.h"

#include "AliMpEncodePair.h"
#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpCDB.h"

#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliCDBManager.h"

#include <Riostream.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TParticle.h>
#include <TIterator.h>
#include <cstdlib>
#include <cassert>

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
fIsValid(kFALSE),
fCurrentIteratorType(kNoIterator),
fCurrentIndex(-1),
fDataX(-1),
fDataY(-1),
fIterator(0x0)
{
  /// ctor
  
  ++fgInstanceCounter;
  
  if ( AliCDBManager::Instance() != NULL &&
       AliCDBManager::Instance()->GetDefaultStorage() == NULL ) {
      AliFatal("CDB default storage not defined.");
  }
  
  Open(filename);

  // Load mapping
  if ( ! AliMpCDB::LoadMpSegmentation() ) {
    AliFatal("Could not access mapping from OCDB !");
  }
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
AliMUONVHitStore* 
AliMUONMCDataInterface::HitStore(Int_t event, Int_t track)
{
  /// Return the hitStore for a given track of one event
  /// Return 0x0 if event and/or track not found
  /// Returned pointer should not be deleted
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if (not IsValid()) return 0x0;

  if (event == fCurrentEvent)
  {
    if (track == fDataX and fHitStore != 0x0)  // using fDataX as track number.
      return fHitStore;
  }
  else
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->LoadHits();
  
  TTree* treeH = fLoader->TreeH();
  if (treeH == 0x0)
  {
    AliError("Could not get treeH");
    return 0x0;
  }
  
  fHitStore = AliMUONVHitStore::Create(*treeH);
  AliDebug(1,"Creating hitStore from treeH");
  if ( fHitStore != 0x0 )
  {
    fHitStore->Connect(*treeH);
    if ( treeH->GetEvent(track) == 0 ) 
    {
      AliError(Form("Could not read track %d",track));
      fHitStore->Clear();
      return 0x0;
    }
    fDataX = track; // using fDataX as track number.
  }
  
  fLoader->UnloadHits();

  return fHitStore;
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONMCDataInterface::SDigitStore(Int_t event)
{
  /// Return the SDigit store for a given event.
  /// Return 0 if event not found
  /// Returned pointer should not be deleted
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if (not IsValid()) return 0x0;
  
  if (event == fCurrentEvent)
  {
    if (fSDigitStore != 0x0)
      return fSDigitStore;
  }
  else
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->LoadSDigits();
  
  TTree* treeS = fLoader->TreeS();
  if (treeS == 0x0)
  {
    AliError("Could not get treeS");
    return 0x0;
  }
  
  fSDigitStore = AliMUONVDigitStore::Create(*treeS);
  if ( fSDigitStore != 0x0 )
  {
    fSDigitStore->Clear();
    fSDigitStore->Connect(*treeS);
    treeS->GetEvent(0);
  }
  
  fLoader->UnloadSDigits();
  
  return fSDigitStore;
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONMCDataInterface::DigitStore(Int_t event)
{
  /// Return a pointer to the digitStore for a given event (or 0 if not found)
  /// Returned pointer should not be deleted
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if (not IsValid()) return 0x0;
  
  if (event == fCurrentEvent)
  {
    if (fDigitStore != 0x0)
      return fDigitStore;
  }
  else
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->LoadDigits();
  
  TTree* treeD = fLoader->TreeD();
  if (treeD == 0x0)
  {
    AliError("Could not get treeD");
    return 0x0;
  }
  
  fDigitStore = AliMUONVDigitStore::Create(*treeD);
  if ( fDigitStore != 0x0 ) 
  {
    fDigitStore->Clear();
    fDigitStore->Connect(*treeD);
    treeD->GetEvent(0);
  }
  
  fLoader->UnloadDigits();
  
  return fDigitStore;
}

//_____________________________________________________________________________
AliStack*
AliMUONMCDataInterface::Stack(Int_t event)
{
  /// Get the Stack (list of generated particles) for one event
  /// Returned pointer should not be deleted
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if ( not IsValid() ) return 0x0;

  if (event != fCurrentEvent)
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->GetRunLoader()->LoadKinematics();
  
  return fLoader->GetRunLoader()->Stack();
}

//_____________________________________________________________________________
TClonesArray*
AliMUONMCDataInterface::TrackRefs(Int_t event, Int_t track)
{
  /// Get the track references for a given (generated) track of one event
  /// Returned pointer should not be deleted
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if ( not IsValid() ) return 0x0;
  
  if (event == fCurrentEvent)
  {
    if (track == fDataX and fTrackRefs != 0x0)  // using fDataX as track number.
      return fTrackRefs;
  }
  else
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->GetRunLoader()->LoadTrackRefs();
  
  TTree* treeTR = fLoader->GetRunLoader()->TreeTR();
  
  if ( fTrackRefs != 0x0 ) fTrackRefs->Clear("C");
  
  if (treeTR != 0x0)
  {
    if ( treeTR->GetEvent(track) > 0 ) 
    {
      TBranch* branch = treeTR->GetBranch("TrackReferences");
      branch->SetAddress(&fTrackRefs);
      branch->GetEvent(track);
      fDataX = track;  // using fDataX as track number.
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
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if (not IsValid()) return 0x0;
  
  if (event == fCurrentEvent)
  {
    if (fTriggerStore != 0x0)
      return fTriggerStore;
  }
  else
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->LoadDigits();
  
  TTree* treeD = fLoader->TreeD();
  if ( treeD == 0x0 ) 
  {
    AliError("Could not get treeD");
    return 0x0;
  }
  
  fTriggerStore = AliMUONVTriggerStore::Create(*treeD);
  if ( fTriggerStore != 0x0 )
  {
    fTriggerStore->Clear();
    fTriggerStore->Connect(*treeD);
    treeD->GetEvent(0);
  }
  
  fLoader->UnloadDigits();
  
  return fTriggerStore;
}

//_____________________________________________________________________________
void
AliMUONMCDataInterface::DumpDigits(Int_t event, Bool_t sorted)
{
  /// Dump the digits for a given event, sorted if requested.
  DigitStore(event);
  
  if ( fDigitStore != 0x0 ) 
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
  
  if ( stack != 0x0 ) 
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
  
  if ( fSDigitStore != 0x0 ) 
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
    if ( fTrackRefs != 0x0 ) 
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

  if ( fTriggerStore != 0x0 ) 
  {
    fTriggerStore->Print();
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONMCDataInterface::LoadEvent(Int_t event)
{
  /// Load event if different from the current one.
  /// Returns kFALSE on error and kTRUE if the event was loaded.
  
  assert( IsValid() );
  
  AliDebug(1,Form("Loading event %d using runLoader %p",event,fLoader->GetRunLoader()));
  if (fLoader->GetRunLoader()->GetEvent(event) == 0)
  {
    fCurrentEvent = event;
    return kTRUE;
  }
  else
    return kFALSE;
}


//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfEvents() const
{
  /// Number of events in the file we're connected to
  if (not IsValid()) return -1;
  return fLoader->GetRunLoader()->GetNumberOfEvents();
}

//_____________________________________________________________________________
Int_t 
AliMUONMCDataInterface::NumberOfTracks(Int_t event)
{
  /// Number of tracks in the event
  if ( not IsValid()) return -1;
  
  if (event != fCurrentEvent)
  {
    ResetStores();
    if ( not LoadEvent(event) ) return -1;
  }
  
  fLoader->LoadHits();
  
  Int_t rv(-1);
  
  TTree* treeH = fLoader->TreeH();
  if (treeH != 0x0)
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
  if ( not IsValid()) return -1;
  
  if (event != fCurrentEvent)
  {
    ResetStores();
    if ( not LoadEvent(event) ) return -1;
  }

  fLoader->GetRunLoader()->LoadTrackRefs();
  
  Int_t rv(-1);
  
  TTree* treeTR = fLoader->GetRunLoader()->TreeTR();
  if (treeTR != 0x0)
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
  
  ResetStores();
  
  fCurrentEvent=-1;

  if ( fLoader != 0x0 )
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
  if (runLoader == 0x0)
  {
    AliError(Form("Cannot open file %s",filename));    
    fIsValid = kFALSE;
  }
  
  // Get run number and set it to CDB manager
  runLoader->LoadHeader();
  if ( ! runLoader->GetHeader() ) {
    AliError("Cannot load header.");    
    fIsValid = kFALSE;
  }
  else {
    Int_t runNumber = runLoader->GetHeader()->GetRun();
    AliCDBManager::Instance()->SetRun(runNumber);
  }  
  runLoader->UnloadHeader(); 

  fLoader = runLoader->GetDetectorLoader("MUON");
  if (fLoader == 0x0)
  {
    AliError("Cannot get AliMUONLoader");
    fIsValid = kFALSE;
  }
  
  if (not IsValid())
  {
    AliError(Form("Could not access %s filename. Object is unuseable",filename));
  }
}

//_____________________________________________________________________________
Bool_t AliMUONMCDataInterface::GetEvent(Int_t event)
{
/// Loads all simulated data for the given event.

  if (HitStore(event, 0) == 0x0) return kFALSE;
  if (SDigitStore(event) == 0x0) return kFALSE;
  if (DigitStore(event) == 0x0) return kFALSE;
  if (TriggerStore(event) == 0x0) return kFALSE;
  if (TrackRefs(event, 0) == 0x0) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfParticles()
{
/// Returns the total number of particles in the kinematics tree.

  AliStack* stack = Stack(fCurrentEvent);
  if ( stack == 0x0 ) return -1;
  return (Int_t) stack->GetNtrack();
}

//_____________________________________________________________________________
TParticle* AliMUONMCDataInterface::Particle(Int_t index)
{
/// Returns the index'th particle in the kinematics tree.
/// @param index  The index number of the particle in the range [0 ... N-1]
///               where N = NumberOfParticles()

  AliStack* stack = Stack(fCurrentEvent);
  if ( stack == 0x0 ) return 0x0;
  return static_cast<TParticle*>( stack->Particle(index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfTracks()
{
/// Returns the number of primary tracks (from primary particles) in the current event.

  return NumberOfTracks(fCurrentEvent);
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfHits(Int_t track)
{
/// Returns the number of hits for a given primary track/particle.
/// @param track  The track number in the range [0 .. N-1]
///               where N = NumberOfTracks()

  TIterator* iter = GetIterator(kHitIterator, track);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONHit* 
AliMUONMCDataInterface::Hit(Int_t track, Int_t index)
{
/// Returns a pointer to the index'th hit object.
/// @param track  The track number in the range [0 .. N-1]
///               where N = NumberOfTracks()
/// @param index  The index number of the hit in the range [0 ... M-1]
///               where M = NumberOfHits(track)

  TIterator* iter = GetIterator(kHitIterator, track);
  return static_cast<AliMUONHit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfSDigits(Int_t detElemId)
{
/// Returns the number of summable digits to be found on a given detector element.
/// @param detElemId  The detector element ID number to search on.

  TIterator* iter = GetIterator(kSDigitIteratorByDetectorElement, detElemId);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVDigit* AliMUONMCDataInterface::SDigit(Int_t detElemId, Int_t index)
{
/// Returns the a pointer to the index'th summable digit on the specified detector element.
/// @param detElemId  The detector element ID number to search on.
/// @param index  The index number of the digit to fetch in the range [0 .. N-1],
///   where N = NumberOfDigits(detElemId)

  TIterator* iter = GetIterator(kSDigitIteratorByDetectorElement, detElemId);
  return static_cast<AliMUONVDigit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfSDigits(Int_t chamber, Int_t cathode)
{
/// Returns the number of summable digits to be found on a specific chamber and cathode.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param cathode  The cathode in the range [0 .. 1], where 0 is the bending and
///   1 is the non-bending plane.

  TIterator* iter = GetIterator(kSDigitIteratorByChamberAndCathode, chamber, cathode);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVDigit* AliMUONMCDataInterface::SDigit(Int_t chamber, Int_t cathode, Int_t index)
{
/// Returns the a pointer to the index'th summable digit on the specified chamber and cathode.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param cathode  The cathode in the range [0 .. 1], where 0 is the bending and
///   1 is the non-bending plane.
/// @param index  The index number of the digit to fetch in the range [0 .. N-1],
///   where N = NumberOfDigits(chamber, cathode)

  TIterator* iter = GetIterator(kSDigitIteratorByChamberAndCathode, chamber, cathode);
  return static_cast<AliMUONVDigit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfDigits(Int_t detElemId)
{
/// Returns the number of simulated digits to be found on a given detector element.
/// @param detElemId  The detector element ID number to search on.

  TIterator* iter = GetIterator(kDigitIteratorByDetectorElement, detElemId);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVDigit* AliMUONMCDataInterface::Digit(Int_t detElemId, Int_t index)
{
/// Returns the a pointer to the index'th simulated digit on the specified detector element.
/// @param detElemId  The detector element ID number to search on.
/// @param index  The index number of the digit to fetch in the range [0 .. N-1],
///   where N = NumberOfDigits(detElemId)

  TIterator* iter = GetIterator(kDigitIteratorByDetectorElement, detElemId);
  return static_cast<AliMUONVDigit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfDigits(Int_t chamber, Int_t cathode)
{
/// Returns the number of simulated digits to be found on a specific chamber and cathode.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param cathode  The cathode in the range [0 .. 1], where 0 is the bending and
///   1 is the non-bending plane.

  TIterator* iter = GetIterator(kDigitIteratorByChamberAndCathode, chamber, cathode);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVDigit* AliMUONMCDataInterface::Digit(Int_t chamber, Int_t cathode, Int_t index)
{
/// Returns the a pointer to the index'th simulated digit on the specified chamber and cathode.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param cathode  The cathode in the range [0 .. 1], where 0 is the bending and
///   1 is the non-bending plane.
/// @param index  The index number of the digit to fetch in the range [0 .. N-1],
///   where N = NumberOfDigits(chamber, cathode)

  TIterator* iter = GetIterator(kDigitIteratorByChamberAndCathode, chamber, cathode);
  return static_cast<AliMUONVDigit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfLocalTriggers()
{
/// Returns the number of simulated local trigger objects.

  TIterator* iter = GetIterator(kLocalTriggerIterator);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONLocalTrigger* AliMUONMCDataInterface::LocalTrigger(Int_t index)
{
/// Returns a pointer to the index'th simulated local trigger object.
/// @param index  The index number of the local trigger object to fetch in the range [0 .. N-1],
///   where N = NumberOfLocalTriggers()

  TIterator* iter = GetIterator(kLocalTriggerIterator);
  return static_cast<AliMUONLocalTrigger*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfRegionalTriggers()
{
/// Returns the number of simulated regional trigger objects.

  TIterator* iter = GetIterator(kRegionalTriggerIterator);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONRegionalTrigger* AliMUONMCDataInterface::RegionalTrigger(Int_t index)
{
/// Returns a pointer to the index'th simulated regional trigger object.
/// @param index  The index number of the regional trigger object to fetch in the range [0 .. N-1],
///   where N = NumberOfRegionalTriggers()

  TIterator* iter = GetIterator(kRegionalTriggerIterator);
  return static_cast<AliMUONRegionalTrigger*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
AliMUONGlobalTrigger* AliMUONMCDataInterface::GlobalTrigger()
{
/// Returns a pointer to the simulated global trigger object for the event.

  AliMUONVTriggerStore* store = TriggerStore(fCurrentEvent);
  if (store == 0x0) return 0x0;
  return store->Global();
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::NumberOfTrackRefs()
{
/// Number of track references in the currently selected event.

  return NumberOfTrackRefs(fCurrentEvent);
}

//_____________________________________________________________________________
TClonesArray* AliMUONMCDataInterface::TrackRefs(Int_t track)
{
/// Returns the track references for a given track in the current event.
/// @param track  The track to returns track references for. In the range [0 .. N-1]
///               where N = NumberOfTrackRefs()

  return TrackRefs(fCurrentEvent, track);
}

//_____________________________________________________________________________
void AliMUONMCDataInterface::ResetStores()
{
/// Deletes all the store objects that have been created and resets the pointers to 0x0.
/// The temporary iterator object is automatically reset. See ResetIterator for more details.

  ResetIterator();
  if (fHitStore != 0x0)
  {
    delete fHitStore;
    fHitStore = 0x0;
  }
  if (fSDigitStore != 0x0)
  {
    delete fSDigitStore;
    fSDigitStore = 0x0;
  }
  if (fDigitStore != 0x0)
  {
    delete fDigitStore;
    fDigitStore = 0x0;
  }
  if (fTrackRefs != 0x0)
  {
    delete fTrackRefs;
    fTrackRefs = 0x0;
  }
  if (fTriggerStore != 0x0)
  {
    delete fTriggerStore;
    fTriggerStore = 0x0;
  }
}

//_____________________________________________________________________________
TIterator* AliMUONMCDataInterface::GetIterator(IteratorType type, Int_t x, Int_t y)
{
/// Creates an appropriate iterator object and returns it.
/// If the iterator has already been created then that one is returned otherwise
/// a new object is created.
/// Depending on the value of 'type' the semantics of parameters x and y can change.
/// @param type  The type of iterator to create.
/// @param x  This is the detector element ID if type equals kDigitIteratorByDetectorElement
///           or kSDigitIteratorByDetectorElement.
///           If type equals kDigitIteratorByChamberAndCathode or
///           kSDigitIteratorByChamberAndCathode then this is the chamber number.
///           For type == kHitIterator the parameter x is the track number.
///           In all other cases this parameter is ignored.
/// @param y  If type equals kDigitIteratorByChamberAndCathode or
///           kSDigitIteratorByChamberAndCathode then this parameter is the cathode
///           number. In all other cases this parameter is ignored.

  if (type == fCurrentIteratorType and fDataX == x and fDataY == y)
  	return fIterator;
  
  if (fCurrentEvent == -1)
  {
    AliError("No event was selected. Try first using GetEvent().");
    return 0x0;
  }
  
  ResetIterator();
  
  switch (type)
  {
  case kHitIterator:
    {
      Int_t track = x;
      AliMUONVHitStore* store = HitStore(fCurrentEvent, track);
      if (store == 0x0) return 0x0;
      fIterator = store->CreateIterator();
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kHitIterator;
      return fIterator;
    }
    
  case kSDigitIteratorByDetectorElement:
    {
      Int_t detElem = x;
      AliMUONVDigitStore* store = SDigitStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      fIterator = store->CreateIterator(detElem, detElem, 2);
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kSDigitIteratorByDetectorElement;
      fDataX = detElem;
      return fIterator;
    }
    
  case kSDigitIteratorByChamberAndCathode:
    {
      Int_t chamber = x;
      Int_t cathode = y;
      if (chamber < 0 or AliMpConstants::NofChambers() <= chamber)
      {
        AliError(Form(
          "Must have give a chamber value in the range [0..%d], but got a value of: %d",
          AliMpConstants::NofChambers() - 1,
          chamber
        ));
        return 0x0;
      }
      if (cathode < 0 or 1 < cathode)
      {
        AliError(Form("Must have give a cathode value in the range [0..1], but got a value of: %d", cathode));
        return 0x0;
      }
      
      AliMUONVDigitStore* store = SDigitStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      MpPair_t pair = AliMpDEManager::GetDetElemIdRange(chamber);
      fIterator = store->CreateIterator(AliMp::PairFirst(pair), AliMp::PairSecond(pair), cathode);
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kSDigitIteratorByChamberAndCathode;
      fDataX = chamber;
      fDataY = cathode;
      return fIterator;
    }
    
  case kDigitIteratorByDetectorElement:
    {
      Int_t detElem = x;
      AliMUONVDigitStore* store = DigitStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      fIterator = store->CreateIterator(detElem, detElem, 2);
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kDigitIteratorByDetectorElement;
      fDataX = detElem;
      return fIterator;
    }
    
  case kDigitIteratorByChamberAndCathode:
    {
      Int_t chamber = x;
      Int_t cathode = y;
      if (chamber < 0 or AliMpConstants::NofChambers() <= chamber)
      {
        AliError(Form(
          "Must have give a chamber value in the range [0..%d], but got a value of: %d",
          AliMpConstants::NofChambers() - 1,
          chamber
        ));
        return 0x0;
      }
      if (cathode < 0 or 1 < cathode)
      {
        AliError(Form("Must have give a cathode value in the range [0..1], but got a value of: %d", cathode));
        return 0x0;
      }
      
      AliMUONVDigitStore* store = DigitStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      MpPair_t pair = AliMpDEManager::GetDetElemIdRange(chamber);
      fIterator = store->CreateIterator(AliMp::PairFirst(pair), AliMp::PairSecond(pair), cathode);
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kDigitIteratorByChamberAndCathode;
      fDataX = chamber;
      fDataY = cathode;
      return fIterator;
    }
    
  case kLocalTriggerIterator:
    {
      AliMUONVTriggerStore* store = TriggerStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      fIterator = store->CreateLocalIterator();
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kLocalTriggerIterator;
      return fIterator;
    }
    
  case kRegionalTriggerIterator:
    {
      AliMUONVTriggerStore* store = TriggerStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      fIterator = store->CreateRegionalIterator();
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kRegionalTriggerIterator;
      return fIterator;
    }
    
  default:
    return 0x0;
  }
}

//_____________________________________________________________________________
void AliMUONMCDataInterface::ResetIterator()
{
/// The temporary iterator object is deleted if it exists and the pointer reset to 0x0.
/// The iterator type and temporary data indicating the state of the iterator are
/// also reset.

  if (fIterator != 0x0) delete fIterator;
  fCurrentIteratorType = kNoIterator;
  fCurrentIndex = fDataX = fDataY = -1;
  fIterator = 0x0;
}

//_____________________________________________________________________________
Int_t AliMUONMCDataInterface::CountObjects(TIterator* iter)
{
/// Counts the number of objects in the iterator and resets it.
/// @return The number of objects in 'iter'.

  if (iter == 0x0) return -1;
  Int_t count = 0;
  iter->Reset();
  while ( iter->Next() != 0x0 ) count++;
  iter->Reset();
  fCurrentIndex = -1;
  return count;
}

//_____________________________________________________________________________
TObject* AliMUONMCDataInterface::FetchObject(TIterator* iter, Int_t index)
{
/// Fetches the index'th object from the iterator counting the first object
/// returned by iterator after it is reset as index == 0. The next object
/// has index == 1 and so on where the last object returned by the iterator
/// has index == N-1 where N = CountObjects(iter)
/// This method will only reset the iterator if index is smaller than
/// fCurrentIndex, which is used to track the iteration progress and is
/// updated when a new object if returned by this method.
/// @param iter  The iterator to fetch an object from.
/// @param index The index number of the object to fetch in the range [0 .. N-1]
///        where N = CountObjects(iter)

  if (index < 0)
  {
    AliError(Form("Index is out of bounds. Got a value of %d.", index));
    return 0x0;
  }

  if (iter == 0x0) return 0x0;
  if (index <= fCurrentIndex)
  {
    iter->Reset();
    fCurrentIndex = -1;
  }
  
  TObject* object = 0x0;
  while (fCurrentIndex < index)
  {
    object = iter->Next();
    if (object == 0x0)
    {
      AliError(Form("Index is out of bounds. Got a value of %d.", index));
      iter->Reset();
      fCurrentIndex = -1;
      return 0x0;
    }
    fCurrentIndex++;
  }
  return object;
}
