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
 
/* $Id$ */

#include "AliMUONDataInterface.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONVDigit.h"
#include "AliMUONVCluster.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMpCDB.h"

#include "AliMpEncodePair.h"
#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpCDB.h"

#include "AliLoader.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliCDBManager.h"

#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TIterator.h>
#include <cstdlib>
#include <cassert>

//-----------------------------------------------------------------------------
/// \class AliMUONDataInterface
///
/// An easy to use interface to the MUON data data stored in
/// TreeS, TreeD and TreeR.
///
/// For MC related information (i.e. TreeH, TreeK, TreeTR), see
/// AliMUONMCDataInterface.
///
///
/// This interface in not necessarily the fastest way to fetch the data but
/// it is the easiest.
///
/// \author Laurent Aphecetche, Subatech & Artur Szostak <artursz@iafrica.com> (University of Cape Town)
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONDataInterface)
/// \endcond


Int_t AliMUONDataInterface::fgInstanceCounter(0);

//______________________________________________________________________________
AliMUONDataInterface::AliMUONDataInterface(const char* filename)
: TObject(), 
fLoader(0x0),
fDigitStore(0x0),
fTriggerStore(0x0),
fClusterStore(0x0),
fCurrentEvent(-1),
fTreeLetter(""),
fIsValid(kFALSE),
fCurrentIteratorType(kNoIterator),
fCurrentIndex(-1),
fDataX(-1),
fDataY(-1),
fIterator(0x0)
{
  /// ctor
  /// @param filename should be the full path to a valid galice.root file
  
  ++fgInstanceCounter;
  
  if ( AliCDBManager::Instance() != NULL &&
       AliCDBManager::Instance()->GetDefaultStorage() == NULL ) {
      AliFatal("CDB default storage not defined.");
  }
  
  Open(filename);

  // Load mapping
  if ( ! AliMpCDB::LoadDDLStore() ) {
    AliFatal("Could not access mapping from OCDB !");
  }
}

//______________________________________________________________________________
AliMUONDataInterface::~AliMUONDataInterface()
{
  /// dtor
  ResetStores();
  if ( fLoader != 0x0 ) 
  {
    delete fLoader->GetRunLoader();
  }
  --fgInstanceCounter;  
}

//______________________________________________________________________________
AliMUONVDigitStore*
AliMUONDataInterface::DigitStore(Int_t event)
{
  /// Return digitStore for a given event.
  /// Return 0x0 if event not found.
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

//______________________________________________________________________________
AliMUONVClusterStore*
AliMUONDataInterface::ClusterStore(Int_t event)
{
  /// Return clusterStore for a given event.
  /// Return 0x0 if event not found.
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
    if (fClusterStore != 0x0)
      return fClusterStore;
  }
  else
  {
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  fLoader->LoadRecPoints();
  
  TTree* treeR = fLoader->TreeR();
  if (treeR == 0x0)
  {
    AliError("Could not get treeR");
    return 0x0;
  }
  
  fClusterStore = AliMUONVClusterStore::Create(*treeR);
  if ( fClusterStore != 0x0 ) 
  {
    fClusterStore->Clear();
    fClusterStore->Connect(*treeR);
    treeR->GetEvent(0);
  }
  
  fLoader->UnloadRecPoints();
  
  return fClusterStore;
}

//_____________________________________________________________________________
AliMUONVTriggerStore*
AliMUONDataInterface::TriggerStore(Int_t event, const char* treeLetter)
{
  /// Return the triggerStore for a given event.
  /// Return 0x0 if event not found.
  /// Returned pointer should not be deleted
  /// treeLetter can be R or D to tell from which tree to read the information
  ///
  /// \note If a previous store has been retrieved by one of the methods of
  /// this class, but for a different event number, then those stores will
  /// be deleted and no longer valid.
  /// If you require access to the data for the earlier retrieved store,
  /// but for different events, then you should deep copy / clone the object.
  
  if (not IsValid()) return 0x0;
  
  if (event == fCurrentEvent)
  {
    if (fTreeLetter == treeLetter)
    {
      if (fTriggerStore != 0x0)
        return fTriggerStore;
    }
    else
    {
      // Reset only the fTriggerStore since the others might still be valid
      // for the same event.
      if (fTriggerStore != 0x0)
      {
        delete fTriggerStore;
        fTriggerStore = 0x0;
      }
    }
  }
  else
  {
    // Event has changed so reset all the stores.
    ResetStores();
    if ( not LoadEvent(event) ) return 0x0;
  }
  
  TTree* tree(0x0);
  
  TString stree(treeLetter);
  stree.ToUpper();
  
  if ( stree == "D" )
  {
    fLoader->LoadDigits();    
    tree = fLoader->TreeD();
  }
  else if ( stree == "R" )
  {
    fLoader->LoadRecPoints();
    tree = fLoader->TreeR();
  }
  
  if ( tree == 0x0 ) 
  {
    AliError(Form("Could not get tree%s",treeLetter));
    return 0x0;
  }
  
  fTriggerStore = AliMUONVTriggerStore::Create(*tree);
  if ( fTriggerStore != 0x0 ) 
  {
    fTriggerStore->Clear();
    fTriggerStore->Connect(*tree);
    tree->GetEvent(0);
  }
  
  if ( stree == "D" )
  {
    fLoader->UnloadDigits();    
  }
  else if ( stree == "R" )
  {
    fLoader->UnloadRecPoints();
  }
  fTreeLetter = stree;
  
  return fTriggerStore;
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpDigits(Int_t event, Bool_t sorted)
{
  /// Dump the digits for a given event, sorted if so required
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

//______________________________________________________________________________
void
AliMUONDataInterface::DumpRecPoints(Int_t event, Bool_t sorted)
{
  /// Dump the recpoints for a given event, sorted if so required
  ClusterStore(event);
  if ( fClusterStore != 0x0 ) 
  {
    if ( sorted ) 
    {
      DumpSorted(*fClusterStore);
    }
    else
    {
      fClusterStore->Print();
    }
  }
}

//_____________________________________________________________________________
void
AliMUONDataInterface::DumpSorted(const AliMUONVStore& store) const
{
  /// Dump the given store, in sorted order
  
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
AliMUONDataInterface::DumpTrigger(Int_t event, const char* treeLetter)
{
  /// Dump trigger for a given event from a given tree (if event>=0)
  /// or loop over all events and build a trigger ntuple if event<0
  /// treeLetter can be R or D to tell from which tree to read the information
  
  if ( event < 0 ) 
  {
    NtupleTrigger(treeLetter);
  }
  else
  {
    TriggerStore(event,treeLetter);
  
    if ( fTriggerStore != 0x0 ) 
    {
      fTriggerStore->Print();
    }
  }
}

//_____________________________________________________________________________
void
AliMUONDataInterface::NtupleTrigger(const char* treeLetter)
{
  //// Loop over events to build trigger ntuples
  ///
  
  TString sTreeLetter(treeLetter);
  sTreeLetter.ToUpper();
  
  if ( sTreeLetter != "R" && sTreeLetter != "D" ) 
  {
    AliError(Form("Cannot handle tree%s. Use D or R",treeLetter));
    return;
  }
  
  // book ntuples
  TNtuple tupleGlo("TgtupleGlo","Global Trigger Ntuple",
                   "ev:slpt:shpt:uplpt:uphpt:lplpt:lplpt");
  TNtuple tupleLoc("TgtupleLoc","Local Trigger Ntuple",
                   "ev:LoCircuit:LoStripX:LoDev:StripY:LoLpt:LoHpt:y11:y21:x11");
  
  // initialize counters
  Int_t sLowpt=0;
  Int_t sHighpt=0;
  Int_t uSLowpt=0;
  Int_t uSHighpt=0;
  Int_t lSLowpt=0;
  Int_t lSHighpt=0;
  
  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData(Form("%s/geometry.root",
                                    gSystem->DirName(fLoader->GetRunLoader()->GetFileName())));
    AliMUONTriggerCircuit triggerCircuit(&transformer);

  // select output file name from selected Tree
  Char_t fileNameOut[30];
  if (sTreeLetter == "D") 
  {
    AliInfo(Form("reading from Digits\n"));
    sprintf(fileNameOut,"TriggerCheckFromDigits.root");
  } 
  else if (sTreeLetter == "R") 
  {
    AliInfo(Form("reading from RecPoints\n"));
    sprintf(fileNameOut,"TriggerCheckFromRP.root");
  }
  
  // loop on events
  Int_t nevents = NumberOfEvents();

  for (Int_t ievent=0; ievent<nevents; ++ievent) 
  {
    if (ievent%100==0) AliInfo(Form("Processing event %d\n",ievent));
    
    AliMUONVTriggerStore* triggerStore = TriggerStore(ievent,treeLetter);
    
    if (!triggerStore)
    {
      AliError(Form("Could not read %s from tree%s","Trigger",treeLetter));
      return;
    }
    
    // get global trigger info
    AliMUONGlobalTrigger* gloTrg = triggerStore->Global();	
    sLowpt+=gloTrg->SingleLpt();
    sHighpt+=gloTrg->SingleHpt();
    uSLowpt+=gloTrg->PairUnlikeLpt(); 
    uSHighpt+=gloTrg->PairUnlikeHpt();
    lSLowpt+=gloTrg->PairLikeLpt(); 
    lSHighpt+=gloTrg->PairLikeHpt();
    
    // loop on local triggers	
    TIter next(triggerStore->CreateIterator());
    AliMUONLocalTrigger* locTrg(0x0);
    while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(next()) ) )
    {
      Bool_t xTrig=locTrg->IsTrigX();
      Bool_t yTrig=locTrg->IsTrigY();
      
      if (xTrig && yTrig) 
      { // fill ntuple if trigger in X and Y		        
        tupleLoc.Fill(ievent,locTrg->LoCircuit(),
                       locTrg->LoStripX(),
                       locTrg->LoDev(),
                       locTrg->LoStripY(),
                       locTrg->LoLpt(),
                       locTrg->LoHpt(),
                       triggerCircuit.GetY11Pos(locTrg->LoCircuit(),locTrg->LoStripX()),
                       triggerCircuit.GetY21Pos(locTrg->LoCircuit(),locTrg->LoStripX()+locTrg->LoDev()+1),
                       triggerCircuit.GetX11Pos(locTrg->LoCircuit(),locTrg->LoStripY()));
      }
      tupleGlo.Fill(ievent,gloTrg->SingleLpt(),gloTrg->SingleHpt(),
                     gloTrg->PairUnlikeLpt(),gloTrg->PairUnlikeHpt(),
                     gloTrg->PairLikeLpt(),gloTrg->PairLikeHpt());
    } // end of loop on local triggers
  } // end of loop on events
  
  // print info and store ntuples
  printf("\n");
  printf("=============================================\n");
  printf("================  SUMMARY  ==================\n");
  printf("\n");
  printf("Total number of events processed %d \n",nevents);
  printf("\n");
  printf(" Global Trigger output       Low pt  High pt\n");
  printf(" number of Single           :\t");
  printf("%i\t%i\t",sLowpt,sHighpt);
  printf("\n");
  printf(" number of UnlikeSign pair  :\t"); 
  printf("%i\t%i\t",uSLowpt,uSHighpt);
  printf("\n");
  printf(" number of LikeSign pair    :\t");  
  printf("%i\t%i\t",lSLowpt,lSHighpt);
  printf("\n");
  printf("=============================================\n");
  fflush(stdout);    
  
  TFile myFile(fileNameOut, "RECREATE");
  tupleGlo.Write();
  tupleLoc.Write();
  myFile.Close();
}

//_____________________________________________________________________________
Bool_t
AliMUONDataInterface::LoadEvent(Int_t event)
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

//______________________________________________________________________________
Int_t
AliMUONDataInterface::NumberOfEvents() const
{
  /// Number of events in the current galice.root file we're attached to 
  if (not IsValid()) return -1;
  return fLoader->GetRunLoader()->GetNumberOfEvents();
}

//_____________________________________________________________________________
void
AliMUONDataInterface::Open(const char* filename)
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
  
  while (AliRunLoader::GetRunLoader(foldername) != 0x0) 
  {
    delete AliRunLoader::GetRunLoader(foldername);
  }
  
  AliRunLoader* runLoader = AliRunLoader::Open(filename,foldername);
  if (runLoader == 0x0) 
  {
    AliError(Form("Cannot open file %s",filename));    
    fIsValid = kFALSE;
  }

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
Bool_t AliMUONDataInterface::GetEvent(Int_t event)
{
/// Loads all reconstructed data for the given event.

  if (DigitStore(event) == 0x0) return kFALSE;
  if (ClusterStore(event) == 0x0) return kFALSE;
  if (TriggerStore(event) == 0x0) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliMUONDataInterface::NumberOfDigits(Int_t detElemId)
{
/// Returns the number of digits to be found on a given detector element.
/// @param detElemId  The detector element ID number to search on.

  TIterator* iter = GetIterator(kDigitIteratorByDetectorElement, detElemId);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVDigit* AliMUONDataInterface::Digit(Int_t detElemId, Int_t index)
{
/// Returns the a pointer to the index'th digit on the specified detector element.
/// @param detElemId  The detector element ID number to search on.
/// @param index  The index number of the digit to fetch in the range [0 .. N-1],
///   where N = NumberOfDigits(detElemId)

  TIterator* iter = GetIterator(kDigitIteratorByDetectorElement, detElemId);
  return static_cast<AliMUONVDigit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONDataInterface::NumberOfDigits(Int_t chamber, Int_t cathode)
{
/// Returns the number of digits to be found on a specific chamber and cathode.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param cathode  The cathode in the range [0 .. 1], where 0 is the bending and
///   1 is the non-bending plane.

  TIterator* iter = GetIterator(kDigitIteratorByChamberAndCathode, chamber, cathode);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVDigit* AliMUONDataInterface::Digit(Int_t chamber, Int_t cathode, Int_t index)
{
/// Returns the a pointer to the index'th digit on the specified chamber and cathode.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param cathode  The cathode in the range [0 .. 1], where 0 is the bending and
///   1 is the non-bending plane.
/// @param index  The index number of the digit to fetch in the range [0 .. N-1],
///   where N = NumberOfDigits(chamber, cathode)

  TIterator* iter = GetIterator(kDigitIteratorByChamberAndCathode, chamber, cathode);
  return static_cast<AliMUONVDigit*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONDataInterface::NumberOfRawClusters(Int_t chamber)
{
/// Returns the number of reconstructed raw clusters on the specified chamber.
/// @param chamber  The chamber number in the range [0 .. 13].

  TIterator* iter = GetIterator(kRawClusterIterator, chamber);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONDataInterface::RawCluster(Int_t chamber, Int_t index)
{
/// Returns a pointer to the index'th raw cluster on the specified chamber.
/// @param chamber  The chamber number in the range [0 .. 13].
/// @param index  The index number of the raw cluster to fetch in the range [0 .. N-1],
///   where N = NumberOfRawClusters(chamber)

  TIterator* iter = GetIterator(kRawClusterIterator, chamber);
  return static_cast<AliMUONVCluster*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONDataInterface::NumberOfLocalTriggers()
{
/// Returns the number of reconstructed local trigger objects.

  TIterator* iter = GetIterator(kLocalTriggerIterator);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONLocalTrigger* AliMUONDataInterface::LocalTrigger(Int_t index)
{
/// Returns a pointer to the index'th local trigger object.
/// @param index  The index number of the local trigger object to fetch in the range [0 .. N-1],
///   where N = NumberOfLocalTriggers()

  TIterator* iter = GetIterator(kLocalTriggerIterator);
  return static_cast<AliMUONLocalTrigger*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
Int_t AliMUONDataInterface::NumberOfRegionalTriggers()
{
/// Returns the number of regional trigger objects reconstructed.

  TIterator* iter = GetIterator(kRegionalTriggerIterator);
  return CountObjects(iter);
}

//_____________________________________________________________________________
AliMUONRegionalTrigger* AliMUONDataInterface::RegionalTrigger(Int_t index)
{
/// Returns a pointer to the index'th regional trigger object.
/// @param index  The index number of the regional trigger object to fetch in the range [0 .. N-1],
///   where N = NumberOfRegionalTriggers()

  TIterator* iter = GetIterator(kRegionalTriggerIterator);
  return static_cast<AliMUONRegionalTrigger*>( FetchObject(iter, index) );
}

//_____________________________________________________________________________
AliMUONGlobalTrigger* AliMUONDataInterface::GlobalTrigger()
{
/// Returns a pointer to the reconstructed global trigger object for the event.

  AliMUONVTriggerStore* store = TriggerStore(fCurrentEvent);
  if (store == 0x0) return 0x0;
  return store->Global();
}

//_____________________________________________________________________________
void AliMUONDataInterface::ResetStores()
{
/// Deletes all the store objects that have been created and resets the pointers to 0x0.
/// The temporary iterator object is automatically reset. See ResetIterator for more details.

  ResetIterator();
  if (fDigitStore != 0x0)
  {
    delete fDigitStore;
    fDigitStore = 0x0;
  }
  if (fTriggerStore != 0x0)
  {
    delete fTriggerStore;
    fTriggerStore = 0x0;
  }
  if (fClusterStore != 0x0)
  {
    delete fClusterStore;
    fClusterStore = 0x0;
  }
}

//_____________________________________________________________________________
TIterator* AliMUONDataInterface::GetIterator(IteratorType type, Int_t x, Int_t y)
{
/// Creates an appropriate iterator object and returns it.
/// If the iterator has already been created then that one is returned otherwise
/// a new object is created.
/// Depending on the value of 'type' the semantics of parameters x and y can change.
/// @param type  The type of iterator to create.
/// @param x  This is the detector element ID if type == kDigitIteratorByDetectorElement
///           If type equals kDigitIteratorByChamberAndCathode or kRawClusterIterator
///           then this is the chamber number. In all other cases this parameter is
///           ignored.
/// @param y  If type == kDigitIteratorByChamberAndCathode then this parameter is the
///           cathode number. In all other cases this parameter is
///           ignored.

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
    
  case kRawClusterIterator:
    {
      Int_t chamber = x;
      AliMUONVClusterStore* store = ClusterStore(fCurrentEvent);
      if (store == 0x0) return 0x0;
      fIterator = store->CreateChamberIterator(chamber, chamber);
      if (fIterator == 0x0) return 0x0;
      fCurrentIteratorType = kRawClusterIterator;
      fDataX = chamber;
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
void AliMUONDataInterface::ResetIterator()
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
Int_t AliMUONDataInterface::CountObjects(TIterator* iter)
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
TObject* AliMUONDataInterface::FetchObject(TIterator* iter, Int_t index)
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
