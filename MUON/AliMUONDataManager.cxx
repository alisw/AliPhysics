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

/// \class AliMUONDataManager
///
/// Utility class to ease access to data stores.
/// This one is not really meant to be used directly, but more as 
/// an implementation helper for other classes (see AliMUONMCDataInterface
/// and AliMUONDataInterface)
/// 
/// Note that in the methods, trees are labelled by a single letter (except
/// for trackrefs) :
/// - K for Kine
/// - H for Hits
/// - S for SDigits
/// - D for Digits / Trigger
/// - R for RecPoints / Trigger
/// - T for Tracks / TriggerTracks
/// - TR for TrackRefs
///
/// \author Laurent Aphecetche

#include "AliMUONDataManager.h"

#include "AliLoader.h"
#include "AliLog.h"
#include "AliMUONVStore.h"
#include "AliRunLoader.h"

#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONDataManager)
/// \endcond

Int_t AliMUONDataManager::fgCount(0);

//_____________________________________________________________________________
AliMUONDataManager::AliMUONDataManager(const char* file)
: TObject(),
  fLoader(0x0),
  fCurrentEvent(-1),
  fIsValid(kTRUE)
{
    /// ctor
    ///

  TString foldername(Form("AliMUONDataManager-%d",fgCount));
  
    ++fgCount;
    
  while (AliRunLoader::GetRunLoader(foldername)) 
  {
    delete AliRunLoader::GetRunLoader(foldername);
  }
    
  AliRunLoader* runLoader = AliRunLoader::Open(file,foldername);
  if (!runLoader) 
  {
    AliError(Form("Cannot open file %s",file));    
    fIsValid = kFALSE;
  }
  fLoader = runLoader->GetDetectorLoader("MUON");
  if (!fLoader) 
  {
    AliError("Cannot get AliMUONLoader");
    fIsValid = kFALSE;
  }
}

//_____________________________________________________________________________
AliMUONDataManager::~AliMUONDataManager()
{
  /// dtor. We kill our runLoader to allow other ones to be created.
  if ( fLoader ) 
  {
    delete fLoader->GetRunLoader();
  }
}

//_____________________________________________________________________________
Int_t
AliMUONDataManager::Load(Int_t event)
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
void 
AliMUONDataManager::Load(const char* tree)
{
  /// Load a tree, labelled using an abbreviation
  
  TString s(tree);
  s.ToUpper();

  const char* opt = "READ";
  
  if ( s == "K" ) fLoader->GetRunLoader()->LoadKinematics(opt);
  if ( s == "H" ) fLoader->LoadHits(opt);
  if ( s == "S" ) fLoader->LoadSDigits(opt);
  if ( s == "T" ) fLoader->LoadTracks(opt);
  if ( s == "D" ) fLoader->LoadDigits(opt);
  if ( s == "R" ) fLoader->LoadRecPoints(opt);
  if ( s == "TR") fLoader->GetRunLoader()->LoadTrackRefs(opt);
}

//_____________________________________________________________________________
Int_t
AliMUONDataManager::NumberOfEvents() const
{
  /// Get the number of events in the file we're connected to.
  if (!IsValid()) return 0;
  return fLoader->GetRunLoader()->GetNumberOfEvents();
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONDataManager::ReadConnectable(Int_t event, const char* tree, const char* what)
{
  /// Read an object from tree, for a given event.
  
  if ( !IsValid() ) 
  {
    AliError("Object is not valid");
    return 0x0;
  }
  
  Int_t rv = Load(event);
  
  if (rv)
  {
    AliError(Form("Cannot read event %d",event));
    return 0x0;
  }
  
  Load(tree);
  
  TTree* theTree = Tree(tree);
  
  if ( theTree ) 
  {
    AliMUONVStore* object = AliMUONVStore::Create(*theTree,what);
    if ( object ) 
    {
      Bool_t ok = object->Connect(*theTree);
      if (ok)
      {
        theTree->GetEvent(0);

      }
      else
      {
        AliError(Form("Could not connect %s to tree%s",what,tree));
        delete object;
        object = 0x0;
      }
    }
    else
    {
      AliError(Form("Could not get %s from tree%s",what,tree));
    }
    Unload(tree);
    return object;
  }
  else
  {
    AliError(Form("Could not get tree%s",tree));
  }
  return 0x0;  
}

//_____________________________________________________________________________
TTree*
AliMUONDataManager::Tree(const char* tree)
{
  /// Return a tree, labelled by a letter.
  TString s(tree);
  s.ToUpper();

  if ( s == "H" ) return fLoader->TreeH();
  if ( s == "T" ) return fLoader->TreeT();
  if ( s == "S" ) return fLoader->TreeS();
  if ( s == "D" ) return fLoader->TreeD();
  if ( s == "R" ) return fLoader->TreeR();
  if ( s == "TR") return fLoader->GetRunLoader()->TreeTR();
  
  return 0x0;
}

//_____________________________________________________________________________
void 
AliMUONDataManager::Unload(const char* tree)
{
  /// Unload a given tree.
  TString s(tree);
  s.ToUpper();

  if ( s == "K" ) fLoader->GetRunLoader()->UnloadKinematics();
  if ( s == "H" ) fLoader->UnloadHits();
  if ( s == "S" ) fLoader->UnloadSDigits();
  if ( s == "T" ) fLoader->UnloadTracks();
  if ( s == "D" ) fLoader->UnloadDigits();
  if ( s == "R" ) fLoader->UnloadRecPoints();
  if ( s == "TR" ) fLoader->GetRunLoader()->UnloadTrackRefs();
}
