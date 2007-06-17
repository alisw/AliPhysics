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

#include <TError.h>
#include <TParticle.h>

#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliMUONDataInterface.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliLog.h"

#include <Riostream.h>

#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONDataManager.h"
#include "AliLog.h"
#include <TList.h>

///
/// \class AliMUONDataInterface
///
/// An easy to use interface to the MUON data data stored in
/// TreeS, TreeD, TreeR and TreeT.
///
/// For MC related information (i.e. TreeH, TreeK, TreeTR), see
/// AliMUONMCDataInterface.
///
///
/// This interface in not necessarily the fastest way to fetch the data but
/// it is the easiest.
///
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONDataInterface)
/// \endcond

//______________________________________________________________________________
AliMUONDataInterface::AliMUONDataInterface(const char* filename)
	: TObject(), 
fDataManager(new AliMUONDataManager(filename))
{
  /// ctor
  /// @param filename should be the full path to a valid galice.root file
  
  if (!IsValid())
  {
    AliError("Improper initialization. Object will be unuseable");
  }
}

//______________________________________________________________________________
AliMUONDataInterface::~AliMUONDataInterface()
{
  /// dtor
  delete fDataManager;
}

//______________________________________________________________________________
AliMUONVClusterStore*
AliMUONDataInterface::ClusterStore(Int_t event) const
{
  /// Return the cluster store for a given event
  if (!IsValid()) return 0x0;
  return dynamic_cast<AliMUONVClusterStore*>(fDataManager->ReadConnectable(event,"R","Cluster"));

}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpRecPoints(Int_t event, Bool_t sorted) const
{
  /// Dump the recpoints for a given event, sorted if so required
  DumpIt("R","Cluster",event,sorted);
}

//______________________________________________________________________________
AliMUONVDigitStore*
AliMUONDataInterface::DigitStore(Int_t event) const
{
  /// Return the digit store for a given event
  if (!IsValid()) return 0x0;
  return dynamic_cast<AliMUONVDigitStore*>(fDataManager->ReadConnectable(event,"D","Digit"));
}

//______________________________________________________________________________
TList* 
AliMUONDataInterface::DigitStoreAsList(Int_t event) const
{
  /// Return the digitStore as a TList
  AliMUONVDigitStore* digitStore = DigitStore(event);
  
  TIter next(digitStore->CreateIterator());

  TList* list = new TList;
  list->SetOwner(kTRUE);
  TObject* object;
  
  while ( ( object = next() ) ) 
  {
    list->Add(object->Clone());
  }
  
  delete digitStore;
  return list;
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpIt(const char* treeLetter, const char* what, 
                             Int_t event, Bool_t sorted) const
{
  /// Generic dump method used by the other DumpXXX methods
  AliMUONVStore* store = fDataManager->ReadConnectable(event,treeLetter,what);
  if (!store)
  {
    AliError(Form("Could not read %s from tree%s",what,treeLetter));
    return;
  }
  
  if ( sorted ) 
  {
    TList list;
    list.SetOwner(kFALSE);
    TIter next(store->CreateIterator());
    TObject* object;
  
    while ( ( object = next() ) ) 
    {
      list.Add(object);
    }

    list.Sort();
  
    list.Print();
  }
  else
  {
    store->Print();
  }
  
  delete store;
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpDigits(Int_t event, Bool_t sorted) const
{
  /// Dump digits of a given event, sorted if so required
  DumpIt("D","Digit",event,sorted);
}

//______________________________________________________________________________
Bool_t
AliMUONDataInterface::IsValid() const
{
  /// Whether we were properly initialized from a valid galice.root file
  return fDataManager->IsValid();
}

//______________________________________________________________________________
Int_t
AliMUONDataInterface::NumberOfEvents() const
{
  /// Number of events in the current galice.root file we're attached to 
  if (!IsValid()) return 0;
  return fDataManager->NumberOfEvents();
}

//______________________________________________________________________________
AliMUONVDigitStore*
AliMUONDataInterface::SDigitStore(Int_t event) const
{
  /// Return the SDigit store for a given event
  if (!IsValid()) return 0x0;
  return dynamic_cast<AliMUONVDigitStore*>(fDataManager->ReadConnectable(event,"S","SDigit"));
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpSDigits(Int_t event, Bool_t sorted) const
{
  /// Dump sdigits for a given event, sorted if so required
  DumpIt("S","Digit",event,sorted);
}


//______________________________________________________________________________
AliMUONVTrackStore* 
AliMUONDataInterface::TrackStore(Int_t event) const
{
  /// Return the track store for a given event
  if (!IsValid()) return 0x0;
  return dynamic_cast<AliMUONVTrackStore*>(fDataManager->ReadConnectable(event,"T","Track"));
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpTracks(Int_t event, Bool_t sorted) const
{
  /// Dump tracks for a given event, sorted if so required
  DumpIt("T","Track",event,sorted);
}

//______________________________________________________________________________
AliMUONVTriggerStore* 
AliMUONDataInterface::TriggerStore(Int_t event, const char* treeLetter) const
{
  /// Return the trigger store for a given event, from a given tree (D or R)
  if (!IsValid()) return 0x0;
  return dynamic_cast<AliMUONVTriggerStore*>(fDataManager->ReadConnectable(event,treeLetter,"Trigger"));  
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpTrigger(Int_t event, const char* treeLetter) const
{
  /// Dump trigger for a given event, from a given tree, sorted if possible
  DumpIt(treeLetter,"Trigger",event,kFALSE);
}

//______________________________________________________________________________
AliMUONVTriggerTrackStore* 
AliMUONDataInterface::TriggerTrackStore(Int_t event) const
{
  /// Return trigger track store for a given event
  if (!IsValid()) return 0x0;
  return dynamic_cast<AliMUONVTriggerTrackStore*>(fDataManager->ReadConnectable(event,"T","TriggerTrack"));  
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpTriggerTracks(Int_t event, Bool_t sorted) const
{
  /// Dump trigger tracks for a given event
  DumpIt("T","Trigger",event,sorted);
}

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

void AliMUONDataInterface::Reset()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
}

Bool_t AliMUONDataInterface::UseCurrentRunLoader()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return kFALSE;
}
  
Int_t AliMUONDataInterface::NumberOfEvents(TString , TString )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfParticles(TString , TString , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


TParticle* AliMUONDataInterface::Particle(
		TString , TString , Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfTracks(TString , TString , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfHits(
		TString , TString , Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONHit* AliMUONDataInterface::Hit(
		TString , TString , Int_t ,
		Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfSDigits(
		TString , TString , Int_t ,
		Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONDigit* AliMUONDataInterface::SDigit(
		TString , TString , Int_t ,
		Int_t , Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfDigits(
		TString , TString , Int_t ,
		Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONDigit* AliMUONDataInterface::Digit(
		TString , TString , Int_t ,
		Int_t , Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfRawClusters(
		TString , TString , Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONRawCluster* AliMUONDataInterface::RawCluster(
		TString , TString , Int_t ,
		Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfLocalTriggers(TString , TString , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONLocalTrigger* AliMUONDataInterface::LocalTrigger(
		TString , TString , Int_t , Int_t 
	)
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}

Bool_t AliMUONDataInterface::SetFile(TString , TString )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Bool_t AliMUONDataInterface::GetEvent(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}

Int_t AliMUONDataInterface::NumberOfParticles()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


TParticle* AliMUONDataInterface::Particle(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfTracks()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfHits(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONHit* 
AliMUONDataInterface::Hit(Int_t , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfSDigits(Int_t , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONDigit* AliMUONDataInterface::SDigit(Int_t , Int_t , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;

}


Int_t AliMUONDataInterface::NumberOfDigits(Int_t , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;

}


AliMUONDigit* AliMUONDataInterface::Digit(Int_t , Int_t , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfRawClusters(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONRawCluster* AliMUONDataInterface::RawCluster(Int_t , Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


Int_t AliMUONDataInterface::NumberOfLocalTriggers()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}


AliMUONLocalTrigger* AliMUONDataInterface::LocalTrigger(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}

Int_t AliMUONDataInterface::NumberOfGlobalTriggers()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}

AliMUONGlobalTrigger* AliMUONDataInterface::GlobalTrigger(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}

Int_t AliMUONDataInterface::NumberOfRecTracks()
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}

AliMUONTrack* AliMUONDataInterface::RecTrack(Int_t )
{
/// \deprecated Method is going to be removed

  AliFatal("Deprecated");
  return 0;
}
