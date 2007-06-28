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
#include "AliMUONDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONHit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"

#include "AliLoader.h"
#include "AliLog.h"
#include "AliLog.h"
#include "AliRunLoader.h"

#include <TError.h>
#include <TParticle.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TNtuple.h>
#include <TSystem.h>

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

//AliLoader* fLoader; //!< Tree accessor
//AliMUONVDigitStore* fDigitStore; //!< current digit store (owner)
//AliMUONVTriggerStore* fTriggerStore; //!< current trigger store (owner)
//AliMUONVClusterStore* fClusterStore; //!< current cluster store (owner)
//AliMUONVTrackStore* fTrackStore; //!< current track store (owner)
//AliMUONVTriggerTrackStore* fTriggerTrackStore; //!< current trigger track store (owner)
//Int_t fCurrentEvent; //!< Current event we've read in
//Bool_t fIsValid; //!< whether we were initialized properly or not

Int_t AliMUONDataInterface::fgInstanceCounter(0);

//______________________________________________________________________________
AliMUONDataInterface::AliMUONDataInterface(const char* filename)
: TObject(), 
fLoader(0x0),
fDigitStore(0x0),
fTriggerStore(0x0),
fClusterStore(0x0),
fTrackStore(0x0),
fTriggerTrackStore(0x0),
fCurrentEvent(-1),
fIsValid(kFALSE)
{
  /// ctor
  /// @param filename should be the full path to a valid galice.root file
  
  ++fgInstanceCounter;
  
  Open(filename);
}

//______________________________________________________________________________
AliMUONDataInterface::~AliMUONDataInterface()
{
  /// dtor
  if ( fLoader ) 
  {
    delete fLoader->GetRunLoader();
  }
  --fgInstanceCounter;  
}

//______________________________________________________________________________
AliMUONVClusterStore*
AliMUONDataInterface::ClusterStore(Int_t event)
{
  /// Return clusterStore for a given event.
  /// Return 0x0 if event not found.
  /// Returned pointer should not be deleted
  
  if ( LoadEvent(event) ) return 0x0;
  
  fLoader->LoadRecPoints();
  
  TTree* treeR = fLoader->TreeR();
  
  if (!treeR)
  {
    AliError("Could not get treeR");
    return 0x0;
  }
  
  if (!fClusterStore)
  {
    fClusterStore = AliMUONVClusterStore::Create(*treeR);
  }
  
  if ( fClusterStore ) 
  {
    fClusterStore->Clear();
    fClusterStore->Connect(*treeR);
    treeR->GetEvent(0);
  }
  
  fLoader->UnloadRecPoints();
  
  return fClusterStore;
}

//______________________________________________________________________________
AliMUONVDigitStore*
AliMUONDataInterface::DigitStore(Int_t event)
{
  /// Return digitStore for a given event.
  /// Return 0x0 if event not found.
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

//______________________________________________________________________________
void
AliMUONDataInterface::DumpDigits(Int_t event, Bool_t sorted)
{
  /// Dump the digits for a given event, sorted if so required
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

//______________________________________________________________________________
void
AliMUONDataInterface::DumpRecPoints(Int_t event, Bool_t sorted)
{
  /// Dump the recpoints for a given event, sorted if so required
  ClusterStore(event);
  if ( fClusterStore ) 
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

//______________________________________________________________________________
void
AliMUONDataInterface::DumpTracks(Int_t event, Bool_t sorted)
{
  /// Dump tracks for a given event, sorted if requested
  
  TrackStore(event);
  
  if ( fTrackStore ) 
  {
    if ( sorted ) 
    {
      DumpSorted(*fTrackStore);
    }
    else
    {
      fTrackStore->Print();
    }
  }
}

//______________________________________________________________________________
void
AliMUONDataInterface::DumpTriggerTracks(Int_t event, Bool_t sorted)
{
  /// Dump trigger tracks for a given event, sorted if requested

  TriggerTrackStore(event);
  
  if ( fTriggerTrackStore ) 
  {
    if ( sorted ) 
    {
      DumpSorted(*fTriggerTrackStore);
    }
    else
    {
      fTriggerTrackStore->Print();
    }
  }
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
  
    if ( fTriggerStore ) 
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
  
  AliMUONGeometryTransformer transformer(kFALSE);
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
    
    AliMUONVTriggerStore* triggerStore = TriggerStore(ievent);
    
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
      Bool_t xTrig=kFALSE;
      Bool_t yTrig=kFALSE;
      
      if ( locTrg->LoSdev()==1 && locTrg->LoDev()==0 && 
           locTrg->LoStripX()==0) xTrig=kFALSE; // no trigger in X
      else xTrig=kTRUE;                         // trigger in X
      if (locTrg->LoTrigY()==1 && 
          locTrg->LoStripY()==15 ) yTrig = kFALSE; // no trigger in Y
      else yTrig = kTRUE;                          // trigger in Y
      
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

//______________________________________________________________________________
Bool_t
AliMUONDataInterface::IsValid() const
{
  /// Whether we were properly initialized from a valid galice.root file
  return fIsValid;
}

//_____________________________________________________________________________
Int_t
AliMUONDataInterface::LoadEvent(Int_t event)
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

//______________________________________________________________________________
Int_t
AliMUONDataInterface::NumberOfEvents() const
{
  /// Number of events in the current galice.root file we're attached to 
  if (!IsValid()) return 0;
  return fLoader->GetRunLoader()->GetNumberOfEvents();
}

//_____________________________________________________________________________
void
AliMUONDataInterface::Open(const char* filename)
{
  /// Connect to a given galice.root file
  
  delete fDigitStore;
  fDigitStore=0x0;
  delete fTriggerStore;
  fTriggerStore=0x0;
  delete fClusterStore;
  fClusterStore=0x0;
  delete fTrackStore;
  fTrackStore=0x0;
  delete fTriggerTrackStore;
  fTriggerTrackStore=0x0;
  
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

//______________________________________________________________________________
AliMUONVTrackStore* 
AliMUONDataInterface::TrackStore(Int_t event)
{
  /// Return the trackStore for a given event.
  /// Return 0x0 if event not found.
  /// Returned pointer should not be deleted
  
  if ( LoadEvent(event) ) return 0x0;
  
  fLoader->LoadTracks();
  
  TTree* treeT = fLoader->TreeT();
  
  if (!treeT)
  {
    AliError("Could not get treeT");
    return 0x0;
  }
  
  if (!fTrackStore)
  {
    fTrackStore = AliMUONVTrackStore::Create(*treeT);
  }
  
  if ( fTrackStore ) 
  {
    fTrackStore->Clear();
    fTrackStore->Connect(*treeT);
    treeT->GetEvent(0);
  }
  
  fLoader->UnloadTracks();
  
  return fTrackStore;
}

//______________________________________________________________________________
AliMUONVTriggerTrackStore* 
AliMUONDataInterface::TriggerTrackStore(Int_t event)
{
  /// Return the triggerTrackStore for a given event.
  /// Return 0x0 if event not found.
  /// Returned pointer should not be deleted
  
  if ( LoadEvent(event) ) return 0x0;
  
  fLoader->LoadTracks();
  
  TTree* treeT = fLoader->TreeT();
  
  if (!treeT)
  {
    AliError("Could not get treeT");
    return 0x0;
  }
  
  if (!fTriggerTrackStore)
  {
    fTriggerTrackStore = AliMUONVTriggerTrackStore::Create(*treeT);
  }
  
  if ( fTriggerTrackStore ) 
  {
    fTriggerTrackStore->Clear();
    fTriggerTrackStore->Connect(*treeT);
    treeT->GetEvent(0);
  }
  
  fLoader->UnloadTracks();
  
  return fTriggerTrackStore;  
}

//_____________________________________________________________________________
AliMUONVTriggerStore*
AliMUONDataInterface::TriggerStore(Int_t event, const char* treeLetter)
{
  /// Return the triggerStore for a given event.
  /// Return 0x0 if event not found.
  /// Returned pointer should not be deleted
  /// treeLetter can be R or D to tell from which tree to read the information
  
  if ( LoadEvent(event) ) return 0x0;
  
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
  
  if ( !tree ) 
  {
    AliError(Form("Could not get tree%s",treeLetter));
    return 0x0;
  }
  
  if (!fTriggerStore)
  {
    fTriggerStore = AliMUONVTriggerStore::Create(*tree);
  }
  
  if ( fTriggerStore ) 
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
  
  return fTriggerStore;
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
