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
/// \class AliMUONTriggerStoreV1
///
/// Implementation of AliMUONVTriggerStore, which is backward compatible,
/// i.e. should be able to read back old TreeR and TreeD files, produced
/// before the introduction of the AliMUONVStore concept.
/// 
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONTriggerStoreV1.h"

#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONTreeManager.h"
#include <Riostream.h>
#include <TClonesArray.h>
#include <TTree.h>

using std::cout;
using std::endl;
/// \cond CLASSIMP
ClassImp(AliMUONTriggerStoreV1)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerStoreV1::AliMUONTriggerStoreV1() : AliMUONVTriggerStore(),
fLocal(new TClonesArray("AliMUONLocalTrigger",234)),
fRegional(new TClonesArray("AliMUONRegionalTrigger",16)),
fGlobal(new TClonesArray("AliMUONGlobalTrigger",1)),
fEmptyLocal(new TClonesArray("AliMUONLocalTrigger",234))
{
  /// ctor
  fLocal->SetOwner(kTRUE);
  fRegional->SetOwner(kTRUE);
  fGlobal->SetOwner(kTRUE);
  fEmptyLocal->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONTriggerStoreV1::~AliMUONTriggerStoreV1()
{
  /// dtor
  delete fLocal;
  delete fRegional;
  delete fGlobal;
  delete fEmptyLocal;
}

//_____________________________________________________________________________
void 
AliMUONTriggerStoreV1::Add(const AliMUONLocalTrigger& localTrigger)
{
  /// Add local information
  /// If the local board has no information (IsNull), we
  /// add it in the fEmpty array
  /// This is really an implementation choice, to store empty boards
  /// in order to be able to return them, if asked for, as is the case
  /// in some client code. Note that only the non-empty boards
  /// are streamed to disk.
  ///
  
  if ( !localTrigger.IsNull() ) 
  {
    new((*fLocal)[fLocal->GetLast()+1]) AliMUONLocalTrigger(localTrigger);

  }
  else
  {
    new((*fEmptyLocal)[fEmptyLocal->GetLast()+1]) AliMUONLocalTrigger(localTrigger);
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONTriggerStoreV1::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this store to the tree
  AliMUONTreeManager tman;
  Bool_t ok(kTRUE);
  
  Bool_t isMaking = ( tree.GetBranch("MUONLocalTrigger") == 0 );
  
  if ( isMaking ) 
  {
    ok = ok && tman.MakeBranch(tree,ClassName(),"TClonesArray",
                               "MUONLocalTrigger",LocalPtr());
    ok = ok && tman.MakeBranch(tree,ClassName(),"TClonesArray",
                               "MUONRegionalTrigger",RegionalPtr());
    ok = ok && tman.MakeBranch(tree,ClassName(),"TClonesArray",
                               "MUONGlobalTrigger",GlobalPtr());
  }
  else
  {
    if ( alone ) tman.UpdateBranchStatuses(tree,"Trigger");
    ok = ok && tman.SetAddress(tree,"MUONLocalTrigger",LocalPtr());
    ok = ok && tman.SetAddress(tree,"MUONRegionalTrigger",RegionalPtr());
    ok = ok && tman.SetAddress(tree,"MUONGlobalTrigger",GlobalPtr());
  }
  return ok;
}

//_____________________________________________________________________________
void 
AliMUONTriggerStoreV1::SetGlobal(const AliMUONGlobalTrigger& globalTrigger)
{
  /// Set the global information
  new((*fGlobal)[0]) AliMUONGlobalTrigger(globalTrigger);
}

//_____________________________________________________________________________
void 
AliMUONTriggerStoreV1::Add(const AliMUONRegionalTrigger& regionalTrigger)
{
  /// Add regional information
  new((*fRegional)[fRegional->GetLast()+1]) AliMUONRegionalTrigger(regionalTrigger);
}

//_____________________________________________________________________________
TIterator* 
AliMUONTriggerStoreV1::CreateLocalIterator() const
{
  /// Return iterator on local cards
  return fLocal->MakeIterator();
}

//_____________________________________________________________________________
TIterator*
AliMUONTriggerStoreV1::CreateRegionalIterator() const
{
  /// Return iterator on regional cards
  return fRegional->MakeIterator();
}

//_____________________________________________________________________________
AliMUONLocalTrigger* 
AliMUONTriggerStoreV1::FindLocal(Int_t boardNumber) const
{
  /// Find a local board, by its *number* (not to be confused with its index,
  /// which used to be the key)
  ///
  
  for ( Int_t i = 0; i <= fLocal->GetLast(); ++i ) 
  {
    AliMUONLocalTrigger* local = static_cast<AliMUONLocalTrigger*>(fLocal->At(i));
    if (local && local->LoCircuit()==boardNumber)
    {
      return local;
    }
  }
  
  for ( Int_t i = 0; i <= fEmptyLocal->GetLast(); ++i ) 
  {
    AliMUONLocalTrigger* local = static_cast<AliMUONLocalTrigger*>(fEmptyLocal->At(i));
    if (local && local->LoCircuit()==boardNumber)
    {
      return local;
    }
  }
  
  if ( boardNumber>=1 && boardNumber<=234 ) 
  {
    AliMUONLocalTrigger empty;
    empty.SetLoCircuit(boardNumber);
    new((*fEmptyLocal)[fEmptyLocal->GetLast()+1]) AliMUONLocalTrigger(empty);
    return FindLocal(boardNumber);
  }
  
  return 0x0;
}

//_____________________________________________________________________________
AliMUONRegionalTrigger* 
AliMUONTriggerStoreV1::FindRegional(Int_t boardNumber) const
{
  /// Return a given regional board
  for ( Int_t i = 0; i <= fRegional->GetLast(); ++i ) 
  {
    AliMUONRegionalTrigger* regional = static_cast<AliMUONRegionalTrigger*>(fRegional->At(i));
    if (regional && regional->GetId()==boardNumber)
    {
      return regional;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONGlobalTrigger*
AliMUONTriggerStoreV1::Global() const
{
  /// Return global trigger
  return static_cast<AliMUONGlobalTrigger*>(fGlobal->At(0));
}

//_____________________________________________________________________________
void
AliMUONTriggerStoreV1::Clear(Option_t*)
{
  /// Reset
  fLocal->Clear("C");
  fRegional->Clear("C");
  fGlobal->Clear("C");
  fEmptyLocal->Clear("C");
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerStoreV1::GetSize() const
{
  /// Number of non-empty local boards we hold
  return fLocal->GetSize();
}

//_____________________________________________________________________________
void
AliMUONTriggerStoreV1::Print(Option_t* what, Option_t* opt) const
{
  /// Printout
  /// \param what used to tell what to print, can be GLOBAL, LOCAL, REGIONAL
  /// or ALL
  /// \param opt is passed to the local, regional, global object
  ///
  
  TString swhat(what);
  swhat.ToUpper();
  
  if ( swhat.Length() == 0 ) swhat = "ALL";
  
  if ( swhat.Contains("GLOBAL") || swhat.Contains("ALL") )
  {
    if ( fGlobal ) 
    {
      cout << "Global:" << endl;
      fGlobal->Print("",opt);
    }
    else 
    {
      cout << "No GlobalTrigger information" << endl;
    }
  }
  
  if ( fLocal && ( swhat.Contains("LOCAL")|| swhat.Contains("ALL") ) ) 
  {    
    // make loops instead of just relying on fLocal
    // to insure backward compatibility with trees where all local boards where
    // stored (even null ones)
        
    TIter next(fLocal);
    AliMUONLocalTrigger* local;
    Int_t n(0);
    
    while ( ( local = static_cast<AliMUONLocalTrigger*>(next()) ) )
    {
      if ( local->IsNull() ) ++n;
    }

    cout << Form("Local: %d cards (and %d null ones)",
                 fLocal->GetLast()+1,n) << endl;
    
    next.Reset();
    
    while ( ( local = static_cast<AliMUONLocalTrigger*>(next()) ) )
    {
      if ( !local->IsNull() ) 
      {
        local->Print(opt);
      }
    }
  }
  
  if ( fRegional && ( swhat.Contains("REGIONAL") || swhat.Contains("ALL") ) )
  {
    fRegional->Print("",opt);
  }
}

