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

#include "AliMUONHitStoreV1.h"

//-----------------------------------------------------------------------------
/// \class AliMUONHitStoreV1
///
/// Implementation of AliMUONVHitStore
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include <TClonesArray.h>
#include <TTree.h>
#include "AliMUONTreeManager.h"
#include "AliMUONHit.h"

/// \cond CLASSIMP
ClassImp(AliMUONHitStoreV1)
/// \endcond

//_____________________________________________________________________________
AliMUONHitStoreV1::AliMUONHitStoreV1(TRootIOCtor* /*dummy*/) : AliMUONVHitStore(),
fHits(0x0)
{
  /// default ctor from file
}

//_____________________________________________________________________________
AliMUONHitStoreV1::AliMUONHitStoreV1() : AliMUONVHitStore(),
 fHits(new TClonesArray("AliMUONHit",10))
{
   /// ctor
   fHits->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONHitStoreV1::~AliMUONHitStoreV1()
{
  /// dtor
  delete fHits;
}

//_____________________________________________________________________________
void 
AliMUONHitStoreV1::Add(const AliMUONHit& hit)
{
  /// add a hit
  new((*fHits)[fHits->GetLast()+1]) AliMUONHit(hit);
}

//_____________________________________________________________________________
TCollection*
AliMUONHitStoreV1::Collection()
{
  return fHits;
}

//_____________________________________________________________________________
Bool_t 
AliMUONHitStoreV1::Connect(TTree& tree, Bool_t /*alone*/) const
{
  /// Connect this to tree.
  AliMUONTreeManager tman;
  Bool_t ok;
  
  if ( tree.GetBranch("MUONHits") )
  {
    ok = tman.SetAddress(tree,"MUONHits",HitsPtr());
  }
  else
  {
    ok = tman.MakeBranch(tree,ClassName(),"TClonesArray","MUONHits",
                         HitsPtr());
  }
  return ok;
}

//_____________________________________________________________________________
TIterator*
AliMUONHitStoreV1::CreateIterator() const
{
  /// create an iterator on hits
  return fHits->MakeIterator();
}

//____________________________________________________________________________
Int_t 
AliMUONHitStoreV1::GetSize() const
{
  return fHits->GetLast()+1;
}

//____________________________________________________________________________
void 
AliMUONHitStoreV1::Clear(Option_t*)
{
  /// reset the internal array
  fHits->Clear("C");
}
