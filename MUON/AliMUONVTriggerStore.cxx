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
/// \class AliMUONVTriggerStore
///
/// Base class of a trigger container, which holds local, regional and
/// global information for one event.
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVTriggerStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONGlobalTrigger.h"

/// \cond CLASSIMP
ClassImp(AliMUONVTriggerStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVTriggerStore::AliMUONVTriggerStore()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVTriggerStore::~AliMUONVTriggerStore()
{
  /// dtor
}

//_____________________________________________________________________________
AliMUONVTriggerStore* 
AliMUONVTriggerStore::Create(TTree& tree)
{
  /// Create a VTriggerStore from the tree (if possible).
  return static_cast<AliMUONVTriggerStore*>(AliMUONVStore::Create(tree,"Trigger"));
}

//_____________________________________________________________________________
TIterator* 
AliMUONVTriggerStore::CreateIterator() const
{
  /// Return local iterator
  return CreateLocalIterator();
}

//_____________________________________________________________________________
Bool_t
AliMUONVTriggerStore::Add(TObject* object)
{
  /// Add an object, if object of type local, regional or global
  if (!object) return kFALSE;
  AliMUONLocalTrigger* local = dynamic_cast<AliMUONLocalTrigger*>(object);
  if (local) 
  {
    Add(*local);
    return kTRUE;
  }  
  AliMUONRegionalTrigger* regional = dynamic_cast<AliMUONRegionalTrigger*>(object);
  if (regional)
  {
    Add(*regional);
    return kTRUE;
  }
  AliMUONGlobalTrigger* global = dynamic_cast<AliMUONGlobalTrigger*>(object);
  if (global)
  {
    SetGlobal(*global);
    return kTRUE;
  }
  return kFALSE;
}
