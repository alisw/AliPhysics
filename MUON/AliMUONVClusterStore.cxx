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
/// \class AliMUONVClusterStore
///
/// An interface of a cluster container
///
/// The object stored are inherited from AliMUONVCluster
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVClusterStore.h"
#include "AliMUONVCluster.h"

/// \cond CLASSIMP
ClassImp(AliMUONVClusterStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVClusterStore::AliMUONVClusterStore()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVClusterStore::~AliMUONVClusterStore()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t
AliMUONVClusterStore::Add(TObject* object)
{
  /// Add an object, if it is of the right class
  AliMUONVCluster* cluster = dynamic_cast<AliMUONVCluster*>(object);
  if (cluster)
  {
    return (Add(*cluster)) ? kTRUE : kFALSE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMUONVClusterStore*
AliMUONVClusterStore::Create(TTree& tree)
{
  /// Create a VClusterStore from the tree
  return static_cast<AliMUONVClusterStore*>(AliMUONVStore::Create(tree,"Cluster"));
}

//______________________________________________________________________________
AliMUONVCluster* AliMUONVClusterStore::FindObject(const TObject *obj) const
{
  /// Find an object, if of AliMUONVCluster type
  const AliMUONVCluster* cluster = dynamic_cast<const AliMUONVCluster*>(obj);
  if (cluster)
    return static_cast<AliMUONVCluster*>(AliMUONVStore::FindObject(obj));
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONVClusterStore::FindObject(UInt_t uniqueID) const
{
  /// Find an object by its uniqueID (default is the same as in AliMUONVStore)
  return static_cast<AliMUONVCluster*>(AliMUONVStore::FindObject(uniqueID));
}

