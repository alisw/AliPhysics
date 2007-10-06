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
/// \class AliMUONClusterStoreV1
///
/// Implementation of VClusterStore.
///
/// This one is a basic implementation, let's say "legacy" one, i.e.
/// compatible with what we stored in MUON.RecPoints.root files before
/// the switch to data stores.
///
/// \author Laurent Aphecetche, Subatech
///
//-----------------------------------------------------------------------------

#include "AliMUONClusterStoreV1.h"

#include "AliLog.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTOTCAStoreIterator.h"
#include "AliMUONTreeManager.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterStoreV1)
/// \endcond

//_____________________________________________________________________________
AliMUONClusterStoreV1::AliMUONClusterStoreV1() 
: AliMUONVClusterStore(), 
fClusters(new TObjArray(AliMpConstants::NofChambers()))
{
  /// ctor. Set correct ownerships
  fClusters->SetOwner(kTRUE);
  for ( Int_t i = 0; i < fClusters->GetSize(); ++i )
  {
    TClonesArray* tca = new TClonesArray("AliMUONRawCluster",100);
    tca->SetOwner(kTRUE);
    fClusters->AddAt(tca,i);
  }
  AliDebug(1,"");
}

//_____________________________________________________________________________
AliMUONClusterStoreV1::AliMUONClusterStoreV1(const AliMUONClusterStoreV1&)
: AliMUONVClusterStore(), 
fClusters(0x0)
{
  /// copy ctor
  AliError("Please implement me");
}

//_____________________________________________________________________________
AliMUONClusterStoreV1& 
AliMUONClusterStoreV1::operator=(const AliMUONClusterStoreV1&)
{
  /// assignment operator
  AliError("Please implement me");
  return *this;
}

//_____________________________________________________________________________
AliMUONClusterStoreV1::~AliMUONClusterStoreV1()
{
  /// dtor
  AliDebug(1,"");
  delete fClusters;
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV1::CreateCluster(Int_t /*chamberId*/, Int_t detElemId, Int_t /*clusterIndex*/) const
{
  /// Create a cluster
  AliMUONVCluster* vCluster = new AliMUONRawCluster();
  (static_cast<AliMUONRawCluster*> (vCluster))->SetDetElemId(detElemId);
  return vCluster;
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterStoreV1::Add(const AliMUONVCluster& vCluster)
{
  /// Add a cluster to this store
  const AliMUONRawCluster* cluster = dynamic_cast<const AliMUONRawCluster*>(&vCluster);
  
  if (!cluster)
  {
    AliError(Form("Cluster is not of the expected type (%s vs AliMUONRawCluster)",
                  vCluster.ClassName()));
    return 0x0;
  }
  
  Int_t iChamber = AliMpDEManager::GetChamberId(cluster->GetDetElemId());
  TClonesArray* array = ChamberClusters(iChamber);
  if (!array) 
  {
    return kFALSE;
  }
  new((*array)[array->GetLast()+1]) AliMUONRawCluster(*cluster);
  return kTRUE;
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV1::Add(Int_t chamberId, Int_t detElemId, Int_t /*clusterIndex*/)
{
  /// Add a cluster to this store
  TClonesArray* array = ChamberClusters(chamberId);
  if (!array) return 0x0;
  
  AliMUONVCluster* vCluster = static_cast<AliMUONVCluster*> (new((*array)[array->GetLast()+1]) AliMUONRawCluster());
  (static_cast<AliMUONRawCluster*> (vCluster))->SetDetElemId(detElemId);
  return vCluster;
}

//_____________________________________________________________________________
TClonesArray*
AliMUONClusterStoreV1::ChamberClusters(Int_t chamberId) const
{
  /// Get the internal array of clusters for a given chamber
  TClonesArray* array = static_cast<TClonesArray*>(fClusters->At(chamberId));
  if (!array) 
  {
    AliError(Form("Cannot get Clusters for chamberId=%d",chamberId));
    return 0x0;
  }
  return array;
}

//_____________________________________________________________________________
TObject**
AliMUONClusterStoreV1::ChamberClustersPtr(Int_t chamberId) const
{
  /// Get the internal array of clusters for a given chamber
  TClonesArray* array = static_cast<TClonesArray*>(fClusters->At(chamberId));
  if (!array) 
  {
    AliError(Form("Cannot get Clusters for chamberId=%d",chamberId));
    return 0x0;
  }
  return fClusters->GetObjectRef(array);
}

//_____________________________________________________________________________
Bool_t
AliMUONClusterStoreV1::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this to the tree, i.e. make the branches or set their addresses.
  
  AliMUONTreeManager tman;
  Bool_t ok(kTRUE);
  
  TBranch* b = tree.GetBranch("MUONRawClusters1");
  
  Bool_t isMaking = (b == 0);
  
  if ( isMaking ) 
  {
    for ( Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i ) 
    {
      TString branchName(Form("MUONRawClusters%d",i+1));
      ok = ok && tman.MakeBranch(tree,ClassName(),"TClonesArray",
                                 branchName.Data(),ChamberClustersPtr(i));
    }
    
  }
  else
  {
    if (alone) tman.UpdateBranchStatuses(tree,"MUONRawClusters");
    for ( Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i ) 
    {
      TString branchName(Form("MUONRawClusters%d",i+1));
      ok = ok && tman.SetAddress(tree,branchName.Data(),
                                   ChamberClustersPtr(i));
    }
  }
  return ok;
}

//_____________________________________________________________________________
AliMUONVCluster*
AliMUONClusterStoreV1::Remove(AliMUONVCluster& cluster)
{
  /// Remove a cluster
  Int_t iChamber = AliMpDEManager::GetChamberId(cluster.GetDetElemId());
  TClonesArray* array = ChamberClusters(iChamber);
  TObject* o = array->Remove(&cluster);
  if (o)
  {
    array->Compress();
  }
  return static_cast<AliMUONVCluster*>(o);
}

//_____________________________________________________________________________
void
AliMUONClusterStoreV1::Clear(Option_t*)
{
  /// Reset internal arrays
  AliDebug(1,"");
  /// Reset the tclonesarray, but keep the tobjarray's size constant.
  for ( Int_t i = 0; i < fClusters->GetSize(); ++i ) 
  {
    ChamberClusters(i)->Clear("C");
  }
}

//_____________________________________________________________________________
TIterator* 
AliMUONClusterStoreV1::CreateIterator() const
{
  /// Return an iterator to loop over our clusters
  return new AliMUONTOTCAStoreIterator(fClusters,0,AliMpConstants::NofTrackingChambers()-1);
}

//_____________________________________________________________________________
TIterator* 
AliMUONClusterStoreV1::CreateChamberIterator(Int_t firstChamber, Int_t lastChamber) const
{
  /// Return an iterator to loop over our clusters
  return new AliMUONTOTCAStoreIterator(fClusters,firstChamber,lastChamber);
}

//_____________________________________________________________________________
Int_t
AliMUONClusterStoreV1::GetSize() const
{
  /// Return the number of clusters we hold
  Int_t n(0);
  for ( Int_t i = 0; i < fClusters->GetSize(); ++i ) 
  {
    n += ChamberClusters(i)->GetLast()+1;
  }
  return n;
}

