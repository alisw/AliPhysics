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
/// \class AliMUONClusterStoreV2
///
/// Implementation of VClusterStore.
///
/// Note that clusters are identified by their UniqueID, so it MUST be correctly set
///
/// \author Philippe Pillot, Subatech
///
//-----------------------------------------------------------------------------

#include "AliMUONClusterStoreV2.h"

#include "AliMUONRawClusterV2.h"
#include "AliMUONClusterStoreV2Iterator.h"
#include "AliMUONTreeManager.h"
#include "AliMpConstants.h"
#include "AliMpExMap.h"

#include "AliLog.h"

#include <TTree.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterStoreV2)
/// \endcond

//_____________________________________________________________________________
AliMUONClusterStoreV2::AliMUONClusterStoreV2() 
: AliMUONVClusterStore(), 
  fClusters(new TClonesArray("AliMUONRawClusterV2",100)),
  fMap(0x0),
  fMapped(kFALSE)
{
  /// Constructor
}

//_____________________________________________________________________________
AliMUONClusterStoreV2::AliMUONClusterStoreV2(const AliMUONClusterStoreV2& store)
: AliMUONVClusterStore(), 
  fClusters(new TClonesArray(*(store.fClusters))),
  fMap(0x0),
  fMapped(kFALSE)
{
  /// Copy constructor
  if (store.fMapped) ReMap();
}

//_____________________________________________________________________________
AliMUONClusterStoreV2& AliMUONClusterStoreV2::operator=(const AliMUONClusterStoreV2& store)
{
  /// Assignment operator
  fClusters = new TClonesArray(*(store.fClusters));
  fMap = 0x0;
  fMapped = kFALSE;
  if (store.fMapped) ReMap();
  return *this;
}

//_____________________________________________________________________________
AliMUONClusterStoreV2::~AliMUONClusterStoreV2()
{
  /// Destructor
  delete fClusters;
  delete fMap;
}

//_____________________________________________________________________________
void AliMUONClusterStoreV2::Clear(Option_t*)
{
  /// Clear the internal cluster array AND the index
  fClusters->Clear("C");
  if (fMap) {
    Int_t nChamber = AliMpConstants::NofTrackingChambers();
    for (Int_t chamber=0; chamber<nChamber; chamber++) {
      AliMpExMap *map = static_cast<AliMpExMap *>(fMap->UncheckedAt(chamber));
      map->Clear("C");
    }
    fMapped = kFALSE;
  }
}

//_____________________________________________________________________________
Bool_t AliMUONClusterStoreV2::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this to the tree, i.e. make the branches or set their addresses.
  
  AliMUONTreeManager tman;
  
  if (tree.GetBranch("MUONRawClusters")) {
    
    if (alone) tman.UpdateBranchStatuses(tree,"MUONRawClusters");
    
    return tman.SetAddress(tree,"MUONRawClusters", 
			 const_cast<TClonesArray**>(&fClusters));
  } else {
    
    return tman.MakeBranch(tree,ClassName(),"TClonesArray", "MUONRawClusters",
			 const_cast<TClonesArray**>(&fClusters));
  }
    
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV2::CreateCluster(Int_t chamberId, Int_t detElemId, Int_t clusterIndex) const
{
  /// Create a cluster
  return new AliMUONRawClusterV2(chamberId, detElemId, clusterIndex);
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV2::Add(const AliMUONVCluster& vCluster)
{
  /// Add a cluster to this store
  const AliMUONRawClusterV2* cluster = dynamic_cast<const AliMUONRawClusterV2*>(&vCluster);
  
  if (!cluster) {
    AliError(Form("Cluster is not of the expected type (%s vs AliMUONRawClusterV2)",
                  vCluster.ClassName()));
    return 0x0;
  }
  
  // check chamberId
  Int_t chamberId = cluster->GetChamberId();
  if (chamberId < 0 || chamberId >= AliMpConstants::NofTrackingChambers()) {
    AliError(Form("ChamberId (%d) out of boundaries [0,%d[",chamberId,AliMpConstants::NofTrackingChambers()));
    return 0x0;
  }
  
  // check that there is no cluster with the same Id
  AliMUONVCluster *c = FindObject(cluster->GetUniqueID());
  if (c) {
    AliError("cluster store already contains a cluster with the same ID --> add() exited:");
    c->Print("FULL");
    return 0x0;
  }
  
  // add new cluster
  c = new((*fClusters)[fClusters->GetLast()+1]) AliMUONRawClusterV2(*cluster);
  
  if (c) UpdateMap(*c);
  
  return c;
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV2::Add(Int_t chamberId, Int_t detElemId, Int_t clusterIndex)
{
  /// Add an empty cluster with an unique ID to this store
  
  // check chamberId
  if (chamberId < 0 || chamberId >= AliMpConstants::NofTrackingChambers()) {
    AliError(Form("ChamberId (%d) out of boundaries [0,%d[",chamberId,AliMpConstants::NofTrackingChambers()));
    return 0x0;
  }
  
  // check that there is no cluster with the same Id
  AliMUONVCluster *c = FindObject(AliMUONVCluster::BuildUniqueID(chamberId, detElemId, clusterIndex));
  if (c) {
    AliError("cluster store already contains a cluster with the same ID --> add() exited:");
    c->Print("FULL");
    return 0x0;
  }
  
  // add new cluster
  c = new((*fClusters)[fClusters->GetLast()+1]) AliMUONRawClusterV2(chamberId, detElemId, clusterIndex);
  
  if (c) UpdateMap(*c);
  
  return c;
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV2::Remove(AliMUONVCluster& cluster)
{
  /// Remove a cluster
  AliMUONVCluster* c = static_cast<AliMUONVCluster*>(fClusters->Remove(&cluster));
  
  if (c) 
  {
    fClusters->Compress();
    fMapped = kFALSE;
  }
  else
  {
    AliError("Could not remove cluster from array");
  }
  
  return c;
}

//_____________________________________________________________________________
void AliMUONClusterStoreV2::ReMap()
{
  /// Recompute the fMap, which map (ch) to an index within the fClusters array
  fMapped = kTRUE;
  
  // Create (or clear) the TClonesArray of map
  Int_t nChamber = AliMpConstants::NofTrackingChambers();
  
  if (!fMap) {
    fMap = new TClonesArray("AliMpExMap",nChamber);
    
    // Create one map per chamber
    AliMpExMap *map;
    for (Int_t chamber=0; chamber<nChamber; chamber++) {
      map = new((*fMap)[chamber]) AliMpExMap;
      map->SetOwner(kFALSE);
    }
  }
  else {
    for (Int_t chamber=0; chamber<nChamber; chamber++) {
      AliMpExMap *map = static_cast<AliMpExMap *>(fMap->UncheckedAt(chamber));
      map->Clear("C");
    }
  }  

  // Fill the maps
  TIter next(fClusters);
  AliMUONVCluster* cluster;
  while ( (cluster = static_cast<AliMUONVCluster*>(next())) ) UpdateMap(*cluster);
}

//_____________________________________________________________________________
void AliMUONClusterStoreV2::UpdateMap(AliMUONVCluster& cluster)
{
  /// Update the internal index given this new cluster
  if (fMapped) static_cast<AliMpExMap*>(fMap->UncheckedAt(cluster.GetChamberId()))->Add(cluster.GetUniqueID(),&cluster);
  else ReMap();
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV2::FindObject(const TObject* object) const
{
  /// Find an object, if of AliMUONVCluster type.
  const AliMUONVCluster* cluster = dynamic_cast<const AliMUONVCluster*>(object);
  if (cluster) return FindObject(cluster->GetUniqueID());
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVCluster* AliMUONClusterStoreV2::FindObject(UInt_t uniqueID) const
{
  /// Find a cluster by its UniqueID
  if (!fMapped) (const_cast<AliMUONClusterStoreV2*>(this))->ReMap();
  AliMpExMap* map = static_cast<AliMpExMap*>(fMap->UncheckedAt(AliMUONVCluster::GetChamberId(uniqueID)));
  return static_cast<AliMUONVCluster*>(map->GetValue(uniqueID));
}

//_____________________________________________________________________________
TIterator* AliMUONClusterStoreV2::CreateIterator() const
{
  /// Return an iterator to loop over all clusters
  return fClusters->MakeIterator();
}

//_____________________________________________________________________________
TIterator* AliMUONClusterStoreV2::CreateChamberIterator(Int_t firstChamber, Int_t lastChamber) const
{
  /// Return an iterator to loop over clusters in the chambers within the given range
  
  // check validity of given chamber IDs
  if (firstChamber < 0 || firstChamber >= AliMpConstants::NofTrackingChambers()) {
    AliError(Form("First chamber out of boundaries [0,%d[", AliMpConstants::NofTrackingChambers()));
    return 0x0;
  }
  if (lastChamber < 0 || lastChamber >= AliMpConstants::NofTrackingChambers()) {
    AliError(Form("Last chamber out of boundaries [0,%d[", AliMpConstants::NofTrackingChambers()));
    return 0x0;
  }
  
  if (!fMapped) (const_cast<AliMUONClusterStoreV2*>(this))->ReMap();
  return new AliMUONClusterStoreV2Iterator(this,firstChamber,lastChamber);
}
