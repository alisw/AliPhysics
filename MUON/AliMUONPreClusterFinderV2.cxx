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

#include "AliMUONPreClusterFinderV2.h"

#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMpVSegmentation.h"
#include "TClonesArray.h"
#include "TVector2.h"
#include "AliMUONPad.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

//-----------------------------------------------------------------------------
/// \class AliMUONPreClusterFinderV2
///
/// Implementation of AliMUONVClusterFinder
///
/// This one ressembles the preclustering stage in the original ClusterFinderAZ
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

ClassImp(AliMUONPreClusterFinderV2)

//_____________________________________________________________________________
AliMUONPreClusterFinderV2::AliMUONPreClusterFinderV2()
: AliMUONVClusterFinder(),
  fClusters(0x0),
  fSegmentations(0x0),
  fDetElemId(0)
{
    /// ctor
  for ( Int_t i = 0; i < 2; ++i )
  {
    fPads[i] = 0x0;
  } 
}

//_____________________________________________________________________________
AliMUONPreClusterFinderV2::~AliMUONPreClusterFinderV2()
{
  /// dtor : note we're owner of the pads and the clusters, but not of
  /// the remaining objects (digits, segmentations)
  delete fClusters;
  for ( Int_t i = 0; i < 2; ++i )
  {
    delete fPads[i];
  }  
}

//_____________________________________________________________________________
Bool_t
AliMUONPreClusterFinderV2::UsePad(const AliMUONPad& pad)
{
  /// Add a pad to the list of pads to be considered
  if ( pad.DetElemId() != fDetElemId )
  {
    AliError(Form("Cannot add pad from DE %d to this cluster finder which is "
                  "currently dealing with DE %d",pad.DetElemId(),fDetElemId));
    return kFALSE;
  }
  
  new ((*fPads[pad.Cathode()])[fPads[pad.Cathode()]->GetLast()+1]) AliMUONPad(pad); 
  // FIXME: should set the ClusterId of that new pad to be -1
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUONPreClusterFinderV2::Prepare(const AliMpVSegmentation* segmentations[2],
                                 const AliMUONVDigitStore& digitStore)
// FIXME : add area on which to look for clusters here.
{
  /// Prepare for clustering, by giving access to segmentations and digit lists
  
  fSegmentations = segmentations;
  
  delete fClusters;
  fClusters = new TClonesArray("AliMUONCluster");
  for ( Int_t i = 0; i < 2; ++i )
  {
    delete fPads[i];
    fPads[i] = new TClonesArray("AliMUONPad");
  }
  
  fDetElemId = -1;
  
  TIter next(digitStore.CreateIterator());
  AliMUONVDigit* d;
  
  while ( ( d = static_cast<AliMUONVDigit*>(next()) ) )
  {
    Int_t ix = d->PadX();
    Int_t iy = d->PadY();
    Int_t cathode = d->Cathode();
    AliMpPad pad = fSegmentations[cathode]->PadByIndices(AliMpIntPair(ix,iy));
    TClonesArray& padArray = *(fPads[cathode]);
    if ( fDetElemId == -1 ) 
    {
      fDetElemId = d->DetElemId();
    }
    else
    {
      if ( d->DetElemId() != fDetElemId ) 
      {
        AliError("Something is seriously wrong with DE. Aborting clustering");
        return kFALSE;
      }
    }
    
    AliMUONPad mpad(fDetElemId,cathode,
                    ix,iy,pad.Position().X(),pad.Position().Y(),
                    pad.Dimensions().X(),pad.Dimensions().Y(),
                    d->Charge());
    if ( d->IsSaturated() ) mpad.SetSaturated(kTRUE); 
    mpad.SetUniqueID(d->GetUniqueID());
    new (padArray[padArray.GetLast()+1]) AliMUONPad(mpad);      
  }
  if ( fPads[0]->GetLast() < 0 && fPads[1]->GetLast() < 0 )
  {
    // no pad at all, nothing to do...
    return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONPreClusterFinderV2::AddPad(AliMUONCluster& cluster, AliMUONPad* pad)
{
  /// Add a pad to a cluster
  
  cluster.AddPad(*pad);
  pad->SetClusterId(cluster.GetUniqueID());
  
  Int_t cathode = pad->Cathode();
  TClonesArray& padArray = *fPads[cathode];
  padArray.Remove(pad);
  TIter next(&padArray);
  
  // Check neighbours
  TObjArray neighbours;
  AliMpPad p = fSegmentations[pad->Cathode()]->PadByIndices(AliMpIntPair(pad->Ix(),pad->Iy()),kTRUE);
  Int_t nn = fSegmentations[pad->Cathode()]->GetNeighbours(p,neighbours);
  for (Int_t in = 0; in < nn; ++in) 
  {
    AliMpPad* p = static_cast<AliMpPad*>(neighbours.At(in));
    
    TIter next(&padArray);
    AliMUONPad* p2;
    
    while ( ( p2 = static_cast<AliMUONPad*>(next()) ) )
    {
        if ( !p2->IsUsed() && p2->Ix()==p->GetIndices().GetFirst() 
             && p2->Iy() == p->GetIndices().GetSecond() &&
             p2->Cathode() == pad->Cathode() )
        {
          AddPad(cluster,p2);
        }
    }
  } // for (Int_t in = 0;
}

namespace
{
//_____________________________________________________________________________
Bool_t
AreOverlapping(const AliMUONPad& pad, const AliMUONCluster& cluster)
{
  /// Whether the pad overlaps with the cluster
  
  static Double_t precision = 1E-4; // cm
  static TVector2 precisionAdjustment(precision,precision);//-precision,-precision);
  for ( Int_t i = 0; i < cluster.Multiplicity(); ++i )
  {
    AliMUONPad* testPad = cluster.Pad(i);
    // Note: we use negative precision numbers, meaning
    // the area of the pads will be *increased* by these small numbers
    // prior to check the overlap by the AreOverlapping method,
    // so pads touching only by the corners will be considered as
    // overlapping.    
    if ( AliMUONPad::AreOverlapping(*testPad,pad,precisionAdjustment) )
    {
      return kTRUE;
    }
  }
  return kFALSE;
}
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONPreClusterFinderV2::NextCluster()
{
  /// Builds the next cluster, and returns it.
  
  // Start a new cluster
  Int_t id = fClusters->GetLast()+1;
  AliMUONCluster* cluster = new ((*fClusters)[id]) AliMUONCluster;
  cluster->SetUniqueID(id);
  
  AliMUONPad* pad;
  TIter next(fPads[0]);
  while (  ( pad = static_cast<AliMUONPad*>(next())) && pad->IsUsed() );

  if (!pad) // protection against no pad in first cathode, which might happen
  {
    // try other cathode
    TIter next(fPads[1]);
    while (  ( pad = static_cast<AliMUONPad*>(next())) && pad->IsUsed() );
    if (!pad) 
    {
      // we are done.
      return 0x0;
    }
    // Builds (recursively) a cluster on second cathode only
    AddPad(*cluster,pad);
  }
  else
  {
    // Builds (recursively) a cluster on first cathode only
      
    AddPad(*cluster,pad);

    // On the 2nd cathode, only add pads overlapping with the current cluster
    TIter next1(fPads[1]);
    AliMUONPad* testPad;
  
    while ( ( testPad = static_cast<AliMUONPad*>(next1())))
    {
      if ( !testPad->IsUsed() && AreOverlapping(*testPad,*cluster) )
      {
        AddPad(*cluster,testPad);
      }
    }
  }
  
  if ( cluster->Multiplicity() <= 1 )
  {
    if ( cluster->Multiplicity() == 0 ) 
    {
      // no pad is suspicious
      AliWarning("Got an empty cluster...");
    }
    // else only 1 pad (not suspicious, but kind of useless, probably noise)
    // so we remove it from our list
    fClusters->Remove(cluster);
    fClusters->Compress();
    // then proceed further
    return NextCluster();
  }
  
  return cluster;
}
