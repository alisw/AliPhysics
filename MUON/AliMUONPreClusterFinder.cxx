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

#include "AliMUONPreClusterFinder.h"

#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMpVSegmentation.h"
#include "TClonesArray.h"
#include "AliMpArea.h"
#include "TVector2.h"
#include "AliMUONPad.h"
#include "AliMUONDigit.h"

/// \class AliMUONPreClusterFinder
///
/// Implementation of AliMUONVClusterFinder
///
/// This class simply find adjacent pads to form clusters
///
/// \author Laurent Aphecetche

ClassImp(AliMUONPreClusterFinder)

//_____________________________________________________________________________
AliMUONPreClusterFinder::AliMUONPreClusterFinder()
: AliMUONVClusterFinder(),
  fClusters(0x0),
  fSegmentations(0x0),
  fDigits(0x0),
  fDetElemId(0)
{
    /// ctor
  for ( Int_t i = 0; i < 2; ++i )
  {
    fPads[i] = 0x0;
  } 
}

//_____________________________________________________________________________
AliMUONPreClusterFinder::~AliMUONPreClusterFinder()
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
AliMUONPreClusterFinder::UsePad(const AliMUONPad& pad)
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
AliMUONPreClusterFinder::Prepare(const AliMpVSegmentation* segmentations[2],
                                 TClonesArray* digits[2]) 
// FIXME : add area on which to look for clusters here.
{
  /// Prepare for clustering, by giving access to segmentations and digit lists
  
  fSegmentations = segmentations;
  fDigits = digits;
  
  delete fClusters;
  fClusters = new TClonesArray("AliMUONCluster");
  for ( Int_t i = 0; i < 2; ++i )
  {
    delete fPads[i];
    fPads[i] = new TClonesArray("AliMUONPad");
  }
  
  fDetElemId = -1;
  
  // Converts digits into pads
  for ( Int_t cathode = 0; cathode < 2; ++cathode )
  {
    if ( !digits[cathode] ) continue;

    AliMUONDigit* d;
    TIter next(digits[cathode]);
    while ( ( d = static_cast<AliMUONDigit*>(next())))
    {
      Int_t ix = d->PadX();
      Int_t iy = d->PadY();
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
                      d->Signal());
      if ( d->IsSaturated() ) mpad.SetSaturated(kTRUE); 
      new (padArray[padArray.GetLast()+1]) AliMUONPad(mpad);      
    }
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
AliMUONPreClusterFinder::AddPad(AliMUONCluster& cluster, AliMUONPad* pad)
{
  /// Add a pad to a cluster
  cluster.AddPad(*pad);
  
  Int_t cathode = pad->Cathode();
  TClonesArray& padArray = *fPads[cathode];
  padArray.Remove(pad);
  padArray.Compress();
  TIter next(&padArray);
  AliMUONPad* testPad;
  
  while ( ( testPad = static_cast<AliMUONPad*>(next())))
  {
    if ( AliMUONPad::AreNeighbours(*testPad,*pad) )
    {
      AddPad(cluster,testPad);
    }
  }
}

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

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONPreClusterFinder::NextCluster()
{
  /// Builds the next cluster, and returns it.
  
  // Start a new cluster
  Int_t id = fClusters->GetLast()+1;
  AliMUONCluster* cluster = new ((*fClusters)[id]) AliMUONCluster;
  cluster->SetUniqueID(id);
  
  AliMUONPad* pad = static_cast<AliMUONPad*>(fPads[0]->First());
  
  if (!pad) // protection against no pad in first cathode, which might happen
  {
    // try other cathode
    pad = static_cast<AliMUONPad*>(fPads[1]->First());
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
    TClonesArray& padArray = *fPads[1];
    TIter next(&padArray);
    AliMUONPad* testPad;
  
    while ( ( testPad = static_cast<AliMUONPad*>(next())))
    {
      if ( AreOverlapping(*testPad,*cluster) )
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
