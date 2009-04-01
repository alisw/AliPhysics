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

#include "AliMUONCluster.h"
#include "AliMUONPad.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

#include "AliMpArea.h"
#include "AliMpConstants.h"
#include "AliMpVSegmentation.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TClonesArray.h>
#include <TVector2.h>

//-----------------------------------------------------------------------------
/// \class AliMUONPreClusterFinder
///
/// Implementation of AliMUONVClusterFinder
///
/// This class simply find adjacent pads to form clusters
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

ClassImp(AliMUONPreClusterFinder)

//_____________________________________________________________________________
AliMUONPreClusterFinder::AliMUONPreClusterFinder()
: AliMUONVClusterFinder(),
  fClusters("AliMUONCluster"),
  fPads(0x0),
  fDetElemId(0),
  fArea(),
  fShouldAbort(kFALSE)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPreClusterFinder::~AliMUONPreClusterFinder()
{
  /// dtor : note we're owner of the clusters, but not of the pads
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
AliMUONPreClusterFinder::Prepare(Int_t detElemId,
                                 TClonesArray* pads[2],
                                 const AliMpArea& area)
{
  /// Prepare for clustering, by giving access to segmentations and digit lists

  fClusters.Clear("C");
  
  fPads = pads;
  fDetElemId = detElemId;
  fArea = area;
  
  fShouldAbort = kFALSE;
  
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONPreClusterFinder::AddPad(AliMUONCluster& cluster, AliMUONPad* pad)
{
  /// Add a pad to a cluster
  
  if ( cluster.Multiplicity() > 500 ) 
  {
    fShouldAbort = kTRUE;
    return;
  }
  
  cluster.AddPad(*pad);
  
  Int_t cathode = pad->Cathode();
  TClonesArray& padArray = *fPads[cathode];
  // WARNING: this Remove method uses the AliMUONPad::IsEqual if that method is
  // present (otherwise just compares pointers) : so that one must be correct
  // if implemented !
  padArray.Remove(pad);
 // TObject* o = padArray.Remove(pad); 
//  if (!o)
//  {
//    AliFatal("Oups. Could not remove pad from pads to consider. Aborting as anyway "
//             " we'll get an infinite loop. Please check the AliMUONPad::IsEqual method"
//             " as the first suspect for failed remove");
//  }  
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
AliMUONPad*
AliMUONPreClusterFinder::GetNextPad(Int_t cathode) const
{
/// Return the next unused pad of given cathode, which is within fArea

  TIter next(fPads[cathode]);
  
  if ( !fArea.IsValid() )
  {
    return static_cast<AliMUONPad*>(next());
  }
  else
  {
    AliMUONPad* pad;
    while ( ( pad = static_cast<AliMUONPad*>(next())) )
    {
      AliMpArea padArea(pad->X(), pad->Y(), pad->DX(), pad->DY());
      
      if (fArea.Overlap(padArea)) return pad;

    }
    return 0x0;
  }
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONPreClusterFinder::NextCluster()
{
  /// Builds the next cluster, and returns it.
//  AliCodeTimerAuto("pre-clustering")
  
  // Start a new cluster
  Int_t id = fClusters.GetLast()+1;
  AliMUONCluster* cluster = new (fClusters[id]) AliMUONCluster;
  cluster->SetUniqueID(id);
  
  AliMUONPad* pad = GetNextPad(0);
  
  if (!pad) // protection against no pad in first cathode, which might happen
  {
    // try other cathode
    pad = GetNextPad(1);
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
    
    if ( !ShouldAbort() ) 
    {
      // On the 2nd cathode, only add pads overlapping with the current cluster
      TClonesArray& padArray = *fPads[1];
      TIter next(&padArray);
      AliMUONPad* testPad;
      
      while ( ( testPad = static_cast<AliMUONPad*>(next())) && !ShouldAbort() )
      {
        if (AreOverlapping(*testPad,*cluster) )
        {
          AddPad(*cluster,testPad);
        }
      }
    }
  }
  
  if ( ShouldAbort() ) 
  {
    AliError("Aborting clustering of that DE because we've got too many pads");
    fClusters.Remove(cluster);
    fClusters.Compress();
    return 0x0;
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
    fClusters.Remove(cluster);
    fClusters.Compress();
    // then proceed further
    return NextCluster();
  }
  
  return cluster;
}
