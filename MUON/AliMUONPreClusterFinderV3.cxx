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

#include "AliMUONPreClusterFinderV3.h"

#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMpVSegmentation.h"
#include "TClonesArray.h"
#include "AliMpArea.h"
#include "TVector2.h"
#include "AliMUONPad.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include <Riostream.h>
//#include "AliCodeTimer.h"

//-----------------------------------------------------------------------------
/// \class AliMUONPreClusterFinderV3
///
/// Implementation of AliMUONVClusterFinder
///
/// This version uses a 2 steps approach :
///
/// we first clusterize each cathode independently to form proto-preclusters, 
/// and then we try to "merge" proto-preclusters from the two cathodes
/// when thoses proto-preclusters overlap, thus ending up with preclusters
/// spanning the two cathodes.
///
/// This implementation, on the contrary to PreClusterFinder or PreClusterFinderV2
/// should not depend on the order of the input digits.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

ClassImp(AliMUONPreClusterFinderV3)

namespace
{
  //___________________________________________________________________________
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
AliMUONPreClusterFinderV3::AliMUONPreClusterFinderV3()
: AliMUONVClusterFinder(),
  fClusters(new TClonesArray("AliMUONCluster",10)),
  fSegmentations(0x0),
  fDetElemId(0),
  fIterator(0x0)
{
    /// ctor
    AliInfo("")
  for ( Int_t i = 0; i < 2; ++i )
  {
    fPads[i] = new TClonesArray("AliMUONPad",100);
    fPreClusters[i] = new TClonesArray("AliMUONCluster",10);
  } 
}

//_____________________________________________________________________________
AliMUONPreClusterFinderV3::~AliMUONPreClusterFinderV3()
{
  /// dtor : note we're owner of the pads and the clusters, but not of
  /// the remaining objects (digits, segmentations)
  delete fClusters;
  for ( Int_t i = 0; i < 2; ++i )
  {
    delete fPads[i];
    delete fPreClusters[i];
  } 
}

//_____________________________________________________________________________
Bool_t
AliMUONPreClusterFinderV3::UsePad(const AliMUONPad& pad)
{
  /// Add a pad to the list of pads to be considered
  if ( pad.DetElemId() != fDetElemId )
  {
    AliError(Form("Cannot add pad from DE %d to this cluster finder which is "
                  "currently dealing with DE %d",pad.DetElemId(),fDetElemId));
    return kFALSE;
  }
  
  AliMUONPad* p = new ((*fPads[pad.Cathode()])[fPads[pad.Cathode()]->GetLast()+1]) AliMUONPad(pad); 
  p->SetClusterId(-1);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUONPreClusterFinderV3::Prepare(const AliMpVSegmentation* segmentations[2],
                                 const AliMUONVDigitStore& digitStore)
{
  /// Prepare for clustering, by giving access to segmentations and digit lists
  // FIXME : add area on which to look for clusters here.
  
  fSegmentations = segmentations;
  
  fClusters->Clear("C");
  for ( Int_t i = 0; i < 2; ++i )
  {
    fPads[i]->Clear("C");
    fPreClusters[i]->Clear("C");
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
  
  MakeCathodePreClusters(0);  
  MakeCathodePreClusters(1);  
  MakeClusters();
  
  delete fIterator;
  fIterator = fClusters->MakeIterator();
  
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONPreClusterFinderV3::DumpPreClusters() const
{
  /// Dump preclusters 
  AliMUONCluster *c;
  TIter next0(fPreClusters[0]);
  TIter next1(fPreClusters[1]);
  cout << "Cath0" << endl;
  while ( ( c = static_cast<AliMUONCluster*>(next0())) ) 
  {
    cout << c->AsString().Data() << endl;
  }
  cout << "Cath1" << endl;
  while ( ( c = static_cast<AliMUONCluster*>(next1())) ) 
  {
    cout << c->AsString().Data() << endl;
  }
}

//_____________________________________________________________________________
void
AliMUONPreClusterFinderV3::AddPreCluster(AliMUONCluster& cluster, AliMUONCluster* preCluster)
{
  /// Add a pad to a cluster

  AliMUONCluster a(*preCluster);

  Int_t cathode = preCluster->Cathode();
  if ( cathode <=1 && !fPreClusters[cathode]->Remove(preCluster) ) 
  {
    AliError(Form("Could not remove %s from preclusters[%d]",
                  preCluster->AsString().Data(),cathode));
    StdoutToAliDebug(1,DumpPreClusters());
    AliFatal("");
  }
             
  cluster.AddCluster(a);
  
  // loop on the *other* cathode
  TIter next(fPreClusters[1-cathode]);
  AliMUONCluster* testCluster;
  
  while ( ( testCluster = static_cast<AliMUONCluster*>(next())))
  {
    if ( AliMUONCluster::AreOverlapping(a,*testCluster) )
    {
      AddPreCluster(cluster,testCluster);
    }
  }
}


//_____________________________________________________________________________
void
AliMUONPreClusterFinderV3::AddPad(AliMUONCluster& cluster, AliMUONPad* pad)
{
  /// Add a pad to a cluster
  AliMUONPad* addedPad = cluster.AddPad(*pad);
  
  Int_t cathode = pad->Cathode();
  TClonesArray& padArray = *fPads[cathode];
  padArray.Remove(pad);
  TIter next(&padArray);
  AliMUONPad* testPad;
  
  while ( ( testPad = static_cast<AliMUONPad*>(next())))
  {
    if ( AliMUONPad::AreNeighbours(*testPad,*addedPad) )
    {
      AddPad(cluster,testPad);
    }
  }
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONPreClusterFinderV3::NextCluster()
{
  /// Returns the next cluster
  
  return static_cast<AliMUONCluster*>(fIterator->Next());
}

//_____________________________________________________________________________
void
AliMUONPreClusterFinderV3::MakeClusters()
{
  /// Associate (proto)preclusters to form (pre)clusters
  
//  AliCodeTimerAuto("")
  
  for ( Int_t cathode = 0; cathode < 2; ++cathode ) 
  {
    TClonesArray& preclusters = *(fPreClusters[cathode]);
    
    TIter next(&preclusters);
    AliMUONCluster* preCluster(0x0);
    
    while ( ( preCluster = static_cast<AliMUONCluster*>(next()) ) )
    {
      Int_t id(fClusters->GetLast()+1);
      AliMUONCluster* cluster = new((*fClusters)[id]) AliMUONCluster;
      cluster->SetUniqueID(id);      
      AddPreCluster(*cluster,preCluster);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONPreClusterFinderV3::MakeCathodePreClusters(Int_t cathode)
{
  /// Build (proto)preclusters from digits on a given cathode
  
//  AliCodeTimerAuto(Form("Cathode %d",cathode))
  
  while ( fPads[cathode]->GetLast() > 0  )
  {  
    TIter next(fPads[cathode]);
    AliMUONPad* pad = static_cast<AliMUONPad*>(next());
  
    if (!pad) AliFatal("");

    Int_t id = fPreClusters[cathode]->GetLast()+1;
    AliMUONCluster* cluster = new ((*fPreClusters[cathode])[id]) AliMUONCluster;
    cluster->SetUniqueID(id);
    
    // Builds (recursively) a cluster on first cathode only
    AddPad(*cluster,pad);
    
    if ( cluster->Multiplicity() <= 1 )
    {
      if ( cluster->Multiplicity() == 0 ) 
      {
        // no pad is suspicious
        AliWarning("Got an empty cluster...");
      }
      // else only 1 pad (not suspicious, but kind of useless, probably noise)
      // so we remove it from our list
      fPreClusters[cathode]->Remove(cluster);
      // then proceed further
    }
  }
  
}
