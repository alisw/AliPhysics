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
#include "TObjArray.h"
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

//_____________________________________________________________________________
AliMUONPreClusterFinderV3::AliMUONPreClusterFinderV3()
: AliMUONVClusterFinder(),
  fClusters(new TClonesArray("AliMUONCluster",10)),
  fkSegmentations(0x0),
  fPads(0x0),
  fDetElemId(0),
  fIterator(0x0)
{
    /// ctor
  AliInfo("");
  for ( Int_t i = 0; i < 2; ++i )
  {
    fPreClusters[i] = new TClonesArray("AliMUONCluster",10);
  } 
}

//_____________________________________________________________________________
AliMUONPreClusterFinderV3::~AliMUONPreClusterFinderV3()
{
  /// dtor
  delete fClusters;
  for ( Int_t i = 0; i < 2; ++i )
  {
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
  
  AliMUONPad* p = new AliMUONPad(pad); 
  p->SetClusterId(-1);
  fPads[pad.Cathode()]->Add(p); 
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUONPreClusterFinderV3::Prepare(Int_t detElemId,
                                   TObjArray* pads[2],
                                   const AliMpArea& area,
                                   const AliMpVSegmentation* seg[2])
{
  /// Prepare for clustering, by giving access to segmentations and digit lists
  
  if ( area.IsValid() ) 
  {
    AliError("Handling of area not yet implemented for this class. Please check.");
  }
  
  fkSegmentations = seg;
  fPads = pads;
  
  fClusters->Clear("C");
  for ( Int_t i = 0; i < 2; ++i )
  {
    fPreClusters[i]->Clear("C");
  }
  
  fDetElemId = detElemId;
  
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
  if ( cathode < 0 ) {
    AliError(Form("Cathod undefined: %d",cathode));
    AliFatal("");
    return;
  }
  
  if ( cathode <=1 && !fPreClusters[cathode]->Remove(preCluster) ) 
  {
    AliError(Form("Could not remove %s from preclusters[%d]",
                  preCluster->AsString().Data(),cathode));
    StdoutToAliDebug(1,DumpPreClusters());
    AliFatal("");
    return;
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
  TObjArray& padArray = *fPads[cathode];
  delete padArray.Remove(pad);
  
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
  
//  AliCodeTimerAuto("",0)
  
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
  
//  AliCodeTimerAuto(Form("Cathode %d",cathode),0)
  
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
