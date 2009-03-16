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

/// \class AliMUONLegacyClusterServer
///
/// Special implementation of AliMUONVClusterServer, which will only return 
/// clusters from a pre-defined cluster store.
///
/// Made to recover the old (i.e. before introduction of VClusterServer) behavior
/// of the MUON recontruction where rec points were always written to TreeR, 
/// and then the tracking picked them from that tree, in order to have the
/// possibility to save full rec points (for debugging the spectro, mainly, should
/// not be an option used during final production).
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONLegacyClusterServer.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerTrackToTrackerClusters.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONRecoParam.h"
#include "AliMpArea.h"
#include <TCollection.h>

/// \cond CLASSIMP
ClassImp(AliMUONLegacyClusterServer)
/// \endcond

//_____________________________________________________________________________
AliMUONLegacyClusterServer::AliMUONLegacyClusterServer(const AliMUONGeometryTransformer& transformer, 
																											 AliMUONVClusterStore* store,
																											 Bool_t bypassSt4, Bool_t bypassSt5)
: AliMUONVClusterServer(), fkTransformer(transformer), fClusterStore(store), fTriggerTrackStore(0x0),
fBypass(0x0),
fBypassSt4(bypassSt4),
fBypassSt5(bypassSt5)
{
  /// ctor. Mode Read : we'll only server clusters from existing store
}

//_____________________________________________________________________________
AliMUONLegacyClusterServer::~AliMUONLegacyClusterServer()
{
  /// dtor
  delete fBypass;
}

//_____________________________________________________________________________
Int_t 
AliMUONLegacyClusterServer::Clusterize(Int_t chamberId, 
                                       AliMUONVClusterStore& clusterStore,
                                       const AliMpArea& /*area*/,
                                       const AliMUONRecoParam* /*recoParam*/)
{
  /// Fills clusterStore with clusters in given chamber
  ///
  /// Return the number of clusters added to clusterStore
  
  AliCodeTimerAuto(Form("Chamber %d",chamberId));

  if ( fBypassSt4 && ( chamberId == 6 || chamberId == 7 ) ) 
  {
    return fBypass->GenerateClusters(chamberId,clusterStore);
  }

	if ( fBypassSt5 && ( chamberId == 8 || chamberId == 9 ) ) 
  {
    return fBypass->GenerateClusters(chamberId,clusterStore);
  }
	
  AliDebug(1,Form("chamberId=%d fClusterStore(%p).GetSize()=%d clusterStore(%p).GetSize()=%d",
                  chamberId,
                  fClusterStore,fClusterStore->GetSize(),
                  &clusterStore,clusterStore.GetSize()));
  
  TIter next(fClusterStore->CreateChamberIterator(chamberId,chamberId));
  AliMUONVCluster* cluster;
  Int_t n(0);
  TObjArray a;
  
  while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
  {
    clusterStore.Add(*cluster);
    a.Add(cluster);
    ++n;
  }
  
  TIter remove(&a);
  while ( ( cluster = static_cast<AliMUONVCluster*>(remove()) ) )
  {
    fClusterStore->Remove(*cluster);
  }
  
  AliDebug(1,Form("n=%d remaining clusters=%d",n,fClusterStore->GetSize()));
  
  return n;
}

//_____________________________________________________________________________
Bool_t 
AliMUONLegacyClusterServer::UseTriggerTrackStore(AliMUONVTriggerTrackStore* trackStore)
{
  /// Tells us to use trigger track store, and thus to bypass St4 and/or 5 clusters
  fTriggerTrackStore = trackStore; // not owner
  delete fBypass;
  fBypass = new AliMUONTriggerTrackToTrackerClusters(fkTransformer,fTriggerTrackStore);
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONLegacyClusterServer::UseDigits(TIter&, AliMUONVDigitStore*)
{
  /// Give the iterator to our delegate if we have one, of issue and error
  
  AliError("Not implemented for this class, as we're not writing clusters, but reading them instead !");
}

