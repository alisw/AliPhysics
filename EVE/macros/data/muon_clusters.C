// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file muon_clusters.C
/// \brief Macro to visualise clusters from MUON spectrometer 
/// (both tracker and trigger).
///
/// Use muon_clusters() in order to run it.
///
/// Needs that alieve_init() is already called.
///
/// \author P. Pillot, L. Aphecetche; Subatech

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEvePointSet.h>

#include <AliMUONVCluster.h>
#include <AliMUONVClusterStore.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif
class TIter;
class TEvePointSet;

//______________________________________________________________________________
void add_muon_clusters(TIter* next, TEvePointSet* clusterList)
{
  // loop over clusters and produce corresponding graphic objects
  AliMUONVCluster* cluster;
  while ( ( cluster = static_cast<AliMUONVCluster*>((*next)()) ) ) 
  {  
    clusterList->SetNextPoint(cluster->GetX(),cluster->GetY(),cluster->GetZ());
  }
}

//______________________________________________________________________________
void muon_clusters()
{
    printf("*** RecPoints MUON ***");
    
  // load clusters
  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("MUON");
  TTree* ct = rl->GetTreeR("MUON",kFALSE);
  if (!ct) return;
  AliMUONVClusterStore* clusterStore = AliMUONVClusterStore::Create(*ct);
  clusterStore->Clear();
  clusterStore->Connect(*ct,kFALSE);
  ct->GetEvent(0);
  rl->UnloadRecPoints("MUON");
  
  if (clusterStore->GetSize() == 0 && !gEve->GetKeepEmptyCont()) {
    delete clusterStore;
    return;
  }
  
  // cluster container
  TEvePointSet* clusterList = new TEvePointSet(10000);
  clusterList->SetName("MUON Clusters");
  clusterList->SetTitle(Form("N=%d",clusterStore->GetSize()));
  clusterList->SetPickable(kFALSE);
  clusterList->SetMarkerStyle(20);
  clusterList->SetMarkerColor(kCyan);
  clusterList->SetMarkerSize(1.);
  
  // add cluster to the container
  TIter next(clusterStore->CreateIterator());
  add_muon_clusters(&next, clusterList);
  delete clusterStore;
  
  // add graphic containers
  gEve->DisableRedraw();
  gEve->AddElement(clusterList);
  gEve->EnableRedraw();
  gEve->Redraw3D();
}
