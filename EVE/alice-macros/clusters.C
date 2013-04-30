#ifndef __CINT__
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TGeoManager.h>
#include <EveBase/AliEveEventManager.h>

#include "AliRunLoader.h"
#include "AliCluster.h"
#include "AliTracker.h"
#include "AliReconstruction.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"

#include "AliITSRecoParam.h"
#endif

void clusters()
{
  AliEveEventManager::AssertGeometry();

  AliRunLoader *rl        = AliEveEventManager::AssertRunLoader();
  AliESDEvent  *esd       = AliEveEventManager::AssertESD();
  AliEveEventManager::AssertESDfriend();
  AliEveEventManager::AssertMagField();

  const char* detNames[] = { "ITS", "TPC", /*"TRD",*/ "TOF", "HMPID" };
  const Int_t detIds[]   = {   0,     1,   /*  2,  */   3,     5   };
  const Int_t detN       = sizeof(detNames)/sizeof(char*);

  // Hack - AliReconstruction does wonders with gGeoManager.
  TGeoManager* xxx = gGeoManager; gGeoManager = 0;
  AliReconstruction* reco = new AliReconstruction;
  gGeoManager = xxx;

  // Hack for ITS - it requires RecoParams outside event-loop!
  reco->SetRecoParam("ITS", AliITSRecoParam::GetLowFluxParam());

  reco->ImportRunLoader(rl);
  {
    TString alldets;
    for (Int_t i = 0; i < detN; ++i) {
      alldets += detNames[i];
      alldets += " ";
    }
    reco->CreateTrackers(alldets);
  }

  TObjArray* clarr = new TObjArray();

  // Load clusters, fill them into clarr.

  for (Int_t i = 0; i < detN; ++i)
  {
    Int_t det = detIds[i];
    rl->LoadRecPoints(detNames[i]);

    TTree *cTree = rl->GetTreeR(detNames[i], false);
    if (cTree == 0)
      continue;

    AliTracker* tracker = reco->GetTracker(det);
    if (tracker == 0) continue;
    tracker->LoadClusters(cTree);
    tracker->FillClusterArray(clarr);
  }

  // Loop over tracks and friends, tag used clusters.

  for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
  {
    AliESDtrack* at = esd->GetTrack(n);

    Int_t idx[200];
    for (Int_t i = 0; i < detN; ++i)
    {
      Int_t det = detIds[i];
      AliTracker* tracker = reco->GetTracker(det);
      if (tracker == 0) continue;
      Int_t nclusters = at->GetClusters(det, idx);
      Int_t p=0;
      for (Int_t c = 0; c < nclusters; )
      {
        Int_t index = idx[p++];
        if (index < 0) continue;
        c++;
        AliCluster* cluster = tracker->GetCluster(index);
        if (cluster) cluster->IncreaseClusterUsage();
        //else printf("Zero cluster pointer for detector: %s\n",detNames[i]);
      }
    }
  }

  for (Int_t i = 0; i < detN; ++i)
    rl->UnloadRecPoints(detNames[i]);

  // Fill visualization structs

  TEveElementList* list = new TEveElementList("Clusters");
  gEve->AddElement(list);

  TEvePointSet* shared = new TEvePointSet("Shared Clusters");
  shared->SetMainColor(2);
  shared->SetMarkerSize(0.4);
  shared->SetMarkerStyle(2);
  list->AddElement(shared);

  TEvePointSet* used = new TEvePointSet("Single-used Clusters");
  used->SetMainColor(3);
  used->SetMarkerSize(0.4);
  used->SetMarkerStyle(2);
  list->AddElement(used);

  TEvePointSet* nonused = new TEvePointSet("Not-used Clusters");
  nonused->SetMainColor(4);
  nonused->SetMarkerSize(0.4);
  nonused->SetMarkerStyle(2);
  list->AddElement(nonused);

  // Loop over all clusters, fill appropriate container based
  // on shared/used information.

  Int_t ncl = clarr->GetEntriesFast();
  for (Int_t i = 0; i < ncl; ++i)
  {
    AliCluster   *cluster = (AliCluster*) clarr->UncheckedAt(i);
    TEvePointSet *dest    = 0;
    if (cluster->IsClusterShared())
      dest = shared;
    else if (cluster->IsClusterUsed())
      dest = used;
    else
      dest = nonused;

    Float_t g[3]; //global coordinates
    cluster->GetGlobalXYZ(g);
    dest->SetNextPoint(g[0], g[1], g[2]);
    dest->SetPointId(cluster);
  }

  delete clarr;

  // ??? What to do with trackers ???
  // I'd propose: have global reconstruction that owns them.
  // AliEveEventManager::AssertAliReconstruction();
  // do we have bit-field for detectors, like
  // enum AliDetectors_e {
  //    kITS = BIT(0),
  //    kTPC = BIT(1),
  //    ...
  //    kCentralTracking = kITS | kTPC | kTRD | kTOF,
  //    ...
  //  };

  gEve->Redraw3D();
}
