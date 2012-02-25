// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObjArray.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliCluster.h>
#include <AliRunLoader.h>
#include <AliTRDcluster.h>
#include <AliEveEventManager.h>
#else
class TEvePointSet;
class TEveElement;
#endif

TEvePointSet* trd_clusters(TEveElement *cont = 0)
{
  const Int_t kMaxClusters = 18 * 6 * 24 *10;

  AliEveEventManager::AssertGeometry();

  AliRunLoader *rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("TRD");

  TTree *recPoints = rl->GetTreeR("TRD", kFALSE);
  if (recPoints == 0)
    return 0;

  TObjArray *TRDcluster = 0x0;
  recPoints->SetBranchAddress("TRDcluster", &TRDcluster);

  TEvePointSet *clusters = new TEvePointSet(kMaxClusters);
  clusters->SetOwnIds(kTRUE);


  Int_t nentr=(Int_t)recPoints->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!recPoints->GetEvent(i)) continue;

    Int_t ncl=TRDcluster->GetEntriesFast();

    while (ncl--) {
      AliTRDcluster *c = (AliTRDcluster*)TRDcluster->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      clusters->SetNextPoint(g[0], g[1], g[2]);
      AliCluster *atp = new AliCluster(*c);
      clusters->SetPointId(atp);
    }
    TRDcluster->Clear();
  }

  rl->UnloadRecPoints("TRD");

  if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("trd_clusters.C", "No TRD clusters");
    delete clusters;
    return 0;
  }

  clusters->SetName("TRD Clusters");

  clusters->SetTitle(Form("N=%d", clusters->Size()));

  const TString viz_tag("REC Clusters TRD");

  clusters->ApplyVizTag(viz_tag, "Clusters");

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}
