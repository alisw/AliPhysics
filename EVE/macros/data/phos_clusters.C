// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObjArray.h>
#include <TBranch.h>
#include <TTree.h>
#include <TString.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliCluster.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

TEvePointSet* phos_clusters(TEveElement* cont=0)
{
  AliEveEventManager::Instance()->AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("PHOS");

  TTree *cTree = rl->GetTreeR("PHOS", false);

	if(!cTree) return 0;

  TEvePointSet* clusters = new TEvePointSet(10000);
  clusters->SetOwnIds(kTRUE);

  TObjArray *arr=NULL;
  TBranch *branch=cTree->GetBranch("PHOSEmcRP");
  branch->SetAddress(&arr);

  Int_t nentr=(Int_t)branch->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!branch->GetEvent(i)) continue;

    Int_t ncl=arr->GetEntriesFast();
    while (ncl--) {
      AliCluster *cl=(AliCluster*)arr->UncheckedAt(ncl);

      Float_t g[3]; //global coordinates
      cl->GetGlobalXYZ(g);

      AliCluster *atp = new AliCluster(*cl);
      clusters->SetNextPoint(g[0], g[1], g[2]);
      clusters->SetPointId(atp);
    }
  }

  Warning("phos_clusters"," %d", clusters->Size());

  if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("phos_clusters", "No PHOS clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.5);
  clusters->SetMarkerColor(4);

  clusters->SetName("PHOS Clusters");

  clusters->SetTitle(Form("N=%d", clusters->Size()));

  const TString viz_tag("REC Clusters PHOS");

  clusters->ApplyVizTag(viz_tag, "Clusters");

  gEve->AddElement(clusters);

  gEve->Redraw3D();

  return clusters;
}
