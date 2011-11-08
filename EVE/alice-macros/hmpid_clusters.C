// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifdef __CINT__

class TEveElement;
class TEvePointSet;

#else

#include <TEveManager.h>
#include <TEvePointSet.h>
#include <EveBase/AliEveEventManager.h>
#include <TBranch.h>
#include <TTree.h>
#include <AliRunLoader.h>
#include <AliCluster3D.h>

#include <TClonesArray.h>

#endif

TEvePointSet* hmpid_clusters(TEveElement* cont=0)
{
  const Int_t nCh=7;
  TClonesArray *cl[nCh] = {0,0,0,0,0,0,0};
  const Char_t *name[nCh]={
    "HMPID0",
    "HMPID1",
    "HMPID2",
    "HMPID3",
    "HMPID4",
    "HMPID5",
    "HMPID6"
  };


  TEvePointSet* clusters = new TEvePointSet(10000);
  clusters->SetOwnIds(kTRUE);

  AliEveEventManager::AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("HMPID");

  TTree *cTree = rl->GetTreeR("HMPID", false);
  if (!cTree) return 0;

  for (Int_t k=0; k<nCh; k++) {
     TBranch *br=cTree->GetBranch(name[k]);
     if (!br) return 0;
     br->SetAddress(&(cl[k]));
  }

  if (!cTree->GetEvent(0)) return 0;


  for (Int_t i=0; i<nCh; i++) {
    TClonesArray *arr=cl[i];
    Int_t ncl=arr->GetEntriesFast();

    while (ncl--) {
      AliCluster3D *c=(AliCluster3D*)arr->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      clusters->SetNextPoint(g[0], g[1], g[2]);
      AliCluster3D *atp = new AliCluster3D(*c);
      clusters->SetPointId(atp);
    }
  }

  rl->UnloadRecPoints("HMPID");

  if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("hmpid_clusters", "No HMPID clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.2);
  clusters->SetMarkerColor(4);

  clusters->SetName("HMPID Clusters");

  clusters->SetTitle(Form("N=%d", clusters->Size()));

  const TString viz_tag("REC Clusters HMPID");

  clusters->ApplyVizTag(viz_tag, "Clusters");

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}
