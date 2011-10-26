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
#include <AliCluster.h>

#include <TClonesArray.h>

#endif

TEvePointSet* its_clusters(TEveElement* cont=0, Float_t maxR=50)
{
  AliEveEventManager::AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("ITS");

  TTree *cTree = rl->GetTreeR("ITS", false);
  if (cTree == 0)
    return 0;

  TClonesArray *cl = NULL;
  TBranch *branch  = cTree->GetBranch("ITSRecPoints");
  branch->SetAddress(&cl);

  TEvePointSet* clusters = new TEvePointSet(10000);
  clusters->SetOwnIds(kTRUE);

  Int_t nentr = (Int_t) cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++)
  {
    if (!cTree->GetEvent(i)) continue;

    Int_t ncl = cl->GetEntriesFast();

    Float_t maxRsqr = maxR*maxR;
    for (Int_t icl = 0; icl < ncl; ++icl)
    {
      AliCluster *c = (AliCluster*) cl->UncheckedAt(icl);
      // This really should not happen, but did in online display once.
      if (c == 0)
      {
	::Warning("its_clusters", "Got NULL AliCluster*, idx=%d, N=%d.",
		  icl, ncl);
	continue;
      }
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      if (g[0]*g[0] + g[1]*g[1] < maxRsqr)
      {
	clusters->SetNextPoint(g[0], g[1], g[2]);
	AliCluster *atp = new AliCluster(*c);
	clusters->SetPointId(atp);
      }
    }
  }

  rl->UnloadRecPoints("ITS");

  if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("its_clusters.C", "No ITS clusters");
    delete clusters;
    return 0;
  }

  clusters->SetName("ITS Clusters");

  clusters->SetTitle(Form("N=%d", clusters->Size()));

  const TString viz_tag("REC Clusters ITS");

  clusters->ApplyVizTag(viz_tag, "Clusters");

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}
