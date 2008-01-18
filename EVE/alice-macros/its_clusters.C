// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/
#ifdef __CINT__

namespace TEveUtil
{
class TEveElement;
class TEvePointSet;
}

#else

#include <TEveManager.h>
#include <TEvePointSet.h>
#include <Alieve/AliEveEventManager.h>

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

  TEvePointSet* clusters = new TEvePointSet(10000);
  clusters->SetOwnIds(kTRUE);

  TClonesArray *cl = NULL;
  TBranch *branch  = cTree->GetBranch("ITSRecPoints");
  branch->SetAddress(&cl);

  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!cTree->GetEvent(i)) continue;

    Int_t ncl=cl->GetEntriesFast();

    Float_t maxRsqr = maxR*maxR;
    while (ncl--) {
      AliCluster *c=(AliCluster*)cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      if (g[0]*g[0]+g[1]*g[1] < maxRsqr)
      {
	clusters->SetNextPoint(g[0], g[1], g[2]);
	AliCluster *atp = new AliCluster(*c);
	clusters->SetPointId(atp);
      }
    }
  }

  if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("its_clusters", "No ITS clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.2);
  clusters->SetMarkerColor(4);

  char form[1000];
  sprintf(form,"ITS Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);
  gEve->AddElement(clusters, cont);
  gEve->Redraw3D();

  return clusters;
}
