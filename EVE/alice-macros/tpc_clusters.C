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

#include <AliRunLoader.h>
#include <AliCluster.h>
#include <TPC/AliTPCClustersRow.h>

#endif

TEvePointSet* tpc_clusters(TEveElement* cont=0, Float_t maxR=270)
{
  const Int_t kMaxCl=100*160;

  AliEveEventManager::AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("TPC");

  TTree *cTree = rl->GetTreeR("TPC", false);
  if (cTree == 0)
    return 0;

  AliTPCClustersRow *clrow = new AliTPCClustersRow();
  clrow->SetClass("AliTPCclusterMI");
  clrow->SetArray(kMaxCl);
  cTree->SetBranchAddress("Segment", &clrow);

  TEvePointSet* clusters = new TEvePointSet(kMaxCl);
  clusters->SetOwnIds(kTRUE);


  Float_t maxRsqr = maxR*maxR;
  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++)
  {
    if (!cTree->GetEvent(i)) continue;

    TClonesArray *cl = clrow->GetArray();
    Int_t ncl = cl->GetEntriesFast();

    while (ncl--)
    {
      AliCluster *c = (AliCluster*) cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      if (g[0]*g[0]+g[1]*g[1] < maxRsqr)
      {
	clusters->SetNextPoint(g[0], g[1], g[2]);
	AliCluster *atp = new AliCluster(*c);
	clusters->SetPointId(atp);
      }
    }
    cl->Clear();
  }

  delete clrow;

  rl->UnloadRecPoints("TPC");

  if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE)
  {
    Warning("tpc_clusters.C", "No TPC clusters");
    delete clusters;
    return 0;
  }

  char form[1000];
  sprintf(form,"TPC Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  const TString viz_tag("TPC Clusters");
  clusters->ApplyVizTag(viz_tag, "Clusters");

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}
