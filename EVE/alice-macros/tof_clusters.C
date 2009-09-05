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

#endif

TEvePointSet* tof_clusters(TEveElement* cont=0, Float_t maxR=390)
{
  AliEveEventManager::AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("TOF");

  TTree *cTree = rl->GetTreeR("TOF", false);
  if (cTree == 0)
    return 0;

  TObjArray *cl = 0x0;
  cTree->SetBranchAddress("TOF", &cl);

  TEvePointSet* clusters = new TEvePointSet(10000);
  clusters->SetOwnIds(kTRUE);

  Float_t maxRsqr = maxR*maxR;

  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!cTree->GetEvent(i)) continue;

    Int_t ncl=cl->GetEntriesFast();
    while (ncl--) {
      AliCluster *c=(AliCluster*)cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);

      if (g[0]*g[0]+g[1]*g[1] < maxRsqr) {
	clusters->SetNextPoint(g[0], g[1], g[2]);
	AliCluster *atp = new AliCluster(*c);
	clusters->SetPointId(atp);
      }

    }
  }

  rl->UnloadRecPoints("TOF");

  if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("tof_clusters.C", "No TOF clusters");
    delete clusters;
    return 0;
  }

  char form[1000];
  sprintf(form,"TOF Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  const TString viz_tag("TOF Clusters");
  clusters->ApplyVizTag(viz_tag, "Clusters");

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}

TEvePointSet* tof_clusters_sec(Int_t selectedSector,
                               TEveElement* cont=0, Float_t maxR=390)
{
  AliEveEventManager::AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("TOF");

  TTree *cTree = rl->GetTreeR("TOF", false);
  if (cTree == 0)
    return 0;

  TObjArray *cl = 0x0;
  cTree->SetBranchAddress("TOF", &cl);

  TEvePointSet* clusters = new TEvePointSet(10000);
  clusters->SetOwnIds(kTRUE);

  Float_t maxRsqr = maxR*maxR;
  Double_t phiAngle = maxR*maxR;

  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!cTree->GetEvent(i)) continue;

    Int_t ncl=cl->GetEntriesFast();
    while (ncl--) {
      AliCluster *c=(AliCluster*)cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);

      phiAngle = 180./TMath::Pi()*(TMath::ATan2(-g[1],-g[0])+TMath::Pi());

      if (g[0]*g[0]+g[1]*g[1] < maxRsqr) {
	if (
	    (selectedSector!=-1 &&
	     phiAngle>=selectedSector*20. && phiAngle<(selectedSector+1)*20.)
	    ||
	    selectedSector==-1)
	  {
	    clusters->SetNextPoint(g[0], g[1], g[2]);
	    AliCluster *atp = new AliCluster(*c);
	    clusters->SetPointId(atp);
	  }
      }

    }
  }

  rl->UnloadRecPoints("TOF");

  if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("tof_clusters.C", "No TOF clusters");
    delete clusters;
    return 0;
  }

  char form[1000];
  sprintf(form,"TOF Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  const TString viz_tag("TOF Clusters");
  // when going to new root call:
  // clusters->ApplyVizTag(viz_tag, "Clusters");
  clusters->ApplyVizTag(viz_tag);

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}
