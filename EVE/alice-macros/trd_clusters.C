// $Id: tpc_clusters.C 23497 2008-01-23 20:43:14Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifdef __CINT__
class TEvePointSet;
class TEveElement;
#else
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <EveBase/AliEveEventManager.h>

#include "AliRunLoader.h"
#include "AliCluster.h"
#include "AliTRDcluster.h"
#endif

TEvePointSet* trd_clusters(TEveElement *cont = 0){
	const Int_t kMaxClusters = 18 * 6 * 24 *10;
	AliEveEventManager::AssertGeometry();

	TEvePointSet *clusters = new TEvePointSet(kMaxClusters);
	clusters->SetOwnIds(kTRUE);

	AliRunLoader *rl = AliEveEventManager::AssertRunLoader();
	rl->LoadRecPoints("TRD");

	TObjArray *TRDcluster = 0x0;
	TTree *recPoints = rl->GetTreeR("TRD", kFALSE);
	recPoints->SetBranchAddress("TRDcluster", &TRDcluster);

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

  if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("trd_clusters", "No TRD clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.2);
  clusters->SetMarkerColor(4);

  char form[1000];
  sprintf(form,"TRD Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);
  gEve->AddElement(clusters, cont);
  gEve->Redraw3D();

  return clusters;
}
