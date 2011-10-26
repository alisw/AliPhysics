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
#include <TColor.h>
#include <TEvePointSet.h>
#include <EveBase/AliEveEventManager.h>
#include <TTree.h>

#include <AliRunLoader.h>
#include <AliCluster.h>
#include <TPC/AliTPCClustersRow.h>
#include <TPC/AliTPCclusterMI.h>

#endif

TEvePointSet* tpc_clusters(TEveElement* cont=0, Float_t maxR=270)
{

  const Int_t kMaxCl=100*160;

  Int_t fNColorBins = 5;

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

  TEvePointSetArray * cc = new TEvePointSetArray("TPC Clusters Colorized");
  cc->SetMainColor(kRed);
  cc->SetMarkerStyle(4);
  cc->SetMarkerSize(0.4);
  cc->InitBins("Cluster Charge", fNColorBins, 0., fNColorBins*60.);
  
  cc->GetBin(0)->SetMainColor(kGray);
  cc->GetBin(0)->SetMarkerSize(0.4);
  cc->GetBin(1)->SetMainColor(kBlue);
  cc->GetBin(1)->SetMarkerSize(0.42);
  cc->GetBin(2)->SetMainColor(kCyan);
  cc->GetBin(2)->SetMarkerSize(0.44);
  cc->GetBin(3)->SetMainColor(kGreen);
  cc->GetBin(3)->SetMarkerSize(0.46);
  cc->GetBin(4)->SetMainColor(kYellow);
  cc->GetBin(4)->SetMarkerSize(0.48);
  cc->GetBin(5)->SetMainColor(kRed);
  cc->GetBin(5)->SetMarkerSize(0.50);
  cc->GetBin(6)->SetMainColor(kMagenta);
  cc->GetBin(6)->SetMarkerSize(0.52);

  Float_t maxRsqr = maxR*maxR;
  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++)
  {
    if (!cTree->GetEvent(i)) continue;

    TClonesArray *cl = clrow->GetArray();
    Int_t ncl = cl->GetEntriesFast();

    while (ncl--)
    {

      AliTPCclusterMI* clusterMi = (AliTPCclusterMI*) cl->At(ncl);

      AliCluster *c = (AliCluster*) cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      if (g[0]*g[0]+g[1]*g[1] < maxRsqr)
      {
        cc->Fill(g[0], g[1], g[2], clusterMi->GetQ());
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

  clusters->SetName("TPC Clusters");

  clusters->SetTitle(Form("N=%d", clusters->Size()));

  const TString viz_tag("REC Clusters TPC");

  clusters->ApplyVizTag(viz_tag, "Clusters");

//  clusters->SetRnrSelf(kFALSE);
//  clusters->SetRnrChildren(kFALSE);    

//  gEve->AddElement(clusters, cont);

  cc->SetRnrSelf(kTRUE);

  gEve->AddElement(cc);

  gEve->Redraw3D();

  return clusters;
}
