// $Id: its_clusters.C 55060 2012-03-09 18:13:17Z quark $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TBranch.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include "AliRunLoader.h"
#include "AliCluster.h"
#include "AliEveEventManager.h"
#include "AliGeomManager.h"
#include "../ITS/UPGRADE/AliITSUGeomTGeo.h"
#include "../ITS/UPGRADE/AliITSUSegmentationPix.h"
#include "../ITS/UPGRADE/AliITSUClusterPix.h"
#include "TGeoManager.h"

#else
class TEveElement;
class TEvePointSet;
class TTree;
class TBranch;
#endif

void itsU_clusters(TEveElement* cont=0, Float_t maxR=50)
{
  AliEveEventManager::Instance()->AssertGeometry();

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadRecPoints("ITS");

  gGeoManager = gEve->GetGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE,kTRUE);
  AliITSUClusterPix::SetGeom(gm);
  TTree *cTree = rl->GetTreeR("ITS", false);
  if (cTree == 0)
    return ;

  TObjArray layerClus;
  int nlr=0;
  while(1) {
    TBranch* br = cTree->GetBranch(Form("ITSRecPoints%d",nlr));
    if (!br) break;
    TClonesArray* clr = 0;
    br->SetAddress(&clr);
    layerClus.AddLast(clr);
    nlr++;
  }

  TEveElementList* evClusters = new TEveElementList("ITS clusters");
  TEvePointSet* layClusters = 0;
  //  layClusters->SetOwnIds(kTRUE);

  Int_t layOld =-1;
  cTree->GetEntry(0);
  for (int ilr=0;ilr<nlr;ilr++) {
    TClonesArray* clr = (TClonesArray*)layerClus.At(ilr);
    int ncl = clr->GetEntries();
    //      AliITSUSegmentationPix* segm = (AliITSUSegmentationPix*)fGM->GetSegmentation(ilr);
    Float_t maxRsqr = maxR*maxR;
    for (Int_t icl = 0; icl < ncl; ++icl) {
      AliITSUClusterPix *c = (AliITSUClusterPix*) clr->UncheckedAt(icl);
      Int_t id = c->GetVolumeId();
      int lay,sta,ssta,mod,chip;
      gm->GetChipId(id, lay,sta,ssta,mod,chip);
      
      if (lay!=layOld) { // assumes order in the digits !
	layOld=lay;
	layClusters = new TEvePointSet(Form("ITS clusters - Layer %d",lay),10000);
	//  layClusters->ApplyVizTag(viz_tag, "Clusters");
	layClusters->SetMarkerColor(kBlue);
	layClusters->SetMarkerStyle(21);//kFullDotMedium);
	layClusters->SetRnrSelf(kTRUE);
	layClusters->SetOwnIds(kTRUE);
	
	evClusters->AddElement(layClusters);   
      }
      
      float g[3]; 
      c->GetGlobalXYZ(g);
      //      gm->LocalToGlobal(mod,l,g);
      if (g[0]*g[0] + g[1]*g[1] < maxRsqr) {
	layClusters->SetNextPoint(g[0], g[1], g[2]);
	AliITSUClusterPix *atp = new AliITSUClusterPix(*c);
	layClusters->SetPointId(atp);
//	printf("%d: mod %d: loc(%.4lf,%.4lf,%.4lf); glob(%.4lf,%.4lf,%.4lf); \n",
//		       icl,atp->GetVolumeId(), l[0],l[1],l[2], g[0],g[1],g[2]);
      }
    }
  }
  
  rl->UnloadRecPoints("ITS");
  
  //  clusters->SetTitle(Form("N=%d", clusters->Size()));
  
  const TString viz_tag("REC Clusters ITS");
  
  evClusters->ApplyVizTag(viz_tag, "Clusters");
  
  gEve->AddElement(evClusters, cont);
  
  gEve->Redraw3D();

}
