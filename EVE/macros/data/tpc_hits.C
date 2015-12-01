// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include <TTree.h>
#include <TEvePointSet.h>
#include <TEveElement.h>
#include <TEveManager.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

TEvePointSet*
tpc_hits(const char  *varexp    = "TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
	 const char  *selection = "TPC2.fArray.fR>80",
         TEveElement *cont      = 0)
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);

  TEvePointSet* points = new TEvePointSet(Form("TPC Hits '%s'", selection));
  points->SetSourceCS(TEvePointSelectorConsumer::kTVT_RPhiZ);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if (points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("tpc_hits", "No hits match '%s'", selection);
    delete points;
    return 0;
  }
  
  points->SetName(Form("TPC Hits"));
  const TString viz_tag("SIM Hits TPC");
  points->ApplyVizTag(viz_tag, "Hits");

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(3);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
