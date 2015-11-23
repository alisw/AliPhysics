// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TTree.h>
#include <TEvePointSet.h>
#include <TEveElement.h>
#include <TEveManager.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

TEvePointSet*
t0_hits(const char *varexp    = "T0.fX:T0.fY:T0.fZ",
	const char *selection = "")
{
  // Extracts  T0 hits.


  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("T0");

  TTree* ht = rl->GetTreeH("T0", false);

  Int_t nTracks = ht->GetEntries();
  // printf("Found %d tracks. \n",nTracks);
  for (Int_t it = 0; it < nTracks; it++) {

    TClonesArray *hits = 0;
    ht->SetBranchAddress("T0",&hits);

    ht->GetEvent(it);
    // Int_t nHits = hits->GetEntriesFast();
    // printf("Found %d hits in track %d.\n", nHits, it);

  }
  TEvePointSet* points = new TEvePointSet(Form("T0 Hits '%s'", selection));
  points->SetSourceCS(TEvePointSelectorConsumer::kTVT_XYZ);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(3);

  points->SetName(Form("T0 Hits"));
  const TString viz_tag("SIM Hits T0");
  points->ApplyVizTag(viz_tag, "Hits");

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
