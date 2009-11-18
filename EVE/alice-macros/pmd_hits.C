// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEvePointSet*
pmd_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "")
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("PMD");

  TTree* ht = rl->GetTreeH("PMD", false);

  TEvePointSet* points = new TEvePointSet(Form("PMD Hits '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  // PD - added tags
  
  points->SetName(Form("PMD Hits"));
  const TString viz_tag("SIM Hits PMD");
  points->ApplyVizTag(viz_tag, "Hits");

  // PD

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
