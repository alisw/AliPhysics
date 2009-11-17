// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEvePointSet*
its_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "",
         TEveElement* cont = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);

  TEvePointSet* points = new TEvePointSet(Form("SIM Hits ITS '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  // ps.SetSubIdExp("fTrack:fStatus");
  ps.Select();

  if(points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("its_hits", Form("No hits match '%s'", selection));
    delete points;
    return 0;
  }

  points->SetTitle(Form("N=%d", points->Size()));

  points->SetName("ITS Hits"));
  const TString viz_tag("SIM Hits ITS");
  points->ApplyVizTag(viz_tag, "Hits");

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
