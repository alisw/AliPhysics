// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEvePointSet*
tof_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "",
	 TEveElement* cont)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("TOF");

  TTree* ht = rl->GetTreeH("TOF", false);
  
  TEvePointSet* points = new TEvePointSet(Form("TOF Hits '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if( points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("tof_hits", Form("No hits match '%s'", selection));
    delete points;
    return 0;
  }

 // PD - added tags
  
  points->SetName(Form("TOF Hits"));
  const TString viz_tag("SIM Hits TOF");
  points->ApplyVizTag(viz_tag, "Hits");

  // PD

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
