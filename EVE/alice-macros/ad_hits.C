// $Id: acorde_hits.C 55060 2012-03-09 18:13:17Z quark $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TString.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliEveEventManager.h>
#include <AliRunLoader.h>
#endif

TEvePointSet*
ad_hits(const char  *varexp    = "AD.fX:AD.fY:AD.fZ",
	    const char  *selection = "",
	    TEveElement *cont      = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("AD");

  TTree* ht = rl->GetTreeH("AD", false);

  TEvePointSet* points = new TEvePointSet(Form("AD Hits '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if(points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE)
  {
    Warning("ad_hits", "No hits match '%s'", selection);
    delete points;
    return 0;
  }

  points->SetName(Form("AD Hits"));
  const TString viz_tag("SIM Hits AD");
  points->ApplyVizTag(viz_tag, "Hits");

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
