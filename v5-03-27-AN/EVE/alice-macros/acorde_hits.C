// $Id$
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
acorde_hits(const char  *varexp    = "ACORDE.fX:ACORDE.fY:ACORDE.fZ",
	    const char  *selection = "",
	    TEveElement *cont      = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("ACORDE");

  TTree* ht = rl->GetTreeH("ACORDE", false);

  TEvePointSet* points = new TEvePointSet(Form("ACORDE Hits '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if(points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE)
  {
    Warning("acorde_hits", "No hits match '%s'", selection);
    delete points;
    return 0;
  }

  points->SetName(Form("ACORDE Hits"));
  const TString viz_tag("SIM Hits ACORDE");
  points->ApplyVizTag(viz_tag, "Hits");

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
