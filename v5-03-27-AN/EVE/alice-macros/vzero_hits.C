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
#include <TEveManager.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

TEvePointSet*
vzero_hits(const char *varexp    = "fX:fY:fZ",
	   const char *selection = "")
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("VZERO");

  TTree* ht = rl->GetTreeH("VZERO", false);

  TEvePointSet* points = new TEvePointSet(Form("SIM Hits VZERO '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  points->SetName(Form("VZERO Hits"));
  const TString viz_tag("SIM Hits VZERO");
  points->ApplyVizTag(viz_tag, "Hits");

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
