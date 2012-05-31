// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TEvePointSet.h>
#include <TEveElement.h>
#include <TEveManager.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

TEvePointSet*
its_hits(const char *varexp    = "fX:fY:fZ",
	 const char *selection = "",
         TEveElement* cont = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);

  TEvePointSet* points = new TEvePointSet(Form("ITS Hits '%s'", selection));

  TEvePointSelector ps(ht, points, varexp, selection);
  // ps.SetSubIdExp("fTrack:fStatus");
  ps.Select();

  if(points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("its_hits", "No hits match '%s'", selection);
    delete points;
    return 0;
  }

  points->SetName(Form("ITS Hits"));
  points->SetTitle(Form("N=%d", points->Size()));

  points->ApplyVizTag("SIM Hits ITS", "Hits");

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
