// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveTreeTools.h>
#include <TEvePointSet.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#include <AliEveFMDLoader.h>
#endif

TEvePointSet*
fmd_hits2(const char *varexp    = "fX:fY:fZ",
	  const char *selection = "")
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("FMD");

  TTree* ht = rl->GetTreeH("FMD", false);

  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  TEvePointSet* points = new TEvePointSet(Form("FMD Hits '%s'", selection));
  char form[1000];
  sprintf(form,"FMD Hits '%s'", selection);
  TEvePointSet* points = new TEvePointSet(form);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor(2);

  gEve->AddElement(points);
  gEve->Redraw3D();

  return points;
}
