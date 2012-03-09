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

void its_hits_layer_split(const char *varexp    = "fX:fY:fZ:GetLayer()",
                          const char *selection = "")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  printf("THIS SCRIPT DOES NOT WORK.\n"
	 "GetLayer() crashes when trying to load ITS geometry.\n"
	 "Needs to be fixed together with ITS experts.\n");
  return;

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("ITS");

  TTree* ht = rl->GetTreeH("ITS", false);

  TEvePointSetArray* l = new TEvePointSetArray("ITS hits - Layer Slices", "");
  l->SetMarkerColor(2);
  l->SetMarkerStyle(2); // cross
  l->SetMarkerSize(.2);

  gEve->AddElement(l);
  l->InitBins("Layer", 6, 0.5, 6.5);

  TEvePointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gEve->Redraw3D();
}
