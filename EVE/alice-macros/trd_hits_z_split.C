// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void trd_hits_z_split(const char *varexp    = "fX:fY:fZ:fZ",
		      const char *selection = "")
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("TRD");

  TTree* ht = rl->GetTreeH("TRD", false);

  TEvePointSetArray* l = new TEvePointSetArray("TRD hits - Z Slices", "");
  l->SetMarkerColor(7);
  l->SetMarkerStyle(20); // full circle
  l->SetMarkerSize(.5);

  gEve->AddElement(l);
  l->InitBins("Z", 20, -360, 360);

  TEvePointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gEve->Redraw3D();
}
