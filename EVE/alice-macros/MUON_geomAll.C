// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void MUON_geomAll()
{
  gGeoManager = gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");

  TEveGeoTopNode* topn_re = new TEveGeoTopNode
    (gGeoManager, gGeoManager->GetTopNode());

  gEve->AddGlobalElement(topn_re);

  gEve->Redraw3D(kTRUE);

}
