// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void geom_ddip()
{
  gGeoManager = gEve->GetDefaultGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("DDIP_1");

  TEveGeoTopNode* re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re);
  gEve->Redraw3D();
}
