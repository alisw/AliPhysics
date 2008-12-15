// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void geom_zdc()
{
  gGeoManager = gEve->GetDefaultGeometry();

  TEveElementList* list = new TEveElementList("ZDC");
  gEve->AddGlobalElement(list);

  TGeoNode* node;
  TEveGeoTopNode* re;

  node = gGeoManager->GetTopVolume()->FindNode("ZDCA_1");
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  node = gGeoManager->GetTopVolume()->FindNode("ZDCC_1");
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  gEve->Redraw3D();
}
