// $Id: geom_acorde.C 55329 2012-03-23 19:29:00Z quark $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoNode.h>

#include <AliEveEventManager.h>
#endif

void geom_ad()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("AD_1");
  if (!node) {
    Warning("geom_ad()", "Node AD_1 not found.");
    return;
  }

  TEveElementList* list = new TEveElementList("AD");
  gEve->AddGlobalElement(list);

  TEveGeoTopNode* re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  gEve->Redraw3D();
}
