// $Id$
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

void geom_vzero()
{
  static const TEveException kEH("geom_vzero() ");

  AliEveEventManager::GetMaster()->AssertGeometry();

  TEveElementList* list = new TEveElementList("VZero");
  gEve->AddGlobalElement(list);

  TGeoNode* node = 0;
  TEveGeoTopNode* re;

  TGeoNode* mnode = gGeoManager->GetTopVolume()->FindNode("VZERO_1");
  if (!mnode) {
    Error(kEH, "mother node not found.");
    return;
  }

  node = mnode->GetVolume()->FindNode("V0RI_1");
  printf("opofoih %p\n", node);
  if (!node) {
    Error(kEH, "V0R not found.");
    return;
  }
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  node = mnode->GetVolume()->FindNode("V0LE_1");
  if (!node) {
    Error(kEH, "V0L not found.");
    return;
  }
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  gEve->Redraw3D();
}
