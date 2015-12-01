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

void geom_t0()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  TEveElementList* list = new TEveElementList("T0");
  gEve->AddGlobalElement(list);

  TGeoNode* node;
  TEveGeoTopNode* re;

  node = gGeoManager->GetTopVolume()->FindNode("0STR_1");
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  node = gGeoManager->GetTopVolume()->FindNode("0STL_1");
  re = new TEveGeoTopNode(gGeoManager, node);
  re->UseNodeTrans();
  gEve->AddGlobalElement(re, list);

  gEve->Redraw3D();
}
