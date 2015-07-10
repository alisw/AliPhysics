// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>

#include <AliEveEventManager.h>
#endif

void geom_its()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("ITSV_1");

  gEve->AddGlobalElement(new TEveGeoTopNode(gGeoManager, node));
  gEve->Redraw3D();
}

void geom_its_spd()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("ITSV_1");
  node = node->GetVolume()->FindNode("ITSSPD_1");

  gEve->AddGlobalElement(new TEveGeoTopNode(gGeoManager, node));

  gEve->Redraw3D();
}

void geom_its_sdd()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  //TEveGeoTopNode *its_re;
  TGeoNode       *n1, *n2;

  n1 = gGeoManager->GetTopVolume()->FindNode("ITSV_1");

  n2 = n1->GetVolume()->FindNode("ITSsddLayer3_1");
  gEve->AddGlobalElement(new TEveGeoTopNode(gGeoManager, n2));

  n2 = n1->GetVolume()->FindNode("ITSsddLayer4_1");
  gEve->AddGlobalElement(new TEveGeoTopNode(gGeoManager, n2));

  gEve->Redraw3D();
}

void geom_its_ssd()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  //TEveGeoTopNode *its_re;
  TGeoNode       *n1, *n2;

  n1 = gGeoManager->GetTopVolume()->FindNode("ITSV_1");

  n2 = n1->GetVolume()->FindNode("ITSssdLayer5_1");
  gEve->AddGlobalElement(new TEveGeoTopNode(gGeoManager, n2));

  n2 = n1->GetVolume()->FindNode("ITSssdLayer6_1");
  gEve->AddGlobalElement(new TEveGeoTopNode(gGeoManager, n2));

  gEve->Redraw3D();
}

void geom_its_dets()
{
  geom_its_spd();
  geom_its_sdd();
  geom_its_ssd();
}
