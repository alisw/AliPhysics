// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file MUON_geom.C
///
/// \author B. Vulpescu, LPC

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoNode.h>
#endif

void MUON_geom()
{
  gGeoManager = gEve->GetGeometry("geometry.root");

  TEveElementList* list = new TEveElementList("DiMUON");
  gEve->AddGlobalElement(list);

  TGeoNode *node1 = gGeoManager->GetTopVolume()->FindNode("DDIP_1");
  TGeoNode *node2 = gGeoManager->GetTopVolume()->FindNode("YOUT1_1");
  TGeoNode *node3 = gGeoManager->GetTopVolume()->FindNode("YOUT2_1");

  TEveGeoTopNode* re1 = new TEveGeoTopNode(gGeoManager,node1);
  re1->UseNodeTrans();
  gEve->AddGlobalElement(re1,list);

  TEveGeoTopNode* re2 = new TEveGeoTopNode(gGeoManager,node2);
  re2->UseNodeTrans();
  gEve->AddGlobalElement(re2,list);

  TEveGeoTopNode* re3 = new TEveGeoTopNode(gGeoManager,node3);
  re3->UseNodeTrans();
  gEve->AddGlobalElement(re3,list);

  gEve->Redraw3D();

  Info("MUON_geom.C", "Done");

}
