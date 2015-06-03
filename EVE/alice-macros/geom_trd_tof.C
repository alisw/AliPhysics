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
#include <TEveGeoNode.h>

#include <AliEveEventManager.h>
#endif

void geom_trd_tof()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("B077_1");

	if(!node)
		return;

  TEveGeoTopNode* its_re = new TEveGeoTopNode(gGeoManager, node);
  gEve->AddGlobalElement(its_re);
  gEve->Redraw3D();
}
