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

void geom_mft()

{
  //  AliEveEventManager::GetMaster()->AssertGeometry();

  gGeoManager = gEve->GetGeometry("geometry.root");

  TEveElementList* list = new TEveElementList("MFT");
  gEve->AddGlobalElement(list);
 
  TGeoNode *node1 = gGeoManager->GetTopVolume()->FindNode("MFT_0");
  if (!node1) {
    Warning("geom_mft()", "Node MFT_0 not found !");
    return;
  }
  
  TEveGeoTopNode* re1 = new TEveGeoTopNode(gGeoManager,node1);
  re1->UseNodeTrans();
  gEve->AddGlobalElement(re1,list);
  
  gEve->Redraw3D();
  
  Info("geom_mft.C", "Done");
}
