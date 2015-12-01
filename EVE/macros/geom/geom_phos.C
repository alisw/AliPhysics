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

void geom_phos()
{
  AliEveEventManager::GetMaster()->AssertGeometry();

  TEveElementList* list = new TEveElementList("PHOS");
  gEve->AddGlobalElement(list);

  for(Int_t i=1; i<=5; ++i) {
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(Form("PHOS_%d", i));
    char form[1000];
    sprintf(form,"PHOS_%d", i);
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(form);

    TEveGeoTopNode* re = new TEveGeoTopNode(gGeoManager, node);
    re->UseNodeTrans();
    gEve->AddGlobalElement(re, list);
  }

  gEve->Redraw3D();
}
