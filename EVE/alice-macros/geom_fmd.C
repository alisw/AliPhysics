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
#include <TEveElement.h>

#include <AliEveEventManager.h>
#endif

void geom_fmd()
{
  gGeoManager = AliEveEventManager::GetMaster()->AssertGeometry();

  TEveElementList* list = new TEveElementList("FMD");
  gEve->AddGlobalElement(list);

  for(Int_t i=1; i<=3; ++i)
  {
    TGeoNode       *node = 0;
    TEveGeoTopNode *re   = 0;

    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH node = gGeoManager->GetTopVolume()->FindNode(Form("F%dMT_%d", i, i));
    char form[1000];
    sprintf(form,"F%dMT_%d", i, i);
    node = gGeoManager->GetTopVolume()->FindNode(form);
    re = new TEveGeoTopNode(gGeoManager, node);
    re->UseNodeTrans();
    gEve->AddGlobalElement(re, list);

    sprintf(form,"F%dMB_%d", i, i);
    node = gGeoManager->GetTopVolume()->FindNode(form);
    re = new TEveGeoTopNode(gGeoManager, node);
    re->UseNodeTrans();
    gEve->AddGlobalElement(re, list);
  }

  gEve->Redraw3D();
}
