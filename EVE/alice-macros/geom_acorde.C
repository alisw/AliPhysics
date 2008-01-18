// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

void geom_acorde()
{
  gGeoManager = gEve->GetGeometry("geometry.root");

  TEveElementList* list = new TEveElementList("ACORDE");
  gEve->AddGlobalElement(list);

  for (Int_t i=1; i<61; ++i)
  {
    char form[1000];
    sprintf(form, "ACORDE1_%d", i);
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode(form);
    TEveGeoTopNode* re = new TEveGeoTopNode(gGeoManager, node);
    re->UseNodeTrans();
    gEve->AddGlobalElement(re, list);
  }

  gEve->Redraw3D();
}
