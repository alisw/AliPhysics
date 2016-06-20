// $Id: geom_gentle.C 30976 2009-02-11 15:55:45Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TEveGeoShapeExtract.h"
#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveElement.h"
#include "TGLViewer.h"

#include "TFile.h"
#include "TStyle.h"
#endif

// -----------------------------------------------------------------
TEveGeoShape* geom_gentle_hlt(Bool_t register_as_global=kTRUE) {
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* elTRD = gsre->FindChild("TRD+TOF");
  elTRD->SetRnrState(kFALSE);

  TEveElement* elHMPID = gsre->FindChild("HMPID");
  elHMPID->SetRnrState(kFALSE);

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kTRUE);
  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
  elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);


  if (register_as_global) {
    gEve->AddGlobalElement(gsre);
  }

  return gsre;
}

// -----------------------------------------------------------------
TEveGeoShape* geom_gentle_rphi() {
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rphi_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kTRUE);
  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
  elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);

  return gsre;
}

// -----------------------------------------------------------------
TEveGeoShape* geom_gentle_rhoz() {
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  return gsre;
}

// -----------------------------------------------------------------
TEveGeoShape* geom_gentle_trd() {
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo_trd.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle TRD");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  gEve->AddGlobalElement(gsre);
  f.Close();

  Int_t sm = 0;

  // Fix visibility, color and transparency
  gsre->SetRnrSelf(kFALSE);
  for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i) {

    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); ++j) {

      TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
      
      if ( sm == 0 || sm == 1 || sm == 7 || sm == 8 || sm == 9 || sm == 10 || sm == 17 )
	lvl2->SetRnrSelf(kTRUE);	
      else 
	lvl2->SetRnrSelf(kFALSE);	

      lvl2->SetMainColor(3);
      lvl2->SetMainTransparency(80);

      ++sm;
    }
    
  }

  return gsre;
}

void DrawDeep(TEveGeoShape *gsre) {
  
  for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i) {
    TEveGeoShape* lvl = (TEveGeoShape*) *i;
    lvl->SetRnrSelf(kFALSE);
    if (!lvl->HasChildren()) {
      lvl->SetRnrSelf(kTRUE);
      lvl->SetMainColor(3);
      lvl->SetMainTransparency(50);
    }
    DrawDeep(lvl);
  }

}

TEveGeoShape* geom_gentle_muon(Bool_t updateScene = kTRUE) {

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo_muon.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  gEve->AddGlobalElement(gsre);
  f.Close();

  gsre->SetRnrSelf(kFALSE);

  DrawDeep(gsre);

  if ( updateScene ) {
    TGLViewer* v = gEve->GetDefaultGLViewer();
    v->UpdateScene();
  }

  return gsre;

}
