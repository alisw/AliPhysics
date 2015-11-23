// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TGLViewer.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoShape.h>
#include <TEveGeoShapeExtract.h>
#endif

/// \ingroup evemacros
/// \file geom_gentle_muon.C
///
/// \author B. Vulpescu, LPC; M. Tadel, CERN

void DrawDeep(TEveGeoShape *gsre) {
  
  if (gsre->HasChildren()) {
    
    gsre->SetRnrSelf(kFALSE);
    for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i) {
      TEveGeoShape* lvl = (TEveGeoShape*) *i;
      DrawDeep(lvl);
    }
    
  } else {
    
    gsre->SetRnrSelf(kTRUE);
    gsre->SetMainColor(kGray);
    gsre->SetMainTransparency(80);
    
  }
  
}

TEveGeoShape* geom_gentle_muon(Bool_t updateScene = kTRUE) {

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo_muon.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  gEve->AddGlobalElement(gsre);
  f.Close();

  DrawDeep(gsre);

  if ( updateScene ) {
    TGLViewer* v = gEve->GetDefaultGLViewer();
    v->UpdateScene();
  }

  return gsre;

}

