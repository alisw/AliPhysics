#include <iostream>
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TEveGeoNode.h"
#include "TEveGeoShape.h"
#include "TEveManager.h"
#include "TFile.h"

#include "TGLViewer.h"

#include "AliEveHltGeometry.h"
ClassImp(AliEveHltGeometry)

using namespace std;

///_______________________________________________________________________
AliEveHltGeometry::AliEveHltGeometry() {

}



///_______________________________________________________________________
AliEveHltGeometry::~AliEveHltGeometry() {
  // see header file for class documentation
  
}

///______________________________________________________________________
TEveGeoShape * AliEveHltGeometry::CreateGentleGeometry( Bool_t register_as_global) {
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

///______________________________________________________________________
TEveGeoTopNode * AliEveHltGeometry::CreateEmcalGeometry(TGeoManager * manager) {
  if(gGeoManager) cout << "have the manager"<<endl;
  else cout << "we don't "<<endl;
  TGeoVolume * volume = manager->GetTopVolume();
  if(!volume) cout << "no volume"<<endl;
  TGeoNode * gEMCALNode = manager->GetTopVolume()->FindNode("XEN1_1");
  
  TEveGeoTopNode* emcal_re = new TEveGeoTopNode(gGeoManager, gEMCALNode);
  emcal_re->SetVisLevel(1);

  for(Int_t i = 4; i < 11; i++) {
    emcal_re->FindChild(Form("SMOD_%d", i))->SetRnrState(kFALSE);
  }
  emcal_re->FindChild("SM10_1")->SetRnrState(kFALSE);
  emcal_re->FindChild("SM10_2")->SetRnrState(kFALSE);


  gEve->AddGlobalElement(emcal_re);

  return emcal_re;
}


// -----------------------------------------------------------------
TEveGeoShape* AliEveHltGeometry::geom_gentle_rphi() {
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
TEveGeoShape* AliEveHltGeometry::geom_gentle_rhoz() {
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  return gsre;
}

// -----------------------------------------------------------------
TEveGeoShape* AliEveHltGeometry::geom_gentle_trd() {
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

///______________________________________________________________________
TEveGeoShape* AliEveHltGeometry::geom_gentle_muon(Bool_t updateScene) {

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
///______________________________________________________________________
void AliEveHltGeometry::DrawDeep(TEveGeoShape *gsre) {
  
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
