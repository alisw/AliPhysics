#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoShape.h>
#include <TEveGeoShapeExtract.h>

#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

void geom_gentle_default(Bool_t register_as_global=kTRUE)
{
  TEveGeoShape* gsre1;
  TEveGeoShape* gsre2;
  TEveGeoShape* gsre3;
  
  // Geometry 3D
  {
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
    gsre1 = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();

    if (register_as_global)
    {
      gEve->AddGlobalElement(gsre1);
    }
  }
  
  // Geometry rphi
  {
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rphi_geo.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
    gsre2 = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();
  }
  
  // Geometry rhoz
  {
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
    gsre3 = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();
  }
  
  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *mv = AliEveMultiView::Instance();
  
  mv->InitGeomGentle(gsre1, gsre2, gsre3, 0);

  gEve->FullRedraw3D(kTRUE, kTRUE);

}
