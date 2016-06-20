#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TGLViewer.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoShape.h>
#include <TEveGeoShapeExtract.h>

#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

void geom_gentle_projector(Bool_t register_as_global=kTRUE)
{
  TEveGeoShape* gsre1;
  TEveGeoShape* gsre2;
  TEveGeoShape* gsre3;
  
  {
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  gsre1 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  if (register_as_global)
  {
    gEve->AddGlobalElement(gsre1);
  }

  // Fix visibility, color and transparency

  gsre1->SetRnrSelf(kFALSE);
  TEveElement::List_i i = gsre1->BeginChildren();

//ITS
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    TEveElement::List_i j = lvl1->BeginChildren();

    TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
    lvl2->SetRnrSelf(kFALSE);
    TEveElement::List_i k = lvl2->BeginChildren();

    TEveGeoShape* its1 = (TEveGeoShape*) *k;
    its1->SetRnrSelf(kTRUE);
    its1->SetMainColor(kBlue-6);
    its1->SetMainTransparency(50);
    k++;

    TEveGeoShape* its2 = (TEveGeoShape*) *k;
    its2->SetRnrSelf(kTRUE);
    its2->SetMainTransparency(50);
    k++;

    TEveGeoShape* its3 = (TEveGeoShape*) *k;
    its3->SetRnrSelf(kTRUE);
    its1->SetMainColor(kGreen-10);
    its3->SetMainTransparency(50);
  }
//TPC

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    TEveElement::List_i j = lvl1->BeginChildren();

    TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
    lvl2->SetRnrSelf(kFALSE);
    TEveElement::List_i k = lvl2->BeginChildren();

    TEveGeoShape* lvl3 = (TEveGeoShape*) *k;
    lvl3->SetRnrSelf(kTRUE);
    lvl3->SetMainTransparency(50);
    TEveElement::List_i l = lvl3->BeginChildren();

    TEveGeoShape* lvl4 = (TEveGeoShape*) *l;
    lvl4->SetRnrSelf(kFALSE);
    TEveElement::List_i m = lvl4->BeginChildren();

    TEveGeoShape* tpc1 = (TEveGeoShape*) *m;
    tpc1->SetRnrSelf(kTRUE);
    tpc1->SetMainTransparency(50);
    m++;

    TEveGeoShape* tpc2 = (TEveGeoShape*) *m;
    tpc2->SetMainColor(kGray);
    tpc2->SetMainTransparency(50);
    m++;

    TEveGeoShape* tpc3 = (TEveGeoShape*) *m;
    tpc3->SetRnrSelf(kTRUE);
    tpc3->SetMainTransparency(50);
    m++;
  }
//TRD+TOF

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    TEveElement::List_i j = lvl1->BeginChildren();

    TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
    lvl2->SetRnrSelf(kTRUE);
    lvl2->SetMainTransparency(50);
    j++;

    TEveGeoShape* lvl3 = (TEveGeoShape*) *j;
    lvl3->SetRnrSelf(kTRUE);
    lvl3->SetMainTransparency(50);
    j++;
  }
//PHOS

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kTRUE);
        lvl2->SetMainTransparency(50);
      }
  }
//HMPID

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kTRUE);
        lvl2->SetMainTransparency(50);
      }
  }
}

{
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rphi_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  gsre2 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  // Fix visibility, color and transparency

  gsre2->SetRnrSelf(kFALSE);
  TEveElement::List_i i = gsre2->BeginChildren();

//ITS
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    TEveElement::List_i j = lvl1->BeginChildren();

    TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
    lvl2->SetRnrSelf(kFALSE);
    TEveElement::List_i k = lvl2->BeginChildren();

    TEveGeoShape* its1 = (TEveGeoShape*) *k;
    its1->SetRnrSelf(kTRUE);
    its1->SetMainColor(kGreen-10);
        its1->SetMainTransparency(70);
    k++;

    TEveGeoShape* its2 = (TEveGeoShape*) *k;
    its2->SetRnrSelf(kTRUE);
    its2->SetMainTransparency(70);
    k++;

    TEveGeoShape* its3 = (TEveGeoShape*) *k;
    its3->SetRnrSelf(kTRUE);
    its3->SetMainColor(kBlue-6);
    its3->SetMainTransparency(70);
  }
}

{
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  gsre3 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  gsre3->SetRnrSelf(kFALSE);
  TEveElement::List_i i = gsre3->BeginChildren();

//ITS
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    TEveElement::List_i j = lvl1->BeginChildren();

    TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
    lvl2->SetRnrSelf(kFALSE);
    TEveElement::List_i k = lvl2->BeginChildren();

    TEveGeoShape* its1 = (TEveGeoShape*) *k;
    its1->SetRnrSelf(kTRUE);
        its1->SetMainTransparency(70);
    its1->SetMainColor(kGreen-10);
    k++;

    TEveGeoShape* its2 = (TEveGeoShape*) *k;
    its2->SetRnrSelf(kTRUE);
        its2->SetMainTransparency(70);
  //  its2->SetMainColor(kRed+2);
    k++;

    TEveGeoShape* its3 = (TEveGeoShape*) *k;
    its3->SetRnrSelf(kTRUE);
    its3->SetMainColor(kBlue-6);
    its3->SetMainTransparency(70);
  }
}
  

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *mv = AliEveMultiView::Instance();

  mv->InitGeomGentle(gsre1, gsre2, gsre3, 0);

  gEve->FullRedraw3D(kTRUE, kTRUE);

}

