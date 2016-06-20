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

void geom_gentle_tpc(Bool_t register_as_global=kTRUE)
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
    its1->SetMainColor(kRed);
    k++;

    TEveGeoShape* its2 = (TEveGeoShape*) *k;
    its2->SetRnrSelf(kTRUE);
    its2->SetMainColor(kRed+2);
    k++;

    TEveGeoShape* its3 = (TEveGeoShape*) *k;
    its3->SetRnrSelf(kTRUE);
    its3->SetMainColor(kRed+4);
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
    lvl3->SetMainColor(kBlue);
    TEveElement::List_i l = lvl3->BeginChildren();

    TEveGeoShape* lvl4 = (TEveGeoShape*) *l;
    lvl4->SetRnrSelf(kFALSE);
    TEveElement::List_i m = lvl4->BeginChildren();

    TEveGeoShape* tpc1 = (TEveGeoShape*) *m;
    tpc1->SetRnrSelf(kTRUE);
    tpc1->SetMainColor(kGray);
    m++;

    TEveGeoShape* tpc2 = (TEveGeoShape*) *m;
    tpc2->SetRnrSelf(kTRUE);
    tpc2->SetMainColor(kGray);
    m++;

    TEveGeoShape* tpc3 = (TEveGeoShape*) *m;
    tpc3->SetRnrSelf(kTRUE);
    tpc3->SetMainColor(kGray);
    m++;
  }
//TRD+TOF

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kFALSE);

      }
  }
//PHOS

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kFALSE);
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
        lvl2->SetRnrSelf(kFALSE);
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
    its1->SetMainColor(kRed);
    k++;

    TEveGeoShape* its2 = (TEveGeoShape*) *k;
    its2->SetRnrSelf(kTRUE);
    its2->SetMainColor(kRed+2);
    k++;

    TEveGeoShape* its3 = (TEveGeoShape*) *k;
    its3->SetRnrSelf(kTRUE);
    its3->SetMainColor(kRed+4);
  }
//TPC

  i++;
  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kFALSE);

      }
  }
//PHOS

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kFALSE);
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
        lvl2->SetRnrSelf(kFALSE);
      }
  }
  
}

{
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  gsre3 = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  // Fix visibility, color and transparency

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
    its1->SetMainColor(kRed);
    k++;

    TEveGeoShape* its2 = (TEveGeoShape*) *k;
    its2->SetRnrSelf(kTRUE);
    its2->SetMainColor(kRed+2);
    k++;

    TEveGeoShape* its3 = (TEveGeoShape*) *k;
    its3->SetRnrSelf(kTRUE);
    its3->SetMainColor(kRed+4);
  }
//TPC

  i++;
  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kFALSE);

      }
  }

//PHOS

  i++;
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);

    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); j++)
      {
        TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
        lvl2->SetRnrSelf(kFALSE);
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
        lvl2->SetRnrSelf(kFALSE);
      }
  }
}

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *mv = AliEveMultiView::Instance();

  mv->InitGeomGentle(gsre1, gsre2, gsre3, 0);

  gEve->FullRedraw3D(kTRUE, kTRUE);   

}
