#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoMatrix.h>
#include <TStyle.h>
#include <TEveRGBAPalette.h>
#include <TEveGeoShape.h>
#include <TEveStraightLineSet.h>
#include <TEveManager.h>
#include <TEveElement.h>

#include <AliESDEvent.h>
#include <AliESDZDC.h>
#include <AliEveEventManager.h>
#endif

TEveRGBAPalette *g_zdc_palette = 0;
Float_t          g_zdc_scale   = 0.1;
Float_t          g_zdc_dist    = 11695;
Float_t          g_zdc_cross   = 20;

TEveGeoShape* zdc_make_shape(const Text_t* name, const Text_t* title_base,
			     Double_t signal,    const Text_t* path)
{
    return 0;
        printf("\n\nesd_zdc_make_shape\n\n");
    
  if ( ! gGeoManager->cd(path))
  {
    Warning("zdc_make_shape", "Module name=%s, path='%s' not found.\n", name, path);
    return 0;
  }

  UChar_t rgb[3];
  g_zdc_palette->ColorFromValue(TMath::Nint(signal), rgb, kFALSE);

  TGeoShape *gs = (TGeoShape*) gGeoManager->GetCurrentVolume()->GetShape()->Clone();

  TEveGeoShape *s = new TEveGeoShape(name, Form("%s %s, E=%.3f", title_base, name, signal));
  s->SetPickable(kTRUE);
  s->SetMainColorRGB(rgb[0], rgb[1], rgb[2]);
  s->SetShape(gs);
  s->RefMainTrans().SetFrom(*gGeoManager->GetCurrentMatrix());

  // Scale z-dictance by 0.1
  Double_t* t = s->RefMainTrans().ArrT();
  t[2] *= g_zdc_scale;

  return s;
}

TEveStraightLineSet*
zdc_make_cross(const Text_t* name, const Text_t* title_base,
               Float_t x, Float_t y, Float_t z, Float_t dx, Float_t dy, Float_t dz)
{
    printf("\n\nesd_zdc_make_cross\n\n");
    
  TEveStraightLineSet* ls = new TEveStraightLineSet(name, Form("%s, x=%.3f, y=%.3f", title_base, x, y));
  ls->SetMainColor(kYellow);
  ls->SetLineWidth(2);
  ls->RefMainTrans().SetPos(x, y, z);
  ls->AddLine(dx, 0,  0, -dx, 0,  0);
  ls->AddLine(0,  dy, 0,  0, -dy, 0);
  ls->AddLine(0,  0,  dz, 0,  0, -dz);

  return ls;
}

// ???? There are 5 towers in ESD, 4 in geom.
// Not sure about assignment A/C <-> 1/2

void esd_zdc() 
{
    printf("\n\nesd_zdc\n\n");
    
  AliEveEventManager::Instance()->AssertGeometry();

  AliESDZDC *esd = AliEveEventManager::Instance()->AssertESD()->GetESDZDC();

  if (g_zdc_palette == 0)
  {
    // Map values from 0, 50 on a spectrum palette.
    g_zdc_palette = new TEveRGBAPalette(0, 50, kTRUE, kFALSE);
    g_zdc_palette->IncRefCount();
    gStyle->SetPalette(1, 0);
    g_zdc_palette->SetupColorArray();
  }

  TEveElementList* l = new TEveElementList("ZDC Data", "");
  gEve->AddElement(l);

  TEveElementList *c  = 0;
  const Double_t  *te = 0;
  const Text_t    *tb = 0;

  // ZNC geometry ------------------------------------
  tb = "ZNC";
  c = new TEveElementList(tb);
  l->AddElement(c);

  te = esd->GetZN1TowerEnergy();

  c->AddElement(zdc_make_shape("Tower 1", tb, te[1],
			       "ALIC_1/ZDCC_1/ZNEU_1/ZNTX_1/ZN1_1"));

  c->AddElement(zdc_make_shape("Tower 2", tb, te[2],
			       "ALIC_1/ZDCC_1/ZNEU_1/ZNTX_1/ZN1_2"));

  c->AddElement(zdc_make_shape("Tower 3", tb, te[3],
			       "ALIC_1/ZDCC_1/ZNEU_1/ZNTX_2/ZN1_1"));

  c->AddElement(zdc_make_shape("Tower 4", tb, te[4],
			       "ALIC_1/ZDCC_1/ZNEU_1/ZNTX_2/ZN1_2"));

  // ZNA geometry
  tb = "ZNA";
  c = new TEveElementList(tb);
  l->AddElement(c);

  te = esd->GetZN2TowerEnergy();

  c->AddElement(zdc_make_shape("Tower 1", tb, te[1],
			       "ALIC_1/ZDCA_1/ZNEU_2/ZNTX_1/ZN1_1"));

  c->AddElement(zdc_make_shape("Tower 2", tb, te[2],
			       "ALIC_1/ZDCA_1/ZNEU_2/ZNTX_1/ZN1_2"));

  c->AddElement(zdc_make_shape("Tower 3", tb, te[3],
			       "ALIC_1/ZDCA_1/ZNEU_2/ZNTX_2/ZN1_1"));

  c->AddElement(zdc_make_shape("Tower 4", tb, te[4],
			       "ALIC_1/ZDCA_1/ZNEU_2/ZNTX_2/ZN1_2"));


  // ZPC geometry ------------------------------------
  tb = "ZPC";
  c = new TEveElementList(tb);
  l->AddElement(c);

  te = esd->GetZP1TowerEnergy();

  c->AddElement(zdc_make_shape("Tower 1", tb, te[1],
			       "ALIC_1/ZDCC_1/ZPRO_1/ZPTX_1/ZP1_1"));

  c->AddElement(zdc_make_shape("Tower 2", tb, te[2],
			       "ALIC_1/ZDCC_1/ZPRO_1/ZPTX_2/ZP1_1"));

  c->AddElement(zdc_make_shape("Tower 3", tb, te[3],
			       "ALIC_1/ZDCC_1/ZPRO_1/ZPTX_3/ZP1_1"));

  c->AddElement(zdc_make_shape("Tower 4", tb, te[4],
			       "ALIC_1/ZDCC_1/ZPRO_1/ZPTX_4/ZP1_1"));

  // ZPA geometry
  tb = "ZPA";
  c = new TEveElementList(tb);
  l->AddElement(c);

  te = esd->GetZP2TowerEnergy();

  c->AddElement(zdc_make_shape("Tower 1", tb, te[1],
			       "ALIC_1/ZDCA_1/ZPRO_2/ZPTX_1/ZP1_1"));

  c->AddElement(zdc_make_shape("Tower 2", tb, te[2],
			       "ALIC_1/ZDCA_1/ZPRO_2/ZPTX_2/ZP1_1"));

  c->AddElement(zdc_make_shape("Tower 3", tb, te[3],
			       "ALIC_1/ZDCA_1/ZPRO_2/ZPTX_3/ZP1_1"));

  c->AddElement(zdc_make_shape("Tower 4", tb, te[4],
			       "ALIC_1/ZDCA_1/ZPRO_2/ZPTX_4/ZP1_1"));


  gEve->Redraw3D();
}
