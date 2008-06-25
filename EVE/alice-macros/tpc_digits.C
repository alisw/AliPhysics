// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class AliEveTPCData;
class AliEveEventManager;

AliEveTPCData      *g_tpc_data       = 0;
AliEveEventManager *g_tpc_last_event = 0;

void tpc_digits(Int_t mode=1)
{
  if (g_tpc_data == 0 || g_tpc_last_event != gAliEveEvent)
  {
    AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
    rl->LoadDigits("TPC");
    TTree* dt = rl->GetTreeD("TPC", false);

    g_tpc_data = new AliEveTPCData;
    g_tpc_data->LoadDigits(dt, kTRUE); // Create all present sectors.

    g_tpc_last_event = gAliEveEvent;
  }

  // Viewport limits.
  /*
  Float_t left, right, top, bottom;
  right  = di->fOut2Seg.fNMaxPads* di->fOut2Seg.fPadWidth;
  left   = -right;
  bottom = di->fInnSeg.fRlow;
  top    = bottom + di->fOut2Seg.fRlow +
    di->fOut2Seg.fNRows*di->fOut2Seg.fPadLength - di->fInnSeg.fRlow;
  */

  gStyle->SetPalette(1, 0);
  Color_t col = 36;

  switch(mode) {

  case 0: { // Display a single sector
    AliEveTPCSector2D* s = new AliEveTPCSector2D();
    s->SetFrameColor(col);
    s->SetDataSource(g_tpc_data);
    gEve->AddElement(s);
    gEve->Redraw3D();

    //TGLViewer* cam = gEve->GetGLViewer();
    //cam->SetCurrentCamera(TGLViewer::kCameraOrthoXOY) ;
    //cam->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 2*left, 2*right, 2*top, bottom);
    //printf("%f %f %f %f\n", left, right, top, bottom);

    break;
  }

  case 1: { // Display all sectors
    gEve->DisableRedraw();
    {
      TEveElementList* l = new TEveElementList("TPC plate 1");
      l->SetTitle("TPC Plate");
      l->SetMainColor(col);
      gEve->AddElement(l);

      for (Int_t i = 0; i<18; i++)
      {
	AliEveTPCSector2D* s = new AliEveTPCSector2D(Form("AliEveTPCSector2D %d", i));
	s->SetSectorID(i);
	s->SetDataSource(g_tpc_data);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	gEve->AddElement(s, l);
      }
    }
    {
      TEveElementList* l = new TEveElementList("TPC plate 2");
      l->SetTitle("TPC Plate");
      l->SetMainColor(col);

      gEve->AddElement(l);
      for (Int_t i = 18; i<36; i++)
      {
	AliEveTPCSector2D* s = new AliEveTPCSector2D(Form("AliEveTPCSector2D %d", i));
	s->SetSectorID(i);
	s->SetDataSource(g_tpc_data);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	gEve->AddElement(s, l);
      }
    }
    gEve->EnableRedraw();

    break;
  }

  case 2 : { // Display a single sector in 3D
    AliEveTPCSector3D* s = new AliEveTPCSector3D();
    s->SetFrameColor(col);
    s->SetDataSource(g_tpc_data);
    gEve->AddElement(s);
    gEve->Redraw3D();
    break;
  }

  case 3: { // Display all sectors in 3D
    gEve->DisableRedraw();
    {
      TEveElementList* l = new TEveElementList("TPC plate 1");
      l->SetTitle("TPC Plate");
      l->SetMainColor(col);
      gEve->AddElement(l);

      for (Int_t i = 0; i<18; i++)
      {
	AliEveTPCSector3D* s = new AliEveTPCSector3D(Form("AliEveTPCSector3D %d", i));
	s->SetSectorID(i);
	s->SetDataSource(g_tpc_data);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	gEve->AddElement(s, l);
      }
    }
    {
      TEveElementList* l = new TEveElementList("TPC plate 2");
      l->SetTitle("TPC Plate");
      l->SetMainColor(col);

      gEve->AddElement(l);
      for (Int_t i = 18; i<36; i++)
      {
	AliEveTPCSector3D* s = new AliEveTPCSector3D(Form("AliEveTPCSector3D %d", i));
	s->SetSectorID(i);
	s->SetDataSource(g_tpc_data);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	gEve->AddElement(s, l);
      }
    }
    gEve->EnableRedraw();

    break;
  }

  } // switch
}


void tpc_digits_2drange(Int_t start, Int_t end)
{
  if (start <  0)  start = 0;
  if (end   > 35)  end   = 35;

  if (g_tpc_data == 0 || g_tpc_last_event != gAliEveEvent) {
    AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
    rl->LoadDigits("TPC");
    TTree* dt = rl->GetTreeD("TPC", false);

    g_tpc_data = new AliEveTPCData;
    g_tpc_data->LoadDigits(dt, kTRUE); // Create all present sectors.

    g_tpc_last_event = gAliEveEvent;
  }

  gStyle->SetPalette(1, 0);
  Color_t col = 36;

  gEve->DisableRedraw();
  {
    TEveElementList* l = new TEveElementList("TPC sectors");
    l->SetMainColor(col);
    gEve->AddElement(l);

    for (Int_t i=start; i<=end; i++)
    {
      AliEveTPCSector2D* s = new AliEveTPCSector2D();
      s->SetSectorID(i);
      s->SetDataSource(g_tpc_data);
      s->SetFrameColor(col);
      s->SetAutoTrans(kTRUE);
      gEve->AddElement(s, l);
    }
  }
  gEve->EnableRedraw();
}
