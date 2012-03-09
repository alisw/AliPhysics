// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveElement.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#include <AliEveTPCData.h>
#include <AliEveTPCSector2D.h>
#include <AliEveTPCSector3D.h>
#endif

void tpc_digits(Int_t mode=1)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("TPC");
  TTree* dt = rl->GetTreeD("TPC", false);
  if (dt == 0)
  {
    throw TEveException("tpc_digits Can not access digits tree.");
  }

  AliEveTPCData *x = new AliEveTPCData;

  x->LoadDigits(dt, kTRUE); // Create all present sectors.

  gStyle->SetPalette(1, 0);
  Color_t col = 36;

  switch(mode) {

  case 0: { // Display a single sector
    AliEveTPCSector2D* s = new AliEveTPCSector2D();
    s->SetFrameColor(col);
    s->SetDataSource(x);
    gEve->AddElement(s);
    gEve->Redraw3D();

    //TGLViewer* cam = gEve->GetDefaultGLViewer();
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
	s->SetDataSource(x);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	l->AddElement(s);
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
	s->SetDataSource(x);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	l->AddElement(s);
      }
    }
    gEve->EnableRedraw();

    break;
  }

  case 2 : { // Display a single sector in 3D
    AliEveTPCSector3D* s = new AliEveTPCSector3D();
    s->SetFrameColor(col);
    s->SetDataSource(x);
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
	s->SetDataSource(x);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	l->AddElement(s);
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
	s->SetDataSource(x);
	s->SetFrameColor(col);
	s->SetAutoTrans(kTRUE);
	l->AddElement(s);
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

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("TPC");
  TTree* dt = rl->GetTreeD("TPC", false);
  if (dt == 0)
  {
    throw TEveException("tpc_digits Can not access digits tree.");
  }

  AliEveTPCData *x = new AliEveTPCData;

  x->LoadDigits(dt, kTRUE); // Create all present sectors.

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
      s->SetDataSource(x);
      s->SetFrameColor(col);
      s->SetAutoTrans(kTRUE);
      gEve->AddElement(s, l);
    }
  }
  gEve->EnableRedraw();
}
