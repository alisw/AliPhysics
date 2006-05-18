// $Id$

#include "TGLViewer.h"

namespace Alieve {
class TPCData;
class Event;
}

Alieve::TPCData* g_tpc_data = 0;
Alieve::Event*   g_tpc_last_event = 0;

void tpc_digits(Int_t mode=0)
{
  if (g_tpc_data == 0 || g_tpc_last_event != Alieve::gEvent) {
    AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
    rl->LoadDigits("TPC");
    TTree* dt = rl->GetTreeD("TPC", false);

    g_tpc_data = new Alieve::TPCData;
    g_tpc_data->LoadDigits(dt, kTRUE); // Create all present sectors.

    g_tpc_last_event = Alieve::gEvent;
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
  UInt_t col = 36;

  switch(mode) {

  case 0: { // Display a single sector
   
    gReve->DisableRedraw();
  
    Alieve::TPCSector2D* s = new Alieve::TPCSector2D();
    s->SetDataSource(g_tpc_data);
    s->SetMainColor(Color_t(col));
    gReve->AddRenderElement(s);
    gReve->DrawRenderElement(s);

    gReve->EnableRedraw();

    TGLViewer* cam = dynamic_cast<TGLViewer*>(gReve->GetCC()->GetViewer3D());
    //cam->SetCurrentCamera(TGLViewer::kCameraOrthoXOY) ;
    //cam->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 2*left, 2*right, 2*top, bottom); 
    //printf("%f %f %f %f\n", left, right, top, bottom);

    break;
  }

  case 1: { // Display all sectors

    gReve->DisableRedraw();
    {
      Reve::RenderElementList* l = new Reve::RenderElementList("TPC plate 1");
      l->SetTitle("TPC Plate");
      l->SetMainColor(Color_t(col));
      TGListTreeItem *ti = gReve->AddRenderElement(l);
      
      for(Int_t i = 0; i<18; i++) {
	Alieve::TPCSector2D* s = new Alieve::TPCSector2D();
	s->SetSegmentID(i);
	s->SetDataSource(g_tpc_data);
	s->SetMainColor(Color_t(col));
	s->SetTrans(true);
	l->AddElement(s);
	gReve->AddRenderElement(ti, s);
      }
      gReve->DrawRenderElement(l);
    }
    {
      Reve::RenderElementList* l = new Reve::RenderElementList("TPC plate 2");
      l->SetTitle("TPC Plate");
      l->SetMainColor(Color_t(col));

      TGListTreeItem *ti = gReve->AddRenderElement(l);
      for(Int_t i = 18; i<36; i++) {
	Alieve::TPCSector2D* s = new Alieve::TPCSector2D();
	s->SetSegmentID(i);
	s->SetDataSource(g_tpc_data);
	s->SetMainColor(Color_t(col));
	s->SetTrans(true);
	l->AddElement(s);
	gReve->AddRenderElement(ti, s);
      }
      gReve->DrawRenderElement(l);
    }
    gReve->EnableRedraw();

    break;
  }

    /* // Almost ready ...
  case 2 : { // Display a single sector in 3D
    Reve::RenderElementList* l = new Reve::RenderElementList("TPC Drift");
    l->SetTitle("TPC Segment Drift");
    l->SetMainColor(Color_t(col));
    TGListTreeItem *ti = gReve->AddRenderElement(l);
  
    Alieve::TPCSector3D* = new Alieve::TPCSector3D(di, 0);
    s->SetMainColor(Color_t(col));
    l->AddElement(s);
    gReve->AddRenderElement(ti, s);
    gReve->DrawRenderElement(l);
    gReve->EnableRedraw();
    break;
  }
    */

  } // switch
}
