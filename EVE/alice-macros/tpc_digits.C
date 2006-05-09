// $Id$

#include "TGLViewer.h"

void tpc_digits(Int_t mode = 0)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("TPC");
  TTree* dt = rl->GetTreeD("TPC", false);

  rl->CdGAFile();
  AliTPCParam* tpcpar = gDirectory->Get("75x40_100x60_150x60");

  Alieve::TPCDigitsInfo* di = new Alieve::TPCDigitsInfo();
  di->SetData(tpcpar, dt);
  // di->Dump();

  // Viewport limits.
  Float_t left, right, top, bottom;
  right  = di->fOut2Seg.fNMaxPads* di->fOut2Seg.fPadWidth;
  left   = -right;
  bottom = di->fInnSeg.fRlow;
  top    = bottom + di->fOut2Seg.fRlow +
    di->fOut2Seg.fNRows*di->fOut2Seg.fPadLength - di->fInnSeg.fRlow;

  gStyle->SetPalette(1, 0);
  UInt_t col = 36;

  if(mode == 0) {
    // create an object, so the camera will not jump
   
    gReve->DisableRedraw();
    Reve::RenderElementList* l = new Reve::RenderElementList("TPC");
    l->SetTitle("TPC Segment");
    l->SetMainColor(Color_t(col));
    TGListTreeItem *ti = gReve->AddRenderElement(l);
  
    Int_t i = 0;
    //  for(Int_t i = 0; i<7; i++) {
    Alieve::TPCSegment* s = new Alieve::TPCSegment();
    s->SetSegmentID(i);
    s->SetInfo(di);
    s->SetMainColor(Color_t(col));
    l->AddElement(s);
    gReve->AddRenderElement(ti, s);
    //  }
    gReve->DrawRenderElement(l);


    /*
    TPolyMarker3D* framebox = new TPolyMarker3D(8);
    const Float_t a = 250.0;
    framebox->SetPoint(0, a, a,a);
    framebox->SetPoint(1, a, -a,a);
    framebox->SetPoint(2, -a, -a,a);
    framebox->SetPoint(3, -a, a,a);

    framebox->SetPoint(4, a, a,-a);
    framebox->SetPoint(5, a, -a,-a);
    framebox->SetPoint(6, -a, a,-a);
    framebox->SetPoint(7, -a, -a,-a);

    framebox->Draw();
    */

    gReve->EnableRedraw();

    TGLViewer* cam = dynamic_cast<TGLViewer*>(gReve->GetCC()->GetViewer3D());
    //cam->SetCurrentCamera(TGLViewer::kCameraOrthoXOY) ;
    //cam->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 2*left, 2*right, 2*top, bottom); 
    //printf("%f %f %f %f\n", left, right, top, bottom);

  } else {

    gReve->DisableRedraw();
    {
      Reve::RenderElementList* l = new Reve::RenderElementList("TPC plate 1");
      l->SetTitle("TPC Plate");
      l->SetMainColor(Color_t(col));
      TGListTreeItem *ti = gReve->AddRenderElement(l);
      
      for(Int_t i = 0; i<18; i++) {
	Alieve::TPCSegment* s = new Alieve::TPCSegment();
	s->SetSegmentID(i);
	s->SetInfo(di);
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
	Alieve::TPCSegment* s = new Alieve::TPCSegment();
	s->SetSegmentID(i);
	s->SetMainColor(Color_t(col));
	s->SetInfo(di);
	s->SetTrans(true);
	l->AddElement(s);
	gReve->AddRenderElement(ti, s);
      }
      gReve->DrawRenderElement(l);
    }
    gReve->EnableRedraw();
  }
}
