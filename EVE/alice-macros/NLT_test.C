
Reve::NLTProjector* NLT_test(Reve::RenderElement* top=0)
{
  using namespace Reve;

  Scene* s = gReve->SpawnNewScene("Projected Event");
  gReve->GetDefViewer()->AddScene(s);

  TGLViewer* v = (TGLViewer *)gReve->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);

  NLTProjector* p = new NLTProjector;
  p->SetProjection(NLTProjection::PT_RhoZ, 0.01);

  gReve->AddToListTree(p, kTRUE);
  gReve->AddRenderElement(p, s);

  top = gReve->GetCurrentEvent();
  if (top)
    p->ImportElements(top);

  gReve->Redraw3D(kTRUE);

  return p;
}
