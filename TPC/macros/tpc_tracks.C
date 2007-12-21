namespace Reve
{
class NLTProjector;
class GeoShapeRnrEl;
class RnrElement;
class RnrElementList;
}

Reve::NLTProjector  * proj = 0;
Reve::GeoShapeRnrEl * geom = 0;

void tpc_tracks()
{
  Reve::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  Reve::LoadMacro("geom_gentle.C");

  Reve::LoadMacro("primary_vertex.C");
  Reve::LoadMacro("esd_tracks.C");
  Reve::LoadMacro("its_clusters.C+");
  Reve::LoadMacro("tpc_clusters.C+");

  Reve::Viewer* nv = gReve->SpawnNewViewer("NLT Projected");
  Reve::Scene*  ns = gReve->SpawnNewScene("NLT"); 
  nv->AddScene(ns);
  TGLViewer* v = nv->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);

  Reve::TrackCounter* g_trkcnt = new Reve::TrackCounter("Primary Counter");
  gReve->AddToListTree(g_trkcnt, kFALSE);

  Reve::NLTProjector* p = new Reve::NLTProjector; proj = p;
  gReve->AddToListTree(p, kTRUE);
  gReve->AddRenderElement(proj, ns);

  // geometry
  Reve::GeoShapeRnrEl* gg = geom_gentle();
  geom = gg;

  // event
  Alieve::gEvent->AddNewEventCommand("on_new_event();");
  Alieve::gEvent->GotoEvent(0);

  gReve->Redraw3D(kTRUE);
}

/**************************************************************************/

void on_new_event()
{
  try {
    //Reve::PointSet* itsc = its_clusters();
    //itsc->SetMarkerColor(5);

    Reve::PointSet* tpcc = tpc_clusters();
    tpcc->SetMarkerColor(4);
  }
  catch(Reve::Exc_t& exc) {
    printf("Exception loading ITS/TPC clusters: %s\n", exc.Data());
  }

  Reve::TrackList* cont = esd_tracks();
  cont->SetLineWidth((Width_t)2);

  Reve::RenderElement* top = gReve->GetCurrentEvent();
  proj->DestroyElements();
  //AliESDEvent* esd = Alieve::Event::AssertESD();

  // geom
  proj->ImportElements(geom);
  // event
  proj->ImportElements(top);
  // top->SetRnrState(kFALSE);
}
