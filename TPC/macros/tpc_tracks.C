
class TEveProjectionManager;
class TEveGeoShape;
class TEveElement;
class TEveElementList;

TEveProjectionManager  * proj = 0;
TEveGeoShape * geom = 0;

void tpc_tracks()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  TEveUtil::LoadMacro("geom_gentle.C");

  TEveUtil::LoadMacro("primary_vertex.C");
  TEveUtil::LoadMacro("esd_tracks.C");
  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");

  TEveViewer* nv = gEve->SpawnNewViewer("NLT Projected");
  TEveScene*  ns = gEve->SpawnNewScene("NLT"); 
  nv->AddScene(ns);
  TGLViewer* v = nv->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);

  TEveTrackCounter* g_trkcnt = new TEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);

  TEveProjectionManager* p = new TEveProjectionManager; proj = p;
  gEve->AddToListTree(p, kTRUE);
  gEve->AddElement(proj, ns);

  // geometry
  TEveGeoShape* gg = geom_gentle();
  geom = gg;

  // event
  gEvent->AddNewEventCommand("on_new_event();");
  gEvent->GotoEvent(0);

  gEve->Redraw3D(kTRUE);
}

/**************************************************************************/

void on_new_event()
{
  try {
    //TEvePointSet* itsc = its_clusters();
    //itsc->SetMarkerColor(5);

    TEvePointSet* tpcc = tpc_clusters();
    tpcc->SetMarkerColor(4);
  }
  catch(TEveException& exc) {
    printf("Exception loading ITS/TPC clusters: %s\n", exc.Data());
  }

  TEveTrackList* cont = esd_tracks();
  cont->SetLineWidth((Width_t)2);

  TEveElement* top = gEve->GetCurrentEvent();
  proj->DestroyElements();
  //AliESDEvent* esd = AliEveEventManager::AssertESD();

  // geom
  proj->ImportElements(geom);
  // event
  proj->ImportElements(top);
  // top->SetRnrState(kFALSE);
}
