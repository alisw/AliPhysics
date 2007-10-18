// $Id$

void NLT_geo_demo(Int_t type = Reve::NLTProjection::CFishEye, Float_t distortion = 0)
{  
   Reve::LoadMacro("region_marker.C");
   region_marker();
   
   Reve::NLTProjector* pr = new Reve::NLTProjector();
   pr->SetProjection(type, distortion);

   make_geo(pr);
   TGLViewer* glv = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
   glv->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
}

/**************************************************************************/
void draw_node(const Text_t* path)
{
   gGeoManager->cd(path);
   TGeoMatrix* mx = gGeoManager->GetCurrentMatrix();
   
   Reve::GeoTopNodeRnrEl* currn_re = new Reve::GeoTopNodeRnrEl
      (gGeoManager, gGeoManager->GetCurrentNode());
   currn_re->SetGlobalTrans(mx);
   gReve->AddGlobalRenderElement(currn_re);
}

/**************************************************************************/
TBuffer3D* make_buffer(const Text_t* path)
{
  TBuffer3D* buff = 0;
  Bool_t valid = gGeoManager->cd(path);
  if(valid)  
  {
    TBuffer3D* buff = gGeoManager->GetCurrentVolume()->GetShape()->MakeBuffer3D();
    TGeoMatrix* mx = gGeoManager->GetCurrentMatrix();

    Int_t N = buff->NbPnts();   
    Double_t* pnts = buff->fPnts;   
    Double_t  master[4];
    for(Int_t i = 0; i<N; i++) 
    {
      mx->LocalToMaster(&pnts[3*i], master);
      pnts[3*i] =   master[0];
      pnts[3*i+1] = master[1];
      pnts[3*i+2] = master[2];
    }
    buff->fColor = gGeoManager->GetCurrentVolume()->GetLineColor();
  }
  return buff;
}

/**************************************************************************/
Reve::NLTPolygonSet* project_node(Reve::NLTProjector* nlt, const Text_t* path, Int_t useBP)
{
  Reve::NLTPolygonSet* ps = 0;
  TBuffer3D* buff = make_buffer(path); 
  if(buff) {
    ps = nlt->ProjectGeoShape(buff, useBP);
    if(ps)
    {
      ps->SetName(Form("NLTPolygonSet %s",path));
      ps->SetFillColor(Color_t(buff->fColor));
      ps->SetLineColor((Color_t)TColor::GetColorBright(buff->fColor));
      //ps->SetLineColor((Color_t)TColor::GetColorDark(buff->fColor));
    }
  }
  return ps;
}

/**************************************************************************/
void project_nodes(Reve::NLTProjector* nlt, const Text_t* parent_path, Int_t useBP)
{
  gGeoManager->cd(parent_path);
  TGeoNode* holder = gGeoManager->GetCurrentNode();
  holder->ls();
  TIter next_node(holder->GetNodes());
  TGeoNode* geon;

  Reve::RenderElementList* el = new Reve::RenderElementList(parent_path);
  gReve->AddGlobalRenderElement(el);
  while((geon = (TGeoNode*)next_node())) 
  {
    TGeoVolume* v = geon->GetVolume();
    if(v) {      
      TString path = Form("%s/%s",parent_path, geon->GetName());
      Reve::NLTPolygonSet* ps = project_node(nlt, path.Data(), useBP);
      if(ps) {
	ps->SetName(geon->GetName());
	gReve->AddGlobalRenderElement(el, ps);
      }
    }
  }
}

/**************************************************************************/
void project_to_pointset(Reve::NLTProjector* nlt, const Text_t* path) 
{  
  TBuffer3D* buff = make_buffer(path);
  if(buff)
  {
    PointSet* ps = new PointSet(buff->NbPnts());
    ps->SetMarkerColor(buff->fColor);
    for(Int_t i = 0; i<buff->NbPnts() ; i++)
      ps->SetPoint(i,  buff->fPnts[3*i], buff->fPnts[3*i+1], buff->fPnts[3*i+2]);
   
    PointSet* pps = nlt->ProjectPointSet(ps);
    gReve->AddGlobalRenderElement(pps);
  }
}
/**************************************************************************/
void make_geo(Reve::NLTProjector* nlt)
{
  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/simple_geo.root");
  Int_t useBuffPols = -1;
  Reve::NLTPolygonSet* ps;

  project_nodes( nlt, "/ALIC_1/ITSV_holder_1/ITSV_1", useBuffPols);
  ps = project_node( nlt, "/ALIC_1/TPC_holder_1/TPC_1/TDGN_1", useBuffPols);
  if(ps) gReve->AddGlobalRenderElement(ps);
  ps = project_node(nlt, "/ALIC_1/TRD TOF_holder_1/B077_1", useBuffPols);
  if(ps) gReve->AddGlobalRenderElement(ps);
  ps = project_node(nlt, "/ALIC_1/TRD TOF_holder_1/BRS4_1", useBuffPols);
  if(ps) gReve->AddGlobalRenderElement(ps);
  ps = project_node(nlt, "/ALIC_1/TRD TOF_holder_1/BRS4_2", useBuffPols);
  if(ps) gReve->AddGlobalRenderElement(ps);
  ps = project_node(nlt, "/ALIC_1/TRD TOF_holder_1/BFMO_1", useBuffPols);
  if(ps) gReve->AddGlobalRenderElement(ps);
  ps = project_node(nlt, "/ALIC_1/TRD TOF_holder_1/BBMO_1", useBuffPols);
  if(ps) gReve->AddGlobalRenderElement(ps);

  for(Int_t i = 1; i<6; i++) {
    ps = project_node(nlt, Form("/ALIC_1/PHOS_holder_1/PHOS_%d", i), useBuffPols);
    ps->SetFillColor((Color_t)(kOrange-4));
    if(ps) gReve->AddGlobalRenderElement(ps);
  }
  project_nodes( nlt, "/ALIC_1/FMD_holder_1", useBuffPols);  
  project_nodes( nlt, "/ALIC_1/HMPID_holder_1", useBuffPols);
}
