// $Id$

Reve::PointSet* clusters_from_index(Int_t index=0, Reve::RenderElement* cont=0)
{
  AliESDEvent* esd = Alieve::Event::AssertESD();

  if (index < 0) {
    Warning("clusters_from_index", "index not set.");
    return 0;
  }

  if (index >= esd->GetNumberOfTracks()) {
    Warning("clusters_from_index", "index out of range");
    return 0;
  }

  Reve::PointSet* clusters = new Reve::PointSet(64);
  clusters->SetOwnIds(kTRUE);

  AliESDtrack* at = esd->GetTrack(index);
  const AliTrackPointArray* pArr = at->GetTrackPointArray();
  if (pArr == 0) {
    Warning("clusters_from_index", "TrackPointArray not stored with ESD track.");
    continue;
  }
  Int_t np =  pArr->GetNPoints();
  const Float_t* x = pArr->GetX();
  const Float_t* y = pArr->GetY();
  const Float_t* z = pArr->GetZ();
  for (Int_t i=0; i<np; ++i) {
    clusters->SetNextPoint(x[i], y[i], z[i]);
    AliTrackPoint *atp = new AliTrackPoint;
    pArr->GetPoint(*atp, i);
    clusters->SetPointId(atp);    }

  
  if(clusters->Size() == 0 && gReve->GetKeepEmptyCont() == kFALSE) {
    Warning("clusters_from_index", Form("No clusters for index '%d'", index));
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.5);
  clusters->SetMarkerColor(4);

  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  clusters->SetName(Form("Clusters idx=%d", index));
  char form[1000];
  sprintf(form,"Clusters idx=%d", index);
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  using namespace Reve;
  gReve->AddRenderElement(clusters);
  gReve->Redraw3D();

  return clusters;
}
