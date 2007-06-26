
Reve::PointSet* phos_clusters(RenderElement* cont=0)
{
  Alieve::Event::AssertGeometry();

  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadRecPoints("PHOS");

  TTree *cTree = rl->GetTreeR("PHOS", false);

  Reve::PointSet* clusters = new Reve::PointSet(10000);
  clusters->SetOwnIds(kTRUE);

  AliPHOSEmcRecPoint *cl=NULL;
  TBranch *branch=cTree->GetBranch("PHOSEmcRP");
  branch->SetAddress(&cl);

  Int_t nentr=(Int_t)branch->GetEntries();
  Warning("phos_clusters",Form(" %d"),nentr);
  for (Int_t i=0; i<nentr; i++) {
    if (!branch->GetEvent(i)) continue;

    Float_t g[3]; //global coordinates
    cl->GetGlobalXYZ(g);

    clusters->SetNextPoint(g[0], g[1], g[2]);
    AliCluster *atp = new AliCluster(*cl);
    clusters->SetPointId(atp);
  }

  if(clusters->Size() == 0 && gReve->GetKeepEmptyCont() == kFALSE) {
    Warning("phos_clusters", "No PHOS clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.5);
  clusters->SetMarkerColor(4);

  char form[1000];
  sprintf(form,"PHOS Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  using namespace Reve;
  gReve->AddRenderElement(clusters);
  gReve->Redraw3D();

  return clusters;
}
