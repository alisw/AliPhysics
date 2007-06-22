
Reve::PointSet* phos_clusters(RenderElement* cont=0)
{
  Alieve::Event::AssertGeometry();

  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadRecPoints("PHOS");

  TTree *cTree = rl->GetTreeR("PHOS", false);

  Reve::PointSet* clusters = new Reve::PointSet(10000);
  clusters->SetOwnIds(kTRUE);

  TClonesArray *cl=NULL;
  TBranch *branch=cTree->GetBranch("PHOSEmcRP");
  branch->SetAddress(&cl);

  Int_t nentr=(Int_t)cTree->GetEntries();
  Warning("phos_clusters",Form(" %d"),nentr);
  for (Int_t i=0; i<nentr; i++) {
    if (!cTree->GetEvent(i)) continue;

    Int_t ncl=cl->GetEntriesFast();

    while (ncl--) {
      AliCluster *c=(AliCluster*)cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);

      clusters->SetNextPoint(g[0], g[1], g[2]);
      AliCluster *atp = new AliCluster(*c);
      clusters->SetPointId(atp);
    }
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
