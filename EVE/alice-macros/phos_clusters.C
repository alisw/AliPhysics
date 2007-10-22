
Reve::PointSet* phos_clusters(Reve::RenderElement* cont=0)
{
  Alieve::Event::AssertGeometry();

  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadRecPoints("PHOS");

  TTree *cTree = rl->GetTreeR("PHOS", false);

  Reve::PointSet* clusters = new Reve::PointSet(10000);
  clusters->SetOwnIds(kTRUE);

  TObjArray *arr=NULL;
  TBranch *branch=cTree->GetBranch("PHOSEmcRP");
  branch->SetAddress(&arr);

  Int_t nentr=(Int_t)branch->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!branch->GetEvent(i)) continue;

    Int_t ncl=arr->GetEntriesFast();
    while (ncl--) {
      AliCluster *cl=(AliCluster*)arr->UncheckedAt(ncl);

      Float_t g[3]; //global coordinates
      cl->GetGlobalXYZ(g);

      AliCluster *atp = new AliCluster(*cl);
      clusters->SetNextPoint(g[0], g[1], g[2]);
      clusters->SetPointId(atp);
    }
  }

  Warning("phos_clusters"," %d",clusters->Size());

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
