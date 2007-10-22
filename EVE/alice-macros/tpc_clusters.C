#ifdef __CINT__

namespace Reve
{
class RenderElement;
class PointSet;
}

#else

#include <Reve/Reve.h>
#include <Reve/ReveManager.h>
#include <Reve/PointSet.h>
#include <Alieve/EventAlieve.h>

#include <AliRunLoader.h>
#include <AliCluster.h>
#include <TPC/AliTPCClustersRow.h>

#endif

Reve::PointSet* tpc_clusters(Reve::RenderElement* cont=0, Float_t maxR=270)
{
  const Int_t kMaxCl=100*160;

  Alieve::Event::AssertGeometry();

  Reve::PointSet* clusters = new Reve::PointSet(kMaxCl);
  clusters->SetOwnIds(kTRUE);

  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadRecPoints("TPC");

  AliTPCClustersRow *clrow=new AliTPCClustersRow();
  clrow->SetClass("AliTPCclusterMI");
  clrow->SetArray(kMaxCl);

  TTree *cTree = rl->GetTreeR("TPC", false);
  cTree->SetBranchAddress("Segment", &clrow);

  Float_t maxRsqr = maxR*maxR;
  TClonesArray *cl=clrow->GetArray(); 
  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!cTree->GetEvent(i)) continue;

    Int_t ncl=cl->GetEntriesFast();

    while (ncl--) {
      AliCluster *c=(AliCluster*)cl->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      if (g[0]*g[0]+g[1]*g[1] < maxRsqr)
      {
	clusters->SetNextPoint(g[0], g[1], g[2]);
	AliCluster *atp = new AliCluster(*c);
	clusters->SetPointId(atp);
      }
    }
    cl->Clear();
  }

  delete clrow;

  if(clusters->Size() == 0 && gReve->GetKeepEmptyCont() == kFALSE) {
    Warning("tpc_clusters", "No TPC clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.2);
  clusters->SetMarkerColor(4);

  char form[1000];
  sprintf(form,"TPC Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  using namespace Reve;
  gReve->AddRenderElement(clusters, cont);
  gReve->Redraw3D();

  return clusters;
}
