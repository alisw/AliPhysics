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
#include <AliCluster3D.h>

#include <TClonesArray.h>

#endif

Reve::PointSet* hmpid_clusters(Reve::RenderElement* cont=0, Float_t maxR=1000)
{
  const Int_t nCh=7;
  TClonesArray *cl[nCh] = {0,0,0,0,0,0,0};
  Char_t *name[nCh]={
    "HMPID0",
    "HMPID1",
    "HMPID2",
    "HMPID3",
    "HMPID4",
    "HMPID5",
    "HMPID6"
  };


  Reve::PointSet* clusters = new Reve::PointSet(10000);
  clusters->SetOwnIds(kTRUE);

  Alieve::Event::AssertGeometry();
  
  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadRecPoints("HMPID");

  TTree *cTree = rl->GetTreeR("HMPID", false);
  if (!cTree) return 0;

  for (Int_t k=0; k<nCh; k++) {
     TBranch *br=cTree->GetBranch(name[k]);
     if (!br) return 0;
     br->SetAddress(&(cl[k]));
  } 

  if (!cTree->GetEvent(0)) return 0;


  for (Int_t i=0; i<nCh; i++) {
    TClonesArray *arr=cl[i];
    Int_t ncl=arr->GetEntriesFast();

    Float_t maxRsqr = maxR*maxR;
    while (ncl--) {
      AliCluster3D *c=(AliCluster3D*)arr->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g);
      if (g[0]*g[0]+g[1]*g[1] < maxRsqr)
      {
	clusters->SetNextPoint(g[0], g[1], g[2]);
	AliCluster3D *atp = new AliCluster3D(*c);
	clusters->SetPointId(atp);
      }
    }
  }

  if(clusters->Size() == 0 && gReve->GetKeepEmptyCont() == kFALSE) {
    Warning("hmpid_clusters", "No HMPID clusters");
    delete clusters;
    return 0;
  }

  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.2);
  clusters->SetMarkerColor(4);

  char form[1000];
  sprintf(form,"HMPID Clusters");
  clusters->SetName(form);

  char tip[1000];
  sprintf(tip,"N=%d", clusters->Size());
  clusters->SetTitle(tip);

  using namespace Reve;
  gReve->AddRenderElement(clusters, cont);
  gReve->Redraw3D();

  return clusters;
}
