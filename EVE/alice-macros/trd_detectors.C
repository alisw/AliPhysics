#ifdef __CINT__
class TEvePointSet;
class TEveElement;
#else
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <EveBase/AliEveEventManager.h>

#include "AliRunLoader.h"
#include "AliCluster.h"
#include "AliTRDcluster.h"
#endif

TEveElementList* trd_detectors(TEveElement *cont = 0)
{
  // Link data containers
	AliEveEventManager::AssertGeometry();
	AliRunLoader *rl = AliEveEventManager::AssertRunLoader();

	// define EVE containers
	Int_t nclusters = 0;
  TEveElementList *clusters = new TEveElementList("TRD clusters"); 	
  TEvePointSet *clustersDet = 0x0;


//   AliTRDgeometry geo;
//   AliEveTRDChamber *chm = 0x0;

  // Fill EVE containers
	TObjArray *TRDcluster = 0x0;
	rl->LoadRecPoints("TRD");
	TTree *recPoints = rl->GetTreeR("TRD", kFALSE);
	recPoints->SetBranchAddress("TRDcluster", &TRDcluster);
	
  Int_t nentr=(Int_t)recPoints->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    if (!recPoints->GetEvent(i)) continue;


    Int_t det = -1;
    Int_t ncl=TRDcluster->GetEntriesFast();
    nclusters+=ncl;
    while (ncl--) {
      AliTRDcluster *c = (AliTRDcluster*)TRDcluster->UncheckedAt(ncl);
      Float_t g[3]; //global coordinates
      c->GetGlobalXYZ(g); 
      if(det<0){
        det = c->GetDetector();
        clustersDet= new TEvePointSet(Form("detector %d", det));
        clustersDet->SetOwnIds(kTRUE);
        clustersDet->SetMarkerStyle(2);
        clustersDet->SetMarkerSize(0.2);
        clustersDet->SetMarkerColor(kBlue);

/*        chm = new AliEveTRDChamber(det);
        chm->SetGeometry(&geo);
        chm->LoadClusters(TRDcluster);
        break;*/
      }
			Int_t id = clustersDet->SetNextPoint(g[0], g[1], g[2]);
      clustersDet->SetPointId(id, new AliTRDcluster(*c));
    }
    clustersDet->SetTitle(Form("Clusters %d", clustersDet->Size()));
    clusters->AddElement(clustersDet);

    //clusters->AddElement(chm);

    TRDcluster->Clear();
  }
  clusters->SetTitle(Form("Clusters %d", nclusters));
  
  gEve->AddElement(clusters, cont);
  gEve->Redraw3D();

  return clusters;
}
