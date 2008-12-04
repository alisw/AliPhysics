//
// How to steer the basic TRD data containers from a macro.
// 
// The loading and looping over events is performed by the
// AliEve mechanism. For a stand alone TRD loop control check 
// the macro "trd_loader.C"
// 
// Usage:
// .x trd_detectors.C(sector)
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
// 
#ifdef __CINT__
class TEvePointSet;
class TEveElement;
#else
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <EveBase/AliEveEventManager.h>
#include <EveDet/AliEveTRDModuleImp.h>

#include "AliRunLoader.h"
#include "AliCluster.h"
#include "TRD/AliTRDcluster.h"
#include "TRD/AliTRDgeometry.h"
#include "TRD/AliTRDdigitsManager.h"
#endif

TEveElementList* trd_detectors(Int_t sector = -1, TEveElement *cont = 0)
{
  // Link data containers
  AliEveEventManager::AssertGeometry();

  AliRunLoader *rl = AliEveEventManager::AssertRunLoader();
  
  // define EVE containers
  TEveElementList *list = new TEveElementList("TRD Detectors");
  
  AliTRDgeometry *geo = new AliTRDgeometry();
  //geo->CreateClusterMatrixArray();
  
  AliEveTRDNode *sm = 0x0, *stk = 0x0; 
  AliEveTRDChamber *chm=0x0;

  // Link TRD containers
  TObjArray *clusters = 0x0;
  rl->LoadRecPoints("TRD");
  TTree *tR = rl->GetTreeR("TRD", kFALSE);
  tR->SetBranchAddress("TRDcluster", &clusters);

  rl->LoadDigits("TRD");
  TTree *tD = rl->GetTreeD("TRD", kFALSE);
  AliTRDdigitsManager dm; dm.ReadDigits(tD);

  for(Int_t i=0; i<tR->GetEntries(); i++) {
    if (!tR->GetEvent(i)) continue;
    if(!clusters->GetEntries()) continue;
    Int_t icl=0; AliTRDcluster *c = 0x0;
    while(!(c = (AliTRDcluster*)clusters->UncheckedAt(icl++))) {;}
    if(!c) continue;

    Int_t idet, ism, istk, ipla; 
    idet = c->GetDetector();
    ism  = geo->GetSector(idet);
    istk = geo->GetStack(idet);
    ipla = geo->GetLayer(idet);
    if(sector>=0 && ism != sector) continue;
    if(!(sm = (AliEveTRDNode*)list->FindChild(Form("SM%03d", ism)))){ 
      list->AddElement(sm = new AliEveTRDNode("SM", ism));
      sm->SetElementTitle(Form("Supermodule %2d", ism));
    }
    if(!(stk=(AliEveTRDNode*)sm->FindChild(Form("Stack%03d", istk)))){
      sm->AddElement(stk = new AliEveTRDNode("Stack", istk));
      stk->SetElementTitle(Form("SM %2d Stack %1d", ism, istk));
    }
    stk->AddElement(chm = new AliEveTRDChamber(idet));
    chm->SetGeometry(geo);
    chm->LoadClusters(clusters);
    chm->LoadDigits(&dm);

    //clusters->Clear();
  }

  rl->UnloadDigits("TRD");
  rl->UnloadRecPoints("TRD");

  gEve->AddElement(list, cont);
  gEve->Redraw3D();

  return list;
}
