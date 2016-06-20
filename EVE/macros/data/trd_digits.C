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

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliCluster.h>
#include <AliRunLoader.h>
#include <AliTRDarrayADC.h>
#include <AliTRDcluster.h>
#include <AliTRDgeometry.h>
#include <AliTRDdigitsManager.h>
#include <AliEveEventManager.h>
#include <AliEveTRDModuleImp.h>
#else
class TEvePointSet;
class TEveElement;
#endif

TEveElementList* trd_digits(Int_t sector = -1, TEveElement *cont = 0)
{
  // Link data containers
  if(!gGeoManager) AliEveEventManager::Instance()->AssertGeometry();

  AliRunLoader *rl = AliEveEventManager::AssertRunLoader();
  
  // define EVE containers
  TEveElementList *list = new TEveElementList("TRD Digits");
  
  AliEveTRDNode *sm(NULL), *stk(NULL); 
  AliEveTRDChamber *chm(NULL);
  // Link TRD digits
  rl->LoadDigits("TRD");
  TTree *tD = rl->GetTreeD("TRD", kFALSE);
  if(!tD){ 
    Error("trd_digits", "Missing digits tree");
    return NULL;
  }
  AliTRDdigitsManager dm;
  dm.ReadDigits(tD);

  AliTRDgeometry *geo = new AliTRDgeometry();
  Int_t sBegin=sector<0?0:sector,
        sEnd  =sector<0?(AliTRDgeometry::kNsector):sector+1;
  Int_t jdet(0);
  for(Int_t isec=sBegin; isec<sEnd; isec++) {
    sm = new AliEveTRDNode("Sector", isec);
    for(Int_t istk(0); istk<AliTRDgeometry::kNstack; istk++) {
      stk = new AliEveTRDNode("Stack", istk);
      stk->SetTitle(Form("Index %d", isec*AliTRDgeometry::kNstack+istk));
      for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++) {
        Int_t det=AliTRDgeometry::GetDetector(ily, istk, isec);
        if(!(dm.GetDigits(det)->GetDim())) continue;
        chm=new AliEveTRDChamber(det);
        chm->SetGeometry(geo);
        chm->LoadDigits(&dm);
        stk->AddElement(chm);
        jdet++;
      }
      if(!stk->HasChildren()){
        delete stk;
        continue;
      }
      sm->AddElement(stk);
    }
    if(!sm->HasChildren()){
      delete sm;
      continue;
    }
    list->AddElement(sm);
  }
  rl->UnloadDigits("TRD");
  gEve->AddElement(list, cont);
  gEve->Redraw3D();

  Info("trd_digits", "TRD chambers with data for current selection %d.", jdet);
  return list;
}
