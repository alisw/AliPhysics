// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTRDModuleImp.h"
#include "AliEveTRDData.h"
#include "EveBase/AliEveEventManager.h"

#include "TTree.h"
#include "TGListTree.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

#include "TEveManager.h"
#include "TEveTrack.h"
#include "TEveGeoNode.h"
#include "TEveTrans.h"


#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDtrackingChamber.h"
#include "AliTRDtrackletMCM.h"

ClassImp(AliEveTRDChamber)
ClassImp(AliEveTRDNode)

//______________________________________________________________________________
AliEveTRDNode::AliEveTRDNode(const char *typ, Int_t det) :
  TEveElement(), AliEveTRDModule(typ, det)
{
  // Xonstructor.
}

//______________________________________________________________________________
void AliEveTRDNode::Paint(Option_t* option)
{
  // Paint object.

  AliEveTRDModule *module(NULL);
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((module = dynamic_cast<AliEveTRDModule*>(*iter))) module->Paint(option);
    iter++;
  }
}

//______________________________________________________________________________
void AliEveTRDNode::Reset()
{
  // Reset.

  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    //(dynamic_cast<AliEveTRDModule*>(*iter))->Reset();
    iter++;
  }
}

//______________________________________________________________________________
void AliEveTRDNode::Collapse()
{
  // Collapse.

  TGListTree *list = gEve->GetListTree();
  AliEveTRDNode *node = 0x0;
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))) node->Collapse();
    list->CloseItem(FindListTreeItem(list));
    iter++;
  }
}

//______________________________________________________________________________
void AliEveTRDNode::Expand()
{
  // Expand.

  TGListTree *list = gEve->GetListTree();
  AliEveTRDNode *node = 0x0;
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))) node->Expand();
    list->OpenItem(FindListTreeItem(list));
    iter++;
  }
}

//______________________________________________________________________________
void AliEveTRDNode::EnableListElements(Bool_t rnr_self, Bool_t rnr_children)
{
  // Enable list elements.

  SetRnrSelf(rnr_self);
  AliEveTRDNode *node(NULL);
  AliEveTRDChamber *chmb(NULL);
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))){
      node->SetRnrSelf(rnr_children);
      node->EnableListElements(rnr_children, rnr_children);
    }
    if((chmb = dynamic_cast<AliEveTRDChamber*>(*iter))) chmb->EnableListElements(rnr_children, rnr_children);
    iter++;
  }
  gEve->Redraw3D();
}

//______________________________________________________________________________
void AliEveTRDNode::DisableListElements(Bool_t rnr_self, Bool_t rnr_children)
{
  // Disable list elements.

  SetRnrSelf(rnr_self);
  AliEveTRDNode *node(NULL);
  AliEveTRDChamber *chmb(NULL);
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))){
      node->SetRnrSelf(rnr_children);
      node->DisableListElements(rnr_children, rnr_children);
    }
    if((chmb = dynamic_cast<AliEveTRDChamber*>(*iter))) chmb->DisableListElements(rnr_children, rnr_children);
    iter++;
  }
  gEve->Redraw3D();
}

//______________________________________________________________________________
void AliEveTRDNode::UpdateLeaves()
{
  // Update leaves.

  AliEveTRDModule *module;
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    module = dynamic_cast<AliEveTRDModule*>(*iter);
    if(!module) continue;

    module->fRnrHits = fRnrHits;
    module->fRnrDigits = fRnrDigits;
    module->fDigitsLog = fDigitsLog;
    module->fDigitsBox = fDigitsBox;
    module->fDigitsThreshold = fDigitsThreshold;
    module->fDigitsNeedRecompute = fDigitsNeedRecompute;
    module->fRnrRecPoints = fRnrRecPoints;
    module->fRnrTracklets = fRnrTracklets;
    iter++;
  }

  AliEveTRDNode *node = 0x0;
  iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))) node->UpdateLeaves();
    iter++;
  }
}


//______________________________________________________________________________
void AliEveTRDNode::UpdateNode()
{
  // Update node.

  // Info("UpdateNode()", Form("%s", GetName()));
  AliEveTRDNode *node = 0x0;
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))) node->UpdateNode();
    iter++;
  }

  Int_t score[11];
  for(int i=0; i<11; i++) score[i] = 0;
  AliEveTRDModule *module;
  iter = fChildren.begin();
  while(iter != fChildren.end()){
    module = dynamic_cast<AliEveTRDModule*>(*iter);
    if(!module) continue;
    score[0] += (module->fLoadHits) ? 1 : 0;
    score[1] += (module->fRnrHits) ? 1 : 0;

    score[2] += (module->fLoadDigits) ? 1 : 0;
    score[3] += (module->fRnrDigits) ? 1 : 0;
    score[4] += (module->fDigitsLog) ? 1 : 0;
    score[5] += (module->fDigitsBox) ? 1 : 0;
    score[6] += (module->fDigitsNeedRecompute) ? 1 : 0;

    score[7] += (module->fLoadRecPoints) ? 1 : 0;
    score[8] += (module->fRnrRecPoints) ? 1 : 0;

    score[9] += (module->fLoadTracklets) ? 1 : 0;
    score[10] += (module->fRnrTracklets) ? 1 : 0;
    iter++;
  }

  Int_t size = fChildren.size();
  fLoadHits      = (score[0] > 0) ? kTRUE : kFALSE;
  fRnrHits       = (score[1] == size) ? kTRUE : kFALSE;

  fLoadDigits    = (score[2] > 0) ? kTRUE : kFALSE;
  fRnrDigits     = (score[3] == size) ? kTRUE : kFALSE;
  fDigitsLog     = (score[4] == size) ? kTRUE : kFALSE;
  fDigitsBox     = (score[5] == size) ? kTRUE : kFALSE;
  fDigitsNeedRecompute = (score[6] == size) ? kTRUE : kFALSE;

  fLoadRecPoints = (score[7] > 0) ? kTRUE : kFALSE;
  fRnrRecPoints  = (score[8] == size) ? kTRUE : kFALSE;

  fLoadTracklets = (score[9] > 0) ? kTRUE : kFALSE;
  fRnrTracklets  = (score[10] == size) ? kTRUE : kFALSE;
}


///////////////////////////////////////////////////////////
////////////      AliEveTRDChamber     ////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDChamber::AliEveTRDChamber(Int_t det) :
  TEveElement()
  ,AliEveTRDModule("Chmb", det)
  ,fDigits(0x0)
  ,fHits(0x0)
  ,fRecPoints(0x0)
  ,fTracklets(0x0)
  ,fGeo(0x0)
  ,fShape(0x0)
  ,fNrows(-1)
  ,fNcols(-1)
  ,fNtime(22)
{
  //
  // Constructor
  //
}


//______________________________________________________________________________
void AliEveTRDChamber::LoadClusters(TObjArray *clusters)
{
  //
  // Draw clusters
  //

  if(!fGeo){
    AliError(Form("Geometry not set for chamber %d. Please call first AliEveTRDChamber::SetGeometry().", fDet));
    return;
  }

  if(!fRecPoints){ 
    AddElement(fRecPoints = new AliEveTRDClusters());
    fRecPoints->SetTitle(Form("Clusters for Det %d", GetID()));
  }
  fRecPoints->Reset();

  Float_t q;
  Float_t g[3]; //global coordinates
  AliTRDcluster *c=0x0;
  Int_t nc = clusters->GetEntriesFast();
  for(int iclus=0; iclus<nc; iclus++){
    c = (AliTRDcluster*)clusters->UncheckedAt(iclus);
    c->GetGlobalXYZ(g); 
    q = c->GetQ();
    Int_t id = fRecPoints->SetNextPoint(g[0], g[1], g[2]);    
    fRecPoints->SetPointId(id, new AliTRDcluster(*c));
  }
  fRecPoints->StampObjProps();
  fLoadRecPoints = kTRUE;
}


//______________________________________________________________________________
void AliEveTRDChamber::LoadClusters(AliTRDtrackingChamber *tc)
{
  if(!fGeo){
    AliError(Form("Geometry not set for chamber %d. Please call first AliEveTRDChamber::SetGeometry().", fDet));
    return;
  }

  if(!fRecPoints){ 
    AddElement(fRecPoints = new AliEveTRDClusters());
    fRecPoints->SetTitle(Form("Clusters for Det %d", GetID()));
  }
  fRecPoints->Reset();

  Float_t g[3]; //global coordinates
  const AliTRDchamberTimeBin *tb = 0x0;
  for(int itb=0; itb<AliTRDseedV1::kNtb; itb++){
    tb = tc->GetTB(itb);
    if(!(Int_t(*tb))) continue;
    const AliTRDcluster *c= 0x0; Int_t ic = 0;
    while((c=tb->GetCluster(ic))){
      c->GetGlobalXYZ(g); 
      Int_t id = fRecPoints->SetNextPoint(g[0], g[1], g[2]);    
      fRecPoints->SetPointId(id, new AliTRDcluster(*c));
      ic++;
    }
  }
  fRecPoints->StampObjProps();
  fLoadRecPoints = kTRUE;
}


//______________________________________________________________________________
void AliEveTRDChamber::LoadDigits(AliTRDdigitsManager *digits)
{
  //
  // Draw digits
  //
  if(!fGeo){
    AliError(Form("Geometry not set for chamber %d. Please call first AliEveTRDChamber::SetGeometry().", fDet));
    return;
  }

  if(!fDigits) AddElement(fDigits = new AliEveTRDDigits(this));

  //fDigits->Reset();
  fDigits->SetData(digits);
  fDigits->StampObjProps();
  fDigitsNeedRecompute = kTRUE;  
  fLoadDigits = kTRUE;
}

//______________________________________________________________________________
void AliEveTRDChamber::LoadHits(TClonesArray *hits, Int_t &idx)
{
  //
  // Draw hits
  //

  if(!fHits){ 
    AddElement(fHits = new AliEveTRDHits());
    fHits->SetTitle(Form("Hits for Det %d", GetID()));
  }
  fLoadHits = kTRUE;
  Int_t nhits = hits->GetEntriesFast();

  AliTRDhit *hit = 0x0;
  while(idx<nhits){
    hit = (AliTRDhit*)hits->UncheckedAt(idx);
    if(hit->GetDetector() != fDet) return;

    Int_t id = fHits->SetNextPoint(hit->X(), hit->Y(), hit->Z());
    fHits->SetPointId(id, new AliTRDhit(*hit));
    fHits->StampObjProps();
    idx++;
  }
  return;
}

//______________________________________________________________________________
void AliEveTRDChamber::LoadTracklets(TTree *trklTree)
{
  //
  // Draw tracks
  //

  if(!fGeo){
    Error("LoadTracklets()", "Geometry not set for chamber %d. Please call first AliEveTRDChamber::SetGeometry().", fDet);
    return;
  }

  if(!fTracklets){
    fTracklets = new TClonesArray("AliEveTRDTrackletOnline",100);
  } else {
    fTracklets->Delete();
    TEveElementList *trklChild = (TEveElementList*) FindChild("Tracklets");
    if (trklChild)
      trklChild->Destroy();
  }


  TBranch *mcmBranch = trklTree->GetBranch("mcmtrklbranch");
  if (!mcmBranch)
    return;

  AliTRDtrackletMCM *trkl = 0x0;
  mcmBranch->SetAddress(&trkl);

  TEveElementList* listOfTracklets = new TEveElementList("Tracklets");
  gEve->AddElement(listOfTracklets, this);

  for(Int_t iTrkl = 0; iTrkl < mcmBranch->GetEntries(); iTrkl++){
    mcmBranch->GetEntry(iTrkl);
    if (trkl->GetDetector() == GetID()) {
      new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliEveTRDTrackletOnline(trkl);
      gEve->AddElement(new AliEveTRDTrackletOnline(trkl), listOfTracklets);
    }
  }
  fLoadTracklets = kTRUE;
}


//______________________________________________________________________________
void AliEveTRDChamber::SetGeometry(AliTRDgeometry *geo)
{
  // Set geometry.

  fGeo = geo;

  Int_t  ism  = geo->GetSector(fDet);
  Int_t  istk = geo->GetStack(fDet);
  Int_t  ilyr = geo->GetLayer(fDet);
  Int_t  icha = istk*6+ilyr;
  Int_t idx(1); 
  if(ism>12&&ism<16) idx=3;
  else if(ism==11||ism==12) idx=2;

  // define pad plane size in pads
  AliTRDpadPlane *pp = fGeo->GetPadPlane(ilyr, istk);
  fNrows   = pp->GetNrows();
  fNcols   = pp->GetNcols();

// this version for setting the rendarable object is not working very nice
// Int_t shape_offset = TEveGeoShape::Class()->GetDataMemberOffset("fShape");
// TEveGeoShape* eg_shape = new TEveGeoShape("geometry");
// eg_shape->RefMainTrans().SetFrom(* gGeoManager->GetCurrentMatrix());
// * (TGeoShape**) (((char*)eg_shape) + shape_offset) = gGeoManager->GetCurrentVolume()->GetShape();
// 
// eg_shape->StampColorSelection();
  if(!(gGeoManager)){
    AliEveEventManager::AssertGeometry();
    if(!(gGeoManager)){
      AliError("Geo manager not available.");
      return;
    }
  }

  // define rendarable volumes
  if(!gGeoManager->cd(Form("/B077_1/BSEGMO%d_1/BTRD%d_1/UTR%d_1/UTS%d_1/UTI%d_1/UT%02d_1", ism, ism, idx, idx, idx, icha))) return;

  fShape = new TEveGeoTopNode(gGeoManager, gGeoManager->GetCurrentNode());
  fShape->RefMainTrans().SetFrom(*gGeoManager->GetCurrentMatrix());
  fShape->DisableListElements();
  fShape->SetRnrSelf(kFALSE);
// try to set the properties but it is crashing !!
//   TEveGeoNode *node = 0x0; 
//   if((node = (TEveGeoNode*)fShape->FindChild(Form("UA%02d_1", icha)))) node->SetRnrState(kTRUE);
//   else AliWarning(Form("Can not retrieve geo node UA%02d_1", icha));

  AddElement(fShape);
}

