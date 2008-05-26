// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTRDModuleImp.h"
#include "AliEveTRDData.h"

#include "TGListTree.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

#include "TEveManager.h"
#include "TEveTrack.h"
#include "TEveGeoNode.h"
#include "TEveTrans.h"


#include "AliLog.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDmcmTracklet.h"

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

  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    (dynamic_cast<AliEveTRDModule*>(*iter))->Paint(option);
    iter++;
  }
}

//______________________________________________________________________________
void AliEveTRDNode::Reset()
{
  // Reset.

  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    (dynamic_cast<AliEveTRDModule*>(*iter))->Reset();
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
void AliEveTRDNode::EnableListElements()
{
  // Enable list elements.

  SetRnrSelf(kTRUE);
  AliEveTRDNode *node = 0x0;
  AliEveTRDChamber *chmb = 0x0;
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))){
      node->SetRnrSelf(kTRUE);
      node->EnableListElements();
    }
    if((chmb = dynamic_cast<AliEveTRDChamber*>(*iter))) chmb->SetRnrSelf(kTRUE);
    iter++;
  }
  gEve->Redraw3D();
}

//______________________________________________________________________________
void AliEveTRDNode::DisableListElements()
{
  // Disable list elements.

  SetRnrSelf(kFALSE);
  AliEveTRDNode *node = 0x0;
  AliEveTRDChamber *chmb = 0x0;
  List_i iter = fChildren.begin();
  while(iter != fChildren.end()){
    if((node = dynamic_cast<AliEveTRDNode*>(*iter))){
      node->SetRnrSelf(kFALSE);
      node->DisableListElements();
    }
    if((chmb = dynamic_cast<AliEveTRDChamber*>(*iter))) chmb->SetRnrSelf(kFALSE);
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
    AddElement(fRecPoints = new AliEveTRDClusters(this));
    fRecPoints->SetMarkerSize(.2);
    fRecPoints->SetMarkerStyle(24);
    fRecPoints->SetMarkerColor(kGreen);
    fRecPoints->SetOwnIds(kTRUE);
  }
  fRecPoints->Reset();

  Float_t q;
  Float_t g[3]; //global coordinates
  AliTRDcluster *c=0x0;
  for(int iclus=0; iclus<clusters->GetEntriesFast(); iclus++){
    c = (AliTRDcluster*)clusters->UncheckedAt(iclus);
    c->GetGlobalXYZ(g); 
    q = c->GetQ();
    Int_t id = fRecPoints->SetNextPoint(g[0], g[1], g[2]);    
    fRecPoints->SetPointId(id, new AliTRDcluster(*c));
  }
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
  // Info("LoadDigits()", Form("digits =0x%x", digits));

  if(!fDigits) AddElement(fDigits = new AliEveTRDDigits(this));

  fDigits->Reset();
  fDigits->SetData(digits);
  fLoadDigits = kTRUE;
}

//______________________________________________________________________________
void AliEveTRDChamber::LoadHits(TClonesArray *hits, Int_t &idx)
{
  //
  // Draw hits
  //
  //Info("AddHit()", Form("%s", GetName()));

  if(!fHits){
    AddElement(fHits = new AliEveTRDHits(this));
    fHits->SetMarkerSize(.1);
    fHits->SetMarkerColor(2);
    fHits->SetOwnIds(kTRUE);
  }
  fLoadHits = kTRUE;
  Int_t nhits = hits->GetEntriesFast();

  AliTRDhit *hit = 0x0;
  while(idx<nhits){
    hit = (AliTRDhit*)hits->UncheckedAt(idx);
    if(hit->GetDetector() != fDet) return;

    Int_t id = fHits->SetNextPoint(hit->X(), hit->Y(), hit->Z());
    fHits->SetPointId(id, new AliTRDhit(*hit));
    idx++;
  }
  return;
}

//______________________________________________________________________________
void AliEveTRDChamber::LoadTracklets(TObjArray *tracks)
{
  //
  // Draw tracks
  //
  if(!fGeo){
    Error("LoadTracklets()", Form("Geometry not set for chamber %d. Please call first AliEveTRDChamber::SetGeometry().", fDet));
    return;
  }
  // Info("LoadTracklets()", Form("tracks = 0x%x", tracks));

  if(!fTracklets){
    fTracklets = new std::vector<TEveTrack*>;
  } else fTracklets->clear();


  AliTRDmcmTracklet *trk = 0x0;
  Double_t cloc[3], cglo[3];
  for(int itrk=0; itrk<tracks->GetEntries();itrk++){
    trk = (AliTRDmcmTracklet*)tracks->At(itrk);
    trk->MakeTrackletGraph(fGeo,.5);
    fTracklets->push_back(new TEveTrack());
    fTracklets->back()->SetLineColor(4);

    cloc[0] = trk->GetTime0(); // x0
    cloc[1] = trk->GetOffset(); // y0
    cloc[2] = trk->GetRowz(); // z
    fGeo->RotateBack(fDet,cloc,cglo);
    fTracklets->back()->SetNextPoint(cglo[0], cglo[1], cglo[2]);

    cloc[0] += 3.7; // x1
    cloc[1] += TMath::Tan(trk->GetSlope()*TMath::Pi()/180.) * 3.7; // y1
    fGeo->RotateBack(fDet,cloc,cglo);
    fTracklets->back()->SetNextPoint(cglo[0], cglo[1], cglo[2]);
  }
  fLoadTracklets = kTRUE;
}

//____________________________________________________
void AliEveTRDChamber::Paint(Option_t* option)
{
  // Paint object.

  //AliInfo(GetName());
  if(!fRnrSelf) return;
  if(fDigits && fRnrDigits){
    if(fDigitsNeedRecompute){
      fDigits->ComputeRepresentation();
      fDigitsNeedRecompute = kFALSE;
    }
    fDigits->Paint(option);
  }

  if(fRecPoints && fRnrRecPoints) fRecPoints->GetObject()->Paint(option);

  if(fHits && fRnrHits) fHits->GetObject()->Paint(option);

  if(fTracklets && fRnrTracklets){
    for(std::vector<TEveTrack*>::iterator i=fTracklets->begin(); i != fTracklets->end(); ++i) (*i)->Paint(option);
  }
}

//______________________________________________________________________________
void AliEveTRDChamber::Reset()
{
  // Reset.

  if(fHits){
    fHits->Reset();
    fLoadHits = kFALSE;
  }
  if(fDigits){
    fDigits->Reset();
    fLoadDigits = kFALSE;
    fDigitsNeedRecompute = kTRUE;
  }
  if(fRecPoints){
    fRecPoints->Reset();
    fLoadRecPoints = kFALSE;
  }
  if(fTracklets){
    fTracklets->clear();
    fLoadTracklets = kFALSE;
  }
}

//______________________________________________________________________________
void AliEveTRDChamber::SetGeometry(AliTRDgeometry *geo)
{
  // Set geometry.

  fGeo = geo;

  Int_t  ism  = geo->GetSector(fDet);
  Int_t  istk = geo->GetChamber(fDet);
  Int_t  ipla = geo->GetPlane(fDet);
  Int_t icha  = istk*6+ipla;

  // define pad plane size in pads
  AliTRDpadPlane *pp = fGeo->GetPadPlane(ipla, istk);
  fNrows   = pp->GetNrows();
  fNcols   = pp->GetNcols();

// this version for setting the rendarable object is not working very nice
// Int_t shape_offset = TEveGeoShape::Class()->GetDataMemberOffset("fShape");
// TEveGeoShape* eg_shape = new TEveGeoShape("geometry");
// eg_shape->RefMainTrans().SetFrom(* gGeoManager->GetCurrentMatrix());
// * (TGeoShape**) (((char*)eg_shape) + shape_offset) = gGeoManager->GetCurrentVolume()->GetShape();
// 
// eg_shape->StampColorSelection();

  // define rendarable volumes
  gGeoManager->cd(Form("/B077_1/BSEGMO%d_1/BTRD%d_1/UTR1_1/UTS1_1/UTI1_1/UT%02d_1", ism, ism, icha));
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

