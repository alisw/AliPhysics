// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTRDModuleImp.h"
#include "AliEveTRDData.h"

#include <TGListTree.h>

#include "TEveManager.h"
#include "TEveTrack.h"

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
  TEveElement(), AliEveTRDModule("Chmb", det),
  fDigits(0), fHits(0), fRecPoints(0), fTracklets(0),
  fRowMax(-1), fColMax(-1), fTimeMax(22), fSamplingFrequency(0),
  fX0(0.), fPla(-1),
  fPadPlane (0),
  fGeo      (0)
{
  //
  // Constructor
  //

  AliTRDCommonParam* parCom = AliTRDCommonParam::Instance();
  fSamplingFrequency = parCom->GetSamplingFrequency();
}

//______________________________________________________________________________
Int_t AliEveTRDChamber::GetSM() const
{
  // Get supermodule.

  if(!fGeo){
    AliWarning("Fail. No TRD geometry defined.");
    return -1;
  }
  return fGeo->GetSector(fDet);
}

//______________________________________________________________________________
Int_t AliEveTRDChamber::GetSTK() const
{
  // Get stack.

  if(!fGeo){
    AliWarning("Fail. No TRD geometry defined.");
    return -1;
  }
  return fGeo->GetChamber(fDet);
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
    fRecPoints = new AliEveTRDClusters(this);
    fRecPoints->SetMarkerSize(1.);
    fRecPoints->SetMarkerStyle(24);
    fRecPoints->SetMarkerColor(6);
    fRecPoints->SetOwnIds(kTRUE);
  } else fRecPoints->Reset();

  Float_t q;
  Double_t cloc[3], cglo[3];

  AliTRDcluster *c=0x0;
  for(int iclus=0; iclus<clusters->GetEntriesFast(); iclus++){
    c = (AliTRDcluster*)clusters->UncheckedAt(iclus);
    cloc[0] = c->GetX();
    cloc[1] = c->GetY();
    cloc[2] = c->GetZ();
    q = c->GetQ();
    fGeo->RotateBack(fDet,cloc,cglo);
    fRecPoints->SetNextPoint(cglo[0], cglo[1], cglo[2]);
    fRecPoints->SetPointId(c);
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

  if(!fDigits) fDigits = new AliEveTRDDigits(this);
  else fDigits->Reset();

  fDigits->SetData(digits);
  fLoadDigits = kTRUE;
}

//______________________________________________________________________________
void AliEveTRDChamber::AddHit(AliTRDhit *hit)
{
  //
  // Draw hits
  //
  // Info("AddHit()", Form("%s", GetName()));

  if(!fHits){
    fHits = new AliEveTRDHits(this);
    fHits->SetMarkerSize(.1);
    fHits->SetMarkerColor(2);
    fHits->SetOwnIds(kTRUE);
  }

  fHits->SetNextPoint(hit->X(), hit->Y(), hit->Z());
  fHits->SetPointId(hit);
  fLoadHits = kTRUE;
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

  fPla = fGeo->GetPlane(fDet);
  fX0 = fGeo->GetTime0(fPla);

  fPadPlane = fGeo->GetPadPlane(fPla,fGeo->GetChamber(fDet));
  fRowMax   = fPadPlane->GetNrows();
  fColMax   = fPadPlane->GetNcols();
}

