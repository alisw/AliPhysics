#include "TRDModuleImp.h"
#include "TRDData.h"

#include "TMath.h"
#include <TGListTree.h>

#include "Reve/RGTopFrame.h"
#include "Reve/Track.h"

#include "AliRun.h"
#include "AliTRDv1.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"
#include "AliTRDdigit.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDmcmTracklet.h"



using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(TRDChamber)
ClassImp(TRDNode)

//________________________________________________________
TRDNode::TRDNode(const char *typ, const Int_t det) :
  Reve::RenderElementListBase(), TRDModule(typ, det)
{
}

//________________________________________________________
void	TRDNode::Paint(Option_t* option)
{
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		(dynamic_cast<TRDModule*>(*iter))->Paint(option);
		iter++;
	}
}

//________________________________________________________
void	TRDNode::Reset()
{
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		(dynamic_cast<TRDModule*>(*iter))->Reset();
		iter++;
	}
}

//________________________________________________________
void TRDNode::Collapse()
{
	TGListTree *list = gReve->GetListTree();
	TRDNode *node = 0x0;
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		if((node = dynamic_cast<TRDNode*>(*iter))) node->Collapse();
		list->CloseItem(FindListTreeItem(list));
		iter++;
	}
}

//________________________________________________________
void TRDNode::Expand()
{
	TGListTree *list = gReve->GetListTree();
	TRDNode *node = 0x0;
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		if((node = dynamic_cast<TRDNode*>(*iter))) node->Expand();
		list->OpenItem(FindListTreeItem(list));
		iter++;
	}
}

//________________________________________________________
void TRDNode::EnableListElements()
{
	SetRnrElement(kTRUE);
	TRDNode *node = 0x0;
	TRDChamber *chmb = 0x0;	
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		if((node = dynamic_cast<TRDNode*>(*iter))){
			node->SetRnrElement(kTRUE);
			node->EnableListElements();
		}
		if((chmb = dynamic_cast<TRDChamber*>(*iter))) chmb->SetRnrElement(kTRUE);
		iter++;
	}
	gReve->Redraw3D();
}

//________________________________________________________
void TRDNode::DisableListElements()
{
	SetRnrElement(kFALSE);
	TRDNode *node = 0x0;
	TRDChamber *chmb = 0x0;	
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		if((node = dynamic_cast<TRDNode*>(*iter))){
			node->SetRnrElement(kFALSE);
			node->DisableListElements();
		}
		if((chmb = dynamic_cast<TRDChamber*>(*iter))) chmb->SetRnrElement(kFALSE);
		iter++;
	}
	gReve->Redraw3D();
}

//________________________________________________________
void TRDNode::UpdateLeaves()
{
	TRDModule *module;
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		module = dynamic_cast<TRDModule*>(*iter);
		if(!module) continue;
		
		module->fRnrHits = fRnrHits;
		module->fRnrDigits = fRnrDigits;
		module->fDigitsLog = fDigitsLog;
		module->fDigitsBox = fDigitsBox;
		module->fDigitsThreshold = fDigitsThreshold;
		module->kDigitsNeedRecompute = kDigitsNeedRecompute;
		module->fRnrRecPoints = fRnrRecPoints;
		module->fRnrTracklets = fRnrTracklets;
		iter++;
	}

	TRDNode *node = 0x0;
	iter = fChildren.begin();
	while(iter != fChildren.end()){
		if((node = dynamic_cast<TRDNode*>(*iter))) node->UpdateLeaves();
		iter++;
	}
}


//________________________________________________________
void TRDNode::UpdateNode()
{
//	Info("UpdateNode()", Form("%s", GetName()));
	TRDNode *node = 0x0;
	List_i iter = fChildren.begin();
	while(iter != fChildren.end()){
		if((node = dynamic_cast<TRDNode*>(*iter))) node->UpdateNode();
		iter++;
	}

	Int_t score[11];
	for(int i=0; i<11; i++) score[i] = 0;
	TRDModule *module;
	iter = fChildren.begin();
	while(iter != fChildren.end()){
		module = dynamic_cast<TRDModule*>(*iter);
		if(!module) continue;
		score[0] += (module->fLoadHits) ? 1 : 0;
		score[1] += (module->fRnrHits) ? 1 : 0;
		
		score[2] += (module->fLoadDigits) ? 1 : 0;
		score[3] += (module->fRnrDigits) ? 1 : 0;
		score[4] += (module->fDigitsLog) ? 1 : 0;
		score[5] += (module->fDigitsBox) ? 1 : 0;
		score[6] += (module->kDigitsNeedRecompute) ? 1 : 0;

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
	kDigitsNeedRecompute = (score[6] == size) ? kTRUE : kFALSE;

	fLoadRecPoints = (score[7] > 0) ? kTRUE : kFALSE;
	fRnrRecPoints  = (score[8] == size) ? kTRUE : kFALSE;

	fLoadTracklets = (score[9] > 0) ? kTRUE : kFALSE;
	fRnrTracklets  = (score[10] == size) ? kTRUE : kFALSE;
}


///////////////////////////////////////////////////////////
/////////////        TRDChamber       /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDChamber::TRDChamber(const Int_t det) :
  Reve::RenderElement(), TRDModule("Chmb", det), rowMax(-1), colMax(-1), timeMax(22), fX0(0.), fPla(-1)
{
  //
  // Constructor
  //
	
	fDigits    = 0x0;
	fHits      = 0x0;
	fRecPoints = 0x0;
	fTracklets = 0x0;
	
	AliTRDCommonParam* parCom = AliTRDCommonParam::Instance();
	samplingFrequency = parCom->GetSamplingFrequency();
	
	fGeo      = 0x0;
	fPadPlane = 0x0;
}

//________________________________________________________
TRDChamber::TRDChamber(const TRDChamber &mod):
  Reve::RenderElement(), TRDModule("Chmb", mod.fDet)
{
  //
  // Copy constructor
  //
	SetMainColor(mod.GetMainColor());

	if(mod.fDigits) {}
	if(mod.fHits) {}
	if(mod.fRecPoints){} 
}

//________________________________________________________
TRDChamber& TRDChamber::operator=(const TRDChamber &mod)
{
  //
  // Assignment operator
  //

  if (this != &mod) {
    fDet    = mod.fDet;
		if(mod.fDigits) {}
		if(mod.fHits) {}
		if(mod.fRecPoints){} 
  }
  return *this;
}

//________________________________________________________
void TRDChamber::LoadClusters(TObjArray *clusters)
{
  //
  // Draw clusters
  //
	
	if(!fGeo){
		Error("LoadClusters()", Form("Geometry not set for chamber %d. Please call first TRDChamber::SetGeometry().", fDet));
		return;
	}
//	Info("LoadClusters()", Form("clusters = 0x%x", clusters));
	if(!fRecPoints){
		fRecPoints = new TRDHits("clusters");
		fRecPoints->SetMarkerSize(1.);
		fRecPoints->SetMarkerStyle(24);
		fRecPoints->SetMarkerColor(6);
	} else fRecPoints->Reset();

	Float_t q, z0;
  Double_t cloc[3], cglo[3];
	
	z0 = fGeo->GetTime0(fPla);	
	AliTRDcluster *c=0x0;
	for(int iclus=0; iclus<clusters->GetEntriesFast(); iclus++){
		c = (AliTRDcluster*)clusters->UncheckedAt(iclus);
		cloc[2] = c->GetZ(); //x
		cloc[1] = c->GetY(); //y
		cloc[0] = z0 - c->GetX(); //z
		q = c->GetQ();
		fGeo->RotateBack(fDet,cloc,cglo);
		fRecPoints->SetNextPoint(cglo[0], cglo[1], cglo[2]);
		fRecPoints->SetPointId(this);
	}
	fLoadRecPoints = kTRUE;
}

//________________________________________________________
void TRDChamber::LoadDigits(AliTRDdigitsManager *digits)
{
  //
  // Draw digits
  //
	if(!fGeo){
		Error("LoadDigits()", Form("Geometry not set for chamber %d. Please call first TRDChamber::SetGeometry().", fDet));
		return;
	}
//	Info("LoadDigits()", Form("digits =0x%x", digits));
	
	if(!fDigits) fDigits = new TRDDigits(this);
	else fDigits->Reset();
	
	fDigits->SetData(digits);
	fLoadDigits = kTRUE;
}

//________________________________________________________
void TRDChamber::AddHit(AliTRDhit *hit)
{
  //
  // Draw hits
  //
//	Info("AddHit()", Form("%s", GetName()));

	if(!fHits){
		fHits = new TRDHits("hits");
		fHits->SetMarkerSize(.1);
		fHits->SetMarkerColor(2);
	}
	
	fHits->SetNextPoint(hit->X(), hit->Y(), hit->Z());
	fHits->SetPointId(this);
	fLoadHits = kTRUE;
}

//________________________________________________________
void TRDChamber::LoadTracklets(TObjArray *tracks)
{
  //
  // Draw tracks
  //
	if(!fGeo){
		Error("LoadTracklets()", Form("Geometry not set for chamber %d. Please call first TRDChamber::SetGeometry().", fDet));
		return;
	}
//	Info("LoadTracklets()", Form("tracks = 0x%x", tracks));
	
	if(!fTracklets){
		fTracklets = new std::vector<Reve::Track*>;
	} else fTracklets->clear();
	
	
	AliTRDmcmTracklet *trk = 0x0;
	Double_t cloc[3], cglo[3];
	for(int itrk=0; itrk<tracks->GetEntries();itrk++){
		trk = (AliTRDmcmTracklet*)tracks->At(itrk);
		trk->MakeTrackletGraph(fGeo,.5);
		fTracklets->push_back(new Reve::Track());
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
void	TRDChamber::Paint(Option_t* option)
{
/*	Info("Paint()", Form("%s", GetName()));*/
	if(!fRnrElement) return;
	if(fDigits && fRnrDigits){
		if(kDigitsNeedRecompute){
			fDigits->ComputeRepresentation();
			kDigitsNeedRecompute = kFALSE;
		}
		fDigits->Paint(option);
	}
	if(fRecPoints && fRnrRecPoints) fRecPoints->GetObject()->Paint(option);
	if(fHits && fRnrHits) fHits->GetObject()->Paint(option);
	if(fTracklets && fRnrTracklets){
		for(vector<Reve::Track*>::iterator i=fTracklets->begin(); i != fTracklets->end(); ++i) (*i)->Paint(option);
	}
}

//________________________________________________________
void	TRDChamber::Reset()
{
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

//________________________________________________________
void TRDChamber::SetGeometry(AliTRDgeometry *geo)
{
	fGeo = geo;
		
	fPla = fGeo->GetPlane(fDet);
	fX0 = fGeo->GetTime0(fPla);
	
	AliTRDCommonParam *parcom = AliTRDCommonParam::Instance();
	fPadPlane = parcom->GetPadPlane(fPla,fGeo->GetChamber(fDet));
	rowMax = fPadPlane->GetNrows();
	colMax = fPadPlane->GetNcols();
}

