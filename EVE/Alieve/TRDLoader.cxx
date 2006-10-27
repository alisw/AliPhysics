#include "TRDLoader.h"
#include "TRDModuleImp.h"

#include <Reve/RGTopFrame.h>
#include <Reve/RGValuators.h>

#include <TSystem.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGFileDialog.h>
#include <TGListTree.h>
#include <TGToolTip.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliCDBManager.h"

#include "AliTRDv1.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"

#include <algorithm>

using namespace Reve;
using namespace Alieve;
using namespace std;
class AliTRDdataArrayI;

ClassImp(Alieve::TRDLoader)
ClassImp(Alieve::TRDLoaderEditor)

///////////////////////////////////////////////////////////
/////////////         TRDLoader       /////////////////////
///////////////////////////////////////////////////////////


//________________________________________________________
TRDLoader::TRDLoader(const Text_t* n, const Text_t* t) : Reve::RenderElementListBase(), TNamed(n,t)
{
//	fChambers = 0x0;
	fRunLoader = 0x0;
//	fTRDLoader = 0x0;
	
	kLoadHits = kFALSE;
	kLoadDigits = kFALSE;
	kLoadClusters = kFALSE;
	kLoadTracks = kFALSE;
	
	AliCDBManager *fCDBManager=AliCDBManager::Instance();
	fCDBManager->SetDefaultStorage("local://$ALICE_ROOT");
	fCDBManager->SetRun(0);
}

//________________________________________________________
TRDLoader::~TRDLoader()
{
//	if(fChambers) {fChambers->clear(); delete fChambers;}
}

//________________________________________________________
template<class T>
class ID
{
public:
	ID( const int value ) : id(value) {}
	bool operator()(const T &t) const {
		return ((dynamic_cast<TRDModule*>(t))->GetID() == id);
	}
private:
	const int id;
};
void	TRDLoader::AddChambers(const int sm, const int stk, const int ly)
{
	Int_t ism_start = (sm == -1) ?  0 : sm;
	Int_t ism_stop  = (sm == -1) ? 18 : sm+1;
	Int_t istk_start= (stk == -1)?  0 : stk;
	Int_t istk_stop = (stk == -1)?  5 : stk+1;
	Int_t ily_start = (ly == -1) ?  0 : ly;
	Int_t ily_stop  = (ly == -1) ?  6 : ly+1;

	lpRE_i ichmb;
	ichmb = fChildren.begin();
	while(ichmb != fChildren.end()){
		(*ichmb)->SetRnrElement(kFALSE);
		ichmb++;
	}

	TRDNode *SM=0x0, *STK=0x0;
	TRDChamber *CHMB = 0x0;
	int det;
	for(int ism=ism_start; ism<ism_stop; ism++){
		ichmb = find_if(fChildren.begin(), fChildren.end(), ID<RenderElement*>(ism));
		if(ichmb != fChildren.end()){
			SM = (TRDNode*)(*ichmb);
			SM->SetRnrElement(kTRUE);
		}else{
			gReve->AddRenderElement(this, SM = new TRDNode("SM", ism));
			SM->FindListTreeItem(gReve->GetListTree())->SetTipText(Form("Supermodule %2d", ism));
		}
		for(int istk=istk_start; istk<istk_stop; istk++){
			ichmb = find_if(SM->begin(), SM->end(), ID<RenderElement*>(istk));
			if(ichmb != SM->end()){
				STK = (TRDNode*)(*ichmb);
				STK->SetRnrElement(kTRUE);
			}else{
				gReve->AddRenderElement(SM, STK = new TRDNode("Stack", istk));
				STK->FindListTreeItem(gReve->GetListTree())->SetTipText(Form("SM %2d Stack %1d", ism, istk));
			}
			for(int ily=ily_start; ily<ily_stop; ily++){
				det = AliTRDgeometry::GetDetector(ily, istk, ism);
				ichmb = find_if(STK->begin(), STK->end(), ID<RenderElement*>(det));
				if(ichmb != STK->end()) (*ichmb)->SetRnrElement(kTRUE);
				else{
					gReve->AddRenderElement(STK, CHMB = new TRDChamber(det));
					CHMB->FindListTreeItem(gReve->GetListTree())->SetTipText(Form("SM %2d Stack %1d Layer %1d", ism, istk, ily));
				}
			}
		}
	}
	gReve->Redraw3D();
}

//________________________________________________________
TRDChamber*	TRDLoader::GetChamber(const int d)
{
	lpRE_i ism, istack, ichmb;
	AliTRDgeometry *geo = fTRD->GetGeometry();
	
	ism = find_if(fChildren.begin(), fChildren.end(), ID<RenderElement*>(geo->GetSector(d)));
	if(ism == fChildren.end()) return 0x0;
	istack = find_if(((TRDNode*)(*ism))->begin(), ((TRDNode*)(*ism))->end(), ID<RenderElement*>(geo->GetChamber(d)));
	if(istack == ((TRDNode*)(*ism))->end()) return 0x0;
	ichmb = find_if(((TRDNode*)(*istack))->begin(), ((TRDNode*)(*istack))->end(), ID<RenderElement*>(d));
	if(ichmb == ((TRDNode*)(*istack))->end()) return 0x0;
	return dynamic_cast<TRDChamber*>(*ichmb);
}

//________________________________________________________
Bool_t	TRDLoader::GoToEvent(const int ev)
{
	Info("GoToEvent", Form("Event = %d", ev));
	if(!fRunLoader) return kFALSE;
	fRunLoader->UnloadAll("TRD");
	Unload();
	
	fRunLoader->GetEvent(ev);
	TTree *t = 0x0;
	if(kLoadHits){
		fRunLoader->LoadHits("TRD", "READ");
		t = fRunLoader->GetTreeH("TRD", kFALSE);
		if(!t) return kFALSE;
		fTRD->SetTreeAddress();
		if(!LoadHits(t)) return kFALSE;
	}
	if(kLoadDigits){
		fRunLoader->LoadDigits("TRD", "READ");
		t = fRunLoader->GetTreeD("TRD", kFALSE);
		if(!t) return kFALSE;
		fTRD->SetTreeAddress();
		if(!LoadDigits(t)) return kFALSE;
	}
	if(kLoadClusters){
		fRunLoader->LoadRecPoints("TRD", "READ");
		t = fRunLoader->GetTreeR("TRD", kFALSE);
		if(!t) return kFALSE;
		if(!LoadClusters(t)) return kFALSE;
	}
	if(kLoadTracks){
		fRunLoader->LoadTracks("TRD", "READ");
		t = fRunLoader->GetTreeT("TRD", kFALSE);
		if(!t) return kFALSE;
		if(!LoadTracklets(t)) return kFALSE;
	}
	gReve->Redraw3D();
	return kTRUE;
}


//________________________________________________________
Bool_t	TRDLoader::LoadClusters(TTree *tC)
{
	Info("LoadClusters()", Form("Clusters tree 0x%x", tC));
	if(!fChildren.size()) return kTRUE;

	TObjArray *clusters = new TObjArray();
	tC->SetBranchAddress("TRDcluster",&clusters);

	TRDChamber *chmb = 0x0;	
	AliTRDcluster *c=0x0;
	for(int idet=0; idet<540; idet++){
		if(!tC->GetEntry(idet)) continue;
		if(clusters->GetEntriesFast()) c = (AliTRDcluster*)clusters->UncheckedAt(0);
		if((chmb = GetChamber(c->GetDetector()))) chmb->LoadClusters(clusters);
	}
	return kTRUE;
}


//________________________________________________________
Bool_t	TRDLoader::LoadDigits(TTree *tD)
{
	Info("LoadDigits()", Form("Digits tree 0x%x", tD));
	
	if(!fChildren.size()) return kTRUE;
	
	TRDChamber *chmb;
	AliTRDdataArrayI *digits;
	AliTRDdigitsManager dm;
	dm.ReadDigits(tD);
	for(int idet=0; idet<540; idet++){
		if(!(chmb=GetChamber(idet))) continue;
		digits = dm.GetDigits(idet);
		if(!digits) continue;
		chmb->LoadDigits(digits);
	}
	return kTRUE;
}


//________________________________________________________
Bool_t	TRDLoader::LoadHits(TTree *tH)
{
	Info("LoadHits()", Form("Hits tree 0x%x", tH));
	if(!fChildren.size()) return kTRUE;
	
	TRDChamber *chmb = 0x0;
	AliTRDhit *hit;
	Int_t d;
	for(int iTrack=0; iTrack<tH->GetEntries(); iTrack++){
		gAlice->ResetHits();
		if(!tH->GetEvent(iTrack)) continue;
		hit = (AliTRDhit*)fTRD->FirstHit(-1);
		if(!hit) continue;
		d = hit->GetDetector();
		chmb = GetChamber(d);
		while(hit){
			if(d != hit->GetDetector()){
				d = hit->GetDetector();
				chmb = GetChamber(d);
			}
			if(chmb) chmb->AddHit(hit);
			hit = (AliTRDhit*)fTRD->NextHit();
		}
	}
	return kTRUE;
}

//________________________________________________________
Bool_t	TRDLoader::LoadTracklets(TTree *tT)
{
	Info("LoadTracklets()", Form("Tracks tree 0x%x", tT));
	if(!fChildren.size()) return kTRUE;

	TObjArray *tracks = new TObjArray();
	tT->SetBranchAddress("TRDmcmTracklet",&tracks);
	
	TRDChamber *chmb = 0x0;
	AliTRDmcmTracklet *trk=0x0;
	for(int idet=0; idet<540; idet++){
		if(!tT->GetEntry(idet)) continue;
		if(tracks->GetEntriesFast()) trk = (AliTRDmcmTracklet*)tracks->UncheckedAt(0);
		if((chmb = GetChamber(trk->GetDetector()))) chmb->LoadTracklets(tracks);
	}
	
	return kTRUE;
}
	
//________________________________________________________
Bool_t	TRDLoader::Open(const char *filename)
{
 	fRunLoader = AliRunLoader::GetRunLoader();
	if(!fRunLoader) fRunLoader = AliRunLoader::Open(filename,
				AliConfig::GetDefaultEventFolderName(),"read");
	if(!fRunLoader){
		Error("Open()", "Couldn't find run loader");
		return kFALSE;
	}
	gAlice = fRunLoader->GetAliRun();
  if(!gAlice) fRunLoader->LoadgAlice();
	if(!gAlice){
		Error("Open()", "Couldn't find gAlice object");
		return kFALSE;
	}
	fTRD = (AliTRDv1*)gAlice->GetDetector("TRD");
	if(!fTRD){
		Error("Open()", "Couldn't find TRD");
		return kFALSE;
	}

	return kTRUE;
}

//________________________________________________________
void TRDLoader::Paint(Option_t *option)
{
	lpRE_i ichmb = fChildren.begin();
	while(ichmb != fChildren.end()){
		(dynamic_cast<TRDModule*>(*ichmb))->Paint(option);
		ichmb++;
	}
}

//________________________________________________________
void TRDLoader::Unload()
{
	lpRE_i ichmb = fChildren.begin();
	while(ichmb != fChildren.end()){
		(dynamic_cast<TRDModule*>(*ichmb))->Reset();
		ichmb++;
	}
}

///////////////////////////////////////////////////////////
/////////////   TRDLoaderEditor       /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDLoaderEditor::TRDLoaderEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
	MakeTitle("TRDLoader");
	
  fFile = 0x0;
	fOpenFile = 0x0;
	Int_t labelW = 42;
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
    TGLabel* l = new TGLabel(g, "File: ");
    g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
    f->AddFrame(g);
    fFile = new TGTextEntry(f);
    fFile->SetWidth(140);
    f->AddFrame(fFile);
/*    
		fFile->Connect("DoubleClicked()",
		   "Alieve::TPCLoaderEditor", this, "FileSelect()");
    fFile->Connect("TextChanged(const char *)",
		   "Alieve::TPCLoaderEditor", this, "FileChanged()");
*/    
		fOpenFile = new TGTextButton(f, "Open");
    f->AddFrame(fOpenFile);
    fOpenFile->Connect("Clicked()",
		       "Alieve::TRDLoaderEditor", this, "FileOpen()");
		AddFrame(f);
  }

  fEvent = new RGValuator(this, "Event:", 110, 0);
  fEvent->SetShowSlider(kFALSE);
  fEvent->SetLabelWidth(labelW);
  fEvent->SetNELength(6);
  fEvent->Build();
  fEvent->SetLimits(-1, 1000);
  fEvent->SetToolTip("Current event number");
/*  fEvent->Connect("ValueSet(Double_t)",
		  "Alieve::TPCLoaderEditor", this, "DoEvent()");
*/  
	AddFrame(fEvent);


	// "Chamber(s) selector" group frame
	TGGroupFrame *fGroupFrame1974 = new TGGroupFrame(this,"Chamber(s) selector");
	
	TGVerticalFrame *fVerticalFrame1974 = new TGVerticalFrame(fGroupFrame1974,128,26,kVerticalFrame);

	// horizontal frame
	TGHorizontalFrame *fHorizontalFrame2020 = new TGHorizontalFrame(fVerticalFrame1974,128,26,kHorizontalFrame);
	TGLabel *fLabel2170 = new TGLabel(fHorizontalFrame2020,"SM   ");
	fLabel2170->SetTextJustify(kTextLeft | kTextCenterY);
	fHorizontalFrame2020->AddFrame(fLabel2170, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	fSMNumber = new TGNumberEntry(fHorizontalFrame2020, (Int_t) -1, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -1, 17);
	fHorizontalFrame2020->AddFrame(fSMNumber, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));

	fVerticalFrame1974->AddFrame(fHorizontalFrame2020, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	// horizontal frame
	TGHorizontalFrame *fHorizontalFrame2025 = new TGHorizontalFrame(fVerticalFrame1974,128,28,kHorizontalFrame);
	TGLabel *fLabel2186 = new TGLabel(fHorizontalFrame2025,"Stack");
	fLabel2186->SetTextJustify(kTextLeft | kTextCenterY);
	fHorizontalFrame2025->AddFrame(fLabel2186, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fStackNumber = new TGNumberEntry(fHorizontalFrame2025, (Int_t) -1, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -1, 4);
	fHorizontalFrame2025->AddFrame(fStackNumber, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));

	fVerticalFrame1974->AddFrame(fHorizontalFrame2025, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	// horizontal frame
	TGHorizontalFrame *fHorizontalFrame2030 = new TGHorizontalFrame(fVerticalFrame1974,128,26,kHorizontalFrame);
	TGLabel *fLabel2208 = new TGLabel(fHorizontalFrame2030,"Plane");
	fLabel2208->SetTextJustify(kTextLeft | kTextCenterY);
	fHorizontalFrame2030->AddFrame(fLabel2208, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fPlaneNumber = new TGNumberEntry(fHorizontalFrame2030, (Int_t) -1, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -1, 5);
	fHorizontalFrame2030->AddFrame(fPlaneNumber, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));

	fVerticalFrame1974->AddFrame(fHorizontalFrame2030, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	
	fGroupFrame1974->AddFrame(fVerticalFrame1974, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY | kLHintsExpandX,2,2,2,2));

	TGTextButton *fTextButton2037 = new TGTextButton(fGroupFrame1974,"Select");
	fTextButton2037->SetTextJustify(36);
	fTextButton2037->Resize(128,22);
	fGroupFrame1974->AddFrame(fTextButton2037, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterX | kLHintsExpandY,2,2,2,2));
	fTextButton2037->Connect("Clicked()",
					"Alieve::TRDLoaderEditor", this, "AddChambers()");

	fGroupFrame1974->SetLayoutManager(new TGHorizontalLayout(fGroupFrame1974));
//	fGroupFrame1974->Resize(164,150);
	AddFrame(fGroupFrame1974, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	
	// "Data selector" group frame
	TGGroupFrame *fGroupFrame1987 = new TGGroupFrame(this,"Data selector");
	fLoadHits = new TGCheckButton(fGroupFrame1987,"  Hits");
	fGroupFrame1987->AddFrame(fLoadHits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fLoadDigits = new TGCheckButton(fGroupFrame1987,"  Digits");
	fGroupFrame1987->AddFrame(fLoadDigits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fLoadClusters = new TGCheckButton(fGroupFrame1987,"  Clusters");
	fGroupFrame1987->AddFrame(fLoadClusters, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fLoadESDs = new TGCheckButton(fGroupFrame1987,"  Tracklets ");
	fGroupFrame1987->AddFrame(fLoadESDs, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	fGroupFrame1987->SetLayoutManager(new TGVerticalLayout(fGroupFrame1987));
//	fGroupFrame1987->Resize(164,116);
	AddFrame(fGroupFrame1987, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	TGTextButton *fTextButton2004 = new TGTextButton(this,"Load");
	fTextButton2004->SetTextJustify(36);
	fTextButton2004->Resize(164,22);
	AddFrame(fTextButton2004, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	fTextButton2004->Connect("Clicked()",
					"Alieve::TRDLoaderEditor", this, "Load()");
}

//________________________________________________________
TRDLoaderEditor::~TRDLoaderEditor()
{}

//_________________________________________________________
void TRDLoaderEditor::SetModel(TObject* obj)
{
	fM = dynamic_cast<TRDLoader*>(obj);
}

//________________________________________________________
void TRDLoaderEditor::AddChambers()
{
	fM->AddChambers((int)fSMNumber->GetNumber(),
		(int)fStackNumber->GetNumber(),
		(int)fPlaneNumber->GetNumber());
}

//________________________________________________________
void TRDLoaderEditor::FileOpen()
{
  TGFileInfo fi;
/*  fi.fIniDir    = StrDup(gSystem->DirName (fM->fFile));
  fi.fFilename  = StrDup(gSystem->BaseName(fM->fFile));
  fi.fFileTypes = tpcfiletypes;
*/
  new TGFileDialog(fClient->GetRoot(), gReve, kFDOpen, &fi);
  if (!fi.fFilename) return;

  fFile->SetToolTipText(gSystem->DirName (fi.fFilename));
  fFile->SetText       (gSystem->BaseName(fi.fFilename));

	fM->Open(gSystem->BaseName(fi.fFilename));
}

//________________________________________________________
void TRDLoaderEditor::Load()
{
	fM->kLoadHits     = kFALSE;
	fM->kLoadDigits   = kFALSE;
	fM->kLoadClusters = kFALSE;
	fM->kLoadTracks   = kFALSE;
	if(fLoadHits->IsDown()) fM->kLoadHits = kTRUE;
	if(fLoadDigits->IsDown()) fM->kLoadDigits = kTRUE;
	if(fLoadClusters->IsDown()) fM->kLoadClusters = kTRUE;
	if(fLoadESDs->IsDown()) fM->kLoadTracks = kTRUE;
	
	fM->GoToEvent((int)fEvent->GetEntry()->GetNumber());
}

