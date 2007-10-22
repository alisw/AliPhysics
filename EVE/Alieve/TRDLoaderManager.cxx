#include "TRDLoaderManager.h"
#include "TRDLoader.h"
#include "TRDLoaderImp.h"

#include <Reve/ReveManager.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <TGListTree.h>
#include <TGString.h>
#include <TGToolTip.h>
#include <TClonesArray.h>

#include "AliLog.h"

using namespace Reve;
using namespace Alieve;
using namespace std;


ClassImp(Alieve::TRDLoaderManager)
ClassImp(Alieve::TRDLoaderManagerEditor)

///////////////////////////////////////////////////////////
/////////////         TRDLoaderManager       //////////////
///////////////////////////////////////////////////////////


//________________________________________________________
TRDLoaderManager::TRDLoaderManager(const Text_t* n, const Text_t* t) : Reve::RenderElementList(n, t)
{

}

//________________________________________________________
TRDLoaderManager::~TRDLoaderManager()
{

}

//________________________________________________________
void	TRDLoaderManager::Add(Int_t type, const Text_t *name, const Text_t *title)
{
	//Info("Add()", Form("type %d, name %s, title %s", type, name, title));
	Alieve::TRDLoader *trdl = 0x0;
	switch(type){
	case 0:
		//fChildren.push_back(new TRDLoaderSim(name, title));
		gReve->AddRenderElement(trdl = new TRDLoaderSim(name, title), this);
		((TRDLoaderSim*)trdl)->FindListTreeItem(gReve->GetListTree())->SetTipText(title);
		break;	
	case 1:
	case 2:
	case 3:
		//fChildren.push_back(new TRDLoader(name, title));
	  gReve->AddRenderElement(trdl = new TRDLoader(name, title), this);
		trdl->FindListTreeItem(gReve->GetListTree())->SetTipText(title);
		trdl->SetDataType((Alieve::TRDDataTypes)type);
		break;
	case 4:
	case 5:
		//fChildren.push_back(new TRDLoaderRaw(name, title));
	  gReve->AddRenderElement(trdl = new TRDLoaderRaw(name, title), this);
		((TRDLoaderRaw*)trdl)->FindListTreeItem(gReve->GetListTree())->SetTipText(title);
		trdl->SetDataType((Alieve::TRDDataTypes)type);
		break;
	}
	
	gReve->Redraw3D();
}



//________________________________________________________
void TRDLoaderManager::Paint(Option_t *option)
{
	List_i ichmb = fChildren.begin();
	while(ichmb != fChildren.end()){
		(dynamic_cast<TRDLoader*>(*ichmb))->Paint(option);
		ichmb++;
	}
}

//________________________________________________________
void	TRDLoaderManager::Remove(Int_t entry)
{
	//printf("TRDLoaderManager::Remove(%d)\n", entry);
	List_i it = fChildren.begin();
	for(int i=0; i<entry; i++) it++;
	gReve->RemoveRenderElement((*it), this);
	fChildren.erase(it);
}

///////////////////////////////////////////////////////////
/////////////   TRDLoaderManagerEditor       //////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDLoaderManagerEditor::TRDLoaderManagerEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
	MakeTitle("TRDLoaderManager");

//	gClient->GetColorByName("#ffffff", bg);
//	ChangeBackground(bg);

	// control frame - always there
	TGHorizontalFrame *fHorizontalFrame539 = new TGHorizontalFrame(this, 300, 26, kHorizontalFrame);//, bg);

	TGLabel *fLabel546 = new TGLabel(fHorizontalFrame539,"Register Loader",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame);//, bg);
	fLabel546->SetTextJustify(36);
	fHorizontalFrame539->AddFrame(fLabel546, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));

	// combo box
	fSelector = new TGComboBox(fHorizontalFrame539,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
	fSelector->AddEntry("MC (gAlice) ",0);
	fSelector->AddEntry("Digits ",1);
	fSelector->AddEntry("Clusters ",2);
	fSelector->AddEntry("Tracklets ",3);
	fSelector->AddEntry("Raw (ROOT) ",4);
	fSelector->AddEntry("Raw (DATE) ",5);
	fSelector->Resize(136,22);
	fHorizontalFrame539->AddFrame(fSelector, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));

	fAdd = new TGTextButton(fHorizontalFrame539, "Add");
	fAdd->SetTextJustify(36);
	fAdd->Resize(31,22);
	fAdd->SetToolTipText("Add selected loader to list");
	fAdd->Connect("Clicked()", "Alieve::TRDLoaderManagerEditor", this, "Add()");
	fHorizontalFrame539->AddFrame(fAdd, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));
	AddFrame(fHorizontalFrame539, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));

	fGroupFrame = 0x0;
	fRemove     = 0x0;
}

//________________________________________________________
TRDLoaderManagerEditor::~TRDLoaderManagerEditor()
{
	
}


//_________________________________________________________
void TRDLoaderManagerEditor::Add()
{
	TGTextLBEntry *entry = (TGTextLBEntry*)fSelector->GetSelectedEntry();
	if(!entry){
		AliWarning("Select first the loader type that you want to use from the drop down list.");
		return;
	}

	if(!fGroupFrame){
		// "TRD Loaders" group frame
		fGroupFrame = new TGGroupFrame(this,"TRD Loaders",kVerticalFrame,TGGroupFrame::GetDefaultGC()(),TGGroupFrame::GetDefaultFontStruct());//, bg);
		fGroupFrame->SetLayoutManager(new TGVerticalLayout(fGroupFrame));
		fGroupFrame->Resize(300,128);
		AddFrame(fGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2,2,2,2));

		fRemove = new TClonesArray("TGTextButton", 3);
	}

	char *title[] = {"MC loader", "Single file loader", "Raw data loader"};
	// char *color[] = {"#ff0000", "#0000ff", "#59d454"};
	int id = fSelector->GetSelected(), type;
	switch(id){
	case 1:
	case 2:
	case 3:
		type = 1;
		break;
	case 4:	
	case 5:
		type = 2;
		break;
	default:
		type = 0;
		break;
	}

	
	// horizontal frame
	TGHorizontalFrame *fHorizontalFrame = new TGHorizontalFrame(fGroupFrame, 264, 26, kHorizontalFrame);//, bg);

// 	TGFont *ufont = gClient->GetFont("-*-helvetica-(null)-*-*-0-*-*-*-*-*-*-*");
// 	TGGC   *uGC;           // will reflect user GC changes
// 	// graphics context changes
// 	GCValues_t vall717;
// 	vall717.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
// 	gClient->GetColorByName(color[type], vall717.fForeground);
// 	gClient->GetColorByName("#c0c0c0", vall717.fBackground);
// 	vall717.fFillStyle = kFillSolid;
// 	vall717.fFont = ufont->GetFontHandle();
// 	vall717.fGraphicsExposures = kFALSE;
// 	uGC = gClient->GetGC(&vall717, kTRUE);
	
	TGLabel *fLabel717 = new TGLabel(fHorizontalFrame, entry->GetText()->GetString()/*, uGC->GetGC(), ufont->GetFontStruct(), kChildFrame*/);//, bg);
	fLabel717->SetTextJustify(36);
	fHorizontalFrame->AddFrame(fLabel717, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	
	Int_t nbutton = fM->fChildren.size();
	fRemoveButton = new((*fRemove)[nbutton]) TGTextButton(fHorizontalFrame, "Remove", nbutton);
	fRemoveButton->SetTextJustify(36);
	fRemoveButton->Resize(53,22);
	fRemoveButton->Connect("Clicked()", "Alieve::TRDLoaderManagerEditor", this, Form("Remove(=%d)", nbutton));
	fRemoveButton->SetToolTipText(Form("Remove %s Loader", entry->GetText()->GetString()));
	fHorizontalFrame->AddFrame(fRemoveButton, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));

	fGroupFrame->AddFrame(fHorizontalFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));


	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();

	fM->Add(id, entry->GetText()->GetString(), title[type]);
}


//_________________________________________________________
void TRDLoaderManagerEditor::Remove(Int_t entry)
{
	TIterator *it = fGroupFrame->GetList()->MakeIterator();
	int ientry = 0;
	while(/*TGFrame *f=(TGFrame*)*/it->Next()){
		//printf("%s\n", f->IsA()->GetName());
		if(entry == ientry){
			//fGroupFrame->RemoveFrame(f);
			break;
		}
		ientry++;
	}


	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();

	//fM->Remove(entry);
}

//_________________________________________________________
void TRDLoaderManagerEditor::SetModel(TObject* obj)
{
	fM = dynamic_cast<TRDLoaderManager*>(obj);
}

