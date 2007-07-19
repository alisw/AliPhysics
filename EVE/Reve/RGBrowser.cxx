#include "RGBrowser.h"
#include "RGTopFrame.h"
#include "Reve.h"
#include "RGEditor.h"
#include "VSDSelector.h"
#include <Reve/PointSet.h>
#include <Reve/Track.h>

#include <Riostream.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TRint.h>
#include <TVirtualX.h>
#include <TEnv.h>

#include <TApplication.h>
#include <TFile.h>
#include <TEventList.h>
#include <TClassMenuItem.h>

#include <TColor.h>
#include <TPolyMarker3D.h>

#include <TGCanvas.h>
#include <TGSplitter.h>
#include <TGStatusBar.h>
#include <TGMenu.h>
#include <TGPicture.h>
#include <TGToolBar.h>
#include <TGLabel.h>
#include <TGXYLayout.h>
#include <TGNumberEntry.h>
#include <KeySymbols.h>

#include <TGLSAViewer.h>
#include <TGLSAFrame.h>
#include <TGTab.h>

#include <TGeoVolume.h>
#include <TGeoNode.h>

using namespace Reve;
using namespace Reve;

/**************************************************************************/

void RGBrowser::SetupCintExport(TClass* cl)
{
 
  TList* l = cl->GetMenuList();
  TClassMenuItem* n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction, cl,
					 "Export to CINT", "ExportToCINT", this, "const char*,TObject*", 1);

  l->AddFirst(n);
}

void RGBrowser::CalculateReparentXY(TGObject* parent, Int_t& x, Int_t& y)
{
  UInt_t   w, h;
  Window_t childdum;
  gVirtualX->GetWindowSize(parent->GetId(), x, y, w, h);
  gVirtualX->TranslateCoordinates(parent->GetId(),
				  gClient->GetDefaultRoot()->GetId(),
				  0, 0, x, y, childdum);
}

/**************************************************************************/

RGBrowser::RGBrowser(const TGWindow *p, UInt_t w, UInt_t h) :
  TGCompositeFrame(p, w, h),
    
  fMainFrame(0), fV1(0), fV2(0),
  fSelectionFrame(0), fTreeView(0),
  fCanvasWindow(0), fDisplayFrame(0),  
  fListTree(0),
  fCtxMenu(0)
{
  fMainFrame = new TGCompositeFrame(this, 100, 10, kHorizontalFrame | kRaisedFrame);
  fMainFrame->SetCleanup(kDeepCleanup);
  fV1 = new TGVerticalFrame(fMainFrame, 250, 10, kSunkenFrame | kFixedWidth);
  fV2 = new TGVerticalFrame(fMainFrame,  50, 10, kSunkenFrame);

  TGLayoutHints *lo;
  lo = new TGLayoutHints(kLHintsLeft | kLHintsExpandY,2,0,2,2);
  fMainFrame->AddFrame(fV1, lo);

  TGVSplitter *splitter = new TGVSplitter(fMainFrame);
  splitter->SetFrame(fV1, kTRUE);
  fMainFrame->AddFrame(splitter,
		       new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 1,1,2,2));
   
  lo = new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY,0,2,2,4);
  fMainFrame->AddFrame(fV2, lo);

  // selection frame
  fSelectionFrame = new TGCompositeFrame(fV1, 250, 10, kVerticalFrame);
  fTreeView = new TGCanvas(fSelectionFrame, 250, 10, kSunkenFrame | kDoubleBorder);
  fListTree = new TGListTree(fTreeView->GetViewPort(), 250, 10, kHorizontalFrame);
  fListTree->SetCanvas(fTreeView);
  fListTree->Associate(this);
  fListTree->SetColorMode(TGListTree::EColorMarkupMode(TGListTree::kColorUnderline | TGListTree::kColorBox));
  fListTree->SetAutoCheckBoxPic(kFALSE);
  fTreeView->SetContainer(fListTree);

  lo= new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
			2, 2, 2, 2);
  fSelectionFrame->AddFrame(fTreeView, lo);

  lo = new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY);
  fV1->AddFrame(fSelectionFrame, lo);
 
  // Classic look vars:
  fCanvasWindow = 0;
  fDisplayFrame = 0;
 
  //display frame

  lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  AddFrame(fMainFrame, lo);


  SetWindowName("Reve List Browser");
  MapSubwindows();
  //Resize(GetDefaultSize()); // this is used here to init layout algoritme

  //MapWindow();

  // popup menu
  
  fCtxMenu = new TContextMenu("", "");


  //-- CINT export now declared in RenderElement with *MENU*
  // SetupCintExport(PointSet::Class());
  // SetupCintExport(Track::Class());
  // SetupCintExport(TrackList::Class());
  
  fListTree->Connect("Clicked(TGListTreeItem*, Int_t, Int_t, Int_t)", "Reve::RGBrowser", 
		     this, "ItemClicked(TGListTreeItem*, Int_t, Int_t, Int_t)");  
  fListTree->Connect("DoubleClicked(TGListTreeItem*, Int_t)", "Reve::RGBrowser", 
		     this, "DbClickListItem(TGListTreeItem*,Int_t )"); 
  //fListTree->Connect("Clicked(TGListTreeItem*, Int_t)", "Reve::RGBrowser", 
  //		     this, "DisplayChildren(TGListTreeItem*, Int_t)");  

  //---------------------------------------------
  // WARNING ... this Connect goes to *gReve*!
  fListTree->Connect("Checked(TObject*,Bool_t)", "Reve::RGTopFrame",
		     gReve, "RenderElementChecked(TObject*, Bool_t)");
}

/**************************************************************************/

void RGBrowser::SetupClassicLook(RGEditor*& editor, TCanvas* glpad)
{
  fCanvasWindow = new TGCanvas(fV2, 25, 250);
  fDisplayFrame = new TGCompositeFrame(fCanvasWindow->GetViewPort(), 0, 0,kVerticalFrame, TGFrame::GetWhitePixel() );
  fCanvasWindow->SetContainer(fDisplayFrame);
  fDisplayFrame->SetCleanup(kDeepCleanup);

  fV2->AddFrame(fCanvasWindow, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 2, 2));
  fV2->MapSubwindows();

  editor = new RGEditor(glpad);
  editor->GetTGCanvas()->ChangeOptions(0);
  editor->SetWindowName("Reve Editor");
}

void RGBrowser::SetupEditorLook(RGEditor*& editor, TCanvas* glpad)
{
  fV2->SetEditDisabled(kEditEnable);
  fV2->SetEditable();
  editor = new RGEditor(glpad);
  editor->GetTGCanvas()->ChangeOptions(0);
  fV2->SetEditable(kEditDisable);
  fV2->SetEditable(kFALSE);

  TGFrameElement *el = 0;
  TIter next(fV2->GetList());
  while ((el = (TGFrameElement *)next())) {
      if (el->fFrame == editor)
         if (el->fLayout) {
            el->fLayout->SetLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY);
            el->fLayout->SetPadLeft(0);
            el->fLayout->SetPadRight(0);
            el->fLayout->SetPadTop(2);
            el->fLayout->SetPadBottom(2);
            break;
         }
  }

  fV2->MapSubwindows();
}

void RGBrowser::SetupGLViewerLook(RGEditor*& editor, TCanvas* glpad)
{

  TGFrameElement *el = 0;
  fV2->SetEditDisabled(kEditEnable);
  fV2->SetEditable();
  TGLSAViewer* v = new TGLSAViewer(fV2, glpad);
  v->GetFrame()->SetMinWidth(200);
  v->GetFrame()->SetCleanup(kNoCleanup);
  fV2->SetEditable(kEditDisable);
  fV2->SetEditable(kFALSE);
  glpad->SetViewer3D(v);
  TIter next2(fV2->GetList());
  while ((el = (TGFrameElement *)next2())) {
     if (el->fFrame == v->GetFrame()) {
         el->fLayout = new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY);
         break;
     }
  }
  fSelectionFrame->Resize(fSelectionFrame->GetWidth(), fSelectionFrame->GetHeight()/2);

  TGHSplitter *splitter = new TGHSplitter(fV1);
  fV1->AddFrame(splitter, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 2, 2, 0));

  fV1->SetEditDisabled(kEditEnable);
  fV1->SetEditable();
  editor = new RGEditor(glpad);
  editor->SetGlobal(kFALSE);
  editor->GetTGCanvas()->ChangeOptions(0);
  editor->ChangeOptions(editor->GetOptions() | kFixedHeight);
  fV1->SetEditable(kEditDisable);
  fV1->SetEditable(kFALSE);

  TIter next(fV1->GetList());
  while ((el = (TGFrameElement *)next())) {
      if (el->fFrame == editor)
         if (el->fLayout) {
            el->fLayout->SetLayoutHints(kLHintsTop | kLHintsExpandX);
            el->fLayout->SetPadLeft(0);
            el->fLayout->SetPadRight(2);
            el->fLayout->SetPadTop(2);
            el->fLayout->SetPadBottom(2);
            break;
         }
  }
  splitter->SetFrame(editor, kFALSE);
  fV1->MapSubwindows();
}


/**************************************************************************/
/**************************************************************************/

void RGBrowser::RedrawListTree()
{
  gClient->NeedRedraw(fListTree);
}

/**************************************************************************/

void RGBrowser::ItemClicked(TGListTreeItem *item, Int_t btn, Int_t x, Int_t y)
{
  //printf("ItemClicked item %s List %d btn=%d, x=%d, y=%d\n",
  //  item->GetText(),fDisplayFrame->GetList()->GetEntries(), btn, x, y);

  RenderElement* re = (RenderElement*)item->GetUserData();
  if(re == 0) return;
  TObject* obj = re->GetObject();

  if(btn == 3) {
    // If control pressed, show menu for renderelement itself.
    // event->fState & kKeyControlMask
    // ??? how do i get current event?
    if (obj) {
      fCtxMenu->Popup(x, y, obj);
    }
    return;
  }

  gReve->EditRenderElement(re);
}

void RGBrowser::DbClickListItem(TGListTreeItem* item, Int_t btn)
{
  static const Exc_t eH("RGBrowser::DbClickListItem ");

  // printf("dbclick item %s\n", item->GetText());
  RenderElement* re = (RenderElement*)item->GetUserData();
  if(re == 0) return;
  TObject* obj = re->GetObject();

  if (obj) {
    //	ListTreeHighlight(item);

    re->ExpandIntoListTree(fListTree, item);
    
    // browse geonodes
    if(obj->IsA()->InheritsFrom("TGeoNode")){
      TGeoNode* n = (TGeoNode*) obj->IsA()->DynamicCast( TGeoNode::Class(), obj );
      // initialization
      if(item->GetFirstChild() == 0 && n->GetNdaughters()){
	UpdateListItems(item, btn);
      }
    }
  }
}

/**************************************************************************/

void RGBrowser::ExportToCINT(Text_t* var_name, TObject* obj)
{
  const char* cname = obj->IsA()->GetName();
  gROOT->ProcessLine(Form("%s* %s = (%s*) %p;", cname, var_name, cname, obj));
}
/**************************************************************************/

void RGBrowser::UpdateListItems(TGListTreeItem* item, Int_t )
{
  if (item->GetUserData()) {
    //	ListTreeHighlight(item);
    RenderElement* re = (RenderElement*)item->GetUserData();
    TObject* obj = re->GetObject();

    // geometry tree
    if(obj->IsA()->InheritsFrom("TGeoNode")){
      // delete exisiting
      fListTree->DeleteChildren(item);
      TGeoNode* n = (TGeoNode*) obj->IsA()->DynamicCast( TGeoNode::Class(), obj );
      //printf("adding items\n");
      if (n->GetNdaughters()) {
	for (Int_t i=0; i< n->GetNdaughters(); i++) { 
	  TString title;
	  title.Form("%d : %s[%d]", i,
		     n->GetDaughter(i)->GetVolume()->GetName(),
		     n->GetDaughter(i)->GetNdaughters());

	  TGListTreeItem* child = fListTree->AddItem( item, title.Data());
	  child->SetUserData( n->GetDaughter(i));
	}
      }
    }
  }
}

/**************************************************************************/
