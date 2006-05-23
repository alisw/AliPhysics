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

RGBrowser::RGBrowser(const TGWindow *p, UInt_t w, UInt_t h)
  : TGCompositeFrame(p, w, h)
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
  
  fCtxMenu = new TContextMenu("Pepe", "Moroder");

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

void RGBrowser::SetupClassicLook()
{
  fCanvasWindow = new TGCanvas(fV2, 25, 250);
  fDisplayFrame = new TGCompositeFrame(fCanvasWindow->GetViewPort(), 0, 0,kVerticalFrame, TGFrame::GetWhitePixel() );
  fCanvasWindow->SetContainer(fDisplayFrame);
  fDisplayFrame->SetCleanup(kDeepCleanup);

  fV2->AddFrame(fCanvasWindow, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 2, 2));
  fV2->MapSubwindows();
}



void RGBrowser::SetupEditorLook(RGEditor* editor)
{
  editor->UnmapWindow();
  fV2->AddFrame(editor, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 2, 2));
  Int_t x, y;
  CalculateReparentXY(fV2, x, y);
  editor->ReparentWindow(fV2, x, y);

  fV2->MapSubwindows();
}

void RGBrowser::SetupGLViewerLook(RGEditor* editor, TVirtualPad* glpad)
{
  TGLayoutHints *lo;

  TGLSAViewer* v = new TGLSAViewer(fV2, glpad);
  v->GetFrame()->SetMinWidth(200);
  lo = new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY);
  fV2->AddFrame(v->GetFrame(), lo);
  glpad->SetViewer3D(v);

  fSelectionFrame->Resize(fSelectionFrame->GetWidth(), fSelectionFrame->GetHeight()/2);

  TGHSplitter *splitter = new TGHSplitter(fV1);
  lo = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 2, 2, 0);
  fSelectionFrame->AddFrame(splitter, lo);
   
  editor->UnmapWindow();
  editor->ChangeOptions(editor->GetOptions() | kFixedHeight);
  lo = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,2,2,2);
  fV1->AddFrame(editor, lo);
  Int_t x, y;
  CalculateReparentXY(fV1, x, y);
  editor->ReparentWindow(fV1, x, y);

  splitter->SetFrame(editor, kFALSE);

  fV1->MapSubwindows();
  v->GetFrame()->MapWindow();
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

  // A pathetic hack to get at least a bit of color coordination
  // for RenderElementObjPtr.
  if(item->GetColor() != re->GetMainColor()) {
    item->SetColor(re->GetMainColor());
    fListTree->GetClient()->NeedRedraw(fListTree);
  }

  if(btn == 3) {
    if (obj) {
      fCtxMenu->Popup(x, y, obj);
    }
    return;
  }

  gReve->EditRenderElement(re);

  // This only available in classic look.
  // Still working but slowly drifting towards obscurity (4.2006).
  DisplayChildren(item, btn);
}

void RGBrowser::DbClickListItem(TGListTreeItem* item, Int_t btn)
{
  static const Exc_t eH("RGBrowser::DbClickListItem ");

  printf("dbclick item %s\n", item->GetText());
  RenderElement* re = (RenderElement*)item->GetUserData();
  if(re == 0) return;
  TObject* obj = re->GetObject();

  if (obj) {
    //	ListTreeHighlight(item);

    {
      RenderElementListBase* rel = dynamic_cast<RenderElementListBase*>(re);
      if(rel != 0) {
	Int_t ni = rel->ExpandIntoListTree(fListTree, item);
	printf("%s expanded by %d\n", eH.Data(), ni);
      }
    }
    
    // browse geonodes
    if(obj->IsA()->InheritsFrom("TGeoNode")){
      TGeoNode* n = (TGeoNode*) obj->IsA()->DynamicCast( TGeoNode::Class(), obj );
      // initialization
      if(item->GetFirstChild() == 0 && n->GetNdaughters()){
	UpdateListItems(item, btn);
      }
    }
  }
  DisplayChildren(item,0);
}

/**************************************************************************/

void RGBrowser::ExportToCINT(Text_t* var_name, TObject* obj)
{
  const char* cname = obj->IsA()->GetName();
  gROOT->ProcessLine(Form("%s* %s = (%s*) %p;", cname, var_name, cname, obj));
}

void RGBrowser::DisplayChildren(TGListTreeItem *item, Int_t btn)
{
  // Only classic mode provides direct children editing.
  if(fDisplayFrame == 0)
    return;

  fDisplayFrame->DestroySubwindows();
  fDisplayFrame->Cleanup();
  printf("DisplayChildren item %s List %d btn=%d\n", item->GetText(),fDisplayFrame->GetList()->GetEntries(), btn);

  if(item->GetFirstChild() == 0) return;

  UInt_t wH = 2;
  UInt_t wW = 7;

  UInt_t fw, fh;
  Int_t nc = 0;  
  TGListTreeItem *child = item->GetFirstChild();
  do {
    child = child->GetNextSibling();
    nc ++;
  } while(child);
  fw = 70;
  fh = UInt_t(nc*2);
  fDisplayFrame->Resize(fw, fh);
  TGXYLayout* xyl = new TGXYLayout(fDisplayFrame);
  fDisplayFrame->SetLayoutManager(xyl);
  xyl->Layout();

  TGXYLayoutHints* lh;
  Float_t x,y;
  y  = 0.;
  nc = 0;
  child = item->GetFirstChild();
  do {
    // generic info 
    wW = 24;
    x = 0.;
    TGTextButton* b1 = new TGTextButton( fDisplayFrame, Form("%s",child->GetText()));
    b1->Resize(wW,wH);
    b1->SetTextJustify(kTextLeft | kTextCenterY);
    lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
    fDisplayFrame->AddFrame(b1,lh);
    x += wW;
    wW = 8;
    TGCheckButton* b2 = new TGCheckButton(fDisplayFrame, "Draw");
    b2->Resize(wW,wH);
    b2->SetTextJustify(kTextLeft | kTextCenterY);
    lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
    fDisplayFrame->AddFrame(b2, lh);
    x += wW;

    RenderElement* re = (RenderElement*)child->GetUserData();
    TObject* obj = re->GetObject();
    if(obj != 0) {
      TGXYLayoutHints* lh;

      Track*          track = dynamic_cast<Track*>(obj); // (Track*) obj->IsA()->DynamicCast(Track::Class(), obj );
      PointSet* hcont = dynamic_cast<PointSet*>(obj);
      TrackList* tcont = dynamic_cast<TrackList*>(obj);
      TGeoNode*          gnode = dynamic_cast<TGeoNode*>(obj);
      
      // Track
      //---------

      if(track) {       
	// printf("display children track \n");
        b2->SetOn(track->GetRnrElement());
        b2->Connect("Toggled(Bool_t)", "Reve::Track", track, "SetRnrElement(Bool_t)"); 
      }

      // PointSet
      //------------------

      if (hcont) { 
	// connect to toggle signal
        wW = 8;
        //printf("add label to %s %d\n", cont->GetName(), cont->GetNPoints());
	TGLabel* b3 = new TGLabel(fDisplayFrame, Form("%d", hcont->GetN()));
	b3->SetTextJustify(kTextLeft | kTextCenterY);
	b3->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(b3, lh);
	x += wW;

	wW = 5; 
	TGColorSelect* b4 = new TGColorSelect(fDisplayFrame, TColor::Number2Pixel(hcont->GetMainColor()));
	b4->Resize();
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(b4,lh);
	b4->Connect("ColorSelected(Pixel_t)",
		    "Reve::PointSet", hcont, "SetMainColor(Pixel_t)");  
	
        x += wW;	
	wW = 8;
	ReveValuator* ne = new ReveValuator(fDisplayFrame, hcont->GetMarkerStyle());
	ne->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	ne->Connect("ValueSet(Long_t)", "Reve::RGBrowser", this, "SetMarkerStyle(Long_t)");
	ne->SetUserData(hcont);
	fDisplayFrame->AddFrame(ne,lh);

	//connect to container
        b2->SetUserData(hcont);
        b2->SetOn(hcont->GetRnrElement());
        b2->Connect("Toggled(Bool_t)", "Reve::PointSet", hcont, "SetRnrElement(Bool_t)");  
      }

      // TrackList
      //------------------

      if (tcont) {
	wW = 8;
        //printf("add label to %s %d\n", cont->GetName(), cont->GetNPoints());
	TGLabel* b3 = new TGLabel(fDisplayFrame, Form("%d", tcont->GetNTracks()));
	b3->SetTextJustify(kTextLeft | kTextCenterY);
	b3->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(b3, lh);
	x += wW;
        // track color
	wW = 5; 
	TGColorSelect* b4 = new TGColorSelect(fDisplayFrame, TColor::Number2Pixel(tcont->GetMainColor()));
	b4->Resize();
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(b4,lh);
	b4->Connect("ColorSelected(Pixel_t)",
		    "Reve::TrackList", tcont, "SetMainColor(Pixel_t)");  
	x += wW;
	wW = 8;
	ReveValuator* ne1 = new ReveValuator(fDisplayFrame, tcont->GetRnrStyle()->fMaxR);
	ne1->SetUserData(tcont);
	ne1->Connect("ValueSet(Long_t)", "Reve::RGBrowser", this, "SetMaxR(Long_t)");
	// ne1->SetToolTipText("Maximum radius [cm]");
	ne1->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(ne1,lh);
 
	x += wW;	
	wW = 8;
	ReveValuator* ne2 = new ReveValuator(fDisplayFrame, tcont->GetRnrStyle()->fMaxZ);
	ne2->SetUserData(tcont);
	ne2->Connect("ValueSet(Long_t)", "Reve::RGBrowser", this, "SetMaxZ(Long_t)");
	// ne2->SetToolTipText("Maximum z [cm]");
	ne2->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(ne2,lh);

	x += wW;	
	wW = 8;
	ReveValuator* ne3 = new ReveValuator(fDisplayFrame, tcont->GetRnrStyle()->fMaxOrbs);
	ne3->SetUserData(tcont);
	ne3->Connect("ValueSet(Long_t)", "Reve::RGBrowser", this, "SetMaxOrbs(Long_t)");
	// ne3->SetToolTipText("Maximum number of orbits");
	ne3->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(ne3,lh);

	x += wW;	
	wW = 8;
	TGCheckButton*  dau  = new TGCheckButton(fDisplayFrame, "Daughters");
        dau->SetOn(tcont->GetRnrStyle()->fFitDaughters);
	dau->SetUserData(tcont);
	dau->Connect("Toggled(Bool_t)", "Reve::TrackList", tcont, "SetFitDaughters(Bool_t)");
	dau->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(dau,lh);

	x += wW;	
	wW = 8;
	TGCheckButton* dec = new TGCheckButton(fDisplayFrame, "Decay");
        dec->SetOn(tcont->GetRnrStyle()->fFitDaughters);
	dec->SetUserData(tcont);
	dec->Connect("Toggled(Bool_t)", "Reve::TrackList", tcont, "SetFitDecay(Bool_t)");
	dec->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(dec,lh);

	//connect to container
        b2->SetUserData(tcont);
        //b2->SetOn(tcont->GetRnrTracks());
        //b2->Connect("Toggled(Bool_t)", "Reve::TrackList", tcont, "SetRnrTracks(Bool_t)");  
        b2->SetOn(tcont->GetRnrElement());
        b2->Connect("Toggled(Bool_t)", "Reve::TrackList", tcont, "SetRnrElement(Bool_t)");  
      }

      // TGeoNode
      //---------

      if(gnode) {
        TGeoVolume* vol = gnode->GetVolume();
	b2->SetOn(gnode->IsVisible());
        b2->Connect("Toggled(Bool_t)", "Reve::RGBrowser", this,"NodeVis(Bool_t)");  
        b2->SetUserData(gnode);

	wW = 11;
	TGCheckButton* b3 = new TGCheckButton(fDisplayFrame, "VisibleDaughters");
	b3->SetTextJustify(kTextLeft | kTextCenterY);
	b3->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(b3,lh);
	b3->SetOn(vol->IsVisibleDaughters());
        b3->Connect("Toggled(Bool_t)", "Reve::RGBrowser", this, "VolumeDaughterVis(Bool_t)");  
	b3->SetUserData(vol);
	x += wW;

	wW = 5; 
	ReveColorSelect* b4 = new ReveColorSelect(fDisplayFrame, TColor::Number2Pixel(vol->GetLineColor()));
	b4->Resize();
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(b4,lh);
	b4->Connect("ColorSelected(Pixel_t)", "Reve::RGBrowser", this, "SetVolumeColor(UInt_t)");  
        b4->SetUserData(vol);
	x += wW;

	wW = 11;
	ReveValuator* ne = new ReveValuator(fDisplayFrame, vol->GetTransparency());
	ne->Connect("ValueSet(Long_t)", "Reve::RGBrowser", this, "SetTransparency(Long_t)");
	ne->SetUserData(vol);

	ne->Resize(wW,wH);
	lh = new TGXYLayoutHints(x, y, wW, wH,0);  lh->SetPadLeft(2); lh->SetPadRight(2);
	fDisplayFrame->AddFrame(ne,lh);
        x += wW;
      }

    }
    y += wH; 
    nc++;
    child = child->GetNextSibling();
  } while(child);
  fDisplayFrame->MapSubwindows();
  fDisplayFrame->MapWindow();
  MapSubwindows();
}

/**************************************************************************/
// Slots
/**************************************************************************/

void RGBrowser::SetTransparency(Long_t )
{
  ReveValuator& rv = *(ReveValuator*)gTQSender;
  // printf("VSet idx=%d, double value=%lf, part=%p\n", val, rv.GetNumber(), rv.GetUserData());
  TGeoVolume* vol = (TGeoVolume*) rv.GetUserData();
  if(vol) {
    //    printf("set volume user data %d \n",val);
    vol->SetTransparency(char(rv.GetNumber()));
  }

  TGFrameElement* fel;
  TList* list = fDisplayFrame->GetList();
  TIter nextin(list);
  ReveValuator* cw;
  while ((fel = (TGFrameElement*)nextin())){
    // printf("RGBrowser::SetTransparency %s  in fDisplayFrame\n", fel->fFrame->GetName());
    cw = dynamic_cast<ReveValuator*>(fel->fFrame);
    if(cw) {

      TGeoVolume* v = dynamic_cast<TGeoVolume*>((RenderElement*)cw->GetUserData());
      if(v) {
	cw->SetNumber(v->GetTransparency());
      }
    }
  }
  gReve->Redraw3D();
}

/**************************************************************************/

void RGBrowser::SetVolumeColor(UInt_t pixel)
{
  Int_t r, g, b;
  TColor::Pixel2RGB(pixel, r, g, b);

  TGColorSelect* w = (TGColorSelect*) gTQSender;
  TGeoVolume* vol = (TGeoVolume*) w->GetUserData();
  Int_t col = TColor::GetColor(pixel);
  vol->SetLineColor(col);

  ReveColorSelect* cw;
  TGFrameElement* fel;
  TList* list = fDisplayFrame->GetList();
  TIter nextin(list);
  while ((fel = (TGFrameElement*)nextin())){
    // printf("%s  in fDisplayFrame\n", fel->fFrame->GetName());
    cw = dynamic_cast<ReveColorSelect*>(fel->fFrame);
    if(cw) {
      TGeoVolume* cv = dynamic_cast<TGeoVolume*>((TObject*)cw->GetUserData());
      if(cv) {
	// printf("TGColorSelect  %d %d\n",pixel, cv->GetLineColor());
	cw->UpdateColor(TColor::Number2Pixel(cv->GetLineColor())); 
      }
    }
  }
  gClient->NeedRedraw(fDisplayFrame);
  gReve->Redraw3D();
}

void RGBrowser::NodeVis(Bool_t vis)
{
  TGCheckButton& rv = *(TGCheckButton*)gTQSender;
  TGeoNode* node = (TGeoNode*) rv.GetUserData();
  if(node) {
    Reve::PadHolder pHolder(false, gReve->GetCC());
    node->SetVisibility(vis);
    gReve->Redraw3D();
  }
}

void RGBrowser::VolumeDaughterVis(Bool_t vis)
{
  TGCheckButton& rv = *(TGCheckButton*)gTQSender;
  // printf("VSet idx=%d, double value=%lf, part=%p\n", val, rv.GetNumber(), rv.GetUserData());
  TGeoVolume* vol = (TGeoVolume*) rv.GetUserData();
  if(vol) {
    Reve::PadHolder pHolder(false, gReve->GetCC());
    vol->VisibleDaughters(vis);
    gReve->Redraw3D();
  }
}

/**************************************************************************/
/**************************************************************************/

void RGBrowser::SetMaxR(Long_t )
{
  ReveValuator*      rv = (ReveValuator*) gTQSender;
  TrackList* tc = (TrackList*) rv->GetUserData();
  if(tc) {
    tc->SetMaxR(rv->GetNumber());
  }
}

void RGBrowser::SetMaxZ(Long_t )
{
  ReveValuator*      rv = (ReveValuator*) gTQSender;
  TrackList* tc = (TrackList*) rv->GetUserData();
  if(tc) {
    tc->SetMaxZ(rv->GetNumber());
  }
}

void RGBrowser::SetMaxOrbs(Long_t )
{
  ReveValuator*      rv = (ReveValuator*) gTQSender;
  TrackList* tc = (TrackList*) rv->GetUserData();
  if(tc) {
    tc->SetMaxOrbs(rv->GetNumber());
  }
}

/**************************************************************************/
/**************************************************************************/

void RGBrowser::SetMarkerStyle(Long_t )
{ 
  ReveValuator*   rv = (ReveValuator*) gTQSender;
  PointSet* pc = (PointSet*) rv->GetUserData();
  if(pc) {
    Reve::PadHolder pHolder(false, gReve->GetCC());
    pc->SetMarkerStyle(short(rv->GetNumber()));
    gReve->Redraw3D();
  }
}

/**************************************************************************/
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
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

/**************************************************************************/
// ReveValuator
/**************************************************************************/

ReveValuator::~ReveValuator()
{}
