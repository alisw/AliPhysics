#include "RGTopFrame.h"

#include "RGBrowser.h"
#include "RGEditor.h"

#include <Reve/EventBase.h>
#include "VSDSelector.h"

#include <TGMenu.h>
#include <TGTab.h>
#include <TGToolBar.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGSplitter.h>
#include <TRootEmbeddedCanvas.h>
#include <TGMimeTypes.h>

#include <TGLSAViewer.h>
#include <TH1F.h>
#include <TView.h>

#include <TROOT.h>
#include <TFile.h>
#include <TMacro.h>
#include <TFolder.h>
#include <TStyle.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TRint.h>
#include <TVirtualX.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TEnv.h>
#include <TStyle.h>
#include <KeySymbols.h>
#include "TVirtualGL.h"
#include "TPluginManager.h"

#include <iostream>

using namespace Reve;
using namespace Reve;

Reve::RGTopFrame* gReve = 0;

/**************************************************************************/

RGTopFrame::RGTopFrame(const TGWindow *p, UInt_t w, UInt_t h, LookType_e look) :
  TGMainFrame(p, w, h),

  fMasterFrame (0),
  fMasterTab   (0),
  fGLCanvas    (0),
  fSelector    (0),
  fBrowser     (0),
  fStatusBar   (0),
  fVSDFile     (""),

  fMacroFolder(0),
  fEditor (0),

  fCurrentEvent   (0),
  fGlobalStore    (0),

  fRedrawDisabled (0),
  fResetCameras   (kFALSE),
  fTimerActive    (kFALSE),
  fRedrawTimer    (),

  fLook           (LT_Editor),
  fGeometries     ()
{
  gReve = this;
  fRedrawTimer.Connect("Timeout()", "Reve::RGTopFrame", this, "DoRedraw3D()");
  fMacroFolder = new TFolder("EVE", "Visualization macros");
  gROOT->GetListOfBrowsables()->Add(fMacroFolder);

  fClient->GetMimeTypeList()->AddType("root/tmacro", "Reve::RMacro",
                                      "tmacro_s.xpm", "tmacro_t.xpm", "");

  // Build GUI

  TGLayoutHints *lay0 = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandY | kLHintsExpandX);
  TGLayoutHints *lay1 = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandY | kLHintsExpandX, 2, 0, 2, 2);
  TGLayoutHints *lay2 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

  fMasterFrame = new TGCompositeFrame(this, w, h, kHorizontalFrame | kRaisedFrame);
  TGVerticalFrame* fV2 = new TGVerticalFrame(fMasterFrame, GetWidth()-40, GetHeight()-40, kSunkenFrame);

  fMasterFrame->AddFrame(fV2, lay1);

  fMasterTab = new TGTab(fV2, GetWidth(), GetHeight());  

  // browser tab
  TGCompositeFrame* tframe1 = fMasterTab->AddTab("Object Browser");
  fBrowser = new RGBrowser(tframe1, w, h);
  tframe1->AddFrame(fBrowser, lay2);

  // tree selection tab
  TGCompositeFrame* tframe2 = fMasterTab->AddTab("Tree Selections");  
  fSelector = new VSDSelector(tframe2);

  // gl-canvas
  Reve::PushPad();
  TGCompositeFrame* tframe3 = fMasterTab->AddTab("GLCanvas");
  TRootEmbeddedCanvas* ecanvas3 = new TRootEmbeddedCanvas("GLCanvas", tframe3, 580, 360);
  tframe3->AddFrame(ecanvas3, lay2);
  fGLCanvas = ecanvas3->GetCanvas();
  fGLCanvas->SetEditable(kFALSE);
  Reve::PopPad();

  fV2->AddFrame(fMasterTab, lay0);
  AddFrame(fMasterFrame, lay0);
   
  // Create status bar
  Int_t parts[] = {45, 45, 10};
  fStatusBar = new TGStatusBar(this, 50, 10, kHorizontalFrame);
  fStatusBar->SetParts(parts, 3);
  TGLayoutHints* lay6 = new TGLayoutHints(kLHintsBottom| kLHintsExpandX, 0, 0, 0, 0);
  AddFrame(fStatusBar, lay6);
  fStatusBar->SetText("GUI created", 0);

  MapSubwindows();

  /**************************************************************************/
  /**************************************************************************/
  
  switch(look) {
  case LT_Classic: {
    fBrowser->SetupClassicLook(fEditor, fGLCanvas);
    fGLCanvas->GetViewer3D("ogl");
    break;
  }

  case LT_Editor: {
    fBrowser->SetupEditorLook(fEditor, fGLCanvas);
    fGLCanvas->GetViewer3D("ogl");
    break;
  }

  case LT_GLViewer: {
    fBrowser->SetupGLViewerLook(fEditor, fGLCanvas);
    break;
  }

  default: {
    printf("RGTopFrame unknown look-type, ignoring.\n");
    break;
  }
  }

  TGLViewer* glv = dynamic_cast<TGLViewer*>(fGLCanvas->GetViewer3D());
  if(glv) {
    glv->SetSmartRefresh(kTRUE);
    glv->SetResetCamerasOnUpdate(kFALSE);
    glv->SetResetCameraOnDoubleClick(kFALSE);
    TGLSAViewer* glsav = dynamic_cast<TGLSAViewer*>(glv);
    if(glsav) {
      TGedEditor* e = glsav->GetGedEditor();
      e->SetModel(e->GetPad(), glsav, kButton1Down);
    }
  }

  /**************************************************************************/
  /**************************************************************************/

  fGlobalStore = new RenderElementList("Geometry", "");
  fGlobalStore->SetDenyDestroy(kTRUE);
  TGListTreeItem* glti = fGlobalStore->AddIntoListTree(GetListTree(), (TGListTreeItem*)0);
  GetListTree()->OpenItem(glti);
  DrawRenderElement(fGlobalStore);

  Resize(GetDefaultSize()); // this is used here to init layout algorithm
  SetWindowName("Reve");
  MapWindow();

  fEditor->DisplayObject(0);

  gSystem->ProcessEvents();
}

/**************************************************************************/

TCanvas* RGTopFrame::AddCanvasTab(const char* name)
{
  TGCompositeFrame    *f  = fMasterTab->AddTab(name);
  TRootEmbeddedCanvas *ec = new TRootEmbeddedCanvas(name, f, 580, 360);
  f->AddFrame(ec, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
				    2, 2, 2, 2));

  f->MapSubwindows();
  f->MapWindow();
  fMasterTab->GetTabTab(name)->MapWindow();
  fMasterTab->Layout();
  
  return ec->GetCanvas();
}

/**************************************************************************/

TGListTree* RGTopFrame::GetListTree() const
{
  return fBrowser->GetListTree();
}

/**************************************************************************/
// Macro management
/**************************************************************************/

TMacro* RGTopFrame::GetMacro(const Text_t* name) const
{
  return dynamic_cast<TMacro*>(fMacroFolder->FindObject(name));
}

/**************************************************************************/
// Editor
/**************************************************************************/

void RGTopFrame::EditRenderElement(RenderElement* rnr_element)
{
  static const Exc_t eH("RGTopFrame::EditRenderElement ");

  fEditor->DisplayRenderElement(rnr_element);
}

/**************************************************************************/
// 3D Pad management
/**************************************************************************/

void RGTopFrame::RegisterRedraw3D()
{
  fRedrawTimer.Start(0, kTRUE);
  fTimerActive = true;
}

void RGTopFrame::DoRedraw3D()
{
  // printf("RGTopFrame::DoRedraw3D redraw triggered\n");
  if (fResetCameras) {
    fGLCanvas->GetViewer3D()->ResetCamerasAfterNextUpdate();
    fResetCameras = kFALSE;
  }

  fGLCanvas->Modified();
  fGLCanvas->Update();
  fTimerActive = kFALSE;
}

/**************************************************************************/

int RGTopFrame::SpawnGuiAndRun(int argc, char **argv)
{
  LookType_e revemode = LT_Editor;
  Int_t w = 540;
  Int_t h = 500;
  if(argc >= 3 && (strcmp(argv[1], "-revemode")==0 || strcmp(argv[1], "-mode")==0)) {
    LookType_e m = LookType_e(atoi(argv[2]));
    if(m >= LT_Classic && m <= LT_GLViewer)
      revemode = m;
    printf("revemode = %d\n", revemode);
    if(revemode == LT_GLViewer) {
      w = 1024; h = 768;
    }
  }

  TRint theApp("App", &argc, argv);

  /* gReve = */ new RGTopFrame(gClient ? gClient->GetRoot() : 0, w, h, revemode);

 run_loop:
  try {
    theApp.Run();
  }
  catch(Exc_t& exc) {
    gReve->GetStatusBar()->SetText(exc.Data());
    fprintf(stderr, "Exception: %s\n", exc.Data());
    goto run_loop;
  }
  return 0;
}

void RGTopFrame::SpawnGui(LookType_e revemode)
{
  Int_t w = 540;
  Int_t h = 500;
  if(revemode == LT_GLViewer) {
    w = 1024; h = 768;
  }

  /* gReve = */ new RGTopFrame(gClient->GetRoot(), w, h, revemode);
}

/**************************************************************************/
/**************************************************************************/

TGListTreeItem* RGTopFrame::AddEvent(EventBase* event)
{
  fCurrentEvent = event;
  fCurrentEvent->SetDenyDestroy(kTRUE);
  TGListTreeItem* elti = event->AddIntoListTree(GetListTree(), (TGListTreeItem*)0);
  GetListTree()->OpenItem(elti);
  DrawRenderElement(event);
  return elti;
}

TGListTreeItem* RGTopFrame::AddRenderElement(RenderElement* rnr_element)
{
  if (fCurrentEvent == 0)
    AddEvent(new EventBase("Event", "Auto-created event directory"));
  return AddRenderElement(fCurrentEvent, rnr_element);
}

TGListTreeItem* RGTopFrame::AddRenderElement(RenderElement* parent,
					     RenderElement* rnr_element)
{
  static const Exc_t eH("RGTopFrame::AddRenderElement ");

  // Here could route rnr-element to several browsers/pads.

  RenderElementListBase* rel = dynamic_cast<RenderElementListBase*>(parent);
  if(rel)
    rel->AddElement(rnr_element);

  TGListTreeItem* newitem =
    rnr_element->AddIntoListTree(GetListTree(), parent);

  return newitem;
}

TGListTreeItem* RGTopFrame::AddGlobalRenderElement(RenderElement* rnr_element)
{
  return AddGlobalRenderElement(fGlobalStore, rnr_element);
}

TGListTreeItem* RGTopFrame::AddGlobalRenderElement(RenderElement* parent,
						   RenderElement* rnr_element)
{
  static const Exc_t eH("RGTopFrame::AddGlobalRenderElement ");

  // Here could route rnr-element to several browsers/pads.

  RenderElementListBase* rel = dynamic_cast<RenderElementListBase*>(parent);
  if(rel)
    rel->AddElement(rnr_element);

  TGListTreeItem* newitem =
    rnr_element->AddIntoListTree(GetListTree(), parent);

  return newitem;
}

/**************************************************************************/

void RGTopFrame::RemoveRenderElement(RenderElement* parent,
				     RenderElement* rnr_element)
{
  rnr_element->RemoveFromListTree(GetListTree());

  RenderElementListBase* rel = dynamic_cast<RenderElementListBase*>(parent);
  if(rel)
    rel->RemoveElement(rnr_element);
}

void RGTopFrame::PreDeleteRenderElement(RenderElement* rnr_element)
{
  if (fEditor->GetRnrElement() == rnr_element)
    fEditor->DisplayObject(0);
}

/**************************************************************************/

void RGTopFrame::DrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad)
{
  if(pad == 0) pad = fGLCanvas;

  { Reve::PadHolder pHolder(false, pad);
    if (pad == fGLCanvas) fGLCanvas->SetEditable(kTRUE);
    rnr_element->GetObject()->Draw();
    if (pad == fGLCanvas) fGLCanvas->SetEditable(kFALSE);
  }
  Redraw3D();
}

void RGTopFrame::UndrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad)
{
  if(pad == 0) pad = fGLCanvas;
  { Reve::PadHolder pHolder(false, pad);
    pad->GetListOfPrimitives()->Remove(rnr_element->GetObject());
  }
  Redraw3D();
}

/**************************************************************************/

void RGTopFrame::RenderElementChecked(TObject* obj, Bool_t state)
{
  // Item's user-data is blindly casted into TObject.
  // We recast it blindly back into the render element.

  RenderElement* rnrEl = (RenderElement*) obj;
  rnrEl->SetRnrElement(state);
  Redraw3D();
}

/**************************************************************************/

void RGTopFrame::NotifyBrowser(TGListTreeItem* parent_lti)
{
  TGListTree* l_tree = GetListTree();
  if(parent_lti)
    l_tree->OpenItem(parent_lti);
  gClient->NeedRedraw(l_tree);
}

void RGTopFrame::NotifyBrowser(RenderElement* parent)
{
  TGListTreeItem* parent_lti = parent ? parent->FindListTreeItem(GetListTree()) : 0;
  NotifyBrowser(parent_lti);
}

/**************************************************************************/
// GeoManager registration
/**************************************************************************/

TGeoManager* RGTopFrame::GetGeometry(const TString& filename)
{
  static const Exc_t eH("RGTopFrame::GetGeometry ");

  TString exp_filename = filename;
  gSystem->ExpandPathName(exp_filename);
  printf("%s loading: '%s' -> '%s'.\n", eH.Data(),
	 filename.Data(), exp_filename.Data());

  std::map<TString, TGeoManager*>::iterator g = fGeometries.find(filename);
  if(g != fGeometries.end()) {
    return g->second;
  } else {
    if(gSystem->AccessPathName(exp_filename, kReadPermission))
      throw(eH + "file '" + exp_filename + "' not readable.");
    gGeoManager = 0;
    TGeoManager::Import(filename); 
    if(gGeoManager == 0)
      throw(eH + "GeoManager import failed.");
    gGeoManager->GetTopVolume()->VisibleDaughters(1);

    // Import colors exported by Gled, if they exist.
    {
      TFile f(exp_filename, "READ");
      TObjArray* collist = (TObjArray*) f.Get("ColorList");
      f.Close();
      if(collist != 0) {
	TSeqCollection* glist = gROOT->GetListOfColors();
	glist->Clear();
	glist->AddAll(collist);
      }
    }

    fGeometries[filename] = gGeoManager;
    return gGeoManager;
  }
}

/**************************************************************************/
// Testing exceptions
/**************************************************************************/

void RGTopFrame::ThrowException(const char* text)
{
  static const Exc_t eH("RGTopFrame::ThrowException ");

  throw(eH + text);
}
