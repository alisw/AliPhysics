#include "RGTopFrame.h"

#include "RGBrowser.h"
#include "RGEditor.h"
#include "VSDSelector.h"

#include <Reve/RenderElement.h>

#include <TGMenu.h>
#include <TGTab.h>
#include <TGToolBar.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGSplitter.h>
#include <TRootEmbeddedCanvas.h>

#include <TGLSAViewer.h>
#include <TH1F.h>
#include <TView.h>

#include <TROOT.h>
#include <TFile.h>
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

namespace {

enum RGBrowserMT {
  M_LAYOUT_1,
  M_LAYOUT_2,
  M_LAYOUT_3
};

const char *xpm_names[] = {
    "lay1.xpm",
    "lay2.xpm",
    "lay3.xpm",
    0
};

ToolBarData_t tb_data[] = {
  { "", "Standard list layout",     kFALSE, M_LAYOUT_1,        NULL },
  { "", "TParticle latout",         kFALSE, M_LAYOUT_2,        NULL },
  { "", "TGeo layout",              kFALSE, M_LAYOUT_3,        NULL },
  { NULL,            NULL,          0,      0,                 NULL }
};

} // unnamed namespace

/**************************************************************************/


void RGTopFrame::Init(){
  fCC          = 0;
  fHistoCanvas = 0;
  fSelector    = 0;
  fBrowser     = 0;
  fStatusBar   = 0;
  fVSDFile     = "";

  fEditor = 0;

  fCurrentEvent    = 0;
  fCurrentEventLTI = 0;
  fGeometryLTI     = 0;

  fRedrawDisabled = false;
  fTimerActive    = false;
  fRedrawTimer.Connect("Timeout()", "Reve::RGTopFrame", this, "DoRedraw3D()");
}


RGTopFrame::RGTopFrame(const TGWindow *p, UInt_t w, UInt_t h, LookType_e look)
  : TGMainFrame(p, w, h)
{
  Init();
  TGLayoutHints *fL0 = new TGLayoutHints(kLHintsCenterX |kLHintsCenterY | kLHintsExpandY|  kLHintsExpandX);
  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsCenterX |kLHintsCenterY | kLHintsExpandY|  kLHintsExpandX,2,0,2,2);
  TGLayoutHints* fL2 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
					 2, 2, 2, 2);
  TGCompositeFrame*  fMainFrame = new TGCompositeFrame(this, w, h,kHorizontalFrame | kRaisedFrame);
  fMainFrame->SetCleanup(kDeepCleanup);
  TGVerticalFrame* fV2 = new TGVerticalFrame(fMainFrame, GetWidth()-40, GetHeight()-40, kSunkenFrame);

  fMainFrame->AddFrame(fV2, fL1);

  // ??? TGCanvas* fCanvasWindow = new TGCanvas(fV2,w,h);
  TGTab*    fDisplayFrame = new TGTab(fV2, GetWidth(), GetHeight());  

  // browser tab
  TGCompositeFrame* tFrame1 = fDisplayFrame->AddTab("Object Browser");
  fBrowser = new RGBrowser(tFrame1, w, h);
  tFrame1->AddFrame(fBrowser, fL2);

  // tree selection tab
  TGCompositeFrame* tFrame2 = fDisplayFrame->AddTab("Tree Selections");  
  fSelector = new VSDSelector(fBrowser->GetListTree(), tFrame2);

  // canvas
  Reve::PushPad();
  TGCompositeFrame* tFrame3 = fDisplayFrame->AddTab("Canvas");
  TRootEmbeddedCanvas* fEmbeddedCanvas3 = new TRootEmbeddedCanvas("fEmbeddedCanvas3", tFrame3, 580, 360);
  tFrame3->AddFrame(fEmbeddedCanvas3, fL2);
  fEmbeddedCanvas3->GetCanvas()->SetBorderMode(0);
  fCC = fEmbeddedCanvas3->GetCanvas();
  // fCC->SetFillColor(1);
  Reve::PopPad();

  // histo canvas
  TGCompositeFrame* frame4 = fDisplayFrame->AddTab("HistoCanvas");
  TRootEmbeddedCanvas* ecanvas4 = new TRootEmbeddedCanvas("HistoCanvas", frame4, 580, 360);
  frame4->AddFrame(ecanvas4, fL2);
  fHistoCanvas =  ecanvas4->GetCanvas();
  fHistoCanvas->SetBorderMode(0);

  fV2->AddFrame(fDisplayFrame, fL0);
  AddFrame(fMainFrame, fL0);
   
  // Create status bar
  Int_t parts[] = {45, 45, 10};
  fStatusBar = new TGStatusBar(this, 50, 10, kHorizontalFrame);
  fStatusBar->SetParts(parts, 3);
  TGLayoutHints* fL6 = new TGLayoutHints(kLHintsBottom| kLHintsExpandX, 0, 0, 0, 0);
  AddFrame(fStatusBar, fL6);
  fStatusBar->SetText("GUI created", 0);

  MapSubwindows();

  /**************************************************************************/
  /**************************************************************************/
  
  fEditor = new RGEditor(fCC);
  fEditor->GetCan()->ChangeOptions(0);

  switch(look) {
  case LT_Classic: {
    fBrowser->SetupClassicLook();
    // Need push/pop pad around here? Not.
    fCC->GetViewer3D("ogl");
    break;
  }

  case LT_Editor: {
    fBrowser->SetupEditorLook(fEditor);
    fCC->GetViewer3D("ogl");
    break;
  }

  case LT_GLViewer: {
    printf("LT_GLViewer this option currently somewhat broken!\n");
    fBrowser->SetupGLViewerLook(fEditor, fCC);
    printf("Crap1 %d %d\n", GetWidth(), GetHeight());
    break;
  }

  default: {
    printf("RGTopFrame unknown look-type, ignoring.\n");
    break;
  }
  }

  TGLViewer* glv = dynamic_cast<TGLViewer*>(fCC->GetViewer3D());
  if(glv) {
    glv->SetSmartRefresh(true);
  }

  /**************************************************************************/
  /**************************************************************************/

  fGeometryLTI = GetListTree()->AddItem(0, "Geometry");
  GetListTree()->OpenItem(fGeometryLTI);

  Resize(GetDefaultSize()); // this is used here to init layout algorithm
  SetWindowName("Reve");
  MapWindow();

  fEditor->DisplayObject(0);

  gSystem->ProcessEvents();
}

/**************************************************************************/

TGListTree* RGTopFrame::GetListTree()
{
  return fBrowser->GetListTree();
}

TGListTreeItem* RGTopFrame::GetEventTreeItem()
{
  // return fBrowser->GetListTree()->FindItemByPathname("Event");
  return fCurrentEventLTI;
}

TGListTreeItem* RGTopFrame::GetGlobalTreeItem()
{
  return fGeometryLTI;
}

/**************************************************************************/
// Editor
/**************************************************************************/

void RGTopFrame::EditRenderElement(RenderElement* rnr_element)
{
  static const Exc_t eH("RGTopFrame::EditRenderElement ");

  TObject* tobj = 0;
  if(rnr_element) tobj = rnr_element->GetObject();
  fEditor->DisplayObject(tobj);
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
  fCC->Modified();
  fCC->Update();
  fTimerActive = false;
}

/**************************************************************************/

int RGTopFrame::SpawnGuiAndRun(int argc, char **argv)
{
  LookType_e revemode = LT_Editor;
  if(argc >= 3 && strcmp(argv[1], "-revemode")==0) {
    LookType_e m = LookType_e(atoi(argv[2]));
    if(m >= LT_Classic && m <= LT_GLViewer)
      revemode = m;
    printf("revemode = %d\n", revemode);
  }

  TRint theApp("App", &argc, argv);
  Int_t w = 800;
  Int_t h = 600;

  gReve = new RGTopFrame(gClient->GetRoot(), w, h, revemode);
 run_loop:
  try {
    theApp.Run();
  }
  catch(std::string exc) {
    gReve->GetStatusBar()->SetText(exc.c_str());
    goto run_loop;
  }
  return 0;
}

/**************************************************************************/
/**************************************************************************/

TGListTreeItem* RGTopFrame::AddEvent(TObject* event)
{
  fCurrentEvent = event;
  RenderElementObjPtr* rnrEv = new RenderElementObjPtr(event);
  fCurrentEventLTI = rnrEv->AddIntoListTree(GetListTree(), 0);
  GetListTree()->OpenItem(fCurrentEventLTI);
  return fCurrentEventLTI;
}

TGListTreeItem* RGTopFrame::AddRenderElement(RenderElement* rnr_element)
{
  return AddRenderElement(GetEventTreeItem(), rnr_element);
}

TGListTreeItem* RGTopFrame::AddRenderElement(TGListTreeItem* parent,
					     RenderElement* rnr_element)
{
  static const Exc_t eH("RGTopFrame::AddRenderElement ");

  // Here could route rnr-element to several browsers/pads.

  TGListTreeItem* newitem =
    rnr_element->AddIntoListTree(GetListTree(), parent);
  NotifyBrowser();

  return newitem;
}

TGListTreeItem* RGTopFrame::AddGlobalRenderElement(RenderElement* rnr_element)
{
  return AddGlobalRenderElement(GetGlobalTreeItem(), rnr_element);
}

TGListTreeItem* RGTopFrame::AddGlobalRenderElement(TGListTreeItem* parent,
						   RenderElement* rnr_element)
{
  static const Exc_t eH("RGTopFrame::AddGlobalRenderElement ");

  // Here could route rnr-element to several browsers/pads.

  TGListTreeItem* newitem =
    rnr_element->AddIntoListTree(GetListTree(), parent);
  NotifyBrowser();

  return newitem;
}

/**************************************************************************/

void RGTopFrame::DrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad)
{
  if(pad == 0) pad = GetCC();
  { Reve::PadHolder pHolder(false, pad);
    rnr_element->GetObject()->Draw();
  }
  Redraw3D();
}

void RGTopFrame::UndrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad)
{
  if(pad == 0) pad = GetCC();
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
}

/**************************************************************************/

void RGTopFrame::NotifyBrowser(TGListTreeItem* parent)
{
  TGListTree* l_tree = GetListTree();
  if(parent)
    l_tree->OpenItem(parent);
  gClient->NeedRedraw(l_tree);
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
