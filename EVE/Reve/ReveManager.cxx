#include "ReveManager.h"

#include <Reve/Viewer.h>
#include <Reve/Scene.h>
#include <Reve/Pad.h>
#include <Reve/EventBase.h>

#include <Reve/RGBrowser.h>
#include <Reve/RGEditor.h>

#include <TGMenu.h>
#include <TGTab.h>
#include <TGToolBar.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGSplitter.h>
#include <TRootEmbeddedCanvas.h>

#include <TGStatusBar.h>

#include <TGLSAViewer.h>

#include <TROOT.h>
#include <TFile.h>
#include <TMacro.h>
#include <TFolder.h>
#include <TStyle.h>
#include <TBrowser.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TRint.h>
#include <TVirtualX.h>
#include <TEnv.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGeoShape.h>
#include <KeySymbols.h>
#include "TVirtualGL.h"
#include "TPluginManager.h"

#include <iostream>

using namespace Reve;
using namespace Reve;

Reve::ReveManager* gReve = 0;

/**************************************************************************/

ReveManager::ReveManager(UInt_t w, UInt_t h) :
  fBrowser     (0),
  fEditor      (0),
  fStatusBar   (0),

  fMacroFolder (0),

  fViewers        (0),
  fScenes         (0),
  fViewer         (0),
  fGlobalScene    (0),
  fEventScene     (0),
  fCurrentEvent   (0),

  fRedrawDisabled (0),
  fResetCameras   (kFALSE),
  fDropLogicals   (kFALSE),
  fKeepEmptyCont  (kFALSE),
  fTimerActive    (kFALSE),
  fRedrawTimer    (),

  fGeometries     ()
{
  // Constructor.

  static const Exc_t eH("ReveManager::ReveManager ");

  if (gReve != 0)
    throw(eH + "There can be only one!");

  gReve = this;

  fRedrawTimer.Connect("Timeout()", "Reve::ReveManager", this, "DoRedraw3D()");
  fMacroFolder = new TFolder("EVE", "Visualization macros");
  gROOT->GetListOfBrowsables()->Add(fMacroFolder);


  // Build GUI
  fBrowser   = new RGBrowser(w, h);
  fStatusBar = fBrowser->GetStatusBar();

  // ListTreeEditor
  fBrowser->StartEmbedding(0);
  fLTEFrame = new RGLTEFrame("REVE");
  fBrowser->StopEmbedding();
  fBrowser->SetTabTitle("Reve", 0);
  fEditor   = fLTEFrame->fEditor;

  // GL viewer
  fBrowser->StartEmbedding(1);
  TGLSAViewer* glv = new TGLSAViewer(gClient->GetRoot(), 0, fEditor);
  //glv->GetFrame()->SetCleanup(kNoCleanup);
  glv->ToggleEditObject();
  fBrowser->StopEmbedding();
  fBrowser->SetTabTitle("GLViewer", 1);

  // Finalize it
  fBrowser->InitPlugins();
  fBrowser->MapWindow();

  // --------------------------------

  fViewers = new ViewerList("Viewers");
  fViewers->IncDenyDestroy();
  AddToListTree(fViewers, kTRUE);

  fViewer  = new Viewer("GL-One");
  fViewer->SetGLViewer(glv);
  fViewer->IncDenyDestroy();
  AddRenderElement(fViewer, fViewers);

  fScenes  = new SceneList ("Scenes");
  fScenes->IncDenyDestroy();
  AddToListTree(fScenes, kTRUE);

  fGlobalScene = new Scene("Geometry scene");
  fGlobalScene->IncDenyDestroy();
  AddRenderElement(fGlobalScene, fScenes);

  fEventScene = new Scene("Event scene");
  fEventScene->IncDenyDestroy();
  AddRenderElement(fEventScene, fScenes);

  fViewer->AddScene(fGlobalScene);
  fViewer->AddScene(fEventScene);

  /**************************************************************************/
  /**************************************************************************/

  fEditor->DisplayObject(GetGLViewer());

  gSystem->ProcessEvents();
}

ReveManager::~ReveManager()
{
  // Destructor.
}

/**************************************************************************/

TCanvas* ReveManager::AddCanvasTab(const char* name)
{
  // Add a new canvas tab.

  fBrowser->StartEmbedding(1, -1);
  TCanvas* c = new TCanvas;
  fBrowser->StopEmbedding();
  fBrowser->SetTabTitle(name, 1, -1);

  return c;
}

TGWindow* ReveManager::GetMainWindow() const
{
  // Get the main window, i.e. the first created reve-browser.

  return fBrowser;
}

TGLViewer* ReveManager::GetGLViewer() const
{
  // Get default TGLViewer.

  return fViewer->GetGLViewer();
}

Viewer* ReveManager::SpawnNewViewer(const Text_t* name, const Text_t* title,
				    Bool_t embed)
{
  // Create a new GL viewer.

  Viewer* v = new Viewer(name, title);

  if (embed)  fBrowser->StartEmbedding(1);
  v->SpawnGLViewer(gClient->GetRoot(), embed ? fEditor : 0);
  v->IncDenyDestroy();
  if (embed)  fBrowser->StopEmbedding(), fBrowser->SetTabTitle(name, 1);
  AddRenderElement(v, fViewers);
  return v;
}

Scene* ReveManager::SpawnNewScene(const Text_t* name, const Text_t* title)
{
  // Create a new scene.

  Scene* s = new Scene(name, title);
  AddRenderElement(s, fScenes);
  return s;
}

/**************************************************************************/
// Macro management
/**************************************************************************/

TMacro* ReveManager::GetMacro(const Text_t* name) const
{
  return dynamic_cast<TMacro*>(fMacroFolder->FindObject(name));
}

/**************************************************************************/
// Editor
/**************************************************************************/

void ReveManager::EditRenderElement(RenderElement* rnr_element)
{
  static const Exc_t eH("ReveManager::EditRenderElement ");

  fEditor->DisplayRenderElement(rnr_element);
}

/**************************************************************************/
// 3D Pad management
/**************************************************************************/

void ReveManager::RegisterRedraw3D()
{
  // Register a request for 3D redraw.

  fRedrawTimer.Start(0, kTRUE);
  fTimerActive = true;
}

void ReveManager::DoRedraw3D()
{
  // Perform 3D redraw of scenes and viewers whose contents has
  // changed.

  // printf("ReveManager::DoRedraw3D redraw triggered\n");

  fScenes ->RepaintChangedScenes();
  fViewers->RepaintChangedViewers(fResetCameras, fDropLogicals);

  fResetCameras = kFALSE;
  fDropLogicals = kFALSE;

  fTimerActive = kFALSE;
}

void ReveManager::FullRedraw3D(Bool_t resetCameras, Bool_t dropLogicals)
{
  // Perform 3D redraw of all scenes and viewers.

  fScenes ->RepaintAllScenes();
  fViewers->RepaintAllViewers(resetCameras, dropLogicals);
}

/**************************************************************************/

void ReveManager::RenderElementChanged(RenderElement* rnr_element)
{
  std::list<RenderElement*> scenes;
  rnr_element->CollectSceneParents(scenes);
  ScenesChanged(scenes);
}

void ReveManager::ScenesChanged(std::list<RenderElement*>& scenes)
{
  for (RenderElement::List_i s=scenes.begin(); s!=scenes.end(); ++s)
    ((Scene*)*s)->Changed();
}

/**************************************************************************/

int ReveManager::SpawnGuiAndRun(int argc, char **argv)
{
  Int_t w = 1024;
  Int_t h =  768;

  TRint theApp("App", &argc, argv);

  Reve::SetupGUI();
  /* gReve = */ new ReveManager(w, h);

 run_loop:
  try {
    theApp.Run();
  }
  catch(Exc_t& exc) {
    gReve->SetStatusLine(exc.Data());
    fprintf(stderr, "Exception: %s\n", exc.Data());
    goto run_loop;
  }
  return 0;
}

void ReveManager::SpawnGui()
{
  Int_t w = 1024;
  Int_t h =  768;

  Reve::SetupGUI();
  /* gReve = */ new ReveManager(w, h);
}

/**************************************************************************/
/**************************************************************************/

TGListTree* ReveManager::GetListTree() const
{
  return fLTEFrame->fListTree;
}

TGListTreeItem*
ReveManager::AddToListTree(RenderElement* re, Bool_t open, TGListTree* lt)
{
  // Add rnr-el as a top-level to a list-tree.
  // Please add a single copy of a render-element as a top level
  // or we will have to check for that, too.

  if (lt == 0) lt = GetListTree();
  TGListTreeItem* lti = re->AddIntoListTree(lt, (TGListTreeItem*)0);
  if (open) lt->OpenItem(lti);
  return lti;
}

void ReveManager::RemoveFromListTree(RenderElement* re, TGListTree* lt, TGListTreeItem* lti)
{
  // Remove top-level rnr-el from list-tree with specified tree-item.

  static const Exc_t eH("ReveManager::RemoveFromListTree ");

  if (lti->GetParent())
    throw(eH + "not a top-level item.");

  re->RemoveFromListTree(lt, 0);
}

/**************************************************************************/

TGListTreeItem* ReveManager::AddEvent(EventBase* event)
{
  fCurrentEvent = event;
  fCurrentEvent->IncDenyDestroy();
  AddRenderElement(fCurrentEvent, fEventScene);
  return AddToListTree(event, kTRUE);
}

TGListTreeItem* ReveManager::AddRenderElement(RenderElement* rnr_element,
					     RenderElement* parent)
{
  if (parent == 0) {
    if (fCurrentEvent == 0)
      AddEvent(new EventBase("Event", "Auto-created event directory"));
    parent = fCurrentEvent;
  }

  return parent->AddElement(rnr_element);
}

TGListTreeItem* ReveManager::AddGlobalRenderElement(RenderElement* rnr_element,
                                                    RenderElement* parent)
{
  if (parent == 0)
    parent = fGlobalScene;

  return parent->AddElement(rnr_element);
}

/**************************************************************************/

void ReveManager::RemoveRenderElement(RenderElement* rnr_element,
				     RenderElement* parent)
{
  parent->RemoveElement(rnr_element);
}

void ReveManager::PreDeleteRenderElement(RenderElement* rnr_element)
{
  if (fEditor->GetRnrElement() == rnr_element)
    fEditor->DisplayObject(0);
}

/**************************************************************************/

void ReveManager::RenderElementSelect(RenderElement* rnr_element)
{
  EditRenderElement(rnr_element);
}

Bool_t ReveManager::RenderElementPaste(RenderElement* rnr_element)
{
  RenderElement* src = fEditor->GetRnrElement();
  if (src)
    return rnr_element->HandleElementPaste(src);
  return kFALSE;
}

void ReveManager::RenderElementChecked(RenderElement* rnrEl, Bool_t state)
{
  rnrEl->SetRnrState(state);

  if (fEditor->GetModel() == rnrEl->GetEditorObject())
    fEditor->DisplayRenderElement(rnrEl);

  rnrEl->ElementChanged();
}

/**************************************************************************/

void ReveManager::NotifyBrowser(TGListTreeItem* parent_lti)
{
  TGListTree* lt = GetListTree();
  if (parent_lti)
    lt->OpenItem(parent_lti);
  lt->ClearViewPort();
}

void ReveManager::NotifyBrowser(RenderElement* parent)
{
  TGListTreeItem* parent_lti = parent ? parent->FindListTreeItem(GetListTree()) : 0;
  NotifyBrowser(parent_lti);
}

/**************************************************************************/
// GeoManager registration
/**************************************************************************/

TGeoManager* ReveManager::GetGeometry(const TString& filename)
{
  static const Exc_t eH("ReveManager::GetGeometry ");

  TString exp_filename = filename;
  gSystem->ExpandPathName(exp_filename);
  printf("%s loading: '%s' -> '%s'.\n", eH.Data(),
	 filename.Data(), exp_filename.Data());

  std::map<TString, TGeoManager*>::iterator g = fGeometries.find(filename);
  if (g != fGeometries.end()) {
    return g->second;
  } else {
    if (gSystem->AccessPathName(exp_filename, kReadPermission))
      throw(eH + "file '" + exp_filename + "' not readable.");
    gGeoManager = 0;
    TGeoManager::Import(filename); 
    if (gGeoManager == 0)
      throw(eH + "GeoManager import failed.");
    gGeoManager->GetTopVolume()->VisibleDaughters(1);

    // Import colors exported by Gled, if they exist.
    {
      TFile f(exp_filename, "READ");
      TObjArray* collist = (TObjArray*) f.Get("ColorList");
      f.Close();
      if (collist != 0) {
	TIter next(gGeoManager->GetListOfVolumes());
	TGeoVolume* vol;
	while ((vol = (TGeoVolume*) next()) != 0)
	{
          Int_t oldID = vol->GetLineColor();
          TColor* col = (TColor*)collist->At(oldID);
          Float_t r, g, b;
	  col->GetRGB(r, g, b);
          Int_t  newID = TColor::GetColor(r,g,b);
          vol->SetLineColor(newID);
	} 
      }
    }

    fGeometries[filename] = gGeoManager;
    return gGeoManager;
  }
}

/**************************************************************************/
// Testing exceptions
/**************************************************************************/

void ReveManager::SetStatusLine(const char* text)
{
  fStatusBar->SetText(text);
}

void ReveManager::ThrowException(const char* text)
{
  static const Exc_t eH("ReveManager::ThrowException ");

  throw(eH + text);
}
