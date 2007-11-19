// $Header$

#ifndef REVE_ReveManager_H
#define REVE_ReveManager_H

#include <TClass.h>
#include <TGeoManager.h>
#include <TROOT.h>
#include <TTimer.h>
#include <TVirtualPad.h>

#include <list>
#include <map>

class TMacro;
class TFolder;
class TCanvas;

class TGTab;
class TGStatusBar;
class TGListTree;
class TGListTreeItem;
class TGStatusBar;
class TGWindow;

class TGLViewer;

namespace Reve {

class RGLTEFrame;
class RGBrowser;
class RGEditor;

class RenderElement;
class PadPrimitive;

class Viewer; class ViewerList;
class Scene;  class SceneList;

class EventBase;


class ReveManager
{
  ReveManager(const ReveManager&);            // Not implemented
  ReveManager& operator=(const ReveManager&); // Not implemented

public:
  class RedrawDisabler
  {
  private:
    RedrawDisabler(const RedrawDisabler&);            // Not implemented
    RedrawDisabler& operator=(const RedrawDisabler&); // Not implemented

    ReveManager* fFrame;
  public:
    RedrawDisabler(ReveManager* f) : fFrame(f)
    { if (fFrame) fFrame->DisableRedraw(); }
    ~RedrawDisabler()
    { if (fFrame) fFrame->EnableRedraw(); }
  };

private:

  RGBrowser           *fBrowser;
  RGLTEFrame          *fLTEFrame;
  RGEditor            *fEditor;
  TGStatusBar         *fStatusBar;

  TFolder             *fMacroFolder;

  ViewerList          *fViewers;
  SceneList           *fScenes;

  Viewer              *fViewer;   // First / default gl-viewer.
  Scene               *fGlobalScene;
  Scene               *fEventScene;
  EventBase           *fCurrentEvent;

  Int_t                fRedrawDisabled;
  Bool_t               fFullRedraw;
  Bool_t               fResetCameras;
  Bool_t               fDropLogicals;
  Bool_t               fKeepEmptyCont;
  Bool_t               fTimerActive;
  TTimer               fRedrawTimer;

protected:
  std::map<TString, TGeoManager*> fGeometries;

public:
  ReveManager(UInt_t w, UInt_t h);
  virtual ~ReveManager();

  RGBrowser*   GetBrowser()   const { return fBrowser;   }
  RGLTEFrame*  GetLTEFrame()  const { return fLTEFrame;  }
  RGEditor*    GetEditor()    const { return fEditor;    }
  TGStatusBar* GetStatusBar() const { return fStatusBar; }

  SceneList*   GetScenes()   const { return fScenes;  }
  ViewerList*  GetViewers()  const { return fViewers; }

  Viewer*      GetDefViewer()    const { return fViewer; }
  Scene*       GetGlobalScene()  const { return fGlobalScene; }
  Scene*       GetEventScene()   const { return fEventScene; }
  EventBase*   GetCurrentEvent() const { return fCurrentEvent; }

  TCanvas*     AddCanvasTab(const char* name);
  TGWindow*    GetMainWindow() const;
  TGLViewer*   GetGLViewer() const;
  Viewer*      SpawnNewViewer(const Text_t* name, const Text_t* title="", Bool_t embed=kTRUE);
  Scene*       SpawnNewScene(const Text_t* name, const Text_t* title="");

  TFolder*  GetMacroFolder() const { return fMacroFolder; }
  TMacro*   GetMacro(const Text_t* name) const;

  void EditRenderElement(RenderElement* rnr_element);

  void DisableRedraw() { ++fRedrawDisabled; }
  void EnableRedraw()  { --fRedrawDisabled; if(fRedrawDisabled <= 0) Redraw3D(); }

  void Redraw3D(Bool_t resetCameras=kFALSE, Bool_t dropLogicals=kFALSE)
  {
    if(fRedrawDisabled <= 0 && !fTimerActive) RegisterRedraw3D();
    if(resetCameras) fResetCameras = kTRUE;
    if(dropLogicals) fDropLogicals = kTRUE;
  }
  void RegisterRedraw3D();
  void DoRedraw3D();
  void FullRedraw3D(Bool_t resetCameras=kFALSE, Bool_t dropLogicals=kFALSE);

  Bool_t GetKeepEmptyCont() const   { return fKeepEmptyCont; }
  void   SetKeepEmptyCont(Bool_t k) { fKeepEmptyCont = k; }


  void RenderElementChanged(RenderElement* rnr_element);
  void ScenesChanged(std::list<Reve::RenderElement*>& scenes);

  static int  SpawnGuiAndRun(int argc, char **argv);
  static void SpawnGui();

  // These are more like ReveManager stuff.
  TGListTree*     GetListTree() const;
  TGListTreeItem* AddToListTree(RenderElement* re, Bool_t open, TGListTree* lt=0);
  void            RemoveFromListTree(RenderElement* re, TGListTree* lt, TGListTreeItem* lti);

  TGListTreeItem* AddEvent(EventBase* event);
  TGListTreeItem* AddRenderElement(RenderElement* rnr_element,
				   RenderElement* parent=0);
  TGListTreeItem* AddGlobalRenderElement(RenderElement* rnr_element,
					 RenderElement* parent=0);

  void RemoveRenderElement(RenderElement* rnr_element, RenderElement* parent);
  void PreDeleteRenderElement(RenderElement* rnr_element);

  void   RenderElementSelect(RenderElement* rnr_element);
  Bool_t RenderElementPaste(RenderElement* rnr_element);
  void   RenderElementChecked(RenderElement* rnrEl, Bool_t state);

  void NotifyBrowser(TGListTreeItem* parent_lti=0);
  void NotifyBrowser(RenderElement* parent);

  // Hmmph ... geometry management?
  TGeoManager* GetGeometry(const TString& filename);

  void SetStatusLine(const char* text);
  void ThrowException(const char* text="foo");

  ClassDef(ReveManager, 0); // Reve application manager.
};

} // namespace Reve

extern Reve::ReveManager* gReve;

#endif
