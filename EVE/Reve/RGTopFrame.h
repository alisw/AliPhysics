// $Header$

#ifndef REVE_RGTopFrame_H
#define REVE_RGTopFrame_H

#include <TGFrame.h>
#include <TGeoManager.h>
#include <TGStatusBar.h>

#include <map>

class TMacro;
class TFolder;
class TCanvas;
class TGTab;
class TGStatusBar;
class TGListTree;
class TGListTreeItem;

namespace Reve {

class VSDSelector;
class RGBrowser;
class RGEditor;

class RenderElement;
class RenderElementList;

class EventBase;


class RGTopFrame : public TGMainFrame
{
  RGTopFrame(const RGTopFrame&);            // Not implemented
  RGTopFrame& operator=(const RGTopFrame&); // Not implemented

public:
  enum LookType_e { LT_Classic, LT_Editor, LT_GLViewer };

  class RedrawDisabler
  {
  private:
    RGTopFrame* fFrame;
  public:
    RedrawDisabler(RGTopFrame* f) : fFrame(f)
    { if (fFrame) fFrame->DisableRedraw(); }
    ~RedrawDisabler()
    { if (fFrame) fFrame->EnableRedraw(); }
  };

private:
  TGCompositeFrame    *fMasterFrame;
  TGTab               *fMasterTab;

  TCanvas             *fGLCanvas;
  VSDSelector         *fSelector;
  RGBrowser           *fBrowser;
  TGStatusBar         *fStatusBar;
  const char          *fVSDFile;

  TFolder             *fMacroFolder;

  RGEditor            *fEditor;

  EventBase           *fCurrentEvent;
  RenderElementList   *fGlobalStore;

  Int_t                fRedrawDisabled;
  Bool_t               fResetCameras;
  Bool_t               fTimerActive;
  TTimer               fRedrawTimer;

protected:
  LookType_e           fLook;

  std::map<TString, TGeoManager*> fGeometries;

public:
  RGTopFrame(const TGWindow *p, UInt_t w, UInt_t h, LookType_e look=LT_Classic);

  TGCompositeFrame* GetMasterFrame() { return fMasterFrame; }
  TGTab*            GetMasterTab()   { return fMasterTab; }

  TCanvas*     GetGLCanvas()   { return fGLCanvas; }
  VSDSelector* GetSelector()   { return fSelector; }
  RGBrowser*   GetBrowser()    { return fBrowser; }
  TGStatusBar* GetStatusBar()  { return fStatusBar; }

  TCanvas*     AddCanvasTab(const char* name);

  TGListTree*        GetListTree() const;
  EventBase*         GetCurrentEvent() const { return fCurrentEvent; }
  RenderElementList* GetGlobalStore()  const { return fGlobalStore; }

  TFolder* GetMacroFolder() const { return fMacroFolder; }
  TMacro*  GetMacro(const Text_t* name) const;

  RGEditor* GetEditor() const { return fEditor; }
  void EditRenderElement(RenderElement* rnr_element);

  void DisableRedraw() { ++fRedrawDisabled; }
  void EnableRedraw()  { --fRedrawDisabled; if(fRedrawDisabled <= 0) Redraw3D(); }

  void Redraw3D(Bool_t resetCameras=kFALSE) {
    if(fRedrawDisabled <= 0 && !fTimerActive) RegisterRedraw3D();
    if(resetCameras) fResetCameras = kTRUE;
  }
  void RegisterRedraw3D();
  void DoRedraw3D();

  static int  SpawnGuiAndRun(int argc, char **argv);
  static void SpawnGui(LookType_e revemode = LT_Editor);

  // These are more like ReveManager stuff.
  TGListTreeItem* AddEvent(EventBase* event); // Could have Reve::Event ...
  TGListTreeItem* AddRenderElement(RenderElement* rnr_element);
  TGListTreeItem* AddRenderElement(RenderElement* parent, RenderElement* rnr_element);
  TGListTreeItem* AddGlobalRenderElement(RenderElement* rnr_element);
  TGListTreeItem* AddGlobalRenderElement(RenderElement* parent, RenderElement* rnr_element);

  void RemoveRenderElement(RenderElement* parent, RenderElement* rnr_element);
  void PreDeleteRenderElement(RenderElement* rnr_element);

  void DrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad=0);
  void UndrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad=0);

  void RenderElementChecked(TObject* obj, Bool_t state);

  void NotifyBrowser(TGListTreeItem* parent_lti=0);
  void NotifyBrowser(RenderElement* parent);

  // Hmmph ... geometry management?
  TGeoManager* GetGeometry(const TString& filename);

  void ThrowException(const char* text="foo");

  ClassDef(RGTopFrame, 0);
};

} // namespace Reve

extern Reve::RGTopFrame* gReve;

#endif
