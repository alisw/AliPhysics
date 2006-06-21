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
class TGStatusBar;
class TGListTree;
class TGListTreeItem;

namespace Reve {

class VSDSelector;
class RGBrowser;
class RGEditor;
class RenderElement;

class RGTopFrame : public TGMainFrame
{
public:
  enum LookType_e { LT_Classic, LT_Editor, LT_GLViewer };

private:
  void                Init();
  TCanvas             *fCC;
  TCanvas             *fHistoCanvas;
  VSDSelector         *fSelector;
  RGBrowser           *fBrowser;
  TGStatusBar         *fStatusBar;
  const char          *fVSDFile;

  TFolder             *fMacroFolder;

  RGEditor            *fEditor;

  TObject             *fCurrentEvent;
  TGListTreeItem      *fCurrentEventLTI;

  TGListTreeItem      *fGeometryLTI;

  Bool_t               fRedrawDisabled;
  Bool_t               fTimerActive;
  TTimer               fRedrawTimer;

protected:
  LookType_e           fLook;

  std::map<TString, TGeoManager*> fGeometries;

public:
  RGTopFrame(const TGWindow *p, UInt_t w, UInt_t h, LookType_e look=LT_Classic);

  TCanvas*     GetCC()         { return fCC; }
  VSDSelector* GetSelector()   { return fSelector; }
  RGBrowser*   GetBrowser()    { return fBrowser; }
  TGStatusBar* GetStatusBar()  { return fStatusBar; }

  TGListTree*     GetListTree();
  TGListTreeItem* GetEventTreeItem();
  TGListTreeItem* GetGlobalTreeItem();

  TFolder* GetMacroFolder() const { return fMacroFolder; }
  TMacro*  GetMacro(const Text_t* name) const;

  RGEditor* GetEditor() const { return fEditor; }
  void EditRenderElement(RenderElement* rnr_element);

  void DisableRedraw() { fRedrawDisabled = true; }
  void EnableRedraw()  { fRedrawDisabled = false; Redraw3D(); }

  void Redraw3D() { if(!fRedrawDisabled && !fTimerActive) RegisterRedraw3D(); }
  void RegisterRedraw3D();
  void DoRedraw3D();

  static int SpawnGuiAndRun(int argc, char **argv);

  // These are more like ReveManager stuff.
  TGListTreeItem* AddEvent(TObject* event); // Could have Reve::Event ...
  TGListTreeItem* AddRenderElement(RenderElement* rnr_element);
  TGListTreeItem* AddRenderElement(TGListTreeItem* parent, RenderElement* rnr_element);
  TGListTreeItem* AddGlobalRenderElement(RenderElement* rnr_element);
  TGListTreeItem* AddGlobalRenderElement(TGListTreeItem* parent, RenderElement* rnr_element);

  void DrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad=0);
  void UndrawRenderElement(RenderElement* rnr_element, TVirtualPad* pad=0);

  void RenderElementChecked(TObject* obj, Bool_t state);

  void NotifyBrowser(TGListTreeItem* parent=0);

  // Hmmph ... geometry management?
  TGeoManager* GetGeometry(const TString& filename);

  ClassDef(RGTopFrame, 0);
};

} // namespace Reve

extern Reve::RGTopFrame* gReve;

#endif
