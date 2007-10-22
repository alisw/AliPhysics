#ifndef REVE_RGBrowser_H
#define REVE_RGBrowser_H

#include "TGNewBrowser.h"
#include <TGListTree.h>

#include <TContextMenu.h>

class TGFileBrowser;
class TGSplitter;

namespace Reve {

class RGEditor;

class RGLTEFrame : public TGMainFrame
{
  RGLTEFrame(const RGLTEFrame&);            // Not implemented
  RGLTEFrame& operator=(const RGLTEFrame&); // Not implemented

  friend class ReveManager;

protected:
  TGCompositeFrame *fFrame;
  TGCompositeFrame *fLTFrame;

  TGCanvas         *fLTCanvas;
  TGListTree       *fListTree;
  TGSplitter       *fSplitter;
  RGEditor         *fEditor;

  TContextMenu     *fCtxMenu;

  TGListTreeItem   *fNewSelected;

  void ResetSelectedTimer(TGListTreeItem* lti);

public:
  RGLTEFrame(const Text_t* name, Int_t width=250, Int_t height=700);
  virtual ~RGLTEFrame();

  void ReconfToHorizontal();
  void ReconfToVertical();

  TGListTree* GetListTree() { return fListTree; }

  void ItemChecked(TObject* obj, Bool_t state);
  void ItemClicked(TGListTreeItem *entry, Int_t btn, Int_t x, Int_t y);
  void ItemDblClicked(TGListTreeItem* item, Int_t btn);
  void ItemKeyPress(TGListTreeItem *entry, UInt_t keysym, UInt_t mask);

  void ResetSelected();

  ClassDef(RGLTEFrame, 0);
};

// ----------------------------------------------------------------

class RGBrowser : public TGNewBrowser
{
  RGBrowser(const RGBrowser&);            // Not implemented
  RGBrowser& operator=(const RGBrowser&); // Not implemented

protected:
  void SetupCintExport(TClass* cl);
  void CalculateReparentXY(TGObject* parent, Int_t& x, Int_t& y);

  TGFileBrowser    *fFileBrowser;
  TGPopupMenu      *fRevePopup;

 public:
  RGBrowser(UInt_t w, UInt_t h);
  virtual ~RGBrowser() {}

  void InitPlugins();

  TGFileBrowser *GetFileBrowser() const { return fFileBrowser; }

  void ReveMenu(Int_t id);

  ClassDef(RGBrowser, 0);
};

} // namespace Reve

#endif
