#ifndef REVE_RGBrowser_H
#define REVE_RGBrowser_H

#include <TGFrame.h>
#include <TGButton.h>
#include <TGListTree.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>

#include <TParticle.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGeoVolume.h>
#include <TContextMenu.h>

namespace Reve {

class RGEditor;

/**************************************************************************/
// TG wrappers used in detail view.
/**************************************************************************/

class ReveValuator : public TGNumberEntry
{
protected:
  void* fUserData;

public:
  ReveValuator(const TGWindow* parent = 0, Double_t val = 0, Int_t digitwidth = 5, Int_t id = -1, TGNumberFormat::EStyle style = kNESReal, TGNumberFormat::EAttribute attr = kNEAAnyNumber, TGNumberFormat::ELimit limits = kNELNoLimits, Double_t min = 0, Double_t max = 1) :
    TGNumberEntry(parent, val, digitwidth, id, style, attr, limits, min, max),
    fUserData(0)
  {}
  virtual ~ReveValuator();

  void* GetUserData() const   { return fUserData; }
  void  SetUserData(void* ud) { fUserData = ud; }

  ClassDef(ReveValuator, 1);
};


class ReveColorSelect: public TGColorSelect
{
 public:
  ReveColorSelect(const TGWindow* p = 0, Pixel_t color = 0, Int_t id = -1) :
    TGColorSelect(p, color,id)
  {}
  virtual ~ReveColorSelect() {}

  void UpdateColor(Pixel_t col){
    fColor=col; 
    fDrawGC.SetForeground(col);  
    gClient->NeedRedraw(this);
  }

  ClassDef(ReveColorSelect, 1);
};


/**************************************************************************/
// RGBrowser
/**************************************************************************/

class RGBrowser : public TGCompositeFrame
{
protected:
  TGCompositeFrame* fMainFrame;
  TGVerticalFrame*  fV1;
  TGVerticalFrame*  fV2;

  TGCompositeFrame* fSelectionFrame; // in fact list-tree frame
  TGCanvas*         fTreeView;

  TGCanvas*         fCanvasWindow;
  TGCompositeFrame* fDisplayFrame;   // detailed-vire frame, used in Classic look
  
  TGListTree*       fListTree;
  TContextMenu*     fCtxMenu;

 protected:
  void SetupCintExport(TClass* cl);

 public:
  RGBrowser(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~RGBrowser() {}

  void SetupClassicLook();
  void SetupEditorLook(RGEditor* editor);
  void SetupGLViewerLook(RGEditor* editor, TVirtualPad* glpad);

  void RedrawListTree();

  void ItemClicked(TGListTreeItem *entry, Int_t btn, Int_t x, Int_t y);
  void ExportToCINT(Text_t* var_name, TObject* obj);
  void DisplayChildren(TGListTreeItem *entry, Int_t btn);

  void SetVolumeColor(UInt_t col);
  void NodeVis(Bool_t vis);
  void VolumeDaughterVis(Bool_t vis);

  void DbClickListItem(TGListTreeItem* item, Int_t btn);
  void UpdateListItems(TGListTreeItem* item, Int_t btn);
  void SetTransparency(Long_t val);

  // TrackRnrStyle
  void SetMaxR(Long_t);
  void SetMaxZ(Long_t);
  void SetMaxOrbs(Long_t);

  // GuiPOintRnrStyle
  void SetMarkerStyle(Long_t style);


  TGListTree* GetListTree() { return fListTree; }

  ClassDef(RGBrowser, 1);
};

} // namespace Reve

#endif
