/***********************************************************************
  This editor appears in the Reve window when cascades are visualize.
It allows to select the cascades as a function of some useful parameters.

Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/


#ifndef REVE_CascadeEditors_H
#define REVE_CascadeEditors_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TRootEmbeddedCanvas;
class TH1F;
class TGCompositeFrame;
class TGTab;

namespace Reve {

class RGValuator;
class RGDoubleValuator;

class CascadeList;

class CascadeListEditor : public TGedFrame
{
  CascadeListEditor(const CascadeListEditor&);            // Not implemented
  CascadeListEditor& operator=(const CascadeListEditor&); // Not implemented

protected:
  CascadeList* fMList; // fModel dynamic-casted to CascadeListEditor

  TGCheckButton*     fRnrV0Daughters;
  TGCheckButton*     fRnrV0path;
  TGCheckButton*     fRnrVtx;
  TGCheckButton*     fRnrBach;
  TGCheckButton*     fRnrCasPath;

  TGTab *fMainTabA;
  TGTab *fMainTabB;
  TGTab *fTabA[3];
  TGTab *fTabB[3];
  static const Int_t fgkNRange = 14;
  RGDoubleValuator    *fRange[fgkNRange];

  static const Int_t fgkNCanvas = 15;
  TRootEmbeddedCanvas *fCanvasA[fgkNCanvas];
  TRootEmbeddedCanvas *fCanvasB[fgkNCanvas];

  TGCompositeFrame*  AddTab(TGTab *tab, Int_t i, Int_t can, char *name);
  TGCompositeFrame** CreateTab(TGTab **pMainTab, TGTab **ptab, Int_t can);

  void UpdateAll(Int_t iCanA);
  void AddSelectTab();
  void AddSeeTab();
  void AddValuator(TGCompositeFrame* frame, char *name,
		   Float_t min, Float_t max, Int_t pres, char *func,
		   Int_t iHist);


public:
  CascadeListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~CascadeListEditor();

  virtual void SetModel(TObject* obj);
  void DoRnrV0Daughters();
  void DoRnrV0path();
  void DoRnrVtx();
  void DoRnrBach();
  void DoRnrCasPath();

  void FillCanvas();
  void UpdateSelectedTab();
  void AdjustHist(Int_t iHist);
  void ResetCuts();

  void MassXiRange();
  void MassOmegaRange();
  void IndexRange();
  void CosPointingRange();
  void BachV0DCARange();
  void RadiusRange();
  void PtRange();
  void PseudoRapRange();
  void NegPtRange();
  void NegEtaRange();
  void PosPtRange();
  void PosEtaRange();
  void BachPtRange();
  void BachEtaRange();

  ClassDef(CascadeListEditor, 1); // Editor for CascadeList
}; // endclass CascadeListEditor

}

#endif
