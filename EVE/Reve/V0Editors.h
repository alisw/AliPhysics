
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/***********************************************************************
  This editor appears in the Reve window when v0 are visualize.
It allows to select the v0 as a function of some useful parameters.

Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/


#ifndef REVE_V0Editors_H
#define REVE_V0Editors_H

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

class V0List;

class V0ListEditor : public TGedFrame
{
  V0ListEditor(const V0ListEditor&);            // Not implemented
  V0ListEditor& operator=(const V0ListEditor&); // Not implemented

protected:
  V0List* fMList; // fModel dynamic-casted to V0ListEditor

  TGCheckButton*     fRnrV0sDaugh;
  TGCheckButton*     fRnrV0vtx;
  TGCheckButton*     fRnrV0path;

  TGTab *fMainTabA;
  TGTab *fMainTabB;
  TGTab *fTabA[3];
  TGTab *fTabB[3];
  static const Int_t fgkNRange = 13;
  RGDoubleValuator    *fRange[fgkNRange];

  static const Int_t fgkNCanvas = 14;
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
  V0ListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~V0ListEditor();

  virtual void SetModel(TObject* obj);
  void DoRnrDaughters();
  void DoRnrV0vtx();
  void DoRnrV0path();
  void FillCanvas();
  void UpdateSelectedTab();
  void AdjustHist(Int_t iHist);
  void ResetCuts();

  void MassK0sRange();
  void MassLamRange();
  void MassAntiLamRange();
  void ESDv0IndexRange();
  void CosPointingRange();
  void DaughterDcaRange();
  void RadiusRange();
  void PtRange();
  void EtaRange();
  void NegPtRange();
  void NegEtaRange();
  void PosPtRange();
  void PosEtaRange();


  ClassDef(V0ListEditor, 1); // Editor for V0List
}; // endclass V0ListEditor

} // end namespace Reve

#endif
