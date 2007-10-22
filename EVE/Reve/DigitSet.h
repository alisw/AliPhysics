// $Header$

#ifndef REVE_DigitSetBase_H
#define REVE_DigitSetBase_H

#include <TNamed.h>
#include <TQObject.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/Reve.h>
#include <Reve/RenderElement.h>
#include <Reve/FrameBox.h>
#include <Reve/RGBAPalette.h>
#include <Reve/Plex.h>
#include <Reve/ZTrans.h>

#include <TObject.h>

namespace Reve {

class DigitSet : public RenderElement,
                 public TNamed, public TQObject,
                 public TAtt3D,
                 public TAttBBox
{
  friend class DigitSetEditor;

  DigitSet(const DigitSet&);            // Not implemented
  DigitSet& operator=(const DigitSet&); // Not implemented

public:
  enum RenderMode_e { RM_AsIs, RM_Line, RM_Fill };

protected:
  struct DigitBase
  {
    Int_t fValue;
    TRef  fId;

    // Here could have additional integer (like time, second threshold).

    DigitBase(Int_t v=0) : fValue(v), fId() {}
  };

  Int_t             fDefaultValue;
  Bool_t            fValueIsColor;
  Bool_t            fOwnIds;       //Flag specifying if id-objects are owned by the DigitSet
  VoidCPlex         fPlex;
  DigitBase*        fLastDigit;    //!

  FrameBox*         fFrame;
  RGBAPalette*      fPalette;
  RenderMode_e      fRenderMode;
  Bool_t            fDisableLigting;
  Bool_t            fEmitSignals;
  Bool_t            fHistoButtons;
  ZTrans            fHMTrans;

  DigitBase* NewDigit();
  void       ReleaseIds();

public:
  DigitSet(const Text_t* n="DigitSet", const Text_t* t="");
  virtual ~DigitSet();

  virtual Bool_t CanEditMainColor() { return kTRUE; }
  virtual void   SetMainColor(Color_t color);

  // virtual void Reset(QuadType_e quadType, Bool_t valIsCol, Int_t chunkSize);

  void RefitPlex();
  void ScanMinMaxValues(Int_t& min, Int_t& max);

  // --------------------------------

  void DigitValue(Int_t value);
  void DigitColor(Color_t ci);
  void DigitColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);
  void DigitColor(UChar_t* rgba);

  void DigitId(TObject* id);

  Bool_t GetOwnIds() const    { return fOwnIds; }
  void   SetOwnIds(Bool_t o)  { fOwnIds = o; }

  DigitBase* GetDigit(Int_t n) { return (DigitBase*) fPlex.Atom(n);   }
  TObject*   GetId(Int_t n)   { return GetDigit(n)->fId.GetObject(); }

  // --------------------------------

  // virtual void ComputeBBox(); // implement in subclass
  virtual void Paint(Option_t* option="");

  virtual void DigitSelected(Int_t idx);
  virtual void CtrlClicked(DigitSet* qs, Int_t idx); // *SIGNAL*

  // --------------------------------

  VoidCPlex* GetPlex() { return &fPlex; }

  FrameBox* GetFrame() const { return fFrame; }
  void      SetFrame(FrameBox* b);
 
  Bool_t GetValueIsColor()  const { return fValueIsColor; }

  RGBAPalette* GetPalette() const { return fPalette; }
  void         SetPalette(RGBAPalette* p);
  RGBAPalette* AssertPalette();

  RenderMode_e  GetRenderMode() const { return fRenderMode; }
  void SetRenderMode(RenderMode_e rm) { fRenderMode = rm; }

  Bool_t GetEmitSignals() const   { return fEmitSignals; }
  void   SetEmitSignals(Bool_t f) { fEmitSignals = f; }

  Bool_t GetHistoButtons() const   { return fHistoButtons; }
  void   SetHistoButtons(Bool_t f) { fHistoButtons = f; }

  ZTrans& RefHMTrans() { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }

  ClassDef(DigitSet, 1);
}; // endclass DigitSet

}

#endif
