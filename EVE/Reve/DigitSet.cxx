// $Header$

#include "DigitSet.h"

#include "ReveManager.h"

#include <TColor.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

using namespace Reve;

//______________________________________________________________________
// DigitSet
//

ClassImp(DigitSet)

DigitSet::DigitSet(const Text_t* n, const Text_t* t) :
  RenderElement   (),
  TNamed          (n, t),

  fDefaultValue   (kMinInt),
  fValueIsColor   (kFALSE),
  fOwnIds         (kFALSE),
  fPlex           (),
  fLastDigit      (0),

  fFrame          (0),
  fPalette        (0),
  fRenderMode     (RM_Fill),
  fDisableLigting (kTRUE),
  fEmitSignals    (kFALSE),
  fHistoButtons   (kTRUE),
  fHMTrans        ()
{

}

DigitSet::~DigitSet()
{
  SetFrame(0);
  SetPalette(0);
  if (fOwnIds)
    ReleaseIds();
}

/**************************************************************************/

DigitSet::DigitBase* DigitSet::NewDigit()
{
  fLastDigit = new (fPlex.NewAtom()) DigitBase(fDefaultValue);
  return fLastDigit;
}

void DigitSet::ReleaseIds()
{
  VoidCPlex::iterator qi(fPlex);
  while (qi.next()) {
    DigitBase& q = * (DigitBase*) qi();
    if (q.fId.GetObject()) {
      delete q.fId.GetObject();
      q.fId = 0;
    }
  }
}

/**************************************************************************/
/**************************************************************************/

void DigitSet::SetMainColor(Color_t color)
{
  // Override from RenderElement, forward to Frame.

  if (fFrame) {
    fFrame->SetFrameColor(color);
    fFrame->UpdateBackPtrItems();
  }
  gReve->Redraw3D();
}

/**************************************************************************/
/**************************************************************************/

void DigitSet::RefitPlex()
{
  // Instruct underlying memory allocator to regroup itself into a
  // contiguous memory chunk.

  fPlex.Refit();
}

/**************************************************************************/

void DigitSet::ScanMinMaxValues(Int_t& min, Int_t& max)
{
  if (fValueIsColor || fPlex.Size() == 0) return;
  min = kMaxInt;
  max = kMinInt;
  for (Int_t c=0; c<fPlex.VecSize(); ++c)
  {
    Char_t* a = fPlex.Chunk(c);
    Int_t   n = fPlex.NAtoms(c);
    while (n--)
    {
      Int_t v = ((DigitBase*)a)->fValue;
      if (v < min) min = v;
      if (v > max) max = v;
      a += fPlex.S();
    }
  }
  if (min == max)
    --min;
}

/**************************************************************************/


void DigitSet::DigitValue(Int_t value)
{
  fLastDigit->fValue = value;
}

void DigitSet::DigitColor(Color_t ci)
{
  ColorFromIdx(ci, (UChar_t*) & fLastDigit->fValue, kTRUE);
}

void DigitSet::DigitColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  UChar_t* x = (UChar_t*) & fLastDigit->fValue;
  x[0] = r; x[1] = g; x[2] = b; x[3] = a;
}

void DigitSet::DigitColor(UChar_t* rgba)
{
  UChar_t* x = (UChar_t*) & fLastDigit->fValue;
  x[0] = rgba[0]; x[1] = rgba[1]; x[2] = rgba[2]; x[3] = rgba[3];
}

/**************************************************************************/

void DigitSet::DigitId(TObject* id)
{
  fLastDigit->fId = id;
}

/**************************************************************************/
/**************************************************************************/

void DigitSet::Paint(Option_t* /*option*/)
{
  static const Exc_t eH("DigitSet::Paint ");

  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = fFrame ? fFrame->GetFrameColor() : 1;
  buff.fTransparency = 0;
  fHMTrans.SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections != TBuffer3D::kNone)
    Error(eH, "only direct GL rendering supported.");
}

void DigitSet::DigitSelected(Int_t idx)
{
  if (fEmitSignals) {
    CtrlClicked(this, idx);
  } else {
    DigitBase* qb = GetDigit(idx);
    TObject* obj = qb->fId.GetObject();
    printf("DigitSet::DigitSelected idx=%d, value=%d, obj=0x%lx\n",
	   idx, qb->fValue, (ULong_t)obj);
    if (obj)
      obj->Print();
  }
}

void DigitSet::CtrlClicked(DigitSet* qs, Int_t idx)
{
  Long_t args[2];
  args[0] = (Long_t) qs;
  args[1] = (Long_t) idx;

  Emit("CtrlClicked(Reve::DigitSet*, Int_t)", args);
}

/**************************************************************************/
// Getters / Setters for Frame, RGBAPalette, ZTrans
/**************************************************************************/

void DigitSet::SetFrame(FrameBox* b)
{
  if (fFrame == b) return;
  if (fFrame) fFrame->DecRefCount(this);
  fFrame = b;
  if (fFrame) {
    fFrame->IncRefCount(this);
    SetMainColorPtr(fFrame->PtrFrameColor());
  } else {
    SetMainColorPtr(0);
  }
}

void DigitSet::SetPalette(RGBAPalette* p)
{
  if (fPalette == p) return;
  if (fPalette) fPalette->DecRefCount();
  fPalette = p;
  if (fPalette) fPalette->IncRefCount();
}

RGBAPalette* DigitSet::AssertPalette()
{
  if (fPalette == 0) {
    fPalette = new RGBAPalette;
    if (!fValueIsColor) {
      Int_t min, max;
      ScanMinMaxValues(min, max);
      fPalette->SetLimits(min, max);
      fPalette->SetMinMax(min, max);
    }
  }
  return fPalette;
}
