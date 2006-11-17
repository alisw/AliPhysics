// $Header$

#include "RGBAPalette.h"

#include <TColor.h>
#include <TStyle.h>

using namespace Reve;

//______________________________________________________________________
// RGBAPalette
//

ClassImp(RGBAPalette)

RGBAPalette::RGBAPalette() :
  TObject(),
  Reve::ReferenceCount(),

  fLowLimit(0), fHighLimit(0), fMinVal(0), fMaxVal(0), fNBins(0),

  fInterpolate     (kFALSE),
  fShowDefValue    (kTRUE),
  fUnderflowAction (LA_Cut),
  fOverflowAction  (LA_Clip),

  fDefaultColor(0),
  fUnderColor  (1),
  fOverColor   (2),
  fColorArray  (0)
{
  SetLimits(0, 1024);
  SetMinMax(0, 100);
}

RGBAPalette::RGBAPalette(Int_t min, Int_t max, Bool_t interp, Bool_t showdef) :
  TObject(),
  Reve::ReferenceCount(),

  fLowLimit(0), fHighLimit(0), fMinVal(0), fMaxVal(0), fNBins(0),

  fInterpolate     (interp),
  fShowDefValue    (showdef),
  fUnderflowAction (LA_Cut),
  fOverflowAction  (LA_Clip),

  fDefaultColor(0),
  fUnderColor  (1),
  fOverColor   (2),
  fColorArray  (0)
{
  SetLimits(0, 1024);
  SetMinMax(min, max);
}

RGBAPalette::~RGBAPalette()
{
  delete [] fColorArray;
}

/**************************************************************************/

void RGBAPalette::SetupColor(Int_t val, UChar_t* pixel) const
{
  using namespace TMath;
  Float_t div  = Max(1, fMaxVal - fMinVal);
  Int_t   nCol = gStyle->GetNumberOfColors();

  Float_t f;
  if      (val >= fMaxVal) f = nCol - 1;
  else if (val <= fMinVal) f = 0;
  else                     f = (val - fMinVal)/div*(nCol - 1);

  if (fInterpolate) {
    Int_t  bin = (Int_t) f;
    Float_t f1 = f - bin, f2 = 1.0f - f1;
    ColorFromIdx(f1, gStyle->GetColorPalette(bin),
		 f2, gStyle->GetColorPalette(Min(bin + 1, nCol - 1)),
		 pixel);
  } else {
    ColorFromIdx(gStyle->GetColorPalette((Int_t) Nint(f)), pixel);    
  }
}

void RGBAPalette::SetupColorArray() const
{
  if(fColorArray) // !!!! should reinit anyway, maybe palette in gstyle changed
    return;

  // !!!! probably should store original palette for editing ...

  fColorArray = new UChar_t [4 * fNBins];
  UChar_t* p = fColorArray;
  for(Int_t v=fMinVal; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);
}

void RGBAPalette::ClearColorArray()
{
  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

/**************************************************************************/

void RGBAPalette::SetLimits(Int_t low, Int_t high)
{
  fLowLimit  = low;
  fHighLimit = high;
  Bool_t changed = kFALSE;
  if (fMaxVal < fLowLimit)  { SetMax(fLowLimit);  changed = kTRUE; }
  if (fMinVal < fLowLimit)  { SetMin(fLowLimit);  changed = kTRUE; }
  if (fMinVal > fHighLimit) { SetMin(fHighLimit); changed = kTRUE; }
  if (fMaxVal > fHighLimit) { SetMax(fHighLimit); changed = kTRUE; }
  if (changed)
    ClearColorArray();
}

void RGBAPalette::SetMin(Int_t min)
{
  fMinVal = TMath::Min(min, fMaxVal);
  fNBins  = fMaxVal - fMinVal + 1;
  ClearColorArray();
}

void RGBAPalette::SetMax(Int_t max)
{
  fMaxVal = TMath::Max(max, fMinVal);
  fNBins  = fMaxVal - fMinVal + 1;
  ClearColorArray();
}

void RGBAPalette::SetMinMax(Int_t min, Int_t max)
{
  fMinVal = min;
  fMaxVal = max;
  fNBins  = fMaxVal - fMinVal + 1;
  ClearColorArray();
}

/**************************************************************************/
/**************************************************************************/

void RGBAPalette::SetInterpolate(Bool_t b)
{
  fInterpolate = b;
  ClearColorArray();
}

/**************************************************************************/
/**************************************************************************/

void RGBAPalette::SetDefaultColor(Color_t ci)
{
  fDefaultColor = ci;
  ColorFromIdx(ci, fDefaultRGBA, kTRUE);
}

void RGBAPalette::SetDefaultColor(Pixel_t pix)
{
  SetDefaultColor(Color_t(TColor::GetColor(pix)));
}

void RGBAPalette::SetDefaultColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  fDefaultColor = Color_t(TColor::GetColor(r, g, b));
  fDefaultRGBA[0] = r;
  fDefaultRGBA[1] = g;
  fDefaultRGBA[2] = b;
  fDefaultRGBA[3] = a;
}

/**************************************************************************/

void RGBAPalette::SetUnderColor(Color_t ci)
{
  fUnderColor = ci;
  ColorFromIdx(ci, fUnderRGBA, kTRUE);
}

void RGBAPalette::SetUnderColor(Pixel_t pix)
{
  SetUnderColor(Color_t(TColor::GetColor(pix)));
}

void RGBAPalette::SetUnderColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  fUnderColor = Color_t(TColor::GetColor(r, g, b));
  fUnderRGBA[0] = r;
  fUnderRGBA[1] = g;
  fUnderRGBA[2] = b;
  fUnderRGBA[3] = a;
}

/**************************************************************************/

void RGBAPalette::SetOverColor(Color_t ci)
{
  fOverColor = ci;
  ColorFromIdx(ci, fOverRGBA, kTRUE);
}

void RGBAPalette::SetOverColor(Pixel_t pix)
{
  SetOverColor(Color_t(TColor::GetColor(pix)));
}

void RGBAPalette::SetOverColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  fOverColor = Color_t(TColor::GetColor(r, g, b));
  fOverRGBA[0] = r;
  fOverRGBA[1] = g;
  fOverRGBA[2] = b;
  fOverRGBA[3] = a;
}

/**************************************************************************/
/**************************************************************************/
