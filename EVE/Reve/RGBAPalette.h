// $Header$

#ifndef REVE_RGBAPalette_H
#define REVE_RGBAPalette_H

#include <Reve/Reve.h>

#include <TObject.h>

namespace Reve {

class RGBAPalette : public TObject, public ReferenceCount
{
  friend class RGBAPaletteEditor;
  friend class RGBAPaletteSubEditor;

private:
  RGBAPalette(const RGBAPalette&);            // Not implemented
  RGBAPalette& operator=(const RGBAPalette&); // Not implemented

protected:
  Int_t     fLowLimit;  // Low  limit for Min/Max values (used by editor)
  Int_t     fHighLimit; // High limit for Min/Max values (used by editor)
  Int_t     fMinVal;
  Int_t     fMaxVal;
  Int_t     fNBins;
  Bool_t    fCutLow;    // Instruct renderer not to display quoads below fMinVal
  Bool_t    fCutHigh;   // Instruct renderer not to display quoads above fMaxVal
  Bool_t    fInterpolate;
  Bool_t    fWrap;
  Color_t   fDefaultColor;
  UChar_t   fDefaultRGBA[4];

  mutable UChar_t* fColorArray; //[4*fNBins]

  void SetupColor(Int_t val, UChar_t* pix) const;

  static RGBAPalette* fgDefaultPalette;

public:
  RGBAPalette();
  RGBAPalette(Int_t min, Int_t max);
  RGBAPalette(Int_t min, Int_t max, Bool_t interp, Bool_t wrap);
  virtual ~RGBAPalette();

  void SetupColorArray() const;
  void ClearColorArray();

  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix, Bool_t alpha=kTRUE) const;

  Bool_t   WithinVisibleRange(Int_t val) const;

  Int_t  GetMinVal() const      { return fMinVal; }
  Int_t  GetMaxVal() const      { return fMaxVal; }
  Bool_t GetInterpolate() const { return fInterpolate; }
  Bool_t GetWrap() const        { return fWrap; }

  void   SetLimits(Int_t low, Int_t high);
  void   SetMinMax(Int_t min, Int_t max);
  void   SetMin(Int_t min);
  void   SetMax(Int_t max);
  void   SetInterpolate(Bool_t b);
  void   SetWrap(Bool_t b);

  Color_t  GetDefaultColor() const { return fDefaultColor; }
  Color_t* PtrDefaultColor() { return &fDefaultColor; }
  UChar_t* GetDefaultRGBA()  { return fDefaultRGBA;  }
  const UChar_t* GetDefaultRGBA() const { return fDefaultRGBA;  }

  void   SetDefaultColor(Color_t ci);
  void   SetDefaultColor(Pixel_t pix);
  void   SetDefaultColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

  // ?? Should we emit some *SIGNALS* ??
  // ?? Should we have a RendererTimeStamp ??

  ClassDef(RGBAPalette, 1);
}; // endclass RGBAPalette


inline UChar_t* RGBAPalette::ColorFromArray(Int_t val) const
{
  if(!fColorArray)  SetupColorArray();
  if(val < fMinVal) val = fWrap ? ((val+1-fMinVal)%fNBins + fMaxVal) : fMinVal;
  if(val > fMaxVal) val = fWrap ? ((val-1-fMaxVal)%fNBins + fMinVal) : fMaxVal;
  return fColorArray + 4 * (val - fMinVal);
}

inline void RGBAPalette::ColorFromArray(Int_t val, UChar_t* pix, Bool_t alpha) const
{
  UChar_t* c = ColorFromArray(val);
  pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2];
  if (alpha) pix[3] = c[3];
}

inline Bool_t RGBAPalette::WithinVisibleRange(Int_t val) const
{
  if ((val < fMinVal && fCutLow) || (val > fMaxVal && fCutHigh))
    return kFALSE;
  else
    return kTRUE;
}

}

#endif
