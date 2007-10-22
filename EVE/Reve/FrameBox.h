// $Header$

#ifndef REVE_FrameBox_H
#define REVE_FrameBox_H

#include <Reve/Reve.h>

#include <TObject.h>

namespace Reve {

class FrameBox : public TObject, public ReferenceBackPtr
{
  friend class FrameBoxGL;

public:
  enum FrameType_e  { FT_None, FT_Quad, FT_Box };

private:
  FrameBox(const FrameBox&);            // Not implemented
  FrameBox& operator=(const FrameBox&); // Not implemented

protected:
  FrameType_e  fFrameType;
  Int_t        fFrameSize;
  Float_t     *fFramePoints;  //[fFrameSize]

  Float_t      fFrameWidth;
  Color_t      fFrameColor;
  Color_t      fBackColor;
  UChar_t      fFrameRGBA[4];
  UChar_t      fBackRGBA[4];
  Bool_t       fFrameFill;
  Bool_t       fDrawBack;

public:
  FrameBox();
  virtual ~FrameBox();

  void SetAAQuadXY(Float_t x, Float_t y, Float_t z, Float_t dx, Float_t dy);
  void SetAAQuadXZ(Float_t x, Float_t y, Float_t z, Float_t dx, Float_t dz);

  void SetAABox(Float_t x,  Float_t y,  Float_t z,
		Float_t dx, Float_t dy, Float_t dz);

  void SetAABoxCenterHalfSize(Float_t x,  Float_t y,  Float_t z,
                              Float_t dx, Float_t dy, Float_t dz);

  // ----------------------------------------------------------------

  FrameType_e  GetFrameType()   const { return fFrameType; }
  Int_t        GetFrameSize()   const { return fFrameSize; }
  Float_t*     GetFramePoints() const { return fFramePoints; }

  Float_t GetFrameWidth() const    { return fFrameWidth; }
  void    SetFrameWidth(Float_t f) { fFrameWidth = f;    }

  Color_t  GetFrameColor() const { return fFrameColor; }
  Color_t* PtrFrameColor() { return &fFrameColor; }
  UChar_t* GetFrameRGBA()  { return fFrameRGBA;  }

  void SetFrameColor(Color_t ci);
  void SetFrameColor(Pixel_t pix);
  void SetFrameColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

  Color_t  GetBackColor() const { return fBackColor; }
  Color_t* PtrBackColor() { return &fBackColor; }
  UChar_t* GetBackRGBA()  { return fBackRGBA;  }

  void SetBackColor(Color_t ci);
  void SetBackColor(Pixel_t pix);
  void SetBackColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

  Bool_t GetFrameFill() const   { return fFrameFill; }
  void   SetFrameFill(Bool_t f) { fFrameFill = f;    }

  Bool_t GetDrawBack() const   { return fDrawBack; }
  void   SetDrawBack(Bool_t f) { fDrawBack = f;    }

  ClassDef(FrameBox, 1);
}; // endclass FrameBox

}

#endif
