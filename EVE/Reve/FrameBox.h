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

  Color_t      fFrameColor;
  UChar_t      fFrameRGBA[4];

  // !!!! Back-refs to RenderElement consumers (as for palette)

public:
  FrameBox();
  virtual ~FrameBox();

  void SetAAQuadXY(Float_t x, Float_t y, Float_t z, Float_t dx, Float_t dy);
  void SetAAQuadXZ(Float_t x, Float_t y, Float_t z, Float_t dx, Float_t dz);

  void SetAABox(Float_t x,  Float_t y,  Float_t z,
		Float_t dx, Float_t dy, Float_t dz);

  // ----------------------------------------------------------------

  Color_t  GetFrameColor() const { return fFrameColor; }
  Color_t* PtrFrameColor() { return &fFrameColor; }
  UChar_t* GetFrameRGBA()  { return fFrameRGBA;  }

  void SetFrameColor(Color_t ci);
  void SetFrameColor(Pixel_t pix);
  void SetFrameColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

  ClassDef(FrameBox, 1);
}; // endclass FrameBox

}

#endif
