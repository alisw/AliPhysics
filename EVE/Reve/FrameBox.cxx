// $Header$

#include <Reve/FrameBox.h>

#include <TColor.h>

using namespace Reve;

//______________________________________________________________________
// FrameBox
//

ClassImp(FrameBox)

FrameBox::FrameBox() :
  fFrameType   (FT_None),
  fFrameSize   (0),
  fFramePoints (0),

  fFrameWidth  (1),
  fFrameColor  (1),
  fBackColor   (0),
  fFrameFill   (kFALSE),
  fDrawBack    (kFALSE)
{
  fFrameRGBA[0] = fFrameRGBA[1] = fFrameRGBA[2] = 0;   fFrameRGBA[3] = 255;
  fBackRGBA [0] = fBackRGBA [1] = fBackRGBA [2] = 255; fBackRGBA [3] = 255;
}

FrameBox::~FrameBox()
{
  delete [] fFramePoints;
}

/**************************************************************************/

void FrameBox::SetAAQuadXY(Float_t x,  Float_t y, Float_t z,
			   Float_t dx, Float_t dy)
{
  fFrameType = FT_Quad;
  fFrameSize = 12;
  delete [] fFramePoints;
  fFramePoints = new Float_t [fFrameSize];
  Float_t* p = fFramePoints;
  p[0] = x;    p[1] = y;    p[2] = z; p += 3;
  p[0] = x+dx; p[1] = y;    p[2] = z; p += 3;
  p[0] = x+dx; p[1] = y+dy; p[2] = z; p += 3;
  p[0] = x ;   p[1] = y+dy; p[2] = z; p += 3;
}

void FrameBox::SetAAQuadXZ(Float_t x,  Float_t y, Float_t z,
			   Float_t dx, Float_t dz)
{
  fFrameType = FT_Quad;
  fFrameSize = 12;
  delete [] fFramePoints;
  fFramePoints = new Float_t [fFrameSize];
  Float_t* p = fFramePoints;
  p[0] = x;    p[1] = y; p[2] = z;    p += 3;
  p[0] = x+dx; p[1] = y; p[2] = z;    p += 3;
  p[0] = x+dx; p[1] = y; p[2] = z+dz; p += 3;
  p[0] = x ;   p[1] = y; p[2] = z+dz; p += 3;
}

void FrameBox::SetAABox(Float_t x,  Float_t y,  Float_t z,
			Float_t dx, Float_t dy, Float_t dz)
{
  fFrameType = FT_Box;
  fFrameSize = 24;
  delete [] fFramePoints;
  fFramePoints = new Float_t [fFrameSize];

  Float_t* p = fFramePoints;
  //bottom
  p[0] = x;       p[1] = y + dy;  p[2] = z;       p += 3;
  p[0] = x + dx;  p[1] = y + dy;  p[2] = z;       p += 3;
  p[0] = x + dx;  p[1] = y;       p[2] = z;       p += 3;
  p[0] = x;       p[1] = y;       p[2] = z;       p += 3;
  //top
  p[0] = x;       p[1] = y + dy;  p[2] = z + dz;  p += 3;
  p[0] = x + dx;  p[1] = y + dy;  p[2] = z + dz;  p += 3;
  p[0] = x + dx;  p[1] = y;       p[2] = z + dz;  p += 3;
  p[0] = x;       p[1] = y;       p[2] = z + dz;
}

void FrameBox::SetAABoxCenterHalfSize(Float_t x,  Float_t y,  Float_t z,
				      Float_t dx, Float_t dy, Float_t dz)
{
  fFrameType = FT_Box;
  fFrameSize = 24;
  delete [] fFramePoints;
  fFramePoints = new Float_t [fFrameSize];

  Float_t* p = fFramePoints;
  //bottom
  p[0] = x - dx;  p[1] = y + dy;  p[2] = z - dz;  p += 3;
  p[0] = x + dx;  p[1] = y + dy;  p[2] = z - dz;  p += 3;
  p[0] = x + dx;  p[1] = y - dy;  p[2] = z - dz;  p += 3;
  p[0] = x - dx;  p[1] = y - dy;  p[2] = z - dz;  p += 3;
  //top
  p[0] = x - dx;  p[1] = y + dy;  p[2] = z + dz;  p += 3;
  p[0] = x + dx;  p[1] = y + dy;  p[2] = z + dz;  p += 3;
  p[0] = x + dx;  p[1] = y - dy;  p[2] = z + dz;  p += 3;
  p[0] = x - dx;  p[1] = y - dy;  p[2] = z + dz;
}

/**************************************************************************/

void FrameBox::SetFrameColor(Color_t ci)
{
  fFrameColor = ci;
  ColorFromIdx(ci, fFrameRGBA, kTRUE);
}

void FrameBox::SetFrameColor(Pixel_t pix)
{
  SetFrameColor(Color_t(TColor::GetColor(pix)));
}

void FrameBox::SetFrameColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  fFrameColor = Color_t(TColor::GetColor(r, g, b));
  fFrameRGBA[0] = r;
  fFrameRGBA[1] = g;
  fFrameRGBA[2] = b;
  fFrameRGBA[3] = a;
}

/**************************************************************************/

void FrameBox::SetBackColor(Color_t ci)
{
  fBackColor = ci;
  ColorFromIdx(ci, fBackRGBA, kTRUE);
}

void FrameBox::SetBackColor(Pixel_t pix)
{
  SetBackColor(Color_t(TColor::GetColor(pix)));
}

void FrameBox::SetBackColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  fBackColor = Color_t(TColor::GetColor(r, g, b));
  fBackRGBA[0] = r;
  fBackRGBA[1] = g;
  fBackRGBA[2] = b;
  fBackRGBA[3] = a;
}
