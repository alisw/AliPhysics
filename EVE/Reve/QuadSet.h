// $Header$

#ifndef REVE_QuadSet_H
#define REVE_QuadSet_H

#include <Gtypes.h>
#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/Reve.h>
#include <Reve/RenderElement.h>
#include <Reve/FrameBox.h>
#include <Reve/RGBAPalette.h>
#include <Reve/Plex.h>
#include <Reve/ZTrans.h>

#include <vector>

class TRandom;

namespace Reve {

struct Quad
{
  Float_t  vertices[12];
  Int_t    color;

  void ColorFromIdx(Color_t ci);

  Quad(Color_t col = 1) : color(0)
  { ColorFromIdx(col); }

  Quad(Color_t col, Float_t* p) : color(0)
  { ColorFromIdx(col); memcpy(vertices, p, 12*sizeof(Float_t)); }

  Quad(TRandom& rnd, Float_t origin, Float_t size);

  Quad(const Quad& org) : color(0) { memcpy(this, &org, sizeof(Quad)); }

  virtual ~Quad() {}

  ClassDef(Quad, 1);
};

class OldQuadSet : public TNamed, public TAtt3D, public TAttBBox
{
  friend class OldQuadSetGL;

protected:
  std::vector<Quad> fQuads;
  Double_t          fMatrix[16];
  Bool_t            fTrans;

public:
  OldQuadSet(const Text_t* n="QuadSet", const Text_t* t="");
  virtual ~OldQuadSet() {}

  Bool_t GetTrans() const { return fTrans; }
  void SetTrans(Bool_t t) { fTrans = t; }

  std::vector<Reve::Quad>& Quads() { return fQuads; }

  void Test(Int_t nquads);

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option = "");

  ClassDef(OldQuadSet, 1);
};

/**************************************************************************/
// To become new implementation of QuadSet ... somewhat finished.
/**************************************************************************/

class QuadSet : public RenderElement,
		public TNamed,
		public TAtt3D,
		public TAttBBox
{
  friend class QuadSetEditor;
  friend class QuadSetGL;

  QuadSet(const QuadSet&);            // Not implemented
  QuadSet& operator=(const QuadSet&); // Not implemented

public:
  enum QuadType_e
  { 
    QT_Undef,                // unknown-ignored
    QT_FreeQuad,             // arbitrary quad: specify 4*(x,y,z) quad corners
    QT_RectangleXY,          // rectangle in x-y plane: specify x, y, z, w, h
    QT_RectangleXYFixedDim,  // rectangle in x-y plane: specify x, y, z; w, h taken from fDefWidth/Height
    QT_RectangleXYFixedZ,    // rectangle in x-y plane: specify x, y, w, h; z taken from fDefCoord
    QT_RectangleXZFixedY,    // rectangle in x-z plane: specify x, z, w, h; y taken from fDefCoord
    QT_RectangleXYFixedDimZ, // rectangle in x-y plane: specify x, y; w, h, z taken from fDefWidth/Height/Coord
    QT_RectangleXZFixedDimY, // rectangle in x-z plane: specify x, z; w, h, y taken from fDefWidth/Height/Coord
    QT_Rectangle_End,
    // line modes (needed for uniform handling of silicon-strip digits)
    QT_LineXYFixedZ,         // line in x-y plane: specify x, y, w(dx), h(dy); z taken from fDefCoord
    QT_LineXZFixedY,         // line in x-z plane: specify x, z, w(dx), h(dz); y taken from fDefCoord
    QT_Line_End,
    // hexagon modes
    QT_HexagonXY,            // horizontal hexagon: specify x, y, z, r
    QT_HexagonYX,            // vertical   hexagon: specify x, y, z, r
    QT_Hexagon_End
    // circle modes:
    // QT_CircleXY,          // specify r, z
    // QT_CircleXYFixedZ,    // specify r
    // QT_CircleXYFixedR,    // specify z
  };

  enum RenderMode_e { RM_AsIs, RM_Line, RM_Fill };

protected:
  struct QuadBase
  {
    Int_t fValue;
    TRef  fId;

    // Here could have additional integer (like time, second threshold).

    QuadBase(Int_t v=0) : fValue(v) {}
  };

  struct QFreeQuad     : public QuadBase      { Float_t fVertices[12]; };

  struct QOrigin       : public QuadBase      { Float_t fX, fY; };

  struct QRectFixDimC  : public QOrigin       { };

  struct QRectFixDim   : public QRectFixDimC  { Float_t fZ; };

  struct QRectFixC     : public QRectFixDimC  { Float_t fW, fH; };

  struct QRect         : public QRectFixDim   { Float_t fW, fH; };

  struct QLineFixC     : public QOrigin       { Float_t fDx, fDy; };

  struct QHex          : public QOrigin       { Float_t fZ, fR; };

protected:
  QuadType_e        fQuadType;
  Int_t             fDefaultValue;
  Bool_t            fValueIsColor;
  Bool_t            fOwnIds;       //Flag specifying if id-objects are owned by the QuadSet
  VoidCPlex         fPlex;
  QuadBase*         fLastQuad;     //!

  Float_t           fDefWidth;
  Float_t           fDefHeight;
  Float_t           fDefCoord;

  FrameBox*         fFrame;
  RGBAPalette*      fPalette;
  RenderMode_e      fRenderMode;
  Bool_t            fDisableLigting;
  ZTrans            fHMTrans;

  static Int_t SizeofAtom(QuadType_e qt);
  QuadBase*    NewQuad();

  void ReleaseIds();

public:
  QuadSet(const Text_t* n="QuadSet", const Text_t* t="");
  QuadSet(QuadType_e quadType, Bool_t valIsCol, Int_t chunkSize,
	  const Text_t* n="QuadSet", const Text_t* t="");
  virtual ~QuadSet();

  virtual Bool_t CanEditMainColor() { return kTRUE; }

  void Reset(QuadType_e quadType, Bool_t valIsCol, Int_t chunkSize);
  void RefitPlex();

  void ScanMinMaxValues(Int_t& min, Int_t& max);

  Float_t GetDefWidth()  const { return fDefWidth;  }
  Float_t GetDefHeight() const { return fDefHeight; }
  Float_t GetDefCoord()  const { return fDefCoord;  }

  void SetDefWidth(Float_t v)  { fDefWidth  = v ; }
  void SetDefHeight(Float_t v) { fDefHeight = v ; }
  void SetDefCoord(Float_t v)  { fDefCoord  = v ; }

  // --------------------------------

  FrameBox* GetFrame() const { return fFrame; }
  void SetFrame(FrameBox* b);

  RGBAPalette* GetPalette() const { return fPalette; }
  void SetPalette(RGBAPalette* p);
  RGBAPalette* AssertPalette();

  RenderMode_e  GetRenderMode() const { return fRenderMode; }
  void SetRenderMode(RenderMode_e rm) { fRenderMode = rm; }

  ZTrans& RefHMTrans() { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }

  // --------------------------------

  void AddQuad(Float_t* verts);

  void AddQuad(Float_t x, Float_t y);
  void AddQuad(Float_t x, Float_t y, Float_t z);
  void AddQuad(Float_t x, Float_t y, Float_t w, Float_t h);
  void AddQuad(Float_t x, Float_t y, Float_t z, Float_t w, Float_t h);

  void AddLine(Float_t x, Float_t y, Float_t w, Float_t h);

  void AddHexagon(Float_t x, Float_t y, Float_t z, Float_t r);

  void QuadValue(Int_t value);
  void QuadColor(Color_t ci);
  void QuadColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

  void QuadId(TObject* id);
  Bool_t GetOwnIds() const    { return fOwnIds; }
  void   SetOwnIds(Bool_t o)  { fOwnIds = o; }

  QuadBase* GetQuad(Int_t n) { return (QuadBase*) fPlex.Atom(n);   }
  TObject*  GetId(Int_t n)   { return GetQuad(n)->fId.GetObject(); }

  virtual void QuadSelected(Int_t idx);

  // --------------------------------

  // void Test(Int_t nquads);

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option="");

  VoidCPlex* GetPlex() { return &fPlex; }

  ClassDef(QuadSet, 1);
};

} // namespace Reve

#endif
