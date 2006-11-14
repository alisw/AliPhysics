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
// To become new implementation of QuadSet ... not finished yet.
/**************************************************************************/

class QuadSet : public RenderElement,
		public TNamed,
		public TAtt3D,
		public TAttBBox
{
  friend class QuadSetGL;

  QuadSet(const QuadSet&);            // Not implemented
  QuadSet& operator=(const QuadSet&); // Not implemented

public:
  enum QuadType_e
  { 
    QT_Undef,
    QT_FreeQuad,
    QT_AxisAligned,
    QT_AxisAlignedFixedDim,
    QT_AxisAlignedFixedZ,
    QT_AxisAlignedFixedY,
    QT_AxisAlignedFixedDimZ,
    QT_AxisAlignedFixedDimY,
    // line modes (needed for uniform handling of silicon-strip digits)
    QT_LineFixedZ,
    QT_LineFixedY
  };

  enum RenderMode_e { RM_AsIs, RM_Line, RM_Fill };

protected:
  struct QuadBase
  {
    Int_t fValue;
    // Here could have additional integer (like time, second threshold).

    QuadBase(Int_t v=0) : fValue(v) {}
  };

  struct FreeQuad : public QuadBase
  {
    Float_t fVertices[12];
  };

  struct AAFixDimZQuad : public QuadBase
  {
    Float_t fX, fY;
  };

  struct AAFixDimQuad : public AAFixDimZQuad
  {
    Float_t fZ;
  };

  struct AAFixZQuad : public AAFixDimZQuad
  {
    Float_t fW, fH;
  };

  struct AAQuad : public AAFixDimQuad
  {
    Float_t fW, fH;
  };

  struct LineFixedZ : public AAFixDimZQuad
  {
    Float_t fDx, fDy;
  };

protected:
  QuadType_e        fQuadType;
  Bool_t            fValueIsColor;
  Int_t             fDefaultValue;
  VoidCPlex         fPlex;
  QuadBase*         fLastQuad;     //!

  Float_t           fDefWidth;
  Float_t           fDefHeight;
  Float_t           fDefCoord;

  FrameBox*         fFrame;
  RGBAPalette*      fPalette;
  RenderMode_e      fRenderMode;
  ZTrans            fHMTrans;

  static Int_t SizeofAtom(QuadType_e qt);
  QuadBase*    NewQuad();

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

  void QuadValue(Int_t value);
  void QuadColor(Color_t ci);
  void QuadColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

  // --------------------------------

  // void Test(Int_t nquads);

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option="");

  ClassDef(QuadSet, 1);
};

} // namespace Reve

#endif
