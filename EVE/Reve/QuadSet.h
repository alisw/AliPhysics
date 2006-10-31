// $Header$

#ifndef REVE_QuadSet_H
#define REVE_QuadSet_H

#include <Gtypes.h>
#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/RenderElement.h>
#include <Reve/Reve.h>
#include "ZTrans.h"

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

public:
  enum QuadType_e
    { QT_FreeQuad,
      QT_AxisAligned,
      QT_AxisAlignedFixedZ,
      QT_AxisAlignedFixedDim,
      QT_AxisAlignedFixedDimZ
    };

protected:
  struct QuadBase {
    UChar_t fColor[4];
  };

  struct FreeQuad : public QuadBase {
    Float_t fVertices[12];
  };

  struct AAFixDimZQuad : public QuadBase {
    Float_t fX, fY;
  };

  struct AAFixDimQuad : public AAFixDimZQuad {
    Float_t fZ;
  };

  struct AAFixZQuad : public AAFixDimZQuad {
    Float_t fW, fH;
  };

  struct AAQuad : public AAFixDimQuad {
    Float_t fW, fH;
  };

protected:
  QuadType_e        fQuadType;
  Int_t             fSizeOf;
  Int_t             fReserveStep;
  Int_t             fLastEntry;
  Int_t             fNumReserved;

  // Missing:
  // * some actual container
  // * RGBAPalette
  // * user specifies a value instead of a color

  Color_t           fDefaultColor;
  Float_t           fDefWidth;
  Float_t           fDefHeight;
  Float_t           fDefZ;

  ZTrans            fHMTrans;

public:
  QuadSet(const Text_t* n="QuadSet", const Text_t* t="");
  virtual ~QuadSet() {}

  ZTrans& RefHMTrans() { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }

  /*
  void Test(Int_t nquads);

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option = "");
  */

  ClassDef(QuadSet, 1);
};

} // namespace Reve

#endif
