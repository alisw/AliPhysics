#ifndef REVE_QuadSet_H
#define REVE_QuadSet_H

#include <Gtypes.h>
#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

#include <vector>

class TRandom;

namespace Reve {

struct Quad
{
  Float_t  vertices[12];
  Int_t    color;

  void ColorFromIdx(Color_t ci);

  Quad(Color_t col = 1)
  { ColorFromIdx(col); }

  Quad(Color_t col, Float_t* p)
  { ColorFromIdx(col); memcpy(vertices, p, 12*sizeof(Float_t)); }

  Quad(TRandom& rnd, Float_t origin, Float_t size);

  Quad(const Quad& org) { memcpy(this, &org, sizeof(Quad)); }

  ClassDef(Quad, 1);
};

/**************************************************************************/

class QuadSet : public TNamed, public TAtt3D, public TAttBBox
{
  friend class QuadSetGL;

  void Init();

protected:
  std::vector<Quad> fQuads;
  Double_t          fMatrix[16];
  Bool_t            fTrans;

public:
  QuadSet(const Text_t* n="QuadSet", const Text_t* t="") : TNamed(n, t)
  { Init(); }

  Bool_t GetTrans() const { return fTrans; }
  void SetTrans(Bool_t t) { fTrans = t; }

  void Test(Int_t nquads);

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option = "");

  ClassDef(QuadSet, 1);
};

} // namespace Reve

#endif
