// $Header$

#ifndef REVE_BoxSet_H
#define REVE_BoxSet_H

#include <Gtypes.h>
#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>
#include <Reve/Reve.h>
#include <vector>

class TRandom;

namespace Reve {

struct Box
{
  Float_t  vertices[24];
  UChar_t  color[4];

  Box(Color_t col = 1)
  { Reve::ColorFromIdx(col, color); }
  Box(Color_t col, Float_t* p)
  { Reve::ColorFromIdx(col, color); memcpy(vertices, p, 24*sizeof(Float_t)); }

  Box(TRandom& rnd, Float_t origin, Float_t size);

  virtual ~Box() {}

  ClassDef(Box, 1);
};

/**************************************************************************/

class BoxSet: public TNamed, public TAtt3D, public TAttBBox
{
  friend class BoxSetGL;

protected:
  Double_t          fMatrix[16];
  Bool_t            fTrans;

public:
  std::vector<Box>  fBoxes;

  BoxSet(const Text_t* n="BoxSet", const Text_t* t="");
  virtual ~BoxSet() {}

  void ClearSet() { fBoxes.clear(); }

  Bool_t GetTrans() const   { return fTrans; }
  void   SetTrans(Bool_t t) { fTrans = t; }
  Double_t* ArrTrans()      { return fMatrix; }

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option = "");

  void Test(Int_t nboxes);

  ClassDef(BoxSet, 1);
}; // endclass BoxSet

}

#endif
