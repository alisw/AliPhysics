// $Header$

#ifndef REVE_BoxSet_H
#define REVE_BoxSet_H

#include <Gtypes.h>
#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>
#include <Reve/RenderElement.h>
#include <Reve/Reve.h>
#include "ZTrans.h"

#include <vector>

class TGeoMatrix;
class TRandom;

namespace Reve {

struct Box
{
  Float_t  vertices[24];
  UChar_t  color[4];

  Box(Color_t col = 1);
  Box(Color_t col, Float_t* p);
  Box(Color_t col, Float_t  x, Float_t  y, Float_t  z,
                   Float_t dx, Float_t dy, Float_t dz);

  Box(TRandom& rnd, Float_t origin, Float_t size);

  virtual ~Box() {}

  void MakeAxisAlignedBox(Float_t  x, Float_t  y, Float_t  z,
			  Float_t dx, Float_t dy, Float_t dz);

  ClassDef(Box, 1);
};

/**************************************************************************/

class BoxSet: public RenderElement,
              public TNamed,
              public TAtt3D,
              public TAttBBox
{
  friend class BoxSetGL;

public:
  enum RenderMode_e { RM_AsIs, RM_Line, RM_Fill };

protected:
  Color_t           fDefaultColor;
  RenderMode_e      fRenderMode;
  ZTrans            fHMTrans;

public:
  std::vector<Box>  fBoxes;

  BoxSet(const Text_t* n="BoxSet", const Text_t* t="");
  virtual ~BoxSet() {}

  void AddBox(const Box& b) { fBoxes.push_back(b); }
  void ClearSet() { fBoxes.clear(); }

  virtual Bool_t CanEditMainColor() { return kTRUE; }

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option = "");

  RenderMode_e  GetRenderMode() const { return fRenderMode; }
  void SetRenderMode(RenderMode_e rm) { fRenderMode = rm; }

  ZTrans& RefHMTrans() { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }

  void Test(Int_t nboxes);

  ClassDef(BoxSet, 1);
}; // endclass BoxSet

}

#endif
