// $Header$

#ifndef REVE_BoxSet_H
#define REVE_BoxSet_H

#include <Reve/DigitSet.h>

class TGeoMatrix;
class TRandom;

namespace Reve {

class BoxSet: public DigitSet
{
  friend class BoxSetGL;

  BoxSet(const BoxSet&);            // Not implemented
  BoxSet& operator=(const BoxSet&); // Not implemented

public:
  enum BoxType_e
  {
    BT_Undef,           // unknown-ignored
    BT_FreeBox,         // arbitrary box: specify 8*(x,y,z) box corners
    BT_AABox,           // axis-aligned box: specify (x,y,z) and (w, h, d)
    BT_AABoxFixedDim    // axis-aligned box w/ fixed dimensions: specify (x,y,z)
  };

protected:

  struct BFreeBox       : public DigitBase { Float_t fVertices[24]; };

  struct BOrigin        : public DigitBase { Float_t fA, fB, fC; };

  struct BAABox         : public BOrigin   { Float_t fW, fH, fD; };

  struct BAABoxFixedDim : public BOrigin {};

protected:
  BoxType_e         fBoxType;      // Type of rendered box.

  Float_t           fDefWidth;     // Breadth assigned to first coordinate  (A).
  Float_t           fDefHeight;    // Breadth assigned to second coordinate (B).
  Float_t           fDefDepth;     // Breadth assigned to third coordinate  (C).

  static Int_t SizeofAtom(BoxType_e bt);

public:
  BoxSet(const Text_t* n="BoxSet", const Text_t* t="");
  virtual ~BoxSet() {}

  void Reset(BoxType_e boxType, Bool_t valIsCol, Int_t chunkSize);
  void Reset();

  void AddBox(const Float_t* verts);
  void AddBox(Float_t a, Float_t b, Float_t c, Float_t w, Float_t h, Float_t d);
  void AddBox(Float_t a, Float_t b, Float_t c);

  virtual void ComputeBBox();
  // virtual void Paint(Option_t* option = "");

  void Test(Int_t nboxes);

  Float_t GetDefWidth()  const { return fDefWidth;  }
  Float_t GetDefHeight() const { return fDefHeight; }
  Float_t GetDefDepth()  const { return fDefDepth;  }

  void SetDefWidth(Float_t v)  { fDefWidth  = v ; }
  void SetDefHeight(Float_t v) { fDefHeight = v ; }
  void SetDefDepth(Float_t v)  { fDefDepth  = v ; }

  ClassDef(BoxSet, 1); // Visual class showing a set of boxes.
}; // endclass BoxSet

}

#endif
