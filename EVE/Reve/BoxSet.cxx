// $Header$

#include "BoxSet.h"
#include <TRandom.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

using namespace Reve;


//______________________________________________________________________
// BoxSet
//

ClassImp(BoxSet)

BoxSet::BoxSet(const Text_t* n, const Text_t* t) :
  DigitSet      (n, t),

  fBoxType      (BT_Undef),
  fDefWidth     (1),
  fDefHeight    (1),
  fDefDepth     (1)
{
  // Override from DigitSet.
  fDisableLigting = kFALSE;
}

/**************************************************************************/

Int_t BoxSet::SizeofAtom(BoxSet::BoxType_e bt)
{
  static const Exc_t eH("BoxSet::SizeofAtom ");

  switch (bt) {
    case BT_Undef:                return 0;
    case BT_FreeBox:              return sizeof(BFreeBox);
    case BT_AABox:                return sizeof(BAABox);
    case BT_AABoxFixedDim:        return sizeof(BAABoxFixedDim);
    default:                      throw(eH + "unexpected atom type.");
  }
  return 0;
}

/**************************************************************************/

void BoxSet::Reset(BoxSet::BoxType_e boxType, Bool_t valIsCol, Int_t chunkSize)
{
  fBoxType      = boxType;
  fValueIsColor = valIsCol;
  fDefaultValue = valIsCol ? 0 : kMinInt;
  if (fOwnIds)
    ReleaseIds();
  fPlex.Reset(SizeofAtom(fBoxType), chunkSize);
}

void BoxSet::Reset()
{
  if (fOwnIds)
    ReleaseIds();
  fPlex.Reset(SizeofAtom(fBoxType), TMath::Max(fPlex.N(), 64));
}

/**************************************************************************/

void BoxSet::AddBox(const Float_t* verts)
{
  // Create a new box from a set of 8 vertices.
  // To be used for box-type BT_FreeBox.

  static const Exc_t eH("BoxSet::AddBox ");

  if (fBoxType != BT_FreeBox)
    throw(eH + "expect free box-type.");

  BFreeBox* b = (BFreeBox*) NewDigit();
  memcpy(b->fVertices, verts, sizeof(b->fVertices));
}

void BoxSet::AddBox(Float_t a, Float_t b, Float_t c, Float_t w, Float_t h, Float_t d)
{
  // Create a new axis-aligned box from at a given position and with
  // specified dimensions.
  // To be used for box-type BT_AABox.

  static const Exc_t eH("BoxSet::AddBox ");

  if (fBoxType != BT_AABox)
    throw(eH + "expect axis-aligned box-type.");

  BAABox* box = (BAABox*) NewDigit();
  box->fA = a; box->fB = b; box->fC = c;
  box->fW = w; box->fH = h; box->fD = d;
}

void BoxSet::AddBox(Float_t a, Float_t b, Float_t c)
{
  // Create a new axis-aligned box from at a given position.
  // To be used for box-type BT_AABoxFixedDim.

  static const Exc_t eH("BoxSet::AddBox ");

  if (fBoxType != BT_AABoxFixedDim)
    throw(eH + "expect axis-aligned fixed-dimension box-type.");

  BAABoxFixedDim* box = (BAABoxFixedDim*) NewDigit();
  box->fA = a; box->fB = b; box->fC = c;
}

/**************************************************************************/

void BoxSet::ComputeBBox()
{
  // Fill bounding-box information of the base-class TAttBBox (virtual method).
  // If member 'FrameBox* fFrame' is set, frame's corners are used as bbox.

  static const Exc_t eH("BoxSet::ComputeBBox ");

  if (fFrame != 0)
  {
    BBoxInit();
    Int_t    n    = fFrame->GetFrameSize() / 3;
    Float_t *bbps = fFrame->GetFramePoints();
    for (int i=0; i<n; ++i, bbps+=3)
      BBoxCheckPoint(bbps);
    return;
  }

  if(fPlex.Size() == 0)
  {
    BBoxZero();
    return;
  }

  BBoxInit();

  VoidCPlex::iterator bi(fPlex);
  switch (fBoxType)
  {

    case BT_FreeBox:
    {
      while (bi.next()) {
        BFreeBox& b = * (BFreeBox*) bi();
        Float_t * p = b.fVertices;
        for(int i=0; i<8; ++i, p+=3)
          BBoxCheckPoint(p);
      }
      break;
    }

    case BT_AABox:
    {
      while (bi.next()) {
        BAABox& b = * (BAABox*) bi();
        BBoxCheckPoint(b.fA, b.fB, b.fC);
        BBoxCheckPoint(b.fA + b.fW, b.fB + b.fH , b.fC + b.fD);
      }
      break;
    }

    case BT_AABoxFixedDim:
    {
      while (bi.next()) {
        BAABoxFixedDim& b = * (BAABoxFixedDim*) bi();
        BBoxCheckPoint(b.fA, b.fB, b.fC);
        BBoxCheckPoint(b.fA + fDefWidth, b.fB + fDefHeight , b.fC + fDefDepth);
      }
      break;
    }

    default:
    {
      throw(eH + "unsupported box-type.");
    }

  } // end switch box-type

  printf("%s BBox is x(%f,%f), y(%f,%f), z(%f,%f)\n", GetName(),
         fBBox[0], fBBox[1], fBBox[2], fBBox[3], fBBox[4], fBBox[5]);
}

/*
void BoxSet::Paint(Option_t* option)
{
  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = 1;
  buff.fTransparency = 0;
  fHMTrans.SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections == TBuffer3D::kNone) {
    // printf("BoxSet::Paint viewer was happy with Core buff3d.\n");
    return;
  }

  if (reqSections & TBuffer3D::kRawSizes) {
    Int_t nbPnts = fBoxes.size()*8;
    Int_t nbSegs = fBoxes.size()*12;
    Int_t nbPoly = fBoxes.size()*6;
    if (!buff.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPoly, nbPoly*6)) {
      return;
    }
    buff.SetSectionsValid(TBuffer3D::kRawSizes);
  }

  if ((reqSections & TBuffer3D::kRaw) && buff.SectionsValid(TBuffer3D::kRawSizes)) {
    // Points
    Int_t pidx = 0;
    for (std::vector<Box>::iterator i=fBoxes.begin(); i!=fBoxes.end(); ++i) {
      for (Int_t k=0; k<24; ++k) {
	buff.fPnts[pidx] = (*i).vertices[k];
	pidx++;
      }
    }
    Int_t c = 2;   // color; !!! wrong
    Int_t eoff;    // polygon or segment offset in fPols and fSegs array
    Int_t voff;    // vertex offset
    Int_t soff;    // offset in counting segments
    Int_t nspp = 4; // number of segments per polygon

    // Segments & Polygons
    for (Int_t q = 0; q < (Int_t)fBoxes.size(); ++q) {
      eoff = q*36;
      soff = q*12;
      voff = q*8;
      //bottom
      buff.fSegs[ 0 + eoff] = c;  buff.fSegs[ 1 + eoff] = 0 + voff;  buff.fSegs[ 2 + eoff] = 1 + voff;
      buff.fSegs[ 3 + eoff] = c;  buff.fSegs[ 4 + eoff] = 1 + voff;  buff.fSegs[ 5 + eoff] = 2 + voff;
      buff.fSegs[ 6 + eoff] = c;  buff.fSegs[ 7 + eoff] = 2 + voff;  buff.fSegs[ 8 + eoff] = 3 + voff;
      buff.fSegs[ 9 + eoff] = c;  buff.fSegs[10 + eoff] = 3 + voff;  buff.fSegs[11 + eoff] = 0 + voff;
      // top
      buff.fSegs[12 + eoff] = c;  buff.fSegs[13 + eoff] = 4 + voff;  buff.fSegs[14 + eoff] = 5 + voff;
      buff.fSegs[15 + eoff] = c;  buff.fSegs[16 + eoff] = 5 + voff;  buff.fSegs[17 + eoff] = 6 + voff;
      buff.fSegs[18 + eoff] = c;  buff.fSegs[19 + eoff] = 6 + voff;  buff.fSegs[20 + eoff] = 7 + voff;
      buff.fSegs[21 + eoff] = c;  buff.fSegs[22 + eoff] = 7 + voff;  buff.fSegs[23 + eoff] = 4 + voff;
      //sides
      buff.fSegs[24 + eoff] = c;  buff.fSegs[25 + eoff] = 0 + voff;  buff.fSegs[26 + eoff] = 4 + voff;
      buff.fSegs[27 + eoff] = c;  buff.fSegs[28 + eoff] = 1 + voff;  buff.fSegs[29 + eoff] = 5 + voff;
      buff.fSegs[30 + eoff] = c;  buff.fSegs[31 + eoff] = 2 + voff;  buff.fSegs[32 + eoff] = 6 + voff;
      buff.fSegs[33 + eoff] = c;  buff.fSegs[34 + eoff] = 3 + voff;  buff.fSegs[35 + eoff] = 7 + voff;

      buff.fPols[ 0 + eoff] = c;         buff.fPols[ 1 + eoff] = nspp;
      buff.fPols[ 2 + eoff] = 0 + soff;  buff.fPols[ 3 + eoff] = 9  + soff;
      buff.fPols[ 4 + eoff] = 4 + soff;  buff.fPols[ 5 + eoff] = 8  + soff;
      buff.fPols[ 6 + eoff] = c;         buff.fPols[ 7 + eoff] = nspp;
      buff.fPols[ 8 + eoff] = 1 + soff;  buff.fPols[ 9 + eoff] = 10 + soff;
      buff.fPols[10 + eoff] = 5 + soff;  buff.fPols[11 + eoff] = 9  + soff;
      buff.fPols[12 + eoff] = c;         buff.fPols[13 + eoff] = nspp;
      buff.fPols[14 + eoff] = 2 + soff;  buff.fPols[15 + eoff] = 11 + soff;
      buff.fPols[16 + eoff] = 6 + soff;  buff.fPols[17 + eoff] = 10 + soff;
      buff.fPols[18 + eoff] = c;         buff.fPols[19 + eoff] = nspp;
      buff.fPols[20 + eoff] = 3 + soff;  buff.fPols[21 + eoff] = 8  + soff;
      buff.fPols[22 + eoff] = 7 + soff;  buff.fPols[23 + eoff] = 11 + soff;
      buff.fPols[24 + eoff] = c;         buff.fPols[25 + eoff] = nspp;
      buff.fPols[26 + eoff] = 0 + soff;  buff.fPols[27 + eoff] = 3  + soff;
      buff.fPols[28 + eoff] = 2 + soff;  buff.fPols[29 + eoff] = 1  + soff;
      buff.fPols[30 + eoff] = c;         buff.fPols[31 + eoff] = nspp;
      buff.fPols[32 + eoff] = 4 + soff;  buff.fPols[33 + eoff] = 5  + soff;
      buff.fPols[34 + eoff] = 6 + soff;  buff.fPols[35 + eoff] = 7  + soff;
    }
    buff.fColor = 2; // colors on polygons are ignored
    buff.SetSectionsValid(TBuffer3D::kRaw);
  }
  gPad->GetViewer3D()->AddObject(buff);
}
*/

/**************************************************************************/

void BoxSet::Test(Int_t nboxes)
{
  Reset(BT_AABox, kTRUE, nboxes);
  TRandom rnd(0);
  const Float_t origin = 10, size = 2;
  Int_t color;
  for(Int_t i=0; i<nboxes; ++i)
  {
    AddBox(origin * rnd.Uniform(-1, 1),
           origin * rnd.Uniform(-1, 1),
           origin * rnd.Uniform(-1, 1),
           size   * rnd.Uniform(0.1, 1),
           size   * rnd.Uniform(0.1, 1),
           size   * rnd.Uniform(0.1, 1));

    Reve::ColorFromIdx(rnd.Integer(256), (UChar_t*)&color);
    DigitValue(color);
  }
}
