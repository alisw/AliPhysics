// $Header$

#include "BoxSet.h"
#include <TRandom.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

using namespace Reve;

/**************************************************************************/
// Box
/**************************************************************************/

Box::Box(Color_t col)
{
  Reve::ColorFromIdx(col, color);
}

Box::Box(Color_t col, Float_t* p)
{
  Reve::ColorFromIdx(col, color);
  memcpy(vertices, p, 24*sizeof(Float_t));
}

Box::Box(Color_t col, Float_t  x, Float_t  y, Float_t  z,
                      Float_t dx, Float_t dy, Float_t dz)
{
  Reve::ColorFromIdx(col, color);
  MakeAxisAlignedBox(x, y, z, dx, dy, dz);
}

Box::Box(TRandom& rnd, Float_t origin, Float_t size)
{
  Reve::ColorFromIdx(Int_t(30*rnd.Rndm()), color);
  Float_t  x = 2*origin*(rnd.Rndm() - 0.5);
  Float_t  y = 2*origin*(rnd.Rndm() - 0.5);
  Float_t  z = 2*origin*(rnd.Rndm() - 0.5);
  Float_t dx = 2*size*(rnd.Rndm() - 0.5);
  Float_t dy = 2*size*(rnd.Rndm() - 0.5);
  Float_t dz = 2*size*(rnd.Rndm() - 0.5);
  MakeAxisAlignedBox(x, y, z, dx, dy, dz);
}

void Box::MakeAxisAlignedBox(Float_t  x, Float_t  y, Float_t  z,
			     Float_t dx, Float_t dy, Float_t dz)
{
  Float_t* p = vertices;
  //bottom
  p[0] = x - dx;  p[1] = y + dy, p[2] = z - dz;  p += 3;
  p[0] = x + dx;  p[1] = y + dy, p[2] = z - dz;  p += 3;
  p[0] = x + dx;  p[1] = y - dy, p[2] = z - dz;  p += 3;
  p[0] = x - dx;  p[1] = y - dy, p[2] = z - dz;  p += 3;
  //top
  p[0] = x - dx;  p[1] = y + dy, p[2] = z + dz;  p += 3;
  p[0] = x + dx;  p[1] = y + dy, p[2] = z + dz;  p += 3;
  p[0] = x + dx;  p[1] = y - dy, p[2] = z + dz;  p += 3;
  p[0] = x - dx;  p[1] = y - dy, p[2] = z + dz;
 
}

//______________________________________________________________________
// BoxSet
//

ClassImp(BoxSet)

BoxSet::BoxSet(const Text_t* n, const Text_t* t) :
  RenderElement(fDefaultColor),
  TNamed(n, t),
  fRenderMode (RM_AsIs)
{}

/**************************************************************************/

void BoxSet::ComputeBBox()
{
  if(fBoxes.empty()) {
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
    bbox_zero();
#else
    BBoxZero();
#endif
    return;
  }
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  bbox_init();
#else
  BBoxInit();
#endif
  for(std::vector<Box>::iterator q=fBoxes.begin(); q!=fBoxes.end(); ++q) {
    Float_t* p = q->vertices;
    for(int i=0; i<8; ++i, p+=3)
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
      bbox_check_point(p);
#else
      BBoxCheckPoint(p);
#endif
  }
  // printf("%s BBox is x(%f,%f), y(%f,%f), z(%f,%f)\n", GetName(),
  //        fBBox[0], fBBox[1], fBBox[2], fBBox[3], fBBox[4], fBBox[5]);
}

void BoxSet::Paint(Option_t* /*option*/)
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

/**************************************************************************/

void BoxSet::Test(Int_t nboxes)
{
  TRandom rnd(0);
  fBoxes.resize(nboxes);
  for(Int_t i=0; i<nboxes; ++i) {
    new (&fBoxes[i]) Box(rnd, 10, 2);
  }
}
