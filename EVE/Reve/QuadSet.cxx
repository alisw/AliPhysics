// $Header$

#include "QuadSet.h"

#include <TColor.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TGeometry.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TROOT.h>
#include <TRandom.h>


using namespace Reve;

/**************************************************************************/
// Quad
/**************************************************************************/
ClassImp(Reve::Quad)

void Quad::ColorFromIdx(Color_t ci)
{
  TColor* c = gROOT->GetColor(ci);
  if(c) {
    UChar_t *x = (UChar_t*) &color;
    x[0] = (UChar_t)(255*c->GetRed());  x[1] = (UChar_t)(255*c->GetGreen());
    x[2] = (UChar_t)(255*c->GetBlue()); x[3] = 255;
  }
}

Quad::Quad(TRandom& rnd, Float_t origin, Float_t size) : color(0)
{
  ColorFromIdx(Int_t(30*rnd.Rndm()));
  Float_t x = 2*origin*(rnd.Rndm() - 0.5);
  Float_t y = 2*origin*(rnd.Rndm() - 0.5);
  Float_t z = 2*origin*(rnd.Rndm() - 0.5);
  Float_t* p = vertices;
  for(int i=0; i<4; ++i) {
    p[0] = x + 2*size*(rnd.Rndm() - 0.5);
    p[1] = y + 2*size*(rnd.Rndm() - 0.5);
    p[2] = z + 2*size*(rnd.Rndm() - 0.5);
    p += 3;
  }
}

/**************************************************************************/
// OldQuadSet
/**************************************************************************/
ClassImp(Reve::OldQuadSet)


OldQuadSet::OldQuadSet(const Text_t* n, const Text_t* t) :
  TNamed(n, t),
  fQuads(),
  fTrans(false)
{}

void OldQuadSet::Test(Int_t nquads)
{
  TRandom rnd(0);
  fQuads.resize(nquads);
  for(Int_t i=0; i<nquads; ++i) {
    new (&fQuads[i]) Quad(rnd, 10, 2);
  }
}

void OldQuadSet::Paint(Option_t* )
{
  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = 1;
  buffer.fTransparency = 0;
  buffer.fLocalFrame   = fTrans; 
  if (fTrans)
    memcpy(buffer.fLocalMaster, fMatrix, 16*sizeof(Double_t));
  buffer.SetSectionsValid(TBuffer3D::kCore);
   
  // We fill kCore on first pass and try with viewer
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    // printf("OldQuadSet::Paint viewer was happy with Core buff3d.\n");
    return;
  }
   
  if (reqSections & TBuffer3D::kRawSizes) {
    Int_t nbPnts = fQuads.size()*4;
    Int_t nbSegs = nbPnts;
    if (!buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, fQuads.size(), fQuads.size()*6)) {
      return;
    }
    buffer.SetSectionsValid(TBuffer3D::kRawSizes); 
  }

  if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
    // Points
    Int_t pidx = 0;
    for (std::vector<Quad>::iterator i=fQuads.begin(); i!=fQuads.end(); ++i) {
      for (Int_t k = 0; k < 12; k++ ){
	buffer.fPnts[pidx] = (*i).vertices[k]; 
	pidx++;
      }
    }

    // Segments
    Int_t sidx = 0;
    for (Int_t q = 0; q < (Int_t)fQuads.size(); ++q) {
      for (Int_t s = 0; s < 4; ++s ) {
	buffer.fSegs[3*sidx ] = 4; 
	buffer.fSegs[3*sidx+1] = sidx;
        if (s == 3)
	  buffer.fSegs[3*sidx+2] = q*4;
	else
	  buffer.fSegs[3*sidx+2] = sidx + 1;
        sidx ++;
      }
    }

    // Polygons
    for (Int_t q = 0; q < (Int_t)fQuads.size(); ++q) {
      buffer.fPols[6*q] = fQuads[q].color;   
      buffer.fPols[6*q +1] = 4;
      buffer.fPols[6*q +2] = 4*q +0;
      buffer.fPols[6*q +3] = 4*q +1;
      buffer.fPols[6*q +4] = 4*q +2;
      buffer.fPols[6*q +5] = 4*q +3;
    }

    buffer.SetSectionsValid(TBuffer3D::kRaw);
    buffer.fColor = 5;
  }
   
  gPad->GetViewer3D()->AddObject(buffer);
}

/**************************************************************************/

void OldQuadSet::ComputeBBox()
{
  if(fQuads.empty()) {
    BBoxZero();
    return;
  }
  BBoxInit();
  for(std::vector<Quad>::iterator q=fQuads.begin(); q!=fQuads.end(); ++q) {
    Float_t* p = q->vertices;
    for(int i=0; i<4; ++i, p+=3)
      BBoxCheckPoint(p);
  }

  // printf("%s BBox is x(%f,%f), y(%f,%f), z(%f,%f)\n", GetName(),
  //        fBBox[0], fBBox[1], fBBox[2], fBBox[3], fBBox[4], fBBox[5]);
}
