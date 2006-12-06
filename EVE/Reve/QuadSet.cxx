// $Header$

#include "QuadSet.h"
#include "RGBAPalette.h"

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
      buffer.fPols[6*q +2] = 4*q + 0;
      buffer.fPols[6*q +3] = 4*q + 1;
      buffer.fPols[6*q +4] = 4*q + 2;
      buffer.fPols[6*q +5] = 4*q + 3;
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

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

//__________________________________________________________________________
// QuadSet
//
// Supports various internal formats that result in rendering of a
// set of rectangular objects.
//
// Names of internal structures and their variables use fixed
// assignment to x, z, y coordinates; the render types can override
// this convention and impose Y as a fixed coordinate.
// For quad modes the deltas are expected to be positive.
// For line modes negative deltas are ok.

ClassImp(Reve::QuadSet)

QuadSet::QuadSet(const Text_t* n, const Text_t* t) :
  RenderElement(),
  TNamed(n, t),

  fQuadType(QT_Undef),
  fValueIsColor(kFALSE),
  fDefaultValue(kMinInt),
  fPlex(),
  fLastQuad(0),

  fDefWidth(1), fDefHeight(1), fDefCoord(0),

  fFrame  (0),
  fPalette(0),
  fRenderMode(RM_Fill),
  fDisableLigting(kTRUE),
  fHMTrans()
{}

QuadSet::QuadSet(QuadType_e quadType, Bool_t valIsCol, Int_t chunkSize,
		 const Text_t* n, const Text_t* t) :
  RenderElement(),
  TNamed(n, t),

  fQuadType(quadType),
  fValueIsColor(valIsCol),
  fDefaultValue(valIsCol ? 0 : kMinInt),
  fPlex(SizeofAtom(quadType), chunkSize),
  fLastQuad(0),

  fDefWidth(1), fDefHeight(1), fDefCoord(0),

  fFrame  (0),
  fPalette(0),
  fRenderMode(RM_Fill),
  fDisableLigting(kTRUE),
  fHMTrans()
{}

QuadSet::~QuadSet()
{
  SetFrame(0);
  SetPalette(0);
}

/**************************************************************************/

Int_t QuadSet::SizeofAtom(QuadSet::QuadType_e qt)
{
  static const Exc_t eH("QuadSet::SizeofAtom ");

  switch (qt) {
    case QT_Undef:                return 0;
    case QT_FreeQuad:             return sizeof(QFreeQuad);
    case QT_RectangleXY:          return sizeof(QRect);
    case QT_RectangleXYFixedDim:  return sizeof(QRectFixDim);
    case QT_RectangleXZFixedY:
    case QT_RectangleXYFixedZ:    return sizeof(QRectFixC);
    case QT_RectangleXZFixedDimY:
    case QT_RectangleXYFixedDimZ: return sizeof(QRectFixDimC);
    case QT_LineXZFixedY:
    case QT_LineXYFixedZ:         return sizeof(QLineFixC);
    case QT_HexagonXY:
    case QT_HexagonYX:            return sizeof(QHex);
    default:                      throw(eH + "unexpected atom type.");
  }
  return 0;
}

/**************************************************************************/

void QuadSet::Reset(QuadSet::QuadType_e quadType, Bool_t valIsCol, Int_t chunkSize)
{
  fQuadType     = quadType;
  fValueIsColor = valIsCol;
  fDefaultValue = valIsCol ? 0 : kMinInt;
  fPlex.Reset(SizeofAtom(fQuadType), chunkSize);
}

void QuadSet::RefitPlex()
{
  // Instruct underlying memory allocator to regroup itself into a
  // contiguous memory chunk.

  fPlex.Refit();
}

/**************************************************************************/

void QuadSet::ScanMinMaxValues(Int_t& min, Int_t& max)
{
  if (fValueIsColor || fPlex.Size() == 0) return;
  min = kMaxInt;
  max = kMinInt;
  for (Int_t c=0; c<fPlex.VecSize(); ++c)
  {
    Char_t* a = fPlex.Chunk(c);
    Int_t   n = fPlex.NAtoms(c);
    while (n--)
    {
      Int_t v = ((QuadBase*)a)->fValue;
      if (v < min) min = v;
      if (v > max) max = v;
      a += fPlex.S();
    }
  }
  if (min == max)
    --min;
}

/**************************************************************************/

void QuadSet::SetFrame(FrameBox* b)
{
  if (fFrame == b) return;
  if (fFrame) fFrame->DecRefCount();
  fFrame = b;
  if (fFrame) {
    fFrame->IncRefCount();
    SetMainColorPtr(fFrame->PtrFrameColor());
  } else {
    SetMainColorPtr(0);
  }
}

void QuadSet::SetPalette(RGBAPalette* p)
{
  if (fPalette == p) return;
  if (fPalette) fPalette->DecRefCount();
  fPalette = p;
  if (fPalette) fPalette->IncRefCount();
}

RGBAPalette* QuadSet::AssertPalette()
{
  if (fPalette == 0) {
    fPalette = new RGBAPalette;
    if (!fValueIsColor) {
      Int_t min, max;
      ScanMinMaxValues(min, max);
      fPalette->SetLimits(min, max);
      fPalette->SetMinMax(min, max);
    }
  }
  return fPalette;
}

/**************************************************************************/

QuadSet::QuadBase* QuadSet::NewQuad()
{
  fLastQuad = new (fPlex.NewAtom()) QuadBase(fDefaultValue);
  return fLastQuad;
}

void QuadSet::AddQuad(Float_t* verts)
{
  static const Exc_t eH("QuadSet::AddQuad ");

  if (fQuadType != QT_FreeQuad)
    throw(eH + "expect free quad-type.");

  QFreeQuad* fq = (QFreeQuad*) NewQuad();
  memcpy(fq->fVertices, verts, sizeof(fq->fVertices));
}

void QuadSet::AddQuad(Float_t x, Float_t y)
{
  AddQuad(x, y, fDefCoord, fDefWidth, fDefHeight);
}

void QuadSet::AddQuad(Float_t x, Float_t y, Float_t z)
{
  AddQuad(x, y, z, fDefWidth, fDefHeight);
}

void QuadSet::AddQuad(Float_t x, Float_t y, Float_t w, Float_t h)
{
  AddQuad(x, y, fDefCoord, w, h);
}

void QuadSet::AddQuad(Float_t x, Float_t y, Float_t z, Float_t w, Float_t h)
{
  static const Exc_t eH("QuadSet::AddAAQuad ");

  QOrigin& fq = * (QOrigin*) NewQuad();
  fq.fX = x; fq.fY = y;
  switch (fQuadType)
  {
    case QT_RectangleXY: {
      QRect& q = (QRect&) fq;
      q.fZ = z; q.fW = w; q.fH = h;
      break;
    }
    case QT_RectangleXYFixedDim: {
      QRectFixDim& q =  (QRectFixDim&) fq;
      q.fZ = z;
      break;
    }
    case QT_RectangleXZFixedY:
    case QT_RectangleXYFixedZ: {
      QRectFixC& q = (QRectFixC&) fq;
      q.fW = w; q.fH = h;
      break;
    }
    case QT_RectangleXZFixedDimY:
    case QT_RectangleXYFixedDimZ: {
      break;
    }
    default:
      throw(eH + "expect axis-aligned quad-type.");
  }
}

void QuadSet::AddLine(Float_t x, Float_t y, Float_t w, Float_t h)
{
  static const Exc_t eH("QuadSet::AddLine ");

  QOrigin& fq = * (QOrigin*) NewQuad();
  fq.fX = x; fq.fY = y;
  switch (fQuadType)
  {
    case QT_LineXZFixedY:
    case QT_LineXYFixedZ: {
      QLineFixC& q = (QLineFixC&) fq;
      q.fDx = w; q.fDy = h;
      break;
    }
    default:
      throw(eH + "expect line quad-type.");
  }
}

void QuadSet::AddHexagon(Float_t x, Float_t y, Float_t z, Float_t r)
{
  static const Exc_t eH("QuadSet::AddHexagon ");

  QOrigin& fq = * (QOrigin*) NewQuad();
  fq.fX = x; fq.fY = y;
  switch (fQuadType)
  {
    case QT_HexagonXY:
    case QT_HexagonYX: {
      QHex& q = (QHex&) fq;
      q.fZ = z; q.fR = r;
      break;
    }
    default:
      throw(eH + "expect line quad-type.");
  }
}

/**************************************************************************/

void QuadSet::QuadValue(Int_t value)
{
  fLastQuad->fValue = value;
}

void QuadSet::QuadColor(Color_t ci)
{
  ColorFromIdx(ci, (UChar_t*) & fLastQuad->fValue, kTRUE);
}

void QuadSet::QuadColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
  UChar_t* x = (UChar_t*) & fLastQuad->fValue;
  x[0] = r; x[1] = g; x[2] = b; x[3] = a;
}

/**************************************************************************/
/**************************************************************************/

void QuadSet::ComputeBBox()
{
  static const Exc_t eH("QuadSet::ComputeBBox ");

  // !!!! Missing handling of FrameBox !!!!
  // It shoud even simpify things ...

  if(fPlex.Size() == 0) {
    BBoxZero();
    return;
  }

  BBoxInit();
  if (fQuadType == QT_RectangleXYFixedZ    ||
      fQuadType == QT_RectangleXYFixedDimZ)
  {
    fBBox[4] = fDefCoord;
    fBBox[5] = fDefCoord;
  }
  else if (fQuadType == QT_RectangleXZFixedY    ||
	   fQuadType == QT_RectangleXZFixedDimY)
  {
    fBBox[2] = fDefCoord;
    fBBox[3] = fDefCoord;
  }

  for (Int_t c=0; c<fPlex.VecSize(); ++c)
  {
    QuadBase* qbp = (QuadBase*) fPlex.Chunk(c);
    Int_t       n = fPlex.NAtoms(c);

    switch (fQuadType)
    {

      case QT_FreeQuad:
      {
	QFreeQuad* qp = (QFreeQuad*) qbp;
	while (n--) {
	  Float_t* p = qp->fVertices;
	  BBoxCheckPoint(p); p += 3;
	  BBoxCheckPoint(p); p += 3;
	  BBoxCheckPoint(p); p += 3;
	  BBoxCheckPoint(p);
	  ++qp;
	}
	break;
      }

      case QT_RectangleXY:
      {
	QRect* qp = (QRect*) qbp;
	while (n--) {
	  QRect& q = * qp;
	  if(q.fX        < fBBox[0]) fBBox[0] = q.fX;
	  if(q.fX + q.fW > fBBox[1]) fBBox[1] = q.fX + q.fW;
	  if(q.fY        < fBBox[2]) fBBox[2] = q.fY;
	  if(q.fY + q.fH > fBBox[3]) fBBox[3] = q.fY + q.fH;
	  if(q.fZ        < fBBox[4]) fBBox[4] = q.fZ;
	  if(q.fZ        > fBBox[5]) fBBox[5] = q.fZ;
	  ++qp;
	}
	break;
      }

      case QT_RectangleXYFixedDim:
      {
	QRectFixDim* qp =  (QRectFixDim*) qbp;
	const Float_t& w = fDefWidth;
	const Float_t& h = fDefHeight;
	while (n--) {
	  QRectFixDim& q = * qp;
	  if(q.fX     < fBBox[0]) fBBox[0] = q.fX;
	  if(q.fX + w > fBBox[1]) fBBox[1] = q.fX + w;
	  if(q.fY     < fBBox[2]) fBBox[2] = q.fY;
	  if(q.fY + h > fBBox[3]) fBBox[3] = q.fY + h;
	  if(q.fZ     < fBBox[4]) fBBox[4] = q.fZ;
	  if(q.fZ     > fBBox[5]) fBBox[5] = q.fZ;
	  ++qp;
	}
	break;
      }

      case QT_RectangleXYFixedZ:
      {
	QRectFixC* qp = (QRectFixC*) qbp;
	while (n--) {
	  QRectFixC& q = * qp;
	  if(q.fX        < fBBox[0]) fBBox[0] = q.fX;
	  if(q.fX + q.fW > fBBox[1]) fBBox[1] = q.fX + q.fW;
	  if(q.fY        < fBBox[2]) fBBox[2] = q.fY;
	  if(q.fY + q.fH > fBBox[3]) fBBox[3] = q.fY + q.fH;
	  ++qp;
	}
	break;
      }

      case QT_RectangleXZFixedY:
      {
	QRectFixC* qp = (QRectFixC*) qbp;
	while (n--) {
	  QRectFixC& q = * qp;
	  if(q.fX        < fBBox[0]) fBBox[0] = q.fX;
	  if(q.fX + q.fW > fBBox[1]) fBBox[1] = q.fX + q.fW;
	  if(q.fY        < fBBox[4]) fBBox[4] = q.fY;
	  if(q.fY + q.fH > fBBox[5]) fBBox[5] = q.fY + q.fH;
	  ++qp;
	}
	break;
      }

      case QT_RectangleXYFixedDimZ:
      {
	QRectFixDimC* qp =  (QRectFixDimC*) qbp;
	const Float_t& w = fDefWidth;
	const Float_t& h = fDefHeight;
	while (n--) {
	  QRectFixDimC& q = * qp;
	  if(q.fX     < fBBox[0]) fBBox[0] = q.fX;
	  if(q.fX + w > fBBox[1]) fBBox[1] = q.fX + w;
	  if(q.fY     < fBBox[2]) fBBox[2] = q.fY;
	  if(q.fY + h > fBBox[3]) fBBox[3] = q.fY + h;
	  ++qp;
	}
	break;
      }

      case QT_RectangleXZFixedDimY:
      {
	QRectFixDimC* qp =  (QRectFixDimC*) qbp;
	const Float_t& w = fDefWidth;
	const Float_t& h = fDefHeight;
	while (n--) {
	  QRectFixDimC& q = * qp;
	  if(q.fX     < fBBox[0]) fBBox[0] = q.fX;
	  if(q.fX + w > fBBox[1]) fBBox[1] = q.fX + w;
	  if(q.fY     < fBBox[4]) fBBox[4] = q.fY;
	  if(q.fY + h > fBBox[5]) fBBox[5] = q.fY + h;
	  ++qp;
	}
	break;
      }

      case QT_LineXYFixedZ:
      {
	QLineFixC* qp =  (QLineFixC*) qbp;
	while (n--) {
	  QLineFixC& q = * qp;
	  BBoxCheckPoint(q.fX,         q.fY,         fDefCoord);
	  BBoxCheckPoint(q.fX + q.fDx, q.fY + q.fDy, fDefCoord);
	  ++qp;
	}
	break;
      }

      case QT_LineXZFixedY:
      {
	QLineFixC* qp =  (QLineFixC*) qbp;
	while (n--) {
	  QLineFixC& q = * qp;
	  BBoxCheckPoint(q.fX,         fDefCoord, q.fY);
	  BBoxCheckPoint(q.fX + q.fDx, fDefCoord, q.fY + q.fDy);
	  ++qp;
	}
	break;
      }

      case QT_HexagonXY:
      case QT_HexagonYX:
      {
	QHex* qp =  (QHex*) qbp;
	while (n--) {
	  QHex& q = * qp;
	  BBoxCheckPoint(q.fX-q.fR, q.fY-q.fR, q.fZ);
	  BBoxCheckPoint(q.fX+q.fR, q.fY+q.fR, q.fZ);
	  ++qp;
	}
	break;
      }

      default: {
	throw(eH + "unsupported quad-type.");
      }

    } // end switch quad-type

  } // end for chunk

  { // Resize bounding box so that it does not have 0 volume.
    // This should be done in TAttBBox (via method AssertMinExtents(epsilon)).
    // Or handled more gracefully in TGLViewer.
    static const Float_t eps = 1e-3;
    for (Int_t i=0; i<6; i+=2) {
      if (fBBox[i+1] - fBBox[i] < eps) {
	Float_t b = 0.5*(fBBox[i] + fBBox[i+1]);
	fBBox[i]   = b - 0.5*eps;
	fBBox[i+1] = b + 0.5*eps;
      }
    }
  }
}

/**************************************************************************/

void QuadSet::Paint(Option_t* /*option*/)
{
  static const Exc_t eH("QuadSet::Paint ");

  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = 1;
  buff.fTransparency = 0;
  fHMTrans.SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections != TBuffer3D::kNone)
    Error(eH, "only direct GL rendering supported.");
}

/**************************************************************************/
