// $Header$

#include "TriangleSet.h"

#include <TVector3.h>
#include <TRandom3.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TGeoMatrix.h>

using namespace Reve;


ClassImp(TriangleSet)

TriangleSet::TriangleSet(Int_t nv, Int_t nt, Bool_t norms, Bool_t cols) :
  RenderElement(fColor),
  TNamed("TriangleSet", 0),
  fNVerts  (nv),
  fNTrings (nt),
  fColor   (2)
{
  fVerts  = new Float_t[3*fNVerts];
  fTrings = new Int_t  [3*fNTrings];
  fTringNorms = (norms) ? new Float_t[3*fNTrings] : 0;
  fTringCols  = (cols)  ? new UChar_t[3*fNTrings] : 0;
}

TriangleSet::~TriangleSet()
{
  delete [] fVerts;
  delete [] fTrings;
  delete [] fTringNorms;
  delete [] fTringCols;
}

/**************************************************************************/

void TriangleSet::GenerateTriangleNormals()
{
  if (fTringNorms == 0)  fTringNorms = new Float_t[3*fNTrings];

  TVector3 e1, e2, n;
  Float_t *N = fTringNorms;
  Int_t   *T = fTrings;
  for(Int_t t=0; t<fNTrings; ++t, N+=3, T+=3)
    {
      Float_t* v0 = Vertex(T[0]);
      Float_t* v1 = Vertex(T[1]);
      Float_t* v2 = Vertex(T[2]);
      e1.SetXYZ(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
      e2.SetXYZ(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
      n = e1.Cross(e2);
      n.SetMag(1);
      n.GetXYZ(N);
    }
}

void TriangleSet::GenerateRandomColors()
{
  if (fTringCols == 0)  fTringCols = new UChar_t[3*fNTrings];

  TRandom r;
  r.SetSeed(0);
  UChar_t *C = fTringCols;
  for(Int_t t=0; t<fNTrings; ++t, C+=3)
    {
      C[0] = (UChar_t) r.Uniform(60, 255);
      C[1] = (UChar_t) r.Uniform(60, 255);
      C[2] = (UChar_t) r.Uniform(60, 255);
    }
}

/**************************************************************************/

void TriangleSet::ComputeBBox()
{
  if (fNVerts <= 0) {
    BBoxZero();
    return;
  }

  BBoxInit();
  Float_t* v = fVerts;
  for (Int_t i=0; i<fNVerts; ++i, v += 3)
    BBoxCheckPoint(v);
}

void TriangleSet::Paint(Option_t* )
{
  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = fColor;
  buffer.fTransparency = 0;
  fHMTrans.SetBuffer3D(buffer);
  buffer.SetSectionsValid(TBuffer3D::kCore);
   
  // We fill kCore on first pass and try with viewer
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    return;
  }

  Error("TriangleSet::Paint", "only direct OpenGL rendering supported.");
}

/**************************************************************************/

#include <stdio.h>

TriangleSet* TriangleSet::ReadTrivialFile(const char* file)
{
  FILE* f = fopen(file, "r");

  Int_t nv, nt;
  fscanf(f, "%d %d", &nv, &nt);

  TriangleSet* ts = new TriangleSet(nv, nt);

  Float_t *V = ts->Vertex(0);
  for (Int_t i=0; i<nv; ++i) {
    fscanf(f, "%f %f %f", V++, V++, V++);
  }

  Int_t *T = ts->Triangle(0);
  for (Int_t i=0; i<nt; ++i) {
    fscanf(f, "%d %d %d", T++, T++, T++);
  }

  fclose(f);

  return ts;
}
