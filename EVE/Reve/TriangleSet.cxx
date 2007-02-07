/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */


#include "TriangleSet.h"
#include "RGBAPalette.h"

#include <TMath.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>

using namespace Reve;

//______________________________________________________________________
// TriangleSet
//
// Made from a list of vertices and a list of triangles (triplets of
// vertex indices).
//
// If input is composed from triangles with direct vertex coordinates
// one should consider finding all occurences of the same vertex
// and specifying it only once.

ClassImp(TriangleSet)

TriangleSet::TriangleSet(Int_t nv, Int_t nt, Bool_t norms, Bool_t cols) :
  RenderElement(fColor),
  TNamed("TriangleSet", 0),
  fNVerts  (nv), fVerts(0),
  fNTrings (nt), fTrings(0), fTringNorms(0), fTringCols(0),
  fColor   (2),
  fHMTrans ()
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

void TriangleSet::GenerateZNormalColors(Float_t fac, Int_t min, Int_t max,
					Bool_t interp, Bool_t wrap)
{
  if (fTringCols  == 0)  fTringCols = new UChar_t[3*fNTrings];
  if (fTringNorms == 0)  GenerateTriangleNormals();

  RGBAPalette pal(min, max, interp, wrap);
  UChar_t *C = fTringCols;
  Float_t *N = fTringNorms;
  for(Int_t t=0; t<fNTrings; ++t, C+=3, N+=3)
    {
      Int_t v = TMath::Nint(fac * N[2]);
      pal.ColorFromValue(v, C, kFALSE);
    }
  gPad->Modified(); gPad->Update();
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
  if (f == 0) {
    ::Error("TriangleSet::ReadTrivialFile", Form("file '%s' not found.", file));
    return 0;
  }

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
