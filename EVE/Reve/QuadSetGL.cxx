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

#include <TMath.h>

#include "QuadSetGL.h"
#include <Reve/FrameBoxGL.h>

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>
#include <TGLIncludes.h>

using namespace Reve;

//______________________________________________________________________
// OldQuadSetGL
//

ClassImp(OldQuadSetGL)

/**************************************************************************/

OldQuadSetGL::OldQuadSetGL() : TGLObject()
{
  // fDLCache = false; // Disable DL.
}

OldQuadSetGL::~OldQuadSetGL()
{}

/**************************************************************************/

Bool_t OldQuadSetGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  return SetModelCheckClass(obj, Reve::OldQuadSet::Class());
}

void OldQuadSetGL::SetBBox()
{
  SetAxisAlignedBBox(((OldQuadSet*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void OldQuadSetGL::DirectDraw(TGLRnrCtx & /*rnrCtx*/) const
{
  // printf("OldQuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", rnrCtx.Style(), rnrCtx.LOD());

  OldQuadSet& Q = * (OldQuadSet*) fExternalObj;

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT, GL_LINE);
  glPolygonMode(GL_BACK,  GL_LINE);
  glDisable(GL_CULL_FACE);

  Float_t c[4]; glGetFloatv(GL_CURRENT_COLOR, c);
  //  UChar_t alpha = (UChar_t)(255*c[3]);

  glBegin(GL_QUADS);
  for(std::vector<Quad>::iterator q=Q.fQuads.begin(); q!=Q.fQuads.end(); ++q) {
    UChar_t* c = (UChar_t*) &q->color;
    //glColor4ub(c[0], c[1], c[2], (c[3]*alpha) >> 8);
    glColor3ub(c[0], c[1], c[2]);
    glVertex3fv(q->vertices);
    glVertex3fv(q->vertices + 3);
    glVertex3fv(q->vertices + 6);
    glVertex3fv(q->vertices + 9);
  }
  glEnd();

  glPopAttrib();
}


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// QuadSetGL
//

ClassImp(QuadSetGL)

/**************************************************************************/

QuadSetGL::QuadSetGL() : TGLObject(), fM(0)
{
  // fDLCache = false; // Disable DL.
}

QuadSetGL::~QuadSetGL()
{}

/**************************************************************************/

Bool_t QuadSetGL::ShouldDLCache(const TGLRnrCtx & rnrCtx) const
{
  if (rnrCtx.DrawPass() == TGLRnrCtx::kPassOutlineLine)
    return kFALSE;
  return TGLObject::ShouldDLCache(rnrCtx);
}

/**************************************************************************/

Bool_t QuadSetGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  Bool_t ok = SetModelCheckClass(obj, Reve::QuadSet::Class());
  fM = ok ? dynamic_cast<Reve::QuadSet*>(obj) : 0;
  return ok;
}

void QuadSetGL::SetBBox()
{
  SetAxisAlignedBBox(fM->AssertBBox());
}

/**************************************************************************/

inline Bool_t QuadSetGL::SetupColor(const DigitSet::DigitBase& q) const
{
  if (fM->fValueIsColor)
  {
    glColor4ubv((UChar_t*) & q.fValue);
    return kTRUE;
  }
  else
  {
    UChar_t c[4];
    Bool_t visible = fM->fPalette->ColorFromValue(q.fValue, fM->fDefaultValue, c);
    if (visible)
      glColor4ubv(c);
    return visible;
  }
}

/**************************************************************************/

void QuadSetGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  static const Exc_t eH("QuadSetGL::DirectDraw ");

  // printf("QuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", rnrCtx.Style(), rnrCtx.LOD());

  if (rnrCtx.DrawPass() == TGLRnrCtx::kPassOutlineLine)
    return;

  QuadSet& mQ = * fM;

  if (mQ.fFrame != 0 && ! rnrCtx.SecSelection())
    FrameBoxGL::Render(mQ.fFrame);

  if (mQ.fPlex.Size() == 0)
    return;
  if ( ! mQ.fValueIsColor && mQ.fPalette == 0)
  {
    mQ.AssertPalette();
  }

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);

  if (mQ.fRenderMode == DigitSet::RM_Fill)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  else if (mQ.fRenderMode == DigitSet::RM_Line)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  if (mQ.fDisableLigting)  glDisable(GL_LIGHTING);

  if (mQ.fQuadType < QuadSet::QT_Rectangle_End)      RenderQuads(rnrCtx);
  else if (mQ.fQuadType < QuadSet::QT_Line_End)      RenderLines(rnrCtx);
  else if (mQ.fQuadType < QuadSet::QT_Hexagon_End)   RenderHexagons(rnrCtx);

  glPopAttrib();

}


void QuadSetGL::RenderQuads(TGLRnrCtx & rnrCtx) const
{
  static const Exc_t eH("QuadSetGL::RenderQuads ");

  QuadSet& mQ = * fM;

  GLenum primitiveType;
  if (mQ.fRenderMode != DigitSet::RM_Line)
  {
    primitiveType = GL_QUADS;
    if (mQ.fQuadType == QuadSet::QT_FreeQuad)
      glEnable(GL_NORMALIZE);
    else
      glNormal3f(0, 0, 1);
  } else {
    primitiveType = GL_LINE_LOOP;
  }

  VoidCPlex::iterator qi(mQ.fPlex);

  if (rnrCtx.SecSelection()) glPushName(0);

  switch (mQ.fQuadType)
  {

    case QuadSet::QT_FreeQuad:
    {
      Float_t e1[3], e2[3], normal[3];
      while (qi.next()) {
	QuadSet::QFreeQuad& q = * (QuadSet::QFreeQuad*) qi();
	if (SetupColor(q))
	{
	  Float_t* p = q.fVertices;
	  e1[0] = p[3] - p[0]; e1[1] = p[4] - p[1]; e1[2] = p[5] - p[2];
	  e2[0] = p[6] - p[0]; e2[1] = p[7] - p[1]; e2[2] = p[8] - p[2];
	  TMath::Cross(e1, e2, normal);
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glNormal3fv(normal);
	  glVertex3fv(p);
	  glVertex3fv(p + 3);
	  glVertex3fv(p + 6);
	  glVertex3fv(p + 9);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXY:
    {
      while (qi.next()) {
	QuadSet::QRect& q = * (QuadSet::QRect*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,        q.fB,        q.fC);
	  glVertex3f(q.fA + q.fW, q.fB,        q.fC);
	  glVertex3f(q.fA + q.fW, q.fB + q.fH, q.fC);
	  glVertex3f(q.fA,        q.fB + q.fH, q.fC);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXZ:
    {
      while (qi.next()) {
	QuadSet::QRect& q = * (QuadSet::QRect*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,        q.fC, q.fB);
	  glVertex3f(q.fA + q.fW, q.fC, q.fB);
	  glVertex3f(q.fA + q.fW, q.fC, q.fB + q.fH);
	  glVertex3f(q.fA,        q.fC, q.fB + q.fH);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleYZ:
    {
      while (qi.next()) {
	QuadSet::QRect& q = * (QuadSet::QRect*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fC, q.fA,        q.fB);
	  glVertex3f(q.fC, q.fA + q.fW, q.fB);
	  glVertex3f(q.fC, q.fA + q.fW, q.fB + q.fH);
	  glVertex3f(q.fC, q.fA,        q.fB + q.fH);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXYFixedDim:
    {
      const Float_t& w = mQ.fDefWidth;
      const Float_t& h = mQ.fDefHeight;
      while (qi.next()) {
	QuadSet::QRectFixDim& q = * (QuadSet::QRectFixDim*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,     q.fB,     q.fC);
	  glVertex3f(q.fA + w, q.fB,     q.fC);
	  glVertex3f(q.fA + w, q.fB + h, q.fC);
	  glVertex3f(q.fA,     q.fB + h, q.fC);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXYFixedZ:
    {
      const Float_t& z = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QRectFixC& q = * (QuadSet::QRectFixC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,        q.fB,        z);
	  glVertex3f(q.fA + q.fW, q.fB,        z);
	  glVertex3f(q.fA + q.fW, q.fB + q.fH, z);
	  glVertex3f(q.fA,        q.fB + q.fH, z);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXZFixedY:
    {
      const Float_t& y = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QRectFixC& q = * (QuadSet::QRectFixC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,        y, q.fB);
	  glVertex3f(q.fA + q.fW, y, q.fB);
	  glVertex3f(q.fA + q.fW, y, q.fB + q.fH);
	  glVertex3f(q.fA,        y, q.fB + q.fH);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleYZFixedX:
    {
      const Float_t& x = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QRectFixC& q = * (QuadSet::QRectFixC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(x, q.fA,        q.fB);
	  glVertex3f(x, q.fA + q.fW, q.fB);
	  glVertex3f(x, q.fA + q.fW, q.fB + q.fH);
	  glVertex3f(x, q.fA,        q.fB + q.fH);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXYFixedDimZ:
    {
      const Float_t& z = mQ.fDefCoord;
      const Float_t& w = mQ.fDefWidth;
      const Float_t& h = mQ.fDefHeight;
      while (qi.next()) {
	QuadSet::QRectFixDimC& q = * (QuadSet::QRectFixDimC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,     q.fB,     z);
	  glVertex3f(q.fA + w, q.fB,     z);
	  glVertex3f(q.fA + w, q.fB + h, z);
	  glVertex3f(q.fA,     q.fB + h, z);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXZFixedDimY:
    {
      const Float_t& y = mQ.fDefCoord;
      const Float_t& w = mQ.fDefWidth;
      const Float_t& h = mQ.fDefHeight;
      while (qi.next()) {
	QuadSet::QRectFixDimC& q = * (QuadSet::QRectFixDimC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fA,     y, q.fB);
	  glVertex3f(q.fA + w, y, q.fB);
	  glVertex3f(q.fA + w, y, q.fB + h);
	  glVertex3f(q.fA,     y, q.fB + h);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleYZFixedDimX:
    {
      const Float_t& x = mQ.fDefCoord;
      const Float_t& w = mQ.fDefWidth;
      const Float_t& h = mQ.fDefHeight;
      while (qi.next()) {
	QuadSet::QRectFixDimC& q = * (QuadSet::QRectFixDimC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(x, q.fA,     q.fB);
	  glVertex3f(x, q.fA + w, q.fB);
	  glVertex3f(x, q.fA + w, q.fB + h);
	  glVertex3f(x, q.fA,     q.fB + h);
	  glEnd();
	}
      }
      break;
    }

    default:
      throw(eH + "unsupported quad-type.");

  } // end switch quad-type

  if (rnrCtx.SecSelection()) glPopName();
}


void QuadSetGL::RenderLines(TGLRnrCtx & rnrCtx) const
{
  static const Exc_t eH("QuadSetGL::RenderLines ");

  QuadSet& mQ = * fM;

  VoidCPlex::iterator qi(mQ.fPlex);

  if (rnrCtx.SecSelection()) glPushName(0);

  switch (mQ.fQuadType)
  {

    case QuadSet::QT_LineXYFixedZ:
    {
      const Float_t& z = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QLineFixC& q = * (QuadSet::QLineFixC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(GL_LINES);
	  glVertex3f(q.fA,         q.fB,         z);
	  glVertex3f(q.fA + q.fDx, q.fB + q.fDy, z);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_LineXZFixedY:
    {
      const Float_t& z = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QLineFixC& q = * (QuadSet::QLineFixC*) qi();
	if (SetupColor(q))
	{
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(GL_LINES);
	  glVertex3f(q.fA,         z, q.fB);
	  glVertex3f(q.fA + q.fDx, z, q.fB + q.fDy);
	  glEnd();
	}
      }
      break;
    }

    default:
      throw(eH + "unsupported quad-type.");

  }

  if (rnrCtx.SecSelection()) glPopName();
}

void QuadSetGL::RenderHexagons(TGLRnrCtx & rnrCtx) const
{
  static const Exc_t eH("QuadSetGL::RenderHexagons ");

  const Float_t sqr3hf = 0.5*TMath::Sqrt(3);

  QuadSet& mQ = * fM;

  GLenum primitveType = (mQ.fRenderMode != DigitSet::RM_Line) ?
    GL_POLYGON : GL_LINE_LOOP;

  glNormal3f(0, 0, 1);

  VoidCPlex::iterator qi(mQ.fPlex);

  if (rnrCtx.SecSelection()) glPushName(0);

  switch (mQ.fQuadType)
  {

    case QuadSet::QT_HexagonXY:
    {
      while (qi.next()) {
	QuadSet::QHex& q = * (QuadSet::QHex*) qi();
	if (SetupColor(q))
	{
	  const Float_t rh = q.fR * 0.5;
	  const Float_t rs = q.fR * sqr3hf;
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitveType);
	  glVertex3f( q.fR + q.fA,       q.fB, q.fC);
	  glVertex3f(   rh + q.fA,  rs + q.fB, q.fC);
	  glVertex3f(  -rh + q.fA,  rs + q.fB, q.fC);
	  glVertex3f(-q.fR + q.fA,       q.fB, q.fC);
	  glVertex3f(  -rh + q.fA, -rs + q.fB, q.fC);
	  glVertex3f(   rh + q.fA, -rs + q.fB, q.fC);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_HexagonYX:
    {
      while (qi.next()) {
	QuadSet::QHex& q = * (QuadSet::QHex*) qi();
	if (SetupColor(q))
	{
	  const Float_t rh = q.fR * 0.5;
	  const Float_t rs = q.fR * sqr3hf;
	  if (rnrCtx.SecSelection()) glLoadName(qi.index());
	  glBegin(primitveType);
	  glVertex3f( rs + q.fA,    rh + q.fB, q.fC);
	  glVertex3f(      q.fA,  q.fR + q.fB, q.fC);
	  glVertex3f(-rs + q.fA,    rh + q.fB, q.fC);
	  glVertex3f(-rs + q.fA,   -rh + q.fB, q.fC);
	  glVertex3f(      q.fA, -q.fR + q.fB, q.fC);
	  glVertex3f( rs + q.fA,   -rh + q.fB, q.fC);
	  glEnd();
	}
      }
      break;
    }

    default:
      throw(eH + "unsupported quad-type.");

  } // end switch quad-type

  if (rnrCtx.SecSelection()) glPopName();
}

/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________________
void QuadSetGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  if (rec.GetN() < 2) return;
  fM->DigitSelected(rec.GetItem(1));
}
