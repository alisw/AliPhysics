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

#include <TGLDrawFlags.h>

#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// OldQuadSetGL
//

ClassImp(OldQuadSetGL)

/**************************************************************************/

OldQuadSetGL::OldQuadSetGL() : TGLObject()
{
  // fCached = false; // Disable DL.
}

OldQuadSetGL::~OldQuadSetGL()
{}

/**************************************************************************/

Bool_t OldQuadSetGL::SetModel(TObject* obj)
{
  return SetModelCheckClass(obj, Reve::OldQuadSet::Class());
}

void OldQuadSetGL::SetBBox()
{
  SetAxisAlignedBBox(((OldQuadSet*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void OldQuadSetGL::DirectDraw(const TGLDrawFlags & ) const
{
  // printf("OldQuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

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
  // fCached = false; // Disable DL.
}

QuadSetGL::~QuadSetGL()
{}

/**************************************************************************/

Bool_t QuadSetGL::SetModel(TObject* obj)
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

inline Bool_t QuadSetGL::SetupColor(const QuadSet::QuadBase& q) const
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

void QuadSetGL::DirectDraw(const TGLDrawFlags & flags) const
{
  static const Exc_t eH("QuadSetGL::DirectDraw ");

  // printf("QuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  QuadSet& mQ = * fM;

  if (mQ.fFrame != 0)
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

  if (mQ.fRenderMode == QuadSet::RM_Fill)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  else if (mQ.fRenderMode == QuadSet::RM_Line)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  if (mQ.fDisableLigting)  glDisable(GL_LIGHTING);

  if (mQ.fQuadType < QuadSet::QT_Rectangle_End)      RenderQuads(flags);
  else if (mQ.fQuadType < QuadSet::QT_Line_End)      RenderLines(flags);
  else if (mQ.fQuadType < QuadSet::QT_Hexagon_End)   RenderHexagons(flags);

  glPopAttrib();

}


void QuadSetGL::RenderQuads(const TGLDrawFlags &flags) const
{
  static const Exc_t eH("QuadSetGL::RenderQuads ");

  QuadSet& mQ = * fM;

  GLenum primitiveType;
  if (mQ.fRenderMode != QuadSet::RM_Line)
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

  if (flags.SecSelection()) glPushName(0);

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
	  if (flags.SecSelection()) glLoadName(qi.index());
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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fX,        q.fY,        q.fZ);
	  glVertex3f(q.fX + q.fW, q.fY,        q.fZ);
	  glVertex3f(q.fX + q.fW, q.fY + q.fH, q.fZ);
	  glVertex3f(q.fX,        q.fY + q.fH, q.fZ);
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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fX,     q.fY,     q.fZ);
	  glVertex3f(q.fX + w, q.fY,     q.fZ);
	  glVertex3f(q.fX + w, q.fY + h, q.fZ);
	  glVertex3f(q.fX,     q.fY + h, q.fZ);
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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fX,        q.fY,        z);
	  glVertex3f(q.fX + q.fW, q.fY,        z);
	  glVertex3f(q.fX + q.fW, q.fY + q.fH, z);
	  glVertex3f(q.fX,        q.fY + q.fH, z);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXZFixedY:
    {
      const Float_t& z = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QRectFixC& q = * (QuadSet::QRectFixC*) qi();
	if (SetupColor(q))
	{
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fX,        z, q.fY);
	  glVertex3f(q.fX + q.fW, z, q.fY);
	  glVertex3f(q.fX + q.fW, z, q.fY + q.fH);
	  glVertex3f(q.fX,        z, q.fY + q.fH);
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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fX,     q.fY,     z);
	  glVertex3f(q.fX + w, q.fY,     z);
	  glVertex3f(q.fX + w, q.fY + h, z);
	  glVertex3f(q.fX,     q.fY + h, z);
	  glEnd();
	}
      }
      break;
    }

    case QuadSet::QT_RectangleXZFixedDimY:
    {
      const Float_t& z = mQ.fDefCoord;
      const Float_t& w = mQ.fDefWidth;
      const Float_t& h = mQ.fDefHeight;
      while (qi.next()) {
	QuadSet::QRectFixDimC& q = * (QuadSet::QRectFixDimC*) qi();
	if (SetupColor(q))
	{
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitiveType);
	  glVertex3f(q.fX,     z, q.fY);
	  glVertex3f(q.fX + w, z, q.fY);
	  glVertex3f(q.fX + w, z, q.fY + h);
	  glVertex3f(q.fX,     z, q.fY + h);
	  glEnd();
	}
      }
      break;
    }

    default:
      throw(eH + "unsupported quad-type.");

  } // end switch quad-type

  if (flags.SecSelection()) glPopName();
}


void QuadSetGL::RenderLines(const TGLDrawFlags &flags) const
{
  static const Exc_t eH("QuadSetGL::RenderLines ");

  QuadSet& mQ = * fM;

  VoidCPlex::iterator qi(mQ.fPlex);

  if (flags.SecSelection()) glPushName(0);

  switch (mQ.fQuadType)
  {

    case QuadSet::QT_LineXYFixedZ:
    {
      const Float_t& z = mQ.fDefCoord;
      while (qi.next()) {
	QuadSet::QLineFixC& q = * (QuadSet::QLineFixC*) qi();
	if (SetupColor(q))
	{
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(GL_LINES);
	  glVertex3f(q.fX,         q.fY,         z);
	  glVertex3f(q.fX + q.fDx, q.fY + q.fDy, z);
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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(GL_LINES);
	  glVertex3f(q.fX,         z, q.fY);
	  glVertex3f(q.fX + q.fDx, z, q.fY + q.fDy);
	  glEnd();
	}
      }
      break;
    }

    default:
      throw(eH + "unsupported quad-type.");

  }

  if (flags.SecSelection()) glPopName();
}

void QuadSetGL::RenderHexagons(const TGLDrawFlags &flags) const
{
  static const Exc_t eH("QuadSetGL::RenderHexagons ");

  const Float_t sqr3hf = 0.5*TMath::Sqrt(3);

  QuadSet& mQ = * fM;

  GLenum primitveType = (mQ.fRenderMode != QuadSet::RM_Line) ?
    GL_POLYGON : GL_LINE_LOOP;

  glNormal3f(0, 0, 1);

  VoidCPlex::iterator qi(mQ.fPlex);

  if (flags.SecSelection()) glPushName(0);

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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitveType);
	  glVertex3f( q.fR + q.fX,       q.fY, q.fZ);
	  glVertex3f(   rh + q.fX,  rs + q.fY, q.fZ);
	  glVertex3f(  -rh + q.fX,  rs + q.fY, q.fZ);
	  glVertex3f(-q.fR + q.fX,       q.fY, q.fZ);
	  glVertex3f(  -rh + q.fX, -rs + q.fY, q.fZ);
	  glVertex3f(   rh + q.fX, -rs + q.fY, q.fZ);
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
	  if (flags.SecSelection()) glLoadName(qi.index());
	  glBegin(primitveType);
	  glVertex3f( rs + q.fX,    rh + q.fY, q.fZ);
	  glVertex3f(      q.fX,  q.fR + q.fY, q.fZ);
	  glVertex3f(-rs + q.fX,    rh + q.fY, q.fZ);
	  glVertex3f(-rs + q.fX,   -rh + q.fY, q.fZ);
	  glVertex3f(      q.fX, -q.fR + q.fY, q.fZ);
	  glVertex3f( rs + q.fX,   -rh + q.fY, q.fZ);
	  glEnd();
	}
      }
      break;
    }

    default:
      throw(eH + "unsupported quad-type.");

  } // end switch quad-type

  if (flags.SecSelection()) glPopName();
}

/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________________
void QuadSetGL::ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*)
{
   // Processes secondary selection from TGLViewer.
   // Calls TPointSet3D::PointSelected(Int_t) with index of selected
   // point as an argument.

   if (ptr[0] < 2) return;
   fM->QuadSelected(ptr[4]);
}
