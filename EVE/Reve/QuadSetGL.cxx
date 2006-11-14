// $Header$

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
  }
  else
  {
    const RGBAPalette& pal = *fM->fPalette;

    if (q.fValue == fM->fDefaultValue)
    {
      glColor4ubv(pal.GetDefaultRGBA());
    }
    else if ( pal.WithinVisibleRange(q.fValue) == kFALSE )
    {
      return kFALSE;
    }
    else
    {
      glColor4ubv(pal.ColorFromArray(q.fValue));
    }
  }
  return kTRUE;
}

/**************************************************************************/

void QuadSetGL::DirectDraw(const TGLDrawFlags & flags) const
{
  static const Exc_t eH("QuadSetGL::DirectDraw ");

  // printf("QuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  QuadSet& mQ = * fM;
  if (mQ.fPlex.Size() == 0)
    return;
  if ( ! mQ.fValueIsColor && mQ.fPalette == 0)
  {
    Int_t min, max;
    mQ.ScanMinMaxValues(min, max);
    mQ.fPalette = new RGBAPalette(min, max, kTRUE, kFALSE);
  }

  if (mQ.fFrame != 0)
    FrameBoxGL::Render(mQ.fFrame);

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);

  if (mQ.fRenderMode == QuadSet::RM_Fill)
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }
  else if (mQ.fRenderMode == QuadSet::RM_Line)
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDisable(GL_LIGHTING);
  }

  if (mQ.fQuadType < QuadSet::QT_LineFixedZ)
    RenderQuads(flags);
  else
    RenderLines(flags);

  glPopAttrib();

}


void QuadSetGL::RenderQuads(const TGLDrawFlags &) const
{
  static const Exc_t eH("QuadSetGL::RenderQuads ");

  QuadSet& mQ = * fM;

  if (mQ.fRenderMode != QuadSet::RM_Line)
  {

    if (mQ.fQuadType == QuadSet::QT_FreeQuad)
      glEnable(GL_NORMALIZE);
    else
      glNormal3f(0, 0, 1);

    glBegin(GL_QUADS);
    for (Int_t c=0; c<mQ.fPlex.VecSize(); ++c)
    {
      QuadSet::QuadBase* qbp = (QuadSet::QuadBase*) mQ.fPlex.Chunk(c);
      Int_t n = mQ.fPlex.NAtoms(c);

      switch (mQ.fQuadType)
      {

	case QuadSet::QT_FreeQuad:
	{
	  QuadSet::FreeQuad* qp = (QuadSet::FreeQuad*) qbp;
	  Float_t e1[3], e2[3], normal[3];
	  while (n--) {
	    if (SetupColor(*qp) == kFALSE) continue;
	    Float_t* p = qp->fVertices;
	    e1[0] = p[3] - p[0]; e1[1] = p[4] - p[1]; e1[2] = p[5] - p[2];
	    e2[0] = p[6] - p[0]; e2[1] = p[7] - p[1]; e2[2] = p[8] - p[2];
	    TMath::Cross(e1, e2, normal);
	    glNormal3fv(normal);
	    glVertex3fv(p);
	    glVertex3fv(p + 3);
	    glVertex3fv(p + 6);
	    glVertex3fv(p + 9);
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAligned:
	{
	  QuadSet::AAQuad* qp = (QuadSet::AAQuad*) qbp;
	  while (n--) {
	    QuadSet::AAQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glVertex3f(q.fX,        q.fY,        q.fZ);
	    glVertex3f(q.fX + q.fW, q.fY,        q.fZ);
	    glVertex3f(q.fX + q.fW, q.fY + q.fH, q.fZ);
	    glVertex3f(q.fX,        q.fY + q.fH, q.fZ);
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedDim:
	{
	  QuadSet::AAFixDimQuad* qp = (QuadSet::AAFixDimQuad*) qbp;
	  const Float_t& w = mQ.fDefWidth;
	  const Float_t& h = mQ.fDefHeight;
	  while (n--) {
	    QuadSet::AAFixDimQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glVertex3f(q.fX,     q.fY,     q.fZ);
	    glVertex3f(q.fX + w, q.fY,     q.fZ);
	    glVertex3f(q.fX + w, q.fY + h, q.fZ);
	    glVertex3f(q.fX,     q.fY + h, q.fZ);
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedZ:
	{
	  QuadSet::AAFixZQuad* qp = (QuadSet::AAFixZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  while (n--) {
	    QuadSet::AAFixZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glVertex3f(q.fX,        q.fY,        z);
	    glVertex3f(q.fX + q.fW, q.fY,        z);
	    glVertex3f(q.fX + q.fW, q.fY + q.fH, z);
	    glVertex3f(q.fX,        q.fY + q.fH, z);
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedY:
	{
	  QuadSet::AAFixZQuad* qp = (QuadSet::AAFixZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  while (n--) {
	    QuadSet::AAFixZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glVertex3f(q.fX,        z, q.fY);
	    glVertex3f(q.fX + q.fW, z, q.fY);
	    glVertex3f(q.fX + q.fW, z, q.fY + q.fH);
	    glVertex3f(q.fX,        z, q.fY + q.fH);
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedDimZ:
	{
	  QuadSet::AAFixDimZQuad* qp = (QuadSet::AAFixDimZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  const Float_t& w = mQ.fDefWidth;
	  const Float_t& h = mQ.fDefHeight;
	  while (n--) {
	    QuadSet::AAFixDimZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glVertex3f(q.fX,     q.fY,     z);
	    glVertex3f(q.fX + w, q.fY,     z);
	    glVertex3f(q.fX + w, q.fY + h, z);
	    glVertex3f(q.fX,     q.fY + h, z);
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedDimY:
	{
	  QuadSet::AAFixDimZQuad* qp = (QuadSet::AAFixDimZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  const Float_t& w = mQ.fDefWidth;
	  const Float_t& h = mQ.fDefHeight;
	  while (n--) {
	    QuadSet::AAFixDimZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glVertex3f(q.fX,     z, q.fY);
	    glVertex3f(q.fX + w, z, q.fY);
	    glVertex3f(q.fX + w, z, q.fY + h);
	    glVertex3f(q.fX,     z, q.fY + h);
	    ++qp;
	  }
	  break;
	}

	default:
	  throw(eH + "unsupported quad-type.");

      } // end switch quad-type

    } // end for chunk
    glEnd();

  }
  else
  {

    for (Int_t c=0; c<mQ.fPlex.VecSize(); ++c)
    {
      QuadSet::QuadBase* qbp = (QuadSet::QuadBase*) mQ.fPlex.Chunk(c);
      Int_t n = mQ.fPlex.NAtoms(c);

      switch (mQ.fQuadType)
      {

	case QuadSet::QT_FreeQuad:
	{
	  QuadSet::FreeQuad* qp = (QuadSet::FreeQuad*) qbp;
	  while (n--) {
	    if (SetupColor(*qp) == kFALSE) continue;
	    SetupColor(*qp);
	    Float_t* p = qp->fVertices;
	    glBegin(GL_LINE_LOOP);
	    glVertex3fv(p);
	    glVertex3fv(p + 3);
	    glVertex3fv(p + 6);
	    glVertex3fv(p + 9);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAligned:
	{
	  QuadSet::AAQuad* qp = (QuadSet::AAQuad*) qbp;
	  while (n--) {
	    QuadSet::AAQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glBegin(GL_LINE_LOOP);
	    glVertex3f(q.fX,        q.fY,        q.fZ);
	    glVertex3f(q.fX + q.fW, q.fY,        q.fZ);
	    glVertex3f(q.fX + q.fW, q.fY + q.fH, q.fZ);
	    glVertex3f(q.fX,        q.fY + q.fH, q.fZ);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedDim:
	{
	  QuadSet::AAFixDimQuad* qp = (QuadSet::AAFixDimQuad*) qbp;
	  const Float_t& w = mQ.fDefWidth;
	  const Float_t& h = mQ.fDefHeight;
	  while (n--) {
	    QuadSet::AAFixDimQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glBegin(GL_LINE_LOOP);
	    glVertex3f(q.fX,     q.fY,     q.fZ);
	    glVertex3f(q.fX + w, q.fY,     q.fZ);
	    glVertex3f(q.fX + w, q.fY + h, q.fZ);
	    glVertex3f(q.fX,     q.fY + h, q.fZ);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedZ:
	{
	  QuadSet::AAFixZQuad* qp = (QuadSet::AAFixZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  while (n--) {
	    QuadSet::AAFixZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glBegin(GL_LINE_LOOP);
	    glVertex3f(q.fX,        q.fY,        z);
	    glVertex3f(q.fX + q.fW, q.fY,        z);
	    glVertex3f(q.fX + q.fW, q.fY + q.fH, z);
	    glVertex3f(q.fX,        q.fY + q.fH, z);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedY:
	{
	  QuadSet::AAFixZQuad* qp = (QuadSet::AAFixZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  while (n--) {
	    QuadSet::AAFixZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glBegin(GL_LINE_LOOP);
	    glVertex3f(q.fX,        z, q.fY);
	    glVertex3f(q.fX + q.fW, z, q.fY);
	    glVertex3f(q.fX + q.fW, z, q.fY + q.fH);
	    glVertex3f(q.fX,        z, q.fY + q.fH);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedDimZ:
	{
	  QuadSet::AAFixDimZQuad* qp = (QuadSet::AAFixDimZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  const Float_t& w = mQ.fDefWidth;
	  const Float_t& h = mQ.fDefHeight;
	  while (n--) {
	    QuadSet::AAFixDimZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glBegin(GL_LINE_LOOP);
	    glVertex3f(q.fX,     q.fY,     z);
	    glVertex3f(q.fX + w, q.fY,     z);
	    glVertex3f(q.fX + w, q.fY + h, z);
	    glVertex3f(q.fX,     q.fY + h, z);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	case QuadSet::QT_AxisAlignedFixedDimY:
	{
	  QuadSet::AAFixDimZQuad* qp = (QuadSet::AAFixDimZQuad*) qbp;
	  const Float_t& z = mQ.fDefCoord;
	  const Float_t& w = mQ.fDefWidth;
	  const Float_t& h = mQ.fDefHeight;
	  while (n--) {
	    QuadSet::AAFixDimZQuad& q = * qp;
	    if (SetupColor(q) == kFALSE) continue;
	    glBegin(GL_LINE_LOOP);
	    glVertex3f(q.fX,     z, q.fY);
	    glVertex3f(q.fX + w, z, q.fY);
	    glVertex3f(q.fX + w, z, q.fY + h);
	    glVertex3f(q.fX,     z, q.fY + h);
	    glEnd();
	    ++qp;
	  }
	  break;
	}

	default:
	  throw(eH + "unsupported quad-type.");

      } // end switch quad-type

    } // end for chunk

  } // end else of RenderMode
}


void QuadSetGL::RenderLines(const TGLDrawFlags &) const
{
  static const Exc_t eH("QuadSetGL::RenderLines ");

  QuadSet& mQ = * fM;

  for (Int_t c=0; c<mQ.fPlex.VecSize(); ++c)
  {
    QuadSet::QuadBase* qbp = (QuadSet::QuadBase*) mQ.fPlex.Chunk(c);
    Int_t n = mQ.fPlex.NAtoms(c);

    switch (mQ.fQuadType)
    {

      case QuadSet::QT_LineFixedZ:
      {
	QuadSet::LineFixedZ* qp = (QuadSet::LineFixedZ*) qbp;
	const Float_t& z = mQ.fDefCoord;
	glBegin(GL_LINES);
	while (n--) {
	  QuadSet::LineFixedZ& q = * qp;
	  SetupColor(q);
	  glVertex3f(q.fX,         q.fY,         z);
	  glVertex3f(q.fX + q.fDx, q.fY + q.fDy, z);
	  ++qp;
	}
	glEnd();
	break;
      }

      case QuadSet::QT_LineFixedY:
      {
	QuadSet::LineFixedZ* qp = (QuadSet::LineFixedZ*) qbp;
	const Float_t& z = mQ.fDefCoord;
	glBegin(GL_LINES);
	while (n--) {
	  QuadSet::LineFixedZ& q = * qp;
	  SetupColor(q);
	  glVertex3f(q.fX,         z, q.fY);
	  glVertex3f(q.fX + q.fDx, z, q.fY + q.fDy);
	  ++qp;
	}
	glEnd();
	break;
      }

      default:
	throw(eH + "unsupported quad-type.");

    }
  }
}
