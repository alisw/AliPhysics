// $Header$

#include "BoxSetGL.h"
#include <Reve/BoxSet.h>

#include <TGLIncludes.h>
#include <TGLRnrCtx.h>
#include <TGLScene.h>
#include <TGLSelectRecord.h>
#include <TGLContext.h>

using namespace Reve;

//______________________________________________________________________________
// BoxSetGL
//
// A GL rendering class for BoxSet.
//

ClassImp(BoxSetGL)

//______________________________________________________________________________
BoxSetGL::BoxSetGL() : fM(0), fBoxDL(0)
{
  // Default constructor.

  // fDLCache = false; // Disable display list.
}

//______________________________________________________________________________
BoxSetGL::~BoxSetGL()
{
  // Destructor. Noop.
}

/**************************************************************************/
// Protected methods
/**************************************************************************/

//______________________________________________________________________________
Int_t BoxSetGL::PrimitiveType() const
{
  // Return GL primitive used to render the boxes, based on the
  // render-mode specified in the model object.

  return (fM->fRenderMode != DigitSet::RM_Line) ? GL_QUADS : GL_LINE_LOOP;
}

//______________________________________________________________________________
inline Bool_t BoxSetGL::SetupColor(const DigitSet::DigitBase& q) const
{
  // Set GL color for given primitive.

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

//______________________________________________________________________________
void BoxSetGL::MakeOriginBox(Float_t p[24], Float_t dx, Float_t dy, Float_t dz) const
{
  // Fill array p to represent a box (0,0,0) - (dx,dy,dz).

  // bottom
  p[0] = 0;  p[1] = dy; p[2] = 0;  p += 3;
  p[0] = dx; p[1] = dy; p[2] = 0;  p += 3;
  p[0] = dx; p[1] = 0;  p[2] = 0;  p += 3;
  p[0] = 0;  p[1] = 0;  p[2] = 0;  p += 3;
  // top
  p[0] = 0;  p[1] = dy; p[2] = dz; p += 3;
  p[0] = dx; p[1] = dy; p[2] = dz; p += 3;
  p[0] = dx; p[1] = 0;  p[2] = dz; p += 3;
  p[0] = 0;  p[1] = 0;  p[2] = dz;
}

//______________________________________________________________________________
inline void BoxSetGL::RenderBox(const Float_t p[24]) const
{
  // Render a box specified by points in array p.

  // bottom: 0123
  glNormal3f(0, 0, -1);
  glVertex3fv(p);      glVertex3fv(p + 3);
  glVertex3fv(p + 6);  glVertex3fv(p + 9);
  // top:    7654
  glNormal3f(0, 0, 1);
  glVertex3fv(p + 21); glVertex3fv(p + 18);
  glVertex3fv(p + 15); glVertex3fv(p + 12);
  // back:  0451
  glNormal3f(0, 1, 0);
  glVertex3fv(p);      glVertex3fv(p + 12);
  glVertex3fv(p + 15); glVertex3fv(p + 3);
  // front:   3267
  glNormal3f(0, -1, 0);
  glVertex3fv(p + 9);   glVertex3fv(p + 6);
  glVertex3fv(p + 18);  glVertex3fv(p + 21);
  // left:    0374
  glNormal3f(-1, 0, 0);
  glVertex3fv(p);       glVertex3fv(p + 9);
  glVertex3fv(p + 21);  glVertex3fv(p + 12);
  // right:   1562
  glNormal3f(1, 0, 0);
  glVertex3fv(p + 3);   glVertex3fv(p + 15);
  glVertex3fv(p + 18);  glVertex3fv(p + 6);
}

//______________________________________________________________________________
void BoxSetGL::MakeDisplayList() const
{
  // Create a display-list for rendering a single box, based on the
  // current box-type.
  // Some box-types don't benefit from the display-list rendering and
  // so display-list is not created.

  if (fM->fBoxType == BoxSet::BT_AABox ||
      fM->fBoxType == BoxSet::BT_AABoxFixedDim)
  {
    if (fBoxDL == 0)
      fBoxDL = glGenLists(1);

    Float_t p[24];
    if (fM->fBoxType == BoxSet::BT_AABox)
      MakeOriginBox(p, 1.0f, 1.0f, 1.0f);
    else
      MakeOriginBox(p, fM->fDefWidth, fM->fDefHeight, fM->fDefDepth);

    glNewList(fBoxDL, GL_COMPILE);
    glBegin(PrimitiveType());
    RenderBox(p);
    glEnd();
    glEndList();
  }
}

/**************************************************************************/
// Virtuals from base-classes
/**************************************************************************/

//______________________________________________________________________________
Bool_t BoxSetGL::ShouldDLCache(const TGLRnrCtx & rnrCtx) const
{
  // Determines if display-list will be used for rendering.
  // Virtual from TGLLogicalShape.

  MakeDisplayList();

  if (rnrCtx.DrawPass() == TGLRnrCtx::kPassOutlineLine)
    return kFALSE;
  return TGLObject::ShouldDLCache(rnrCtx);
}

//______________________________________________________________________________
void BoxSetGL::DLCacheDrop()
{
  // Called when display lists have been destroyed externally and the
  // internal display-list data needs to be cleare.
  // Virtual from TGLLogicalShape.

  fBoxDL = 0;
  TGLObject::DLCacheDrop();
}

//______________________________________________________________________________
void BoxSetGL::DLCachePurge()
{
  // Called when display-lists need to be returned to the system.
  // Virtual from TGLLogicalShape.

  static const Exc_t eH("BoxSetGL::DLCachePurge ");

  if (fBoxDL == 0) return;
  if (fScene)
  {
    fScene->GetGLCtxIdentity()->RegisterDLNameRangeToWipe(fBoxDL, 1);
  }
  else
  {
    Warning(eH, "Scene unknown, attempting direct deletion.");
    glDeleteLists(fBoxDL, 1);
  }
  TGLObject::DLCachePurge();
}

/**************************************************************************/

//______________________________________________________________________________
Bool_t BoxSetGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  // Set model object.
  // Virtual from TGLObject.

  Bool_t isok = SetModelCheckClass(obj, Reve::BoxSet::Class());
  fM = isok ? dynamic_cast<Reve::BoxSet*>(obj) : 0;
  return isok;
}

//______________________________________________________________________________
void BoxSetGL::SetBBox()
{
  // Fill the bounding-box data of the logical-shape.
  // Virtual from TGLObject.

  SetAxisAlignedBBox(fM->AssertBBox());
}

//______________________________________________________________________________
void BoxSetGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // Actual rendering code.
  // Virtual from TGLLogicalShape.

  static const Exc_t eH("BoxSetGL::DirectDraw ");

  if (rnrCtx.DrawPass() == TGLRnrCtx::kPassOutlineLine)
    return;

  BoxSet& mB = * fM;
  // printf("BoxSetGL::DirectDraw N boxes %d\n", mB.fPlex.Size());
  if(mB.fPlex.Size() == 0)
    return;

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  if (mB.fRenderMode == DigitSet::RM_Fill)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  else if (mB.fRenderMode == DigitSet::RM_Line)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  if (mB.fDisableLigting) glDisable(GL_LIGHTING);

  if (rnrCtx.SecSelection()) glPushName(0);

  Int_t boxSkip = 0;
  if (rnrCtx.ShapeLOD() < 50)
    boxSkip = 6 - (rnrCtx.ShapeLOD()+1)/10;

  VoidCPlex::iterator bi(mB.fPlex);

  switch (mB.fBoxType)
  {

    case BoxSet::BT_FreeBox:
    {
      GLenum primitiveType = PrimitiveType();
      while (bi.next())
      {
        BoxSet::BFreeBox& b = * (BoxSet::BFreeBox*) bi();
        if (SetupColor(b))
	{
          if (rnrCtx.SecSelection()) glLoadName(bi.index());
          glBegin(primitiveType);
          RenderBox(b.fVertices);
          glEnd();
        }
        if (boxSkip) { Int_t s = boxSkip; while (s--) bi.next(); }
      }
      break;
    } // end case free-box

    case BoxSet::BT_AABox:
    {
      glEnable(GL_NORMALIZE);
      while (bi.next())
      {
        BoxSet::BAABox& b = * (BoxSet::BAABox*) bi();
        if (SetupColor(b))
	{
          if (rnrCtx.SecSelection()) glLoadName(bi.index());
          glPushMatrix();
          glTranslatef(b.fA, b.fB, b.fC);
          glScalef    (b.fW, b.fH, b.fD);
          glCallList(fBoxDL);
          glPopMatrix();
        }
        if (boxSkip) { Int_t s = boxSkip; while (s--) bi.next(); }
      }
      break;
    }

    case BoxSet::BT_AABoxFixedDim:
    {
      while (bi.next())
      {
        BoxSet::BAABoxFixedDim& b = * (BoxSet::BAABoxFixedDim*) bi();
        if (SetupColor(b))
	{
          if (rnrCtx.SecSelection()) glLoadName(bi.index());
          glTranslatef(b.fA, b.fB, b.fC);
          glCallList(fBoxDL);
          glTranslatef(-b.fA, -b.fB, -b.fC);
        }
        if (boxSkip) { Int_t s = boxSkip; while (s--) bi.next(); }
      }
      break;
    }

    default:
    {
      throw(eH + "unsupported box-type.");
    }

  } // end switch box-type

  if (rnrCtx.SecSelection()) glPopName();

  glPopAttrib();
}

/**************************************************************************/

//______________________________________________________________________________
void BoxSetGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  if (rec.GetN() < 2) return;
  fM->DigitSelected(rec.GetItem(1));
}

//______________________________________________________________________________
void BoxSetGL::Render(TGLRnrCtx & rnrCtx)
{
  // Interface for direct rendering from classes that include BoxSet
  // as a member.

  MakeDisplayList();
  DirectDraw(rnrCtx);
  glDeleteLists(fBoxDL, 1);
  fBoxDL = 0;
}
