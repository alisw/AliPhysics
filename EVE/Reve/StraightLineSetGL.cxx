// $Header$

#include "StraightLineSetGL.h"
#include <Reve/StraightLineSet.h>
#include <Reve/GLUtilNS.h>

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>

#include <TGLIncludes.h>

using namespace Reve;

//______________________________________________________________________
// StraightLineSetGL
//

ClassImp(StraightLineSetGL)

StraightLineSetGL::StraightLineSetGL() : TGLObject(), fM(0)
{
  // fDLCache = false; // Disable display list.
}

StraightLineSetGL::~StraightLineSetGL()
{}

/**************************************************************************/

Bool_t StraightLineSetGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  if(SetModelCheckClass(obj, StraightLineSet::Class())) {
    fM = dynamic_cast<StraightLineSet*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

void StraightLineSetGL::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  SetAxisAlignedBBox(((StraightLineSet*)fExternalObj)->AssertBBox());
}

//______________________________________________________________________________
Bool_t StraightLineSetGL::ShouldCache(TGLRnrCtx & rnrCtx) const
{
   // Override from TGLDrawable.
   // To account for large point-sizes we modify the projection matrix
   // during selection and thus we need a direct draw.

   if (rnrCtx.Selection()) return kFALSE;
   return fDLCache;
}

/**************************************************************************/

void StraightLineSetGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // printf("StraightLineSetGL::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  StraightLineSet& mL = * fM;

  glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_ENABLE_BIT);
  GLUtilNS::GL_Capability_Switch lights_off(GL_LIGHTING, false);
 
  if(mL.fRnrLines && mL.fLinePlex.Size() > 0)
  {
    UChar_t color[4];
    ColorFromIdx(mL.GetMainColor(), color);
    glColor4ubv(color);

    VoidCPlex::iterator li(mL.fLinePlex);
    if(rnrCtx.SecSelection()) 
    {  
      GLuint name = 0;
      glPushName(1);
      glPushName(0);  
      while (li.next()) 
      {
	StraightLineSet::Line& l = * (StraightLineSet::Line*) li();
	glLoadName(name);
	{
	  glBegin(GL_LINES);  
	  glVertex3f(l.fV1[0], l.fV1[1], l.fV1[2]);
	  glVertex3f(l.fV2[0], l.fV2[1], l.fV2[2]);
	  glEnd();
	}
	name ++;
      }  
      glPopName(); 
      glPopName();
    } 
    else 
    {
      glBegin(GL_LINES);    
      while (li.next()) 
      {
	StraightLineSet::Line& l = * (StraightLineSet::Line*) li();
	glVertex3f(l.fV1[0], l.fV1[1], l.fV1[2]);
	glVertex3f(l.fV2[0], l.fV2[1], l.fV2[2]);
      }    
      glEnd();
    }
  }

  if(mL.fRnrMarkers && mL.fMarkerPlex.Size() > 0)
  {
    UChar_t color[4];
    ColorFromIdx(mL.GetMarkerColor(), color);
    glColor4ubv(color);

    VoidCPlex::iterator mi(mL.fMarkerPlex);
    Float_t* pnts = new Float_t[mL.fMarkerPlex.Size()*3];
    Float_t* pnt  = pnts;
    Int_t lidx = -1; 
    while (mi.next()) 
    {
      StraightLineSet::Marker& m = * (StraightLineSet::Marker*) mi();
      lidx = m.fLineID;
      StraightLineSet::Line& l = * (StraightLineSet::Line*) mL.fLinePlex.Atom(lidx);
      pnt[0] = l.fV1[0] + (l.fV2[0] - l.fV1[0])*m.fPos;
      pnt[1] = l.fV1[1] + (l.fV2[1] - l.fV1[1])*m.fPos;
      pnt[2] = l.fV1[2] + (l.fV2[2] - l.fV1[2])*m.fPos;;
      pnt   += 3;
    }
    if(rnrCtx.SecSelection()) glPushName(2);
    GLUtilNS::RenderPolyMarkers((TAttMarker&)mL, pnts, mL.fMarkerPlex.Size(),
				rnrCtx.Selection(), rnrCtx.SecSelection());
    if(rnrCtx.SecSelection()) glPopName();
    delete [] pnts;
  }

  glPopAttrib();
}

/**************************************************************************/

void StraightLineSetGL::ProcessSelection(TGLRnrCtx       & /*rnrCtx*/,
					 TGLSelectRecord & rec)
{ 
  if (rec.GetN() != 3) return;
  if(rec.GetItem(1) == 1)
  {
    printf("selected line %d\n", rec.GetItem(2));
  }
  else 
  {
    StraightLineSet::Marker& m = * (StraightLineSet::Marker*) fM->fMarkerPlex.Atom(rec.GetItem(2));
    printf("Selected point %d on line %d\n", rec.GetItem(2), m.fLineID);
  }
}
