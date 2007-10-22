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

  // lines
  GLUtilNS::GL_Capability_Switch lights_off(GL_LIGHTING, false);
  if(mL.GetRnrLines() && mL.GetLinePlex().Size() > 0)
  {
    glDisable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    UChar_t color[4];
    Reve::ColorFromIdx(mL.GetLineColor(), color);
    glColor4ubv(color);  
    glLineWidth(mL.GetLineWidth());
    if (mL.GetLineStyle() > 1) {
      Int_t    fac = 1;
      UShort_t pat = 0xffff;
      switch (mL.GetLineStyle()) {
	case 2:  pat = 0x3333; break;
	case 3:  pat = 0x5555; break;
	case 4:  pat = 0xf040; break;
	case 5:  pat = 0xf4f4; break;
	case 6:  pat = 0xf111; break;
	case 7:  pat = 0xf0f0; break;
	case 8:  pat = 0xff11; break;
	case 9:  pat = 0x3fff; break;
	case 10: pat = 0x08ff; fac = 2; break;
      }
      glLineStipple(1, pat);
      glEnable(GL_LINE_STIPPLE);
    }

    VoidCPlex::iterator li(mL.GetLinePlex());
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
  glPopAttrib();
 

  // markers
  if(mL.GetRnrMarkers() && mL.GetMarkerPlex().Size() > 0)
  {
    VoidCPlex::iterator mi(mL.GetMarkerPlex());
    Float_t* pnts = new Float_t[mL.GetMarkerPlex().Size()*3];
    Float_t* pnt  = pnts;
    Int_t lidx = -1; 
    while (mi.next()) 
    {
      StraightLineSet::Marker& m = * (StraightLineSet::Marker*) mi();
      lidx = m.fLineID;
      StraightLineSet::Line& l = * (StraightLineSet::Line*) mL.GetLinePlex().Atom(lidx);
      pnt[0] = l.fV1[0] + (l.fV2[0] - l.fV1[0])*m.fPos;
      pnt[1] = l.fV1[1] + (l.fV2[1] - l.fV1[1])*m.fPos;
      pnt[2] = l.fV1[2] + (l.fV2[2] - l.fV1[2])*m.fPos;;
      pnt   += 3;
    }
    if(rnrCtx.SecSelection()) glPushName(2);
    GLUtilNS::RenderPolyMarkers((TAttMarker&)mL, pnts, mL.GetMarkerPlex().Size(),
				rnrCtx.Selection(), rnrCtx.SecSelection());
    if(rnrCtx.SecSelection()) glPopName();
    delete [] pnts;
  }

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
    StraightLineSet::Marker& m = * (StraightLineSet::Marker*) fM->GetMarkerPlex().Atom(rec.GetItem(2));
    printf("Selected point %d on line %d\n", rec.GetItem(2), m.fLineID);
  }
}
