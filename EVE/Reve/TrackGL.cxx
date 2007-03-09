// $Header$

#include "TrackGL.h"
#include <Reve/Track.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// TrackGL
//

ClassImp(TrackGL)

TrackGL::TrackGL() : LineGL()
{
  // fCached = false; // Disable display list.
}

TrackGL::~TrackGL()
{}

/**************************************************************************/

Bool_t TrackGL::SetModel(TObject* obj)
{
  if(LineGL::SetModel(obj) == kFALSE) return kFALSE;
  if(SetModelCheckClass(obj, Track::Class())) {
    fTrack = dynamic_cast<Track*>(obj);
    return kTRUE;
  }
  return kFALSE;
}
/**************************************************************************/

void TrackGL::ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  Int_t n = ptr[0];
  printf("TrackGL::ProcessSelection %d names on the stack (z1=%g, z2=%g).\n",
	 n, Float_t(ptr[1])/0x7fffffff, Float_t(ptr[2])/0x7fffffff);
  ptr += 3;
  printf("  Names: ");
  for (Int_t j=0; j<n; ++j, ++ptr) printf ("%d ", *ptr);
  printf("\n");

  ((Track*)fM)->CtrlClicked((Track*)fM);
}

/**************************************************************************/
void TrackGL::DirectDraw(const TGLDrawFlags & flags) const
{
  // Render line and path marks

  LineGL::DirectDraw(flags);

  if ( ! fTrack->fPathMarks.empty()){

    TrackRnrStyle* rs = fTrack->GetRnrStyle();
    Int_t  style = rs->fPMStyle;

    UChar_t color[4];
    ColorFromIdx(rs->fPMColor, color);
    glColor4ubv(color);
    
    glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);

    Int_t ms = rs->fPMStyle;
    // points
    if (ms != 2 && ms != 3 && ms != 5 && ms != 28) 
    {
      Float_t size = 5*rs->fPMSize;
      if (style == 4 || style == 20 || style == 24) 
      {
	if (style == 4 || style == 24)
	  glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
      } 
      else 
      {
	glDisable(GL_POINT_SMOOTH);
	if      (style == 1) size = 1;
	else if (style == 6) size = 2;
	else if (style == 7) size = 3;
      }
      glPointSize(size);

      glBegin(GL_POINTS);
      Bool_t accept;
      std::vector<PathMark*>& pm = fTrack->fPathMarks;
      for(std::vector<PathMark*>::iterator i=pm.begin(); i!=pm.end(); ++i) 
      {
	accept = kFALSE;
	switch((*i)->type)
	{
	  case(PathMark::Daughter):
	    if(rs->fRnrDaughters) accept = kTRUE;
	    break;
	  case(PathMark::Reference):
	    if(rs->fRnrReferences) accept = kTRUE;
	    break;
	  case(PathMark::Decay):
	    if(rs->fRnrDecay) accept = kTRUE;
	    break;
	} 
	if(accept)
	{
	  if((TMath::Abs((*i)->V.z) < rs->fMaxZ) && ((*i)->V.Perp() < rs->fMaxR))
	    glVertex3f((*i)->V.x, (*i)->V.y,(*i)->V.z);
	}
      } 
      glEnd();
    } // end render points
    else 
    {
      // crosses
      if ( style== 28) 
      {
	glEnable(GL_BLEND);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(2);
      } 
      else 
      {
	glDisable(GL_LINE_SMOOTH);
      }

      glBegin(GL_LINES);
      Bool_t accept;
      Float_t d = 2* rs->fPMSize;
      Float_t p[3];
      std::vector<PathMark*>& pm = fTrack->fPathMarks;     
      for(std::vector<PathMark*>::iterator i=pm.begin(); i!=pm.end(); ++i) 
      {
	accept = kFALSE;
	switch((*i)->type)
	{
	  case(PathMark::Daughter):
	    if(rs->fRnrDaughters) accept = kTRUE;
	    break;
	  case(PathMark::Reference):
	    if(rs->fRnrReferences) accept = kTRUE;
	    break;
	  case(PathMark::Decay):
	    if(rs->fRnrDecay) accept = kTRUE;
	    break;
	} 
	if(accept)
	{
	  if((TMath::Abs((*i)->V.z) < rs->fMaxZ) && ((*i)->V.Perp() < rs->fMaxR))
	  {
	    p[0] = (*i)->V.x; p[1] = (*i)->V.y; p[2] = (*i)->V.z;
	    glVertex3f(p[0]-d, p[1], p[2]); glVertex3f(p[0]+d, p[1], p[2]);
	    glVertex3f(p[0], p[1]-d, p[2]); glVertex3f(p[0], p[1]+d, p[2]);
	    glVertex3f(p[0], p[1], p[2]-d); glVertex3f(p[0], p[1], p[2]+d);
	  }
	}
      } 
      glEnd();
    } // end render corsses
    glPopAttrib();
  } //if PM not empty
}
