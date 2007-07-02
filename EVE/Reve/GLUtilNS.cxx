#include "GLUtilNS.h"

#include <TAttMarker.h>
#include <TGLIncludes.h>

namespace GLUtilNS {

void RenderPolyMarkers(TAttMarker& marker, Float_t* p, Int_t n,
		       Bool_t selection, Bool_t sec_selection)
{
  // Store attributes GL_POINT_BIT and GL_LINE_BIT before call this function !

  Int_t s = marker.GetMarkerStyle(); 
  if (s == 2 || s == 3 || s == 5 || s == 28)
    RenderCrosses(marker, p, n, sec_selection); 
  else
    RenderPoints(marker, p, n, selection, sec_selection); 
}

//______________________________________________________________________________
void RenderPoints(TAttMarker& marker, Float_t* p, Int_t n,
		  Bool_t selection, Bool_t sec_selection)
{
  // Render markers as circular or square points.
  {
    Int_t ms = marker.GetMarkerStyle();
    Float_t size = 5*marker.GetMarkerSize();
    if (ms == 4 || ms == 20 || ms == 24) 
    {
      if (ms == 4 || ms == 24)
	glEnable(GL_BLEND);
      glEnable(GL_POINT_SMOOTH);
    } else 
    {
      glDisable(GL_POINT_SMOOTH);
      if      (ms == 1) size = 1;
      else if (ms == 6) size = 2;
      else if (ms == 7) size = 3;
    }
    glPointSize(size);
   
    // During selection extend picking region for large point-sizes.
    static const Int_t sPickRadius = 3; // Hardcoded also in TGLViewer::RequestSelect()
    Bool_t changePM = kFALSE;
    if (selection && size > sPickRadius) 
    {
      changePM = kTRUE;
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      Float_t pm[16];
      glGetFloatv(GL_PROJECTION_MATRIX, pm);
      Float_t scale = (Float_t) sPickRadius / size;
      for (Int_t i=0; i<=12; i+=4) {
	pm[i] *= scale; pm[i+1] *= scale;
      }
      glLoadMatrixf(pm);
    }

    if (sec_selection) 
    {
      glPushName(0);
      for (Int_t i=0; i<n; ++i, p+=3) 
      {
	glLoadName(i);
	glBegin(GL_POINTS);
	glVertex3fv(p);
	glEnd();
      }
      glPopName();
    } 
    else 
    {
      glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
      glVertexPointer(3, GL_FLOAT, 0, p);
      glEnableClientState(GL_VERTEX_ARRAY);
      { // Circumvent bug in ATI's linux drivers.
	Int_t nleft = n;
	Int_t ndone = 0;
	const Int_t maxChunk = 8192;
	while (nleft > maxChunk) 
	{
	  glDrawArrays(GL_POINTS, ndone, maxChunk);
	  nleft -= maxChunk;
	  ndone += maxChunk;
	}
	glDrawArrays(GL_POINTS, ndone, nleft);
      }
      glPopClientAttrib();
    }

    if (changePM) 
    {
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
    }
  }
}

//______________________________________________________________________________
void RenderCrosses(TAttMarker& marker, Float_t* p, Int_t n,
		   Bool_t sec_selection)
{
  // Render markers as crosses.
  //
  if (marker.GetMarkerStyle() == 28) 
  {
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2);
  } 
  else 
  {
    glDisable(GL_LINE_SMOOTH);
  }

  // cross dim
  const Float_t  d = marker.GetMarkerSize();
  if (sec_selection) 
  {
    glPushName(0);
    for (Int_t i=0; i<n; ++i, p+=3) 
    {
      glLoadName(i);
      glBegin(GL_LINES);
      glVertex3f(p[0]-d, p[1], p[2]); glVertex3f(p[0]+d, p[1], p[2]);
      glVertex3f(p[0], p[1]-d, p[2]); glVertex3f(p[0], p[1]+d, p[2]);
      glVertex3f(p[0], p[1], p[2]-d); glVertex3f(p[0], p[1], p[2]+d);
      glEnd();
    }
    glPopName();
  } 
  else 
  {
    glBegin(GL_LINES);
    for (Int_t i=0; i<n; ++i, p+=3) 
    {
      glVertex3f(p[0]-d, p[1], p[2]); glVertex3f(p[0]+d, p[1], p[2]);
      glVertex3f(p[0], p[1]-d, p[2]); glVertex3f(p[0], p[1]+d, p[2]);
      glVertex3f(p[0], p[1], p[2]-d); glVertex3f(p[0], p[1], p[2]+d);
    }
    glEnd();
  }
}

/**************************************************************************/
} // end namespace GLUtilNS
