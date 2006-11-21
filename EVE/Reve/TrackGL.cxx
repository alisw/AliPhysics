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

void TrackGL::ProcessSelection(UInt_t* ptr, TGLViewer* v, TGLScene* s)
{
   // Processes secondary selection from TGLViewer.
   // Calls TPointSet3D::PointSelected(Int_t) with index of selected
   // point as an argument.
  TGLLogicalShape::ProcessSelection(ptr, v, s);
  /*
  if (ptr[0] < 2) return;
  TPointSet3D& q = * (TPointSet3D*) fExternalObj;
  q.PointSelected(ptr[4]);
  */
}
