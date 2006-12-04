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

  /*
    if (ptr[0] < 2) return;
    TPointSet3D& q = * (TPointSet3D*) fExternalObj;
    q.PointSelected(ptr[4]);
  */
}
