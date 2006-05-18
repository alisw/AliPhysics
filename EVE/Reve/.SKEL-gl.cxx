// $Header$

#include "CLASS.h"
#include <Reve/STEM.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// CLASS
//

ClassImp(CLASS)

CLASS::CLASS()
{
  // fCached = false; // Disable display list.
}

CLASS::~CLASS()
{}

/**************************************************************************/

Bool_t CLASS::SetModel(TObject* obj)
{
  return set_model(obj, "Reve::STEM");
}

void CLASS::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  set_axis_aligned_bbox(((STEM*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void CLASS::DirectDraw(const TGLDrawFlags & flags) const
{
  // printf("CLASS::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());
}
