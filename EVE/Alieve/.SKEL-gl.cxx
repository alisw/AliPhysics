// $Header$

#include "CLASS.h"
#include <Alieve/STEM.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;
using namespace Alieve;

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
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  if(set_model(obj, "Alieve::STEM")) {
#else
  if(SetModelCheckClass(obj, "Alieve::STEM")) {
#endif
    fSector = (TPCSector3D*) fExternalObj;
    return kTRUE;
  }
  return kFALSE;
}

void CLASS::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  #if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  set_axis_aligned_bbox(((STEM*)fExternalObj)->AssertBBox());
#else
  SetAxisAlignedBBox(((STEM*)fExternalObj)->AssertBBox());
#endif

}

/**************************************************************************/

void CLASS::DirectDraw(const TGLDrawFlags & flags) const
{
  // printf("CLASS::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());
}
