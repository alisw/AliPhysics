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

CLASS::CLASS() : TGLObject(), fM(0)
{
  // fCached = false; // Disable display list.
}

CLASS::~CLASS()
{}

/**************************************************************************/

Bool_t CLASS::SetModel(TObject* obj)
{
  if(SetModelCheckClass(obj, STEM::Class())) {
    fM = dynamic_cast<STEM*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

void CLASS::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  SetAxisAlignedBBox(((STEM*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void CLASS::DirectDraw(const TGLDrawFlags & flags) const
{
  // printf("CLASS::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());
}
