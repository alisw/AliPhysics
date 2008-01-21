// $Header$

#include "CLASS.h"
#include <Alieve/STEM.h>

#include <TGLRnrCtx.h>
#include <TGLIncludes.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// CLASS
//

ClassImp(CLASS)

CLASS::CLASS() : TGLObject(), fM(0)
{
  // fDLCache = kFALSE; // Disable display list.
}

CLASS::~CLASS()
{}

/**************************************************************************/

Bool_t CLASS::SetModel(TObject* obj, const Option_t* /*opt*/)
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

void CLASS::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // printf("CLASS::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());
}
