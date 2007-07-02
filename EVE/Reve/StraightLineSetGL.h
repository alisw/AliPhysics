// $Header$

#ifndef REVE_StraightLineSetGL_H
#define REVE_StraightLineSetGL_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class StraightLineSet;

class StraightLineSetGL : public TGLObject
{
private:
  StraightLineSetGL(const StraightLineSetGL&);            // Not implemented
  StraightLineSetGL& operator=(const StraightLineSetGL&); // Not implemented

protected:
  StraightLineSet* fM; // fModel dynamic-casted to StraightLineSetGL

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  StraightLineSetGL();
  virtual ~StraightLineSetGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  // To support two-level selectionvirtual 
  Bool_t IgnoreSizeForOfInterest() const { return kTRUE; }

  virtual Bool_t ShouldCache(TGLRnrCtx & rnrCtx) const;
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  ClassDef(StraightLineSetGL, 0);
}; // endclass StraightLineSetGL

}

#endif
