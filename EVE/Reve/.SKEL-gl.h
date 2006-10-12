// $Header$

#ifndef REVE_CLASS_H
#define REVE_CLASS_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class STEM;

class CLASS : public TGLObject
{
private:
  CLASS(const CLASS&);            // Not implemented
  CLASS& operator=(const CLASS&); // Not implemented

protected:
  STEM* fM; // fModel dynamic-casted to CLASS

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  CLASS();
  virtual ~CLASS();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  // virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

  ClassDef(CLASS, 0);
}; // endclass CLASS

}

#endif
