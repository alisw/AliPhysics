// $Header$

#ifndef AliEveCLASS_H
#define AliEveCLASS_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

namespace Alieve {

class STEM;

class CLASS : public TGLObject
{
private:
  CLASS(const CLASS&);            // Not implemented
  CLASS& operator=(const CLASS&); // Not implemented

protected:
  STEM* fM; // fModel dynamic-casted to CLASS

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  CLASS();
  virtual ~CLASS();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  // virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  ClassDef(CLASS, 0);
}; // endclass CLASS

}

#endif
