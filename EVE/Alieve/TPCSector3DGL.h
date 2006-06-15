// $Header$

#ifndef ALIEVE_TPCSector3DGL_H
#define ALIEVE_TPCSector3DGL_H

#include <TGLObject.h>

#include <Alieve/TPCSectorData.h>

namespace Reve {
class BoxSetGL;
}

namespace Alieve {

class TPCSector3D;

class TPCSector3DGL : public TGLObject
{
protected:
  TPCSector3D*    fSector; // fModel dynamic-casted to TPCSector3DGL
  Reve::BoxSetGL* fBoxRnr;

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  TPCSector3DGL();
  virtual ~TPCSector3DGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  void DrawSegmentFrame(const TPCSectorData::SegmentInfo& s,
                        Int_t botExtraPads=0, Int_t topExtraPads=0) const;

  ClassDef(TPCSector3DGL, 0);
}; // endclass TPCSector3DGL

}

#endif
