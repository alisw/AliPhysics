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
  TPCSector3DGL(const TPCSector3DGL&);            // Not implemented
  TPCSector3DGL& operator=(const TPCSector3DGL&); // Not implemented

protected:
  TPCSector3D*    fSector; // fModel dynamic-casted to TPCSector3DGL
  Reve::BoxSetGL* fBoxRnr;

  mutable UInt_t  fRTS;

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  TPCSector3DGL();
  virtual ~TPCSector3DGL();

  virtual Bool_t   ShouldDLCache(const TGLRnrCtx&) const { return kFALSE; }
  virtual ELODAxes SupportedLODAxes()              const { return kLODAxesAll; }
  virtual Short_t  QuantizeShapeLOD(Short_t shapeLOD, Short_t combiLOD) const;

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  void DrawSegmentFrame(const TPCSectorData::SegmentInfo& s,
                        Int_t botExtraPads=0, Int_t topExtraPads=0) const;

  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void   ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  ClassDef(TPCSector3DGL, 0);
}; // endclass TPCSector3DGL

}

#endif
