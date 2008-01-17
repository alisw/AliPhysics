#ifndef ALIEVE_MUONChamberGL_H
#define ALIEVE_MUONChamberGL_H

#include <TGLObject.h>
#include <TEveQuadSetGL.h>

class TEveQuadSetGL;

namespace Alieve {

class MUONChamber;

class MUONChamberGL : public TGLObject
{

  MUONChamberGL(const MUONChamberGL&);            // Not implemented
  MUONChamberGL& operator=(const MUONChamberGL&); // Not implemented

 protected:

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;
  void DrawChamberFrame() const;
  void DrawQuads(TGLRnrCtx& rnrCtx) const;
  void DrawPoints() const;

  MUONChamber*             fChamber; // fModel dynamic-casted to MUONChamberGL
  TEveQuadSetGL            fQS1;
  TEveQuadSetGL            fQS2;

  mutable UInt_t           fRTS;     // render time stamp

 public:

  MUONChamberGL();
  virtual ~MUONChamberGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  ClassDef(MUONChamberGL,1);   // the GL drawing class of one chamber

};

}


#endif
