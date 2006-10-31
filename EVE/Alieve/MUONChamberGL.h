#ifndef ALIEVE_MUONChamberGL_H
#define ALIEVE_MUONChamberGL_H

#include <TGLObject.h>

namespace Reve {
class QuadSetGL;
}

namespace Alieve {

class MUONChamber;

class MUONChamberGL : public TGLObject
{

  MUONChamberGL(const MUONChamberGL&);            // Not implemented
  MUONChamberGL& operator=(const MUONChamberGL&); // Not implemented

 protected:

  virtual void DirectDraw(const TGLDrawFlags & flags) const;
  void DrawChamberFrame() const;
  void DrawQuads() const;

  MUONChamber*             fChamber; // fModel dynamic-casted to MUONChamberGL
  mutable UInt_t           fRTS;     // render time stamp

 public:

  MUONChamberGL();
  virtual ~MUONChamberGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  ClassDef(MUONChamberGL,1);   // the GL drawing class of one chamber

};

}


#endif
