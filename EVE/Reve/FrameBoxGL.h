// $Header$

#ifndef REVE_FrameBoxGL_H
#define REVE_FrameBoxGL_H

#include <Reve/Reve.h>

namespace Reve {

class FrameBox;

class FrameBoxGL
{
private:
  FrameBoxGL();                             // Not implemented
  FrameBoxGL(const FrameBoxGL&);            // Not implemented
  FrameBoxGL& operator=(const FrameBoxGL&); // Not implemented

  static void RenderFrame(const FrameBox& b, Bool_t fillp);

public:
  virtual ~FrameBoxGL() {}

  static void Render(const FrameBox* box);

  ClassDef(FrameBoxGL, 0);
}; // endclass FrameBoxGL

}

#endif
