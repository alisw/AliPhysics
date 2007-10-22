// $Header$

#ifndef REVE_NLTPointSetGL_H
#define REVE_NLTPointSetGL_H

#include <TPointSet3DGL.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class NLTPointSet;

class NLTPointSetGL : public TPointSet3DGL
{
private:
  NLTPointSetGL(const NLTPointSetGL&);            // Not implemented
  NLTPointSetGL& operator=(const NLTPointSetGL&); // Not implemented

protected:

public:
  NLTPointSetGL();
  virtual ~NLTPointSetGL();

  ClassDef(NLTPointSetGL, 0);
}; // endclass NLTPointSetGL

}

#endif
