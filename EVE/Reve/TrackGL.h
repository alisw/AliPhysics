// $Header$

#ifndef REVE_TrackGL_H
#define REVE_TrackGL_H

#include <Reve/LineGL.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class Track;

class TrackGL : public LineGL
{
private:
  TrackGL(const TrackGL&);            // Not implemented
  TrackGL& operator=(const TrackGL&); // Not implemented

protected:

public:
  TrackGL();
  virtual ~TrackGL();

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

  ClassDef(TrackGL, 0);
}; // endclass TrackGL

}

#endif
