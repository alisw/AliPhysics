// $Header$

#ifndef REVE_NLTTrackGL_H
#define REVE_NLTTrackGL_H

#include <Reve/TrackGL.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class NLTTrack;

class NLTTrackGL : public TrackGL
{
private:
  NLTTrackGL(const NLTTrackGL&);            // Not implemented
  NLTTrackGL& operator=(const NLTTrackGL&); // Not implemented

protected:
  NLTTrack* fM; // fModel dynamic-casted to NLTTrackGL

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  NLTTrackGL();
  virtual ~NLTTrackGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);

  ClassDef(NLTTrackGL, 0);
}; // endclass NLTTrackGL

}

#endif
