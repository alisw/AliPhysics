// $Header$

#ifndef REVE_SceneInfo_H
#define REVE_SceneInfo_H

#include <Reve/RenderElement.h>

class TGLSceneBase;
class TGLSceneInfo;

namespace Reve {

class Viewer;
class Scene;

class SceneInfo : public RenderElement,
		  public TNamed
{
private:
  SceneInfo(const SceneInfo&);            // Not implemented
  SceneInfo& operator=(const SceneInfo&); // Not implemented

protected:
  Viewer       *fViewer;
  Scene        *fScene;
  TGLSceneInfo *fGLSceneInfo;

public:
  SceneInfo(Viewer* viewer, Scene* scene, TGLSceneInfo* sinfo);
  virtual ~SceneInfo();

  Viewer       * GetViewer()      const { return fViewer; }
  Scene        * GetScene()       const { return fScene;  }
  TGLSceneInfo * GetGLSceneInfo() const { return fGLSceneInfo; }
  TGLSceneBase * GetGLScene()     const;

  virtual void   SetRnrSelf(Bool_t rnr);
  virtual void   SetRnrState(Bool_t rnr);

  virtual Bool_t AcceptRenderElement(RenderElement* el);
  virtual Bool_t HandleElementPaste(RenderElement* el);

  ClassDef(SceneInfo, 0); // Reve representation of TGLSceneInfo.
}; // endclass SceneInfo

}

#endif
