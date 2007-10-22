// $Header$

#ifndef REVE_Scene_H
#define REVE_Scene_H

#include <Reve/RenderElement.h>
#include <Reve/Pad.h>

class TGLScenePad;

namespace Reve {

/**************************************************************************/
// Scene
/**************************************************************************/

class Scene : public RenderElementList
{
private:
  Scene(const Scene&);            // Not implemented
  Scene& operator=(const Scene&); // Not implemented

protected:
  Pad         *fPad;
  TGLScenePad *fGLScene;

  Bool_t       fChanged;
  Bool_t       fSmartRefresh;

public:
  Scene(const Text_t* n="Scene", const Text_t* t="");
  virtual ~Scene();

  virtual void CollectSceneParents(List_t& scenes);

  void   Changed()    { fChanged = kTRUE; }
  Bool_t IsChanged() const { return fChanged;  }
  void   Repaint();

  TGLScenePad* GetGLScene() const { return fGLScene; }
  void SetGLScene(TGLScenePad* s) { fGLScene = s; }

  virtual void SetName(const Text_t* n);
  virtual void Paint(Option_t* option = "");

  virtual const TGPicture* GetListTreeIcon() { return RenderElement::fgListTreeIcons[2]; }
  ClassDef(Scene, 0);
}; // endclass Scene


/**************************************************************************/
// SceneList
/**************************************************************************/

class SceneList : public RenderElementList
{
private:
  SceneList(const SceneList&);            // Not implemented
  SceneList& operator=(const SceneList&); // Not implemented

protected:

public:
  SceneList(const Text_t* n="SceneList", const Text_t* t="");
  virtual ~SceneList();

  void RepaintChangedScenes();
  void RepaintAllScenes();

  ClassDef(SceneList, 0);
}; // endclass SceneList

}

#endif
