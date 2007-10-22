// $Header$

#ifndef REVE_Viewer_H
#define REVE_Viewer_H

#include <Reve/RenderElement.h>

class TGWindow;
class TGedEditor;
class TGLViewer;

namespace Reve {

class Scene;

/**************************************************************************/
// Viewer
/**************************************************************************/

class Viewer : public RenderElementList
{
private:
  Viewer(const Viewer&);            // Not implemented
  Viewer& operator=(const Viewer&); // Not implemented

protected:
  TGLViewer *fGLViewer;

public:
  Viewer(const Text_t* n="Viewer", const Text_t* t="");
  virtual ~Viewer();

  TGLViewer* GetGLViewer() const { return fGLViewer; }
  void SetGLViewer(TGLViewer* s);
  void SpawnGLViewer(const TGWindow* parent, TGedEditor* ged);

  virtual void AddScene(Scene* scene);

  virtual void RemoveElementLocal(RenderElement* el);
  virtual void RemoveElementsLocal();

  virtual TObject* GetEditorObject() const;

  virtual Bool_t HandleElementPaste(RenderElement* el);

  virtual const TGPicture* GetListTreeIcon() { return RenderElement::fgListTreeIcons[1]; }

  ClassDef(Viewer, 0);
}; // endclass Viewer


/**************************************************************************/
// ViewerList
/**************************************************************************/

class ViewerList : public RenderElementList
{
private:
  ViewerList(const ViewerList&);            // Not implemented
  ViewerList& operator=(const ViewerList&); // Not implemented

protected:

public:
  ViewerList(const Text_t* n="ViewerList", const Text_t* t="");
  virtual ~ViewerList();

  void RepaintChangedViewers(Bool_t resetCameras, Bool_t dropLogicals);
  void RepaintAllViewers(Bool_t resetCameras, Bool_t dropLogicals);

  void SceneDestructing(Scene* scene);

  ClassDef(ViewerList, 0);
}; // endclass ViewerList

}

#endif
