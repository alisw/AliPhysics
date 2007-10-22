// $Header$

#include "Scene.h"
#include "Viewer.h"
#include <Reve/ReveManager.h>

#include <TList.h>
#include <TGLScenePad.h>

using namespace Reve;

//______________________________________________________________________
// Reve::Scene
//

ClassImp(Scene)

Scene::Scene(const Text_t* n, const Text_t* t) :
  RenderElementList(n, t),
  fPad    (0),
  fGLScene(0),
  fChanged      (kFALSE),
  fSmartRefresh (kTRUE)
{
  fPad = new Pad;
  fPad->GetListOfPrimitives()->Add(this);
  fGLScene = new TGLScenePad(fPad);
  fGLScene->SetName(n);
  fGLScene->SetAutoDestruct(kFALSE);
}

Scene::~Scene()
{
  gReve->GetViewers()->SceneDestructing(this);
}

/**************************************************************************/

void Scene::CollectSceneParents(List_t& scenes)
{
  scenes.push_back(this);
}

/**************************************************************************/

void Scene::Repaint()
{
  fGLScene->PadPaint(fPad);
  fChanged = kFALSE;
}

/**************************************************************************/

void Scene::SetName(const Text_t* n)
{
  RenderElementList::SetName(n);
  fGLScene->SetName(n);
}

void Scene::Paint(Option_t* option)
{
  if (fRnrChildren)
  {
    for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
      (*i)->PadPaint(option);
  }
}


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// Reve::SceneList
//

ClassImp(SceneList)

SceneList::SceneList(const Text_t* n, const Text_t* t) :
  RenderElementList(n, t)
{
  SetChildClass(Scene::Class());
}

SceneList::~SceneList()
{}

/**************************************************************************/

void SceneList::RepaintChangedScenes()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
  {
    Scene* s = (Scene*) *i;
    if (s->IsChanged())
    {
      // printf(" Scene '%s' changed ... repainting.\n", s->GetName());
      s->Repaint();
    }
  }
}

void SceneList::RepaintAllScenes()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
  {
    Scene* s = (Scene*) *i;
    // printf(" Scene '%s' repainting.\n", s->GetName());
    s->Repaint();
  }
}
