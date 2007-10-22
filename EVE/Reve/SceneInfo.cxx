// $Header$

#include "SceneInfo.h"
#include "Scene.h"
#include "ReveManager.h"

#include <TGLSceneInfo.h>

using namespace Reve;

//______________________________________________________________________
// SceneInfo
//

ClassImp(SceneInfo)

SceneInfo::SceneInfo(Viewer* viewer, Scene* scene, TGLSceneInfo* sinfo) :
  RenderElement (),
  TNamed        (Form("SI - %s", scene->GetName()),
		 Form("SceneInfo of scene '%s'", scene->GetName())),
  fViewer       (viewer),
  fScene        (scene),
  fGLSceneInfo  (sinfo)
{}

SceneInfo::~SceneInfo()
{}

/**************************************************************************/

TGLSceneBase* SceneInfo::GetGLScene() const
{
  return fGLSceneInfo->GetScene();
}

/**************************************************************************/

void SceneInfo::SetRnrSelf(Bool_t rnr)
{
  RenderElement::SetRnrSelf(rnr);
  fGLSceneInfo->SetActive(fRnrSelf);
}

void SceneInfo::SetRnrState(Bool_t rnr)
{
  RenderElement::SetRnrState(rnr);
  fGLSceneInfo->SetActive(fRnrSelf);
}

/**************************************************************************/

Bool_t SceneInfo::AcceptRenderElement(RenderElement* /*el*/)
{
  static const Exc_t eH("SceneInfo::AcceptRenderElement ");

  gReve->SetStatusLine(eH + "this class does not accept children.");
  return kFALSE;
}

Bool_t SceneInfo::HandleElementPaste(RenderElement* /*el*/)
{
  static const Exc_t eH("SceneInfo::HandleElementPaste ");

  gReve->SetStatusLine(eH + "this class does not accept children.");
  return kFALSE;
}
