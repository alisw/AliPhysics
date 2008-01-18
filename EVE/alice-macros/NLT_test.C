// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

TEveProjectionManager* NLT_test(TEveElement* top=0)
{
  TEveScene* s = gEve->SpawnNewScene("Projected AliEveEventManager");
  gEve->GetDefViewer()->AddScene(s);

  TGLViewer* v = (TGLViewer *)gEve->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);

  TEveProjectionManager* p = new TEveProjectionManager;
  p->SetProjection(TEveProjection::kPT_RhoZ, 0.01);

  gEve->AddToListTree(p, kTRUE);
  gEve->AddElement(p, s);

  top = gEve->GetCurrentEvent();
  if (top)
    p->ImportElements(top);

  gEve->Redraw3D(kTRUE);

  return p;
}
