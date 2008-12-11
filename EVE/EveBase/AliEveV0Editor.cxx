// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveV0Editor.h"
#include "AliEveV0.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"


//______________________________________________________________________________
// GUI editor for AliEveV0.
//

ClassImp(AliEveV0Editor)

//______________________________________________________________________________
AliEveV0Editor::AliEveV0Editor(const TGWindow *p, Int_t width, Int_t height,
                               UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fXButton(0)
  // Initialize widget pointers to 0
{
  // Constructor.

  MakeTitle("AliEveV0");

  fInfoLabel0 = new TGLabel(this);
  fInfoLabel0->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fInfoLabel1 = new TGLabel(this);
  fInfoLabel1->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fXButton = new TGTextButton(this, "Detailed View");
  AddFrame(fXButton, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
  fXButton->Connect("Clicked()", "AliEveV0Editor", this, "DisplayDetailed()");
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0Editor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveV0*>(obj);

  // Set values of widgets
  fInfoLabel0->SetText(Form("Radius = %f, DCA = %f", fM->GetRadius(), fM->GetDaughterDCA()));
  fInfoLabel1->SetText(Form("Pt = %f", fM->GetPt()));
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveV0Editor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

#include <TEveManager.h>
#include <TEveWindow.h>
#include <TEveViewer.h>
#include <TEveScene.h>

#include <TGLCamera.h>
#include <TGLViewer.h>

#include <TLatex.h>
#include <TRootEmbeddedCanvas.h>

void AliEveV0Editor::DisplayDetailed()
{
  printf("\n Checking invariant mass calculations: K0s %.3f lambda %.3f antilambda %.3f \n",fM->GetK0sInvMass(),fM->GetLambdaInvMass(),fM->GetAntiLambdaInvMass());

  TEveWindowSlot *slot = TEveWindow::CreateWindowMainFrame();
  TEveWindowPack *pack = slot->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();

  TEveViewer *viewer = gEve->SpawnNewViewer("V0 Detailed View");
  TEveScene  *scene  = gEve->SpawnNewScene("V0 Detail Scene");
  viewer->AddScene(scene);
  scene->AddElement(fM);

  viewer->GetGLViewer()->ResetCamerasAfterNextUpdate();

  TGLCamera& cam = viewer->GetGLViewer()->CurrentCamera();
  cam.SetExternalCenter(kTRUE);
  cam.SetCenterVec(fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);

  slot = pack->NewSlot();
  
  TEveWindowFrame *frame = slot->MakeFrame(new TRootEmbeddedCanvas());
  frame->SetElementName("Details");

  TLatex* ltx = new TLatex(0.2, 0.2, "#eta = #frac{1}{2} #times Ln(#frac{E+p_{z}}{E-p_{z}} )");
  ltx->SetTextSize(0.08);
  ltx->Draw();

  gEve->Redraw3D();
}
