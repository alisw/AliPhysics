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
  fInfoLabelNegDaughter(0),
  fInfoLabelPosDaughter(0),
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

  fInfoLabelNegDaughter = new TGLabel(this);
  fInfoLabelNegDaughter->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabelNegDaughter, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
						    8, 0, 2, 0));

  fInfoLabelPosDaughter = new TGLabel(this);
  fInfoLabelPosDaughter->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabelPosDaughter, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
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
  fInfoLabelNegDaughter->SetText(Form("Neg. Daughter Prob= %.2f for Pdg= %d", fM->GetNegMaxProbPid(), fM->GetNegMaxProbPdg()));
  fInfoLabelPosDaughter->SetText(Form("Pos. Daughter Prob= %.2f for Pdg= %d", fM->GetPosMaxProbPid(), fM->GetPosMaxProbPdg()));

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
  TEveWindowSlot *slot = TEveWindow::CreateWindowMainFrame();
  TEveWindowPack *pack = slot->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->SetHorizontal();

  //
  // This part is for the bending plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *bpViewer = gEve->SpawnNewViewer("V0 bending plane View");
  TEveScene  *bpScene  = gEve->SpawnNewScene("V0 bending plane Scene");
  bpViewer->AddScene(bpScene);
  bpScene->AddElement(fM);
  // This is the to-do list for the bending plane:
  // 1. fix the view to orthographic XOY (no rotation allowed but moving the center ok) --tbu 
  // 2. show axis and tickles along X and Y
  // 3. show the center, the main vertex and the detectors for this view;
  // 4. show V0 direction in the bending plane with arrow length proportional to pT;
  // 5. show angles with respect to axis (phi angle).
  bpViewer->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  //  bpViewer->GetGLViewer()->GetOrthoXOYCamera()->SetEnableRotate(kFALSE);
  bpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  //  bpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  TGLCamera& bpCam = bpViewer->GetGLViewer()->CurrentCamera();
  bpCam.SetExternalCenter(kTRUE);
  //  bpCam.SetEnableRotate(kFALSE);
  bpCam.SetCenterVec(fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);
  // end of the bending plane part

  //
  // This part is for the decay plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *dpViewer = gEve->SpawnNewViewer("V0 decay plane View");
  TEveScene  *dpScene  = gEve->SpawnNewScene("V0 decay plane Scene");
  dpViewer->AddScene(dpScene);
  dpScene->AddElement(fM);
  // This is the to-do list for the decay plane:
  // 1. fix the view to decay plane (no rotation allowed but moving the center ok)
  // 2. show V0 direction with a vertical arrow length proportional to pT;
  // 3. show the center, the main vertex and the detectors for this view;
  // 4. show the x,y and z axis and the different angles;
  // 5. draw the dca between daughters and the extrapolation to the main vertex.
  dpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  TGLCamera& dpCam = dpViewer->GetGLViewer()->CurrentCamera();
  dpCam.SetExternalCenter(kTRUE);
  dpCam.SetCenterVec(fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);
  // end of the decay plane part


  // This part is for displaying the information
  slot = pack->NewSlot();
  
  TEveWindowFrame *frame = slot->MakeFrame(new TRootEmbeddedCanvas());
  frame->SetElementName("Details");

  // Print and show detailed information about the V0
  // Calculation of the invariant mass with the max prob PID hypothesis first
  // pseudorapidity, phi angle, pt, radius, dcas
  char info[100] = {0};
  sprintf(info,"#phi = %.3f",fM->GetPhi());
  TLatex* ltx = new TLatex(0.1, 0.9, info);
  ltx->SetTextSize(0.08);
  ltx->Draw();

  sprintf(info,"radius = %.3f [cm]",fM->GetRadius());
  ltx->DrawLatex(0.1, 0.8, info);

  sprintf(info,"daughters dca = %.3f [cm]",fM->GetDaughterDCA());
  ltx->DrawLatex(0.1, 0.7, info);

  sprintf(info,"#eta = #frac{1}{2} #times Ln(#frac{E+p_{z}}{E-p_{z}} ) = %.3f",fM->GetEta());
  ltx->DrawLatex(0.1, 0.6, info);

  sprintf(info,"mass_{K^{0}_{s}} = %.3f [GeV/c^{2}]",fM->GetK0sInvMass());
  ltx->DrawLatex(0.1, 0.3, info);

  sprintf(info,"mass_{#Lambda} = %.3f [GeV/c^{2}]",fM->GetLambdaInvMass());
  ltx->DrawLatex(0.1, 0.2, info);

  sprintf(info,"mass_{#bar{#Lambda}} = %.3f [GeV/c^{2}]",fM->GetAntiLambdaInvMass());
  ltx->DrawLatex(0.1, 0.1, info);

  gEve->Redraw3D();
}
