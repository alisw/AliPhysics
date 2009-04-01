// $Id$
// Author: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveKinkEditor.h"
#include "AliEveKink.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"
#include "TMath.h"

//______________________________________________________________________________
// GUI editor for AliEveKink.
//

ClassImp(AliEveKinkEditor)

//______________________________________________________________________________
AliEveKinkEditor::AliEveKinkEditor(const TGWindow *p, Int_t width, Int_t height,
                               UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fInfoLabelDaughter(0),
  fXButton(0)
  // Initialize widget pointers to 0
{
  // Constructor.

  MakeTitle("AliEveKink");

  fInfoLabel0 = new TGLabel(this);
  fInfoLabel0->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fInfoLabel1 = new TGLabel(this);
  fInfoLabel1->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));
					  
  fInfoLabelDaughter = new TGLabel(this);
  fInfoLabelDaughter->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabelDaughter, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
						    8, 0, 2, 0));					  

  fXButton = new TGTextButton(this, "Detailed View");
  AddFrame(fXButton, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
  fXButton->Connect("Clicked()", "AliEveKinkEditor", this, "DisplayDetailed()");
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveKinkEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveKink*>(obj);

  // Set values of widgets
  fInfoLabel0->SetText(Form("Radius = %f, Kink Angle = %f", fM->GetKinkRadius(), (fM->GetKinkAngle(2))*TMath::RadToDeg() ));
  fInfoLabel1->SetText(Form("Mother Pt = %f", fM->GetKinkPMotherPerp()));
  fInfoLabelDaughter->SetText(Form("Daughter Prob= %f for Pdg= %d", fM->GetDaugMaxProbPid(), fM->GetDaugMaxProbPdg()));
 
}

/******************************************************************************/

#include <TEveManager.h>
#include <TEveWindow.h>
#include <TEveViewer.h>
#include <TEveScene.h>

#include <TGLCamera.h>
#include <TGLViewer.h>

#include <TLatex.h>
#include <TRootEmbeddedCanvas.h>

void AliEveKinkEditor::DisplayDetailed()
{
  TEveWindowSlot *slot = TEveWindow::CreateWindowMainFrame();
  TEveWindowPack *pack = slot->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->SetHorizontal();

  //
  // This part is for the bending plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *bpViewer = gEve->SpawnNewViewer("Kink bending plane View");
  TEveScene  *bpScene  = gEve->SpawnNewScene("Kink bending plane Scene");
  bpViewer->AddScene(bpScene);
  bpScene->AddElement(fM);
  bpViewer->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  bpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  TGLCamera& bpCam = bpViewer->GetGLViewer()->CurrentCamera();
  bpCam.SetExternalCenter(kTRUE);
  bpCam.SetCenterVec(fM->fRecKinkPosition.fX, fM->fRecKinkPosition.fY, fM->fRecKinkPosition.fZ);

  //
  // This part is for the decay plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *dpViewer = gEve->SpawnNewViewer("Kink decay plane View");
  TEveScene  *dpScene  = gEve->SpawnNewScene("Kink decay plane Scene");
  dpViewer->AddScene(dpScene);
  dpScene->AddElement(fM);
  dpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  TGLCamera& dpCam = dpViewer->GetGLViewer()->CurrentCamera();
  dpCam.SetExternalCenter(kTRUE);
  dpCam.SetCenterVec(fM->fRecKinkPosition.fX, fM->fRecKinkPosition.fY, fM->fRecKinkPosition.fZ);

  // This part is for displaying the information
  slot = pack->NewSlot();
  
  TEveWindowFrame *frame = slot->MakeFrame(new TRootEmbeddedCanvas());
  frame->SetElementName("Details");
  
  char info[100] = {0};
  TLatex* ltx = new TLatex(0.1, 0.9, info);
  ltx->SetTextSize(0.08);

  sprintf(info,"radius = %.3f [cm]",fM->GetKinkRadius());
  ltx->DrawLatex(0.1, 0.8, info);
  sprintf(info,"mass_{K^{+, -}} = %.3f [GeV/c^{2}]",fM->GetInvMass(13));
  ltx->DrawLatex(0.1, 0.3, info);

  gEve->Redraw3D();
}
