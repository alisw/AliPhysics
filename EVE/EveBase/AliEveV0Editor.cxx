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
#include <TEveGeoNode.h>
#include <TEveProjectionManager.h>

#include <TGLCamera.h>
#include <TGLViewer.h>
#include "TGLCameraOverlay.h"

#include <TLatex.h>
#include <TRootEmbeddedCanvas.h>
#include <TInterpreter.h>

void AliEveV0Editor::DisplayDetailed()
{
  TEveWindowSlot *slot = TEveWindow::CreateWindowMainFrame();
  TEveWindowPack *pack = slot->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->SetHorizontal();

  //
  // This part is for getting the different objects to display
  //
  char displayInfo[100] = {0};
  sprintf(displayInfo,"pt = %.3f",fM->GetPt());
  TEveLine *lv0TransverseMomentumDirection = new TEveLine(displayInfo);
  lv0TransverseMomentumDirection->SetLineColor(kOrange+8);
  lv0TransverseMomentumDirection->SetLineWidth(2);
  lv0TransverseMomentumDirection->SetLineStyle(2);
  lv0TransverseMomentumDirection->SetLineWidth(2);
  Float_t scalePt = 100.; // this needs to be available as a ruler
  lv0TransverseMomentumDirection->SetPoint(0,fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);
  lv0TransverseMomentumDirection->SetPoint(1,scalePt*fM->fRecDecayP.fX, scalePt*fM->fRecDecayP.fY,0);

  TEvePointSet *pvlocation = new TEvePointSet("PV location");
  pvlocation->SetNextPoint(fM->fRecBirthV.fX, fM->fRecBirthV.fY, fM->fRecBirthV.fZ);
  pvlocation->SetTitle("pv location");
  pvlocation->SetMarkerStyle(4);
  pvlocation->SetMarkerSize(2.5);
  pvlocation->SetMarkerColor(7);

  TEvePointSet *v0location = new TEvePointSet("V0 location");
  v0location->SetNextPoint(fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);
  v0location->SetTitle("v0 location");
  v0location->SetMarkerStyle(4);
  v0location->SetMarkerSize(2.5);
  v0location->SetMarkerColor(kOrange+8);

  TEveUtil::LoadMacro("clusters_from_index.C");

  TEveTrack *negTrack = fM->GetNegTrack();
  TEveTrack *posTrack = fM->GetPosTrack();


  char macroWithIndex[100] = {0};
  Int_t daughterIndex = 0;
  TEvePointSet *negDaughterCluster = 0;
  TEvePointSet *posDaughterCluster = 0;
  
  daughterIndex = negTrack->GetIndex();
  sprintf(macroWithIndex,"clusters_from_index(%d)",daughterIndex);
  Long_t negResult = gInterpreter->ProcessLine(macroWithIndex);
  if (negResult) {
    negDaughterCluster = reinterpret_cast<TEvePointSet*>(negResult);
    if (negDaughterCluster){
      negDaughterCluster->SetMarkerStyle(4);	
      negDaughterCluster->SetMarkerSize(1.5);	
      negDaughterCluster->SetMarkerColor(kBlue+3);
    }
  }
  else
  {
    Warning("DisplayDetailed", "Import of negative daughter's clusters failed.");
  }

  daughterIndex = posTrack->GetIndex();
  sprintf(macroWithIndex,"clusters_from_index(%d)",daughterIndex);
  Long_t posResult = gInterpreter->ProcessLine(macroWithIndex);
  if (posResult) {
    posDaughterCluster = reinterpret_cast<TEvePointSet*>(posResult);
    if (posDaughterCluster){
      posDaughterCluster->SetMarkerStyle(4);	
      posDaughterCluster->SetMarkerSize(1.5);	
      posDaughterCluster->SetMarkerColor(kRed+3);
    }
  }
  else
  {
    Warning("DisplayDetailed", "Import of positive daughter's clusters failed.");
  }

  //
  // This part is for the bending plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *bpViewer = gEve->SpawnNewViewer("V0 bending plane View");
  TEveScene  *bpScene  = gEve->SpawnNewScene("V0 bending plane Scene");

  TEveUtil::LoadMacro("geom_gentle.C");
  Long_t result = gInterpreter->ProcessLine("geom_gentle_rphi()");
  if (result)
  {
    TEveGeoShape *geomRPhi = reinterpret_cast<TEveGeoShape*>(result);
    geomRPhi->IncDenyDestroy();
    TEveProjectionManager *projMgr = new TEveProjectionManager();
    projMgr->ImportElements(geomRPhi, bpScene);
  }
  else
  {
    Warning("DisplayDetailed", "Import of R-Phi geometry failed.");
  }

  bpViewer->AddScene(bpScene);
  bpScene->AddElement(fM);
  bpScene->AddElement(lv0TransverseMomentumDirection);
  bpScene->AddElement(pvlocation);
  bpScene->AddElement(v0location);

  if (negDaughterCluster) bpScene->AddElement(negDaughterCluster);
  if (posDaughterCluster) bpScene->AddElement(posDaughterCluster);

  // This is the to-do list for the bending plane:
  // 1. fix the view to orthographic XOY (no rotation allowed but moving the center ok) ->done! 
  // 2. show axis and tickles along X and Y ->done!
  //       -> note for the projection the cartesian scales are not very useful
  //       -> propose a phi and R scale which rotates with a reset at 0;
  //       -> propose a transformation for an eta scale (keep the z one);
  // 3. show the center, the main vertex and the detectors for this view ->done!
  // 4. show V0 direction in the bending plane with arrow length proportional to pT ->done!
  // 5. show angles with respect to axis (phi angle) ->almost.
  // 6. show clusters in the ITS and in the TPC associated with the daughter tracks
  //       -> include a radius cut for plotting only ITS and TPC ->done!
  bpViewer->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  bpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  TGLViewer *lbpGLViewer = bpViewer->GetGLViewer();
  TGLCameraOverlay* co = lbpGLViewer->GetCameraOverlay();
  co->SetShowOrthographic(true); //(false);
  co->SetOrthographicMode(TGLCameraOverlay::kAxis); // ::kPlaneIntersect or ::kBar
  // end of the bending plane part

  //
  // This part is for the decay plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *dpViewer = gEve->SpawnNewViewer("V0 decay plane View");
  TEveScene  *dpScene  = gEve->SpawnNewScene("V0 decay plane Scene");

  dpViewer->AddScene(dpScene);

  result = gInterpreter->ProcessLine("geom_gentle(kFALSE)");
  if (result)
  {
    TEveGeoShape *geom = reinterpret_cast<TEveGeoShape*>(result);
    geom->IncDenyDestroy();
    geom->FindChild("TRD+TOF")->SetRnrState(kFALSE);
    geom->FindChild("PHOS")   ->SetRnrState(kFALSE);
    geom->FindChild("HMPID")  ->SetRnrState(kFALSE);
    dpScene->AddElement(geom);
  }
  else
  {
    Warning("DisplayDetailed", "Import of 3D geometry failed.");
  }

  dpScene->AddElement(fM);
  dpScene->AddElement(lv0TransverseMomentumDirection);
  dpScene->AddElement(pvlocation);
  dpScene->AddElement(v0location);
  if (negDaughterCluster) dpScene->AddElement(negDaughterCluster);
  if (posDaughterCluster) dpScene->AddElement(posDaughterCluster);

  // This is the to-do list for the decay plane:
  // 1. fix the view to decay plane (no rotation allowed but moving the center ok)
  // 2. show V0 direction with a vertical arrow length proportional to pT -> done!
  // 3. show the center, the main vertex and the detectors for this view -> done!
  // 4. show the x,y and z axis and the different angles
  //       -> this needs a referential object that we can move around
  //          or fix to a selected point (origin being the default)
  // 5. draw the dca between daughters and the extrapolation to the main vertex.
  //       -> this is an issue since we only store the distance: check with J.Belikov
  // 6. show clusters in the ITS and in the TPC associated with the daughter tracks
  //       -> include a radius cut for plotting only ITS and TPC ->done!
  dpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
  TGLCamera& dpCam = dpViewer->GetGLViewer()->CurrentCamera();
  dpCam.SetExternalCenter(kTRUE);
  dpCam.SetCenterVec(fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);
  dpCam.RotateRad(0,-TMath::Pi()/2.); // RotateRad rotates in radians (hRotate,vRotate)
  //  Here rotate the _view_ (not the camera) by (fM->GetPhi() - TMath::Pi()/2.)

  // In the end maybe truck and rotate properly...
  //  dpCam.Truck(0,200);// Truck center wrt the view panel: (x=0 pixels, y pixels)
  //  dpCam.Rotate(0,50,0,0); // Rotate in pixels (xDelta,yDelta)

  // end of the decay plane part

  //
  // This part is for displaying the information
  //
  slot = pack->NewSlot();
  
  TEveWindowFrame *frame = slot->MakeFrame(new TRootEmbeddedCanvas());
  frame->SetElementName("Details");

  // Print and show detailed information about the V0
  // Calculation of the invariant mass with the max prob PID hypothesis first
  // pseudorapidity, phi angle, pt, radius, dcas
  char info[100] = {0};
  sprintf(info,"#phi = %.3frad = %.1f°",fM->GetPhi(),(180./TMath::Pi())*fM->GetPhi());
  TLatex* ltx = new TLatex(0.05, 0.9, info);
  ltx->SetTextSize(0.08);
  ltx->Draw();

  sprintf(info,"radius = %.3f [cm]",fM->GetRadius());
  ltx->DrawLatex(0.05, 0.8, info);

  sprintf(info,"p_{T} = %.3f [GeV/c]",fM->GetPt());
  ltx->DrawLatex(0.05, 0.7, info);

  sprintf(info,"daughters dca = %.3f [cm]",fM->GetDaughterDCA());
  ltx->DrawLatex(0.05, 0.6, info);

  sprintf(info,"#eta = - ln( tan(#theta/2) ) = %.3f",fM->GetEta());
  ltx->DrawLatex(0.05, 0.5, info);

  sprintf(info,"mass_{K^{0}_{s}} = %.3f [GeV/c^{2}]",fM->GetK0sInvMass());
  ltx->DrawLatex(0.05, 0.3, info);

  sprintf(info,"mass_{#Lambda} = %.3f [GeV/c^{2}]",fM->GetLambdaInvMass());
  ltx->DrawLatex(0.05, 0.2, info);

  sprintf(info,"mass_{#bar{#Lambda}} = %.3f [GeV/c^{2}]",fM->GetAntiLambdaInvMass());
  ltx->DrawLatex(0.05, 0.1, info);

  gEve->Redraw3D();
}
