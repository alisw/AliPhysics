// @(#)root/eve:$Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHFEditor.h"
#include "AliEveHF.h"

#include "TVirtualPad.h"
#include "TColor.h"
#include <TDatabasePDG.h>

#include <TEveTrack.h>

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"


//______________________________________________________________________________
// GUI editor for AliEveHF.
//

ClassImp(AliEveHFEditor)

//______________________________________________________________________________
AliEveHFEditor::AliEveHFEditor(const TGWindow *p, Int_t width, Int_t height,
                               UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fnProng(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fInfoLabel2(0),
  fXButton(0)
  // Initialize widget pointers to 0
{
  // Constructor.

  MakeTitle("AliEveHF");

  fInfoLabel0 = new TGLabel(this);
  fInfoLabel0->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fInfoLabel1 = new TGLabel(this);
  fInfoLabel1->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fInfoLabel2 = new TGLabel(this);
  fInfoLabel2->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel2, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fXButton = new TGTextButton(this, "Detailed View");
  AddFrame(fXButton, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
  fXButton->Connect("Clicked()", "AliEveHFEditor", this, "DisplayDetailed()");
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveHFEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = static_cast<AliEveHF*>(obj);

  // Set values of widgets
  fInfoLabel0->SetText(Form("CosPointingAngle = %f",  fM->GetCosPointingAngle()));
  fInfoLabel1->SetText(Form("Pt = %f", fM->GetPt()));

  fInfoLabel2->SetText(Form("PDG Invariant Mass = %f", TDatabasePDG::Instance()->GetParticle(421)->Mass()));

}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveHFEditor::DoXYZZ()
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

void AliEveHFEditor::DisplayDetailed()
{
  TEveWindowSlot *slot = TEveWindow::CreateWindowMainFrame();
  TEveWindowPack *pack = slot->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->SetHorizontal();

  //
  // This part is for getting the different objects to display
  //
  char displayInfo[100] = {0};
  snprintf(displayInfo,100,"pt = %.3f",fM->GetPt());
  TEveLine *lhfTransverseMomentumDirection = new TEveLine(displayInfo);
  lhfTransverseMomentumDirection->SetLineColor(kOrange+8);
  lhfTransverseMomentumDirection->SetLineWidth(2);
  lhfTransverseMomentumDirection->SetLineStyle(2);
  lhfTransverseMomentumDirection->SetLineWidth(2);
  Float_t scalePt = 100.; // this needs to be available as a ruler
  lhfTransverseMomentumDirection->SetPoint(0,fM->fRecDecayHF.fX, fM->fRecDecayHF.fY, fM->fRecDecayHF.fZ);
  lhfTransverseMomentumDirection->SetPoint(1,scalePt*fM->fRecDecayMomHF.fX, scalePt*fM->fRecDecayMomHF.fY,0);

  TEvePointSet *pvlocation = new TEvePointSet("PV location");
  pvlocation->SetNextPoint(fM->fRecBirthHF.fX, fM->fRecBirthHF.fY, fM->fRecBirthHF.fZ);
  pvlocation->SetTitle("pv location");
  pvlocation->SetMarkerStyle(4);
  pvlocation->SetMarkerSize(2.5);
  pvlocation->SetMarkerColor(7);

  TEvePointSet *hflocation = new TEvePointSet("HF location");
  hflocation->SetNextPoint(fM->fRecDecayHF.fX, fM->fRecDecayHF.fY, fM->fRecDecayHF.fZ);
  hflocation->SetTitle("HF location");
  hflocation->SetMarkerStyle(4);
  hflocation->SetMarkerSize(2.5);
  hflocation->SetMarkerColor(kOrange+8);

  TEveTrack *negTrack =  fM->GetNegTrack();
  TEveTrack *posTrack =  fM->GetPosTrack();


  char macroWithIndex[100] = {0};
  Int_t daughterIndex = 0;
  TEvePointSet *negDaughterCluster = 0;
  TEvePointSet *posDaughterCluster = 0;

  daughterIndex = negTrack->GetIndex();
    
    AliEveTrack *eveTrack = new AliEveTrack();
    
    Long_t negResult = eveTrack->ImportClustersFromIndex(daughterIndex);
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

    Long_t posResult = eveTrack->ImportClustersFromIndex(daughterIndex);
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
    
    delete eveTrack;eveTrack=0;
    

  //
  // This part is for the bending plane view
  //
  pack->NewSlot()->MakeCurrent();
  TEveViewer *bpViewer = gEve->SpawnNewViewer("HF bending plane View");
  TEveScene  *bpScene  = gEve->SpawnNewScene("HF bending plane Scene");

  TEveUtil::LoadMacro("geom_gentle.C");
  Long_t result = gInterpreter->ProcessLine("geom_gentle_rphi()");
  if (result)
  {
    TEveGeoShape *geomRPhi = reinterpret_cast<TEveGeoShape*>(result);
    geomRPhi->IncDenyDestroy();
    TEveProjectionManager *projMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    projMgr->ImportElements(geomRPhi, bpScene);
  }
  else
  {
    Warning("DisplayDetailed", "Import of R-Phi geometry failed.");
  }

  bpViewer->AddScene(bpScene);
  bpScene->AddElement(fM);
  bpScene->AddElement(lhfTransverseMomentumDirection);
  bpScene->AddElement(pvlocation);
  bpScene->AddElement(hflocation);

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
  TEveViewer *dpViewer = gEve->SpawnNewViewer("HF decay plane View");
  TEveScene  *dpScene  = gEve->SpawnNewScene("HF decay plane Scene");

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
  dpScene->AddElement(lhfTransverseMomentumDirection);
  dpScene->AddElement(pvlocation);
  dpScene->AddElement(hflocation);
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
  dpCam.SetCenterVec(fM->fRecDecayHF.fX, fM->fRecDecayHF.fY, fM->fRecDecayHF.fZ);
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

  // Print and show detailed information about the HF
  // Calculation of the invariant mass with the max prob PID hypothesis first
  // pseudorapidity, phi angle, pt, radius, dcas
  char info[100] = {0};

  snprintf(info,100,"Cos Pointing angle = %.3f ",fM->GetCosPointingAngle());
  TLatex* ltx = new TLatex(0.05, 0.9, info);
  ltx->SetTextSize(0.08);
  ltx->DrawLatex(0.05, 0.8, info);

  snprintf(info,100,"p_{T} = %.3f [GeV/c]",fM->GetPt());
  ltx->DrawLatex(0.05, 0.7, info);

  snprintf(info,100,"#eta = - ln( tan(#theta/2) ) = %.3f",fM->GetEta());
  ltx->DrawLatex(0.05, 0.6, info);

  snprintf(info,100,"D^{0} inv. mass = %.3f", fM->GetInvariantMassPart());
  ltx->DrawLatex(0.05, 0.5, info);

  snprintf(info,100,"D^{0}bar inv. mass = %.3f", fM->GetInvariantMassAntiPart());
  ltx->DrawLatex(0.05, 0.4, info);

  gEve->Redraw3D();
}
