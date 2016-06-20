// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveCascadeEditor.h"
#include "AliEveCascade.h"
#include "AliEveTrack.h"

//#include "TVirtualPad.h"
//#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
//#include "TGNumberEntry.h"
//#include "TGColorSelect.h"
//#include "TGDoubleSlider.h"


//______________________________________________________________________________
// GUI editor for AliEveCascade.
//

ClassImp(AliEveCascadeEditor)

//______________________________________________________________________________
AliEveCascadeEditor::AliEveCascadeEditor(const TGWindow *p, Int_t width, Int_t height,
                                         UInt_t options, Pixel_t back) :
TGedFrame(p, width, height, options | kVerticalFrame, back),
fM(0),
fInfoLabelRadius(0),
fInfoLabelDCA(0),
fInfoLabelCharge(0),
fInfoLabelPhi(0),
fInfoLabelTheta(0),
fInfoLabelPtot(0),
fInfoLabelPt(0),
fInfoLabelEta(0),
fXButtonDetailedView(0),
fXButtonMassHyp(0)
// Initialize widget pointers to 0
{
    // Constructor.
    
    MakeTitle("AliEveCascade");
    
    fInfoLabelRadius = new TGLabel(this);
    fInfoLabelRadius->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelRadius, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelDCA = new TGLabel(this);
    fInfoLabelDCA->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelDCA,    new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelCharge = new TGLabel(this);
    fInfoLabelCharge->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelCharge, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelPhi = new TGLabel(this);
    fInfoLabelPhi->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelPhi,    new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelTheta = new TGLabel(this);
    fInfoLabelTheta->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelTheta,  new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelPtot = new TGLabel(this);
    fInfoLabelPtot->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelPtot,   new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelPt = new TGLabel(this);
    fInfoLabelPt->SetTextJustify(kTextLeft);
    AddFrame(fInfoLabelPt,     new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    fInfoLabelEta = new TGLabel(this);
    fInfoLabelEta->SetTextJustify(kTextLeft);
    AddFrame( fInfoLabelEta,   new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 8, 0, 2, 0));
    
    
    fXButtonDetailedView = new TGTextButton(this, "Detailed View");
    AddFrame(fXButtonDetailedView, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    fXButtonDetailedView->Connect("Clicked()", "AliEveCascadeEditor", this, "DisplayDetailed()");
    
    fXButtonMassHyp = new TGTextButton(this, "Mass Hypotheses according to charge");
    AddFrame(fXButtonMassHyp, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    fXButtonMassHyp->Connect("Clicked()", "AliEveCascadeEditor", this, "DisplayMassHyp()");
    
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeEditor::SetModel(TObject* obj)
{
    // Set model object.
    
    fM = static_cast<AliEveCascade*>(obj);
    
    // Set values of widgets
    fInfoLabelRadius->SetText(Form("Radius = %f cm",     fM->GetRadius() ));
    fInfoLabelDCA   ->SetText(Form("DCA (Xi dghters) = %f cm", fM->GetDaughterDCA()));
    fInfoLabelCharge->SetText(Form("Charge = %d",        fM->GetCharge() ));
    fInfoLabelPhi   ->SetText(Form("Phi     = %f deg",   fM->GetPhi()   * 180.0/TMath::Pi()    ));
    fInfoLabelTheta ->SetText(Form("Theta  = %f deg",    fM->GetTheta() * 180.0/TMath::Pi()    ));
    fInfoLabelPtot  ->SetText(Form("Ptot    = %f GeV/c", fM->GetPtot()   ));
    fInfoLabelPt    ->SetText(Form("Pt      = %f GeV/c", fM->GetPt()     ));
    fInfoLabelEta   ->SetText(Form("Eta    = %f",        fM->GetEta()    ));
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveCascadeEditor::DoXYZZ()
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
#include <AliEveTrack.h>

#include <TGLCamera.h>
#include <TGLViewer.h>
#include "TGLCameraOverlay.h"

#include <TLatex.h>
#include <TRootEmbeddedCanvas.h>
#include <TInterpreter.h>


void AliEveCascadeEditor::DisplayDetailed()
{
    // Display a detailed view (bending plane + transverse plane)
    
    printf("\n--> Detailed View :\n");
    TEveWindowSlot *slot = TEveWindow::CreateWindowMainFrame();
    TEveWindowPack *pack = slot->MakePack();
    pack->SetShowTitleBar(kFALSE);
    pack->SetHorizontal();
    
    //----------
    // Part 0 : get the different objects to display
    //
    
    TEvePointSet *lPrimVtxLocation = new TEvePointSet("Prim Vtx location");
    lPrimVtxLocation->SetNextPoint(fM->fRecBirthV.fX, fM->fRecBirthV.fY, fM->fRecBirthV.fZ);
    lPrimVtxLocation->SetTitle("prim vtx location");
    lPrimVtxLocation->SetMarkerStyle(4);
    lPrimVtxLocation->SetMarkerSize(2.5);
    lPrimVtxLocation->SetMarkerColor(kCyan);
    
    TEvePointSet *lCascadeLocation = new TEvePointSet("Cascade decay location");
    lCascadeLocation->SetNextPoint(fM->fRecDecayV.fX, fM->fRecDecayV.fY, fM->fRecDecayV.fZ);
    lCascadeLocation->SetTitle("Cascade decay location");
    lCascadeLocation->SetMarkerStyle(4);
    lCascadeLocation->SetMarkerSize(2.5);
    lCascadeLocation->SetMarkerColor(kMagenta-9);
    
    AliEveTrack *bacTrack = fM->GetBacTrack();
    AliEveTrack *negTrack = fM->GetNegTrack();
    AliEveTrack *posTrack = fM->GetPosTrack();
    
    char macroWithIndex[100] = {0};
    Int_t daughterIndex = 0;
    
    TEvePointSet *bacDaughterCluster = 0;
    TEvePointSet *negDaughterCluster = 0;
    TEvePointSet *posDaughterCluster = 0;
    
    // Clusters linked with the bachelor track
    daughterIndex = bacTrack->GetIndex();
    
    AliEveTrack *eveTrack = new AliEveTrack();
    
    bacDaughterCluster = eveTrack->ImportClustersFromIndex(daughterIndex);
    if (bacDaughterCluster) {
        bacDaughterCluster->SetMarkerStyle(4);
        bacDaughterCluster->SetMarkerSize(1.5);
        bacDaughterCluster->SetMarkerColor(kMagenta);
    }
    else
    {
        Warning("DisplayDetailed", "Cascade : Import of bachelor clusters failed.");
    }
    
    // Clusters linked with the negative daughter track (V0)
    daughterIndex = negTrack->GetIndex();
    
    negDaughterCluster = eveTrack->ImportClustersFromIndex(daughterIndex);
    if (negDaughterCluster) {
        negDaughterCluster->SetMarkerStyle(4);
        negDaughterCluster->SetMarkerSize(1.5);
        negDaughterCluster->SetMarkerColor(kBlue+3);
    }
    else
    {
        Warning("DisplayDetailed", "Cascade : Import of negative daughter's clusters failed.");
    }
    
    // Clusters linked with the positive daughter track (V0)
    daughterIndex = posTrack->GetIndex();
    
    posDaughterCluster = eveTrack->ImportClustersFromIndex(daughterIndex);
    if (posDaughterCluster) {
        posDaughterCluster->SetMarkerStyle(4);
        posDaughterCluster->SetMarkerSize(1.5);
        posDaughterCluster->SetMarkerColor(kRed+3);
    }
    else
    {
        Warning("DisplayDetailed", "Cascade : Import of positive daughter's clusters failed.");
    }
    
    delete eveTrack;eveTrack=0;
    
    
    //----------
    // Part 1 : bending plane view
    //
    pack->NewSlot()->MakeCurrent();
    TEveViewer *bpViewer = gEve->SpawnNewViewer("Cascade bending plane");
    TEveScene  *bpScene  = gEve->SpawnNewScene("Cascade bending plane Scene");
    
    TEveProjectionManager *projMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    bpScene->AddElement(projMgr);
    
    TEveUtil::LoadMacro("geom_gentle.C");
    Long_t result = gInterpreter->ProcessLine("geom_gentle_rphi()");
    if (result)
    {
        TEveGeoShape *geomRPhi = reinterpret_cast<TEveGeoShape*>(result);
        geomRPhi->IncDenyDestroy();
        projMgr->SetCurrentDepth(-10);
        projMgr->ImportElements(geomRPhi);
        projMgr->SetCurrentDepth(0);
    }
    else
    {
        Warning("DisplayDetailed", "Cascade (bending plane view) : Import of R-Phi geometry failed.");
    }
    
    projMgr->ImportElements(fM);
    projMgr->ImportElements(lPrimVtxLocation);
    projMgr->ImportElements(lCascadeLocation);
    bpViewer->AddScene(bpScene);
    
    if (bacDaughterCluster) projMgr->ImportElements(bacDaughterCluster);
    if (negDaughterCluster) projMgr->ImportElements(negDaughterCluster);
    if (posDaughterCluster) projMgr->ImportElements(posDaughterCluster);
    
    // This is the to-do list for the bending plane:
    // 1. show the V0 daughters track + corresponding clusters                        -> to do ...
    // 2. show axis and tickles along X and Y                                         -> done!
    //       -> note for the projection the cartesian scales are not very useful
    //       -> propose a phi and R scale which rotates with a reset at 0;
    //       -> propose a transformation for an eta scale (keep the z one);
    // 3. show the center, the main vertex and the detectors for this view            -> done!
    // 4. show V0 direction in the bending plane with arrow length proportional to pT -> done!
    // 5. show angles with respect to axis (phi angle)                                -> almost.
    // 6. show clusters in the ITS and in the TPC associated with the daughter tracks
    //       -> include a radius cut for plotting only ITS and TPC                    -> done!
    bpViewer->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    bpViewer->GetGLViewer()->ResetCamerasAfterNextUpdate();
    TGLViewer *lbpGLViewer = bpViewer->GetGLViewer();
    TGLCameraOverlay* co = lbpGLViewer->GetCameraOverlay();
    co->SetShowOrthographic(true); //(false);
    co->SetOrthographicMode(TGLCameraOverlay::kAxis); // ::kPlaneIntersect or ::kBar
    // lbpGLViewer->RequestDraw();
    // end of the bending plane part
    
    
    
    
    //----------
    // Part 2 : decay plane view
    //
    pack->NewSlot()->MakeCurrent();
    TEveViewer *dpViewer = gEve->SpawnNewViewer("Cascade decay plane");
    TEveScene  *dpScene  = gEve->SpawnNewScene("Cascade decay plane Scene");
    
    
    
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
        Warning("DisplayDetailed", "Cascade (decay plane view) : Import of 3D geometry  failed.");
    }
    
    dpViewer->AddScene(dpScene);
    dpScene->AddElement(fM);
    dpScene->AddElement(lPrimVtxLocation);
    dpScene->AddElement(lCascadeLocation);
    if (bacDaughterCluster) dpScene->AddElement(bacDaughterCluster);
    if (negDaughterCluster) dpScene->AddElement(negDaughterCluster);
    if (posDaughterCluster) dpScene->AddElement(posDaughterCluster);
    
    
    
    // This is the to-do list for the decay plane:
    // 1. fix the view to decay plane (no rotation allowed but moving the center ok)
    // 3. show the center, the main vertex and the detectors for this view -> done!
    // 4. show the x,y and z axis and the different angles
    //       -> this needs a referential object that we can move around
    //          or fix to a selected point (origin being the default)
    // 5. draw the dca between Xi daughters and the extrapolation to the main vertex.
    //       -> this is an issue since we only store the distance: check with J.Belikov
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
    
    
    
    
    //----------
    // Part 3 : displaying extra information
    //
    slot = pack->NewSlot();
    
    TEveWindowFrame *frame = slot->MakeFrame(new TRootEmbeddedCanvas());
    frame->SetElementName("Details");
    
    // Print and show detailed information about the V0
    // Calculation of the invariant mass with the max prob PID hypothesis first
    // pseudorapidity, phi angle, pt, radius, dcas
    char info[100] = {0};
    snprintf(info,100,"#phi = %.3f rad = %.1f deg",fM->GetPhi(),(180./TMath::Pi())*fM->GetPhi());
    TLatex* ltx = new TLatex(0.05, 0.9, info);
    ltx->SetTextSize(0.08);
    ltx->Draw();
    
    snprintf(info,100,"radius = %.3f [cm]",fM->GetRadius());
    ltx->DrawLatex(0.05, 0.8, info);
    
    snprintf(info,100,"p_{T} = %.3f [GeV/c]",fM->GetPt());
    ltx->DrawLatex(0.05, 0.7, info);
    
    snprintf(info,100,"Xi dghtrs dca = %.4f [cm]",fM->GetDaughterDCA());
    ltx->DrawLatex(0.05, 0.6, info);
    
    snprintf(info,100,"#eta = - ln( tan(#theta/2) ) = %.3f",fM->GetEta());
    ltx->DrawLatex(0.05, 0.5, info);
    
    
    if(fM->GetCharge() < 0){
        snprintf(info,100,"mass_{#Xi^{-}} : %.3f [GeV/c^{2}]", fM->GetXiMinusInvMass()    );
        ltx->DrawLatex(0.05, 0.3, info);
        snprintf(info,100,"mass_{#Omega^{-}} : %.3f [GeV/c^{2}]", fM->GetOmegaMinusInvMass() );
        ltx->DrawLatex(0.05, 0.2, info);
    }
    else {
        snprintf(info,100,"mass_{#Xi^{+}} : %.3f [GeV/c^{2}]", fM->GetXiPlusInvMass()     );
        ltx->DrawLatex(0.05, 0.3, info);
        snprintf(info,100,"mass_{#Omega^{+}} : %.3f [GeV/c^{2}]", fM->GetOmegaPlusInvMass()  );
        ltx->DrawLatex(0.05, 0.2, info);
    }
    
    gEve->Redraw3D();
}

void AliEveCascadeEditor::DisplayMassHyp()
{	
    // Printf the proper mass hypotheses, according to the charge of the Bachelor
    printf("\n--> Mass Hyp:");
    
    if(fM->GetCharge() < 0){
        printf("Xi-    mass hyp : %f GeV/c^2}\n", fM->GetXiMinusInvMass()    );
        printf("Omega- mass hyp : %f GeV/c^2\n\n", fM->GetOmegaMinusInvMass() );
    }
    else{
        printf("Xi+    mass hyp : %f GeV/c^2\n", fM->GetXiPlusInvMass()     );
        printf("Omega+ mass hyp : %f GeV/c^2\n\n", fM->GetOmegaPlusInvMass()  );
    }
}
