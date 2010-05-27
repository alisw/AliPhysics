// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliESDEvent.h"
#include "AliEveLego.h"
#include "AliEveEventManager.h"
#include "AliEveMultiView.h"
#include "AliPhysicsSelection.h"
#include "AliEveEventSelector.h"

#include "TH2F.h"
#include "TMath.h"
#include "TGLViewer.h"
#include "TEveWindow.h"
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TEveCaloLegoOverlay.h"
#include "TEveCalo.h"
#include "TEveCaloData.h"
#include "TEveLegoEventHandler.h"
#include "TEveTrans.h"
#include "TEveProjectionManager.h"
#include "TEveProjectionAxes.h"
#include "TGLWidget.h"
#include "TGLOverlayButton.h"

//______________________________________________________________________________
// Allow histograms visualization in 2D and 3D.
//

ClassImp(AliEveLego)
Double_t pi = TMath::Pi();

//______________________________________________________________________________
AliEveLego::AliEveLego(const char* name) :
  TEveElementList(name),
  fChargeId(1),
  fTracksId(1),
  fEventsId(1),
  fMaxPt(10000),
  fChargeIdAE(0),
  fTracksIdAE(0),
  fMaxPtAE(0),
  fEsd(0),
  fPhysicsSelection(0),
  fHistopos(0),
  fHistoposclone(0),
  fHistopos_all_events(0),
  fHistoneg(0),
  fHistonegclone(0),
  fHistoneg_all_events(0),
  fData(0),
  fData_all_events(0),
  fLego(0),
  fLego_all_events(0),
  fCalo3d(0),
  fCalo3d_all_events(0),
  fGlv(0),
  fHisto2d_v(0),
  fHisto2d_s(0),
  fHisto2d_s2(0),
  fHisto2d_all_events_v0(0),
  fHisto2d_all_events_v1(0),
  fHisto2d_all_events_v2(0),
  fHisto2d_all_events_v3(0),
  fHisto2d_all_events_s0(0),
  fHisto2d_all_events_s1(0),
  fHisto2d_all_events_s2(0),
  fHisto2d_all_events_s3(0),
  fAl(0),
  fHisto2d_lego_overlay(0),
  fHisto2d_all_events_lego_overlay(0),
  fHisto2d_all_events_slot(0),
  fEventSelector(0),
  fShowEventsInfo(0),
  fGButton(0),
  fB1(0),
  fB2(0)
{
  // Constructor.
  gEve->AddToListTree(this,0);

  // Get Current ESD event
  fEsd = AliEveEventManager::AssertESD();

  fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection->Initialize(fEsd->GetRunNumber());

  fEventSelector = AliEveEventManager::GetMaster()->GetEventSelector();

  fHistopos = new TH2F("histopos","Histo 2d positive", 100, -1.5, 1.5, 80, -pi, pi);
  fHistoneg = new TH2F("histoneg","Histo 2d negative", 100, -1.5, 1.5, 80, -pi, pi);
//  fHistoposclone = new TH2F("histoposclone","Histo 2d positive", 100, -1.5, 1.5, 80, -pi, pi);
//  fHistonegclone = new TH2F("histonegclone","Histo 2d positive", 100, -1.5, 1.5, 80, -pi, pi);

  fData = new TEveCaloDataHist();
  fData->AddHistogram(fHistoneg);
  fData->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
  fData->AddHistogram(fHistopos);
  fData->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
  fData->GetEtaBins()->SetTitleFont(120);
  fData->GetEtaBins()->SetTitle("h");
  fData->GetPhiBins()->SetTitleFont(120);
  fData->GetPhiBins()->SetTitle("f");
  fData->IncDenyDestroy();

  fCalo3d = new TEveCalo3D(fData);
  fCalo3d->SetBarrelRadius(550);
  fCalo3d->SetEndCapPos(550);

  // plotting histo
  fLego = new TEveCaloLego(fData);

  // projections
  fAl = AliEveMultiView::Instance();
  fAl->ImportEventRPhi(fCalo3d);
  fAl->ImportEventRhoZ(fCalo3d);

  // glbutton
  fGButton = new TGLOverlayButton(0, "", 10.0, -10.0, 190.0, 20.0);
  fGButton->SetAlphaValues(1.5,1.5);
  fB1 = new TGLOverlayButton(0, "", 10.0, -30.0, 190.0, 20.0);
  fB1->SetAlphaValues(1.5,1.5);
  fB2 = new TGLOverlayButton(0, "", 10.0, -50.0, 190.0, 20.0);
  fB2->SetAlphaValues(1.5,1.5);

  // Update
  Update();
}

//______________________________________________________________________________
AliEveLego::~AliEveLego()
{
   delete fEsd;
   delete fPhysicsSelection;
   delete fHistopos;
   delete fHistopos_all_events;
   delete fHistoneg;
   delete fHistoneg_all_events;

   delete fData;
   delete fData_all_events;
   delete fLego;
   delete fLego_all_events;
   delete fCalo3d;
   delete fCalo3d_all_events;
   delete fGlv;

   delete fHisto2d_v;
   delete fHisto2d_s;
   delete fHisto2d_s2;
   delete fHisto2d_all_events_v0;
   delete fHisto2d_all_events_v1;
   delete fHisto2d_all_events_v2;
   delete fHisto2d_all_events_v3;
   delete fHisto2d_all_events_s0;
   delete fHisto2d_all_events_s1;
   delete fHisto2d_all_events_s2;
   delete fHisto2d_all_events_s3;

   delete fAl;
   delete fHisto2d_lego_overlay;
   delete fHisto2d_all_events_lego_overlay;
   delete fHisto2d_all_events_slot;

   delete fEventSelector;
   delete fGButton;
   delete fB1;
   delete fB2;
}

//______________________________________________________________________________
Double_t getphi(Double_t phi)
{
   Double_t pi = TMath::Pi();

   if (phi > pi) {
      phi -= 2*pi;
   }
   return phi;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadData()
{

   fHistopos->Reset();
   fHistoneg->Reset();

   // Getting current tracks, filling histograms
   for (int n = 0; n < fEsd->GetNumberOfTracks(); ++n) {

      if (fEsd->GetTrack(n)->GetSign() > 0) {
         fHistopos->Fill(fEsd->GetTrack(n)->Eta(),
                         getphi(fEsd->GetTrack(n)->Phi()),
                         fabs(fEsd->GetTrack(n)->Pt()));
      }

      if (fEsd->GetTrack(n)->GetSign() < 0) {
         fHistoneg->Fill(fEsd->GetTrack(n)->Eta(),
                         getphi(fEsd->GetTrack(n)->Phi()),
                         fabs(fEsd->GetTrack(n)->Pt()));
      }
   }

//   fHistoposclone = (TH2F*)fHistopos->Clone("histoposclone");
//   fHistonegclone = (TH2F*)fHistoneg->Clone("histonegclone");
//   fHistoposclone->SetName("histoposclone");
//   fHistonegclone->SetName("histonegclone");

   fData->DataChanged();

   FilterData();

   return fData;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadAllData()
{

   fHistopos_all_events->Reset();
   fHistoneg_all_events->Reset();

   TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

   // Getting current tracks for each event, filling histograms
   for (int event = 0; event < t->GetEntries(); event++) {
      t->GetEntry(event);
         for (int n = 0; n < fEsd->GetNumberOfTracks(); ++n) {

            if (fEsd->GetTrack(n)->GetSign() > 0) {
               fHistopos_all_events->Fill(fEsd->GetTrack(n)->Eta(),
                                          getphi(fEsd->GetTrack(n)->Phi()),
                                          fabs(fEsd->GetTrack(n)->Pt()));
            } else {
               fHistoneg_all_events->Fill(fEsd->GetTrack(n)->Eta(),
                                          getphi(fEsd->GetTrack(n)->Phi()),
                                          fabs(fEsd->GetTrack(n)->Pt()));
            }
         }
   }

   fData_all_events->DataChanged();

   return fData_all_events;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::FilterData()
{
   // Tracks selection
   if ( fTracksId == 2 )
   {
      fHistopos->Reset();
      fHistoneg->Reset();

      const AliESDVertex *pv  = fEsd->GetPrimaryVertex();

      for (Int_t n = 0; n < pv->GetNIndices(); n++ )
      {
         AliESDtrack *at = fEsd->GetTrack(pv->GetIndices()[n]);
         if (at->GetSign() > 0) {
            fHistopos->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
         }
         if (at->GetSign() < 0) {
            fHistoneg->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
         }
      }
   }

   fData->DataChanged();

   // Max Pt threshold
   if (GetPtMax() >= fMaxPt){
      for (Int_t binx = 1; binx <= 100; binx++) {
         for (Int_t biny = 1; biny <= 80; biny++) {
            if (fHistopos->GetBinContent(binx, biny) >= fMaxPt)
            {
               fHistopos->SetBinContent(binx, biny, fMaxPt);
            }
            if (fHistoneg->GetBinContent(binx, biny) >= fMaxPt)
            {
               fHistoneg->SetBinContent(binx, biny, fMaxPt);
            }           
         }
      }
   }

   // Positive only
   if ( fChargeId == 2 ) fHistoneg->Reset();

   // Negative only
   if ( fChargeId == 3 ) fHistopos->Reset();

   fData->DataChanged();

   return fData;
}


//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::FilterAllData()
{
   // Tracks selection
   if ( fTracksIdAE == 2 )
   {
      fHistopos_all_events->Reset();
      fHistoneg_all_events->Reset();

      TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

      // Getting current tracks for each event, filling histograms
      for (int event = 0; event < t->GetEntries(); event++) {
         t->GetEntry(event);

      const AliESDVertex *pv  = fEsd->GetPrimaryVertex();

      for (Int_t n = 0; n < pv->GetNIndices(); n++ )
      {
         AliESDtrack *at = fEsd->GetTrack(pv->GetIndices()[n]);
         if (at->GetSign() > 0) {
            fHistopos_all_events->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
         }
         if (at->GetSign() < 0) {
            fHistoneg_all_events->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
         }
      }
   }
   } else {
      LoadAllData();
   }
   
   fData_all_events->DataChanged();
   
   // Max Pt threshold
   if (GetPtMaxAE() >= fMaxPtAE){
      for (Int_t binx = 1; binx <= 100; binx++) {
         for (Int_t biny = 1; biny <= 80; biny++) {
            if (fHistopos_all_events->GetBinContent(binx, biny) >= fMaxPtAE)
            {
               fHistopos_all_events->SetBinContent(binx, biny, fMaxPtAE);
            }
            if (fHistoneg_all_events->GetBinContent(binx, biny) >= fMaxPtAE)
            {
               fHistoneg_all_events->SetBinContent(binx, biny, fMaxPtAE);
            }           
         }
      }
   }

   // Positive only
   if ( fChargeIdAE == 2 ) fHistoneg_all_events->Reset();

   // Negative only
   if ( fChargeIdAE == 3 ) fHistopos_all_events->Reset();

   fData_all_events->DataChanged();

   gEve->Redraw3D(kTRUE);

   return fData_all_events;
}

//______________________________________________________________________________
void AliEveLego::Update()
{
   // Load/Reload data
   LoadData();

   // Create new histo2d
   CreateHistoLego();

   // Create 3d view
   Create3DView();

   // Show information about event;
   ShowEventSeletion(fShowEventsInfo, kTRUE);

   // Update the viewers
   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
TEveCaloLego* AliEveLego::CreateHistoLego()
{
   // Viewer initialization, tab creation
   if (fHisto2d_v == 0) {
      TEveWindowSlot *fslot    = 0;
      TEveBrowser    *fbrowser = gEve->GetBrowser();

      fslot = TEveWindow::CreateWindowInTab(fbrowser->GetTabRight());
      fslot->MakeCurrent();
      fHisto2d_v = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      fHisto2d_s = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      fHisto2d_v->AddScene(fHisto2d_s);
      fHisto2d_v->SetElementName("2D Lego Viewer");
      fHisto2d_s->SetElementName("2D Lego Scene");

      fGlv = fHisto2d_v->GetGLViewer();
      fHisto2d_lego_overlay = new TEveCaloLegoOverlay();
      fGlv->AddOverlayElement(fHisto2d_lego_overlay);
      fGlv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);

      fHisto2d_s->AddElement(fLego);

      // move to real world coordinates
      fLego->InitMainTrans();
      Float_t sc = TMath::Min(fLego->GetEtaRng(), fLego->GetPhiRng());
      fLego->RefMainTrans().SetScale(sc, sc, sc);

      // set event handler to move from perspective to orthographic view.
      fGlv->SetEventHandler(new TEveLegoEventHandler(fGlv->GetGLWidget(), fGlv, fLego));

      fHisto2d_lego_overlay->SetCaloLego(fLego);
   }

   return fLego;
}

//______________________________________________________________________________
TEveCaloLego* AliEveLego::CreateHistoLego(TEveWindowSlot *slot)
{
   // Viewer initialization, tab creation
   if (fHisto2d_all_events_v0 == 0) {

      slot->MakeCurrent();
      fHisto2d_all_events_v0 = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      fHisto2d_all_events_s0 = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      fHisto2d_all_events_v0->AddScene(fHisto2d_all_events_s0);
      fHisto2d_all_events_v0->SetElementName("2D Lego Viewer");
      fHisto2d_all_events_s0->SetElementName("2D Lego Scene");

      TGLViewer* glv = fHisto2d_all_events_v0->GetGLViewer();
      fHisto2d_all_events_lego_overlay = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(fHisto2d_all_events_lego_overlay);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);

      // Plotting histogram lego
      fLego_all_events = new TEveCaloLego(fData_all_events);
      fHisto2d_all_events_s0->AddElement(fLego_all_events);

      // Move to real world coordinates
      fLego_all_events->InitMainTrans();
      Float_t sc = TMath::Min(fLego_all_events->GetEtaRng(), fLego_all_events->GetPhiRng());
      fLego_all_events->RefMainTrans().SetScale(sc, sc, sc);

      // Set event handler to move from perspective to orthographic view.
      glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, fLego_all_events));

      fHisto2d_all_events_lego_overlay->SetCaloLego(fLego_all_events);
   }

   return fLego_all_events;
}

//______________________________________________________________________________
TEveCalo3D* AliEveLego::Create3DView()
{
   //initialization
   if (fHisto2d_s2 == 0) {
      fHisto2d_s2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      gEve->GetDefaultViewer()->AddScene(fHisto2d_s2);
      fHisto2d_s2->SetElementName("3D Histogram Scene");
      fHisto2d_s2->AddElement(fCalo3d);
   }

   return fCalo3d;
}

//______________________________________________________________________________
TEveCalo3D* AliEveLego::Create3DView(TEveWindowSlot *slot)
{
   if ( fHisto2d_all_events_v1 == 0 ) {

      slot->MakeCurrent();
      fHisto2d_all_events_v1 = gEve->SpawnNewViewer("3D Histogram", "3D Histogram");
      fHisto2d_all_events_s1 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      fHisto2d_all_events_v1->AddScene(fHisto2d_all_events_s1);
      fHisto2d_all_events_v1->SetElementName("3D Histogram Viewer");
      fHisto2d_all_events_s1->SetElementName("3D Histogram Scene");

      fCalo3d_all_events = new TEveCalo3D(fData_all_events);

      fCalo3d_all_events->SetBarrelRadius(550);
      fCalo3d_all_events->SetEndCapPos(550);
      fHisto2d_all_events_s1->AddElement(fCalo3d_all_events);
   }

   return fCalo3d_all_events;
}

//______________________________________________________________________________
void AliEveLego::CreateProjections(TEveWindowSlot* slot1, TEveWindowSlot* slot2){

   if (fHisto2d_all_events_v2 == 0) {

      slot1->MakeCurrent();
      fHisto2d_all_events_v2 = gEve->SpawnNewViewer("RPhi projection", "RPhi projection");
      fHisto2d_all_events_s2 = gEve->SpawnNewScene("RPhi projection", "RPhi projection");
      fHisto2d_all_events_v2->AddScene(fHisto2d_all_events_s2);
      fHisto2d_all_events_v2->SetElementName("RPhi Projection Viewer");
      fHisto2d_all_events_s2->SetElementName("RPhi Projection Scene");

      TEveProjectionManager* mng1 = new TEveProjectionManager();
      mng1->SetProjection(TEveProjection::kPT_RPhi);

      TEveProjectionAxes* axeg_histo2d_all_events_s1 = new TEveProjectionAxes(mng1);
      fHisto2d_all_events_s2->AddElement(axeg_histo2d_all_events_s1);
      TEveCalo2D* fcalo2d1 = (TEveCalo2D*) mng1->ImportElements(fCalo3d_all_events);
      fHisto2d_all_events_s2->AddElement(fcalo2d1);

      fHisto2d_all_events_v2->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   }

   if (fHisto2d_all_events_v3 == 0) {

      slot2->MakeCurrent();
      fHisto2d_all_events_v3 = gEve->SpawnNewViewer("RhoZ projection", "RhoZ projection");
      fHisto2d_all_events_s3 = gEve->SpawnNewScene("RhoZ projection", "RhoZ projection");
      fHisto2d_all_events_v3->AddScene(fHisto2d_all_events_s3);
      fHisto2d_all_events_v3->SetElementName("RhoZ Projection Viewer");
      fHisto2d_all_events_s3->SetElementName("RhoZ Projection Viewer");

      TEveProjectionManager* mng2 = new TEveProjectionManager();
      mng2->SetProjection(TEveProjection::kPT_RhoZ);

      TEveProjectionAxes* axeg_histo2d_all_events_s2 = new TEveProjectionAxes(mng2);
      fHisto2d_all_events_s3->AddElement(axeg_histo2d_all_events_s2);
      TEveCalo2D* fcalo2d2 = (TEveCalo2D*) mng2->ImportElements(fCalo3d_all_events);
      fHisto2d_all_events_s3->AddElement(fcalo2d2);

      fHisto2d_all_events_v3->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   }

   return;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadAllEvents()
{
   if ( fHisto2d_all_events_slot == 0 ) {

      printf("Filling histogram...\n");

      // Creating 2D histograms
      fHistopos_all_events = new TH2F("fHistopos_all_events","Histo 2d positive",
                                 100,-1.5,1.5,80,-pi,pi);
      fHistoneg_all_events = new TH2F("fHistoneg_all_events","Histo 2d negative",
                                 100,-1.5,1.5,80,-pi,pi);

      fData_all_events = new TEveCaloDataHist();
      fData_all_events->AddHistogram(fHistoneg_all_events);
      fData_all_events->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
      fData_all_events->AddHistogram(fHistopos_all_events);
      fData_all_events->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
      fData_all_events->GetEtaBins()->SetTitleFont(120);
      fData_all_events->GetEtaBins()->SetTitle("h");
      fData_all_events->GetPhiBins()->SetTitleFont(120);
      fData_all_events->GetPhiBins()->SetTitle("f");
      fData_all_events->IncDenyDestroy();

      // Creating frames
      fHisto2d_all_events_slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
      TEveWindowPack* packH = fHisto2d_all_events_slot->MakePack();
      packH->SetElementName("Projections");
      packH->SetHorizontal();
      packH->SetShowTitleBar(kFALSE);

      fHisto2d_all_events_slot = packH->NewSlot();
      TEveWindowPack* pack0 = fHisto2d_all_events_slot->MakePack();
      pack0->SetShowTitleBar(kFALSE);
      TEveWindowSlot*  slotLeftTop   = pack0->NewSlot();
      TEveWindowSlot* slotLeftBottom = pack0->NewSlot();

      fHisto2d_all_events_slot = packH->NewSlot();
      TEveWindowPack* pack1 = fHisto2d_all_events_slot->MakePack();
      pack1->SetShowTitleBar(kFALSE);
      TEveWindowSlot* slotRightTop    = pack1->NewSlot();
      TEveWindowSlot* slotRightBottom = pack1->NewSlot();

      // Creating viewers and scenes
      Create3DView(slotLeftTop);
      CreateHistoLego(slotLeftBottom);
      CreateProjections(slotRightTop, slotRightBottom);

      LoadAllData();

      gEve->Redraw3D(kTRUE);

      printf("Filling histogram... Finished\n");
   }

   return fData_all_events;
}

//______________________________________________________________________________
Float_t AliEveLego::GetPtMax()
{
   return fData->GetMaxVal(fLego->GetPlotEt());
}

//______________________________________________________________________________
Float_t AliEveLego::GetPtMaxAE()
{
   return fData_all_events->GetMaxVal(fLego_all_events->GetPlotEt());
}

//______________________________________________________________________________
void AliEveLego::SetMaxPt(Double_t val)
{
   // Add new maximum
   fMaxPt = val;   
   Update();
}

//______________________________________________________________________________
void AliEveLego::SetMaxPtAE(Double_t val)
{
   // Add new maximum
   fMaxPtAE = val;
   FilterAllData();
}

//______________________________________________________________________________
void AliEveLego::SetThreshold(Double_t val)
{  
   // Setting up the new threshold for all histograms
   fData->SetSliceThreshold(0,val);
   fData->SetSliceThreshold(1,val);
   fData->DataChanged();

   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
void AliEveLego::SetThresholdAE(Double_t val)
{
   // Setting up the new threshold for all histograms
   fData_all_events->SetSliceThreshold(0,val);
   fData_all_events->SetSliceThreshold(1,val);
   fData_all_events->DataChanged();

   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
void AliEveLego::SetEventSelection()
{
   if (fShowEventsInfo == 0)
   {
      fShowEventsInfo = 1;
   } else {
      fShowEventsInfo = 0;
   }

   ShowEventSeletion(fShowEventsInfo);
}

//______________________________________________________________________________
void AliEveLego::ShowEventSeletion(Bool_t show, Bool_t updateonly)
{
   if (show == 0)
   {
      gEve->GetDefaultGLViewer()->RemoveOverlayElement(fGButton);
      fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fGButton);
      fHisto2d_v->GetGLViewer()->RemoveOverlayElement(fGButton);

      gEve->GetDefaultGLViewer()->RemoveOverlayElement(fB1);
      fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fB1);
      fHisto2d_v->GetGLViewer()->RemoveOverlayElement(fB1);

      gEve->GetDefaultGLViewer()->RemoveOverlayElement(fB2);
      fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fB2);
      fHisto2d_v->GetGLViewer()->RemoveOverlayElement(fB2);

   } else {

      //Collision candidate
      if (updateonly == kFALSE) {
         gEve->GetDefaultGLViewer()->AddOverlayElement(fGButton);
         fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fGButton);
         fHisto2d_v->GetGLViewer()->AddOverlayElement(fGButton);
      }

      Bool_t ev = fPhysicsSelection->IsCollisionCandidate(fEsd);

      if (ev == 1)
      {
         fGButton->SetText("Collision candidate: YES");
      } else {
         fGButton->SetText("Collision candidate: NO ");
      }

      // Beam 1 & 2 setup: method 1
      if (updateonly == kFALSE) {
         gEve->GetDefaultGLViewer()->AddOverlayElement(fB1);
         fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fB1);
         fHisto2d_v->GetGLViewer()->AddOverlayElement(fB1);

         gEve->GetDefaultGLViewer()->AddOverlayElement(fB2);
         fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fB2);
         fHisto2d_v->GetGLViewer()->AddOverlayElement(fB2);
      }

      Bool_t b1  = fEsd->IsTriggerClassFired("CINT1A-ABCE-NOPF-ALL");
      Bool_t b2  = fEsd->IsTriggerClassFired("CINT1C-ABCE-NOPF-ALL");
      Bool_t b12 = fEsd->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL");

      if (b1 == 1 || b12 == 1)
      {
         fB1->SetText("Beam 1: YES");
         fB1->SetBackColor(0x00ff00);
      } else {
         fB1->SetText("Beam 1: NO");
         fB1->SetBackColor(0xff0000);
      }

      if (b2 == 1 || b12 == 1)
      {
         fB2->SetText("Beam 2: YES");
         fB2->SetBackColor(0x00ff00);
      } else {
         fB2->SetText("Beam 2: NO");
         fB2->SetBackColor(0xff0000);
      }
   }

   gEve->Redraw3D(kTRUE);

}

//______________________________________________________________________________
void AliEveLego::SelectEventSelection(Int_t id)
{
   if (id == 0)
   {
      fEventSelector->SetSelectOnTriggerType(kFALSE);
   } else {
      if (id == 1) fEventSelector->SetTriggerType("CINT1A-ABCE-NOPF-ALL");
      if (id == 2) fEventSelector->SetTriggerType("CINT1C-ABCE-NOPF-ALL");
      if (id == 3) fEventSelector->SetTriggerType("CINT1B-ABCE-NOPF-ALL");
      fEventSelector->SetSelectOnTriggerType(kTRUE);
   }
}

//______________________________________________________________________________
void AliEveLego::ShowPrevEvent()
{
   AliEveEventManager::GetMaster()->PrevEvent();
}

//______________________________________________________________________________
void AliEveLego::ShowNextEvent()
{
   AliEveEventManager::GetMaster()->NextEvent();
}
/******************************************************************************/


