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
#include "TStopwatch.h"

//______________________________________________________________________________
// Allow histograms visualization in 2D and 3D.
//

ClassImp(AliEveLego)
Double_t kPi = TMath::Pi();

//______________________________________________________________________________
AliEveLego::AliEveLego(const char* name) :
  TEveElementList(name),
  fChargeId(1),
  fTracksId(1),
  fEventsId(1),
  fMaxPt(10000),
  fChargeIdAE(0),
  fTracksIdAE(0),
  fMaxPtAE(10000),
  fEsd(0),
  fPhysicsSelection(0),
  fHistopos(0),
  fHistoposAllEvents(0),
  fHistoneg(0),
  fHistonegAllEvents(0),
  fData(0),
  fDataAllEvents(0),
  fLego(0),
  fLegoAllEvents(0),
  fCalo3d(0),
  fCalo3dAllEvents(0),
  fGlv(0),
  fHisto2dv(0),
  fHisto2ds(0),
  fHisto2ds2(0),
  fHisto2dAllEventsv0(0),
  fHisto2dAllEventsv1(0),
  fHisto2dAllEventsv2(0),
  fHisto2dAllEventsv3(0),
  fHisto2dAllEventss0(0),
  fHisto2dAllEventss1(0),
  fHisto2dAllEventss2(0),
  fHisto2dAllEventss3(0),
  fAl(0),
  fHisto2dLegoOverlay(0),
  fHisto2dAllEventsLegoOverlay(0),
  fHisto2dAllEventsSlot(0),
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

  fHistopos = new TH2F("histopos","Histo 2d positive", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoneg = new TH2F("histoneg","Histo 2d negative", 100, -1.5, 1.5, 80, -kPi, kPi);

  fHistopos->SetDirectory(0);
  fHistoneg->SetDirectory(0);

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
   // deleting variables
   delete fEsd;
   delete fPhysicsSelection;
   delete fHistopos;
   delete fHistoposAllEvents;
   delete fHistoneg;
   delete fHistonegAllEvents;

   delete fData;
   delete fDataAllEvents;
   delete fLego;
   delete fLegoAllEvents;
   delete fCalo3d;
   delete fCalo3dAllEvents;
   delete fGlv;

   delete fHisto2dv;
   delete fHisto2ds;
   delete fHisto2ds2;
   delete fHisto2dAllEventsv0;
   delete fHisto2dAllEventsv1;
   delete fHisto2dAllEventsv2;
   delete fHisto2dAllEventsv3;
   delete fHisto2dAllEventss0;
   delete fHisto2dAllEventss1;
   delete fHisto2dAllEventss2;
   delete fHisto2dAllEventss3;

   delete fAl;
   delete fHisto2dLegoOverlay;
   delete fHisto2dAllEventsLegoOverlay;
   delete fHisto2dAllEventsSlot;

   delete fEventSelector;
   delete fGButton;
   delete fB1;
   delete fB2;
}

namespace
{
  //____________________________________________________________________________
  Double_t getphi(Double_t phi)
  {
    // phi correction for alice

    if (phi > TMath::Pi()) {
      phi -= TMath::TwoPi();
    }
    return phi;
  }
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadData()
{
   // Load data from ESD tree
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

   fData->DataChanged();

   FilterData();

   return fData;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadAllData()
{
   // load data from all events ESD
   fHistoposAllEvents->Reset();
   fHistonegAllEvents->Reset();

   TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

   // Getting current tracks for each event, filling histograms
   for (int event = 0; event < t->GetEntries(); event++) {
      t->GetEntry(event);
         for (int n = 0; n < fEsd->GetNumberOfTracks(); ++n) {

            if (fEsd->GetTrack(n)->GetSign() > 0) {
               fHistoposAllEvents->Fill(fEsd->GetTrack(n)->Eta(),
                                          getphi(fEsd->GetTrack(n)->Phi()),
                                          fabs(fEsd->GetTrack(n)->Pt()));
            } else {
               fHistonegAllEvents->Fill(fEsd->GetTrack(n)->Eta(),
                                          getphi(fEsd->GetTrack(n)->Phi()),
                                          fabs(fEsd->GetTrack(n)->Pt()));
            }
         }
   }

   fDataAllEvents->DataChanged();

   return fDataAllEvents;
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
      fHistoposAllEvents->Reset();
      fHistonegAllEvents->Reset();

      TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

      // Getting current tracks for each event, filling histograms
      for (int event = 0; event < t->GetEntries(); event++) {

         t->GetEntry(event);

      const AliESDVertex *pv  = fEsd->GetPrimaryVertex();

      for (Int_t n = 0; n < pv->GetNIndices(); n++ )
      {
         AliESDtrack *at = fEsd->GetTrack(pv->GetIndices()[n]);

         if (at->GetSign() > 0) {
            fHistoposAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
         }
         if (at->GetSign() < 0) {
            fHistonegAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
         }
      }
      }
   } else {
      LoadAllData();
   }
   
   fDataAllEvents->DataChanged();

   // Max Pt threshold
   if (GetPtMaxAE() >= fMaxPtAE){
      for (Int_t binx = 1; binx <= 100; binx++) {
         for (Int_t biny = 1; biny <= 80; biny++) {
            if (fHistoposAllEvents->GetBinContent(binx, biny) >= fMaxPtAE)
            {
               fHistoposAllEvents->SetBinContent(binx, biny, fMaxPtAE);
            }
            if (fHistonegAllEvents->GetBinContent(binx, biny) >= fMaxPtAE)
            {
               fHistonegAllEvents->SetBinContent(binx, biny, fMaxPtAE);
            }           
         }
      }
   }

   // Positive only
   if ( fChargeIdAE == 2 ) fHistonegAllEvents->Reset();

   // Negative only
   if ( fChargeIdAE == 3 ) fHistoposAllEvents->Reset();

   fDataAllEvents->DataChanged();

   gEve->Redraw3D(kTRUE);

   return fDataAllEvents;
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
   if (fHisto2dv == 0) {
      TEveWindowSlot *fslot    = 0;
      TEveBrowser    *fbrowser = gEve->GetBrowser();

      fslot = TEveWindow::CreateWindowInTab(fbrowser->GetTabRight());
      fslot->MakeCurrent();
      fHisto2dv = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      fHisto2ds = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      fHisto2dv->AddScene(fHisto2ds);
      fHisto2dv->SetElementName("2D Lego Viewer");
      fHisto2ds->SetElementName("2D Lego Scene");

      fGlv = fHisto2dv->GetGLViewer();
      fHisto2dLegoOverlay = new TEveCaloLegoOverlay();
      fGlv->AddOverlayElement(fHisto2dLegoOverlay);
      fGlv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);

      fHisto2ds->AddElement(fLego);

      // move to real world coordinates
      fLego->InitMainTrans();
      Float_t sc = TMath::Min(fLego->GetEtaRng(), fLego->GetPhiRng());
      fLego->RefMainTrans().SetScale(sc, sc, sc);

      // set event handler to move from perspective to orthographic view.
      fGlv->SetEventHandler(new TEveLegoEventHandler(fGlv->GetGLWidget(), fGlv, fLego));

      fHisto2dLegoOverlay->SetCaloLego(fLego);
   }

   return fLego;
}

//______________________________________________________________________________
TEveCaloLego* AliEveLego::CreateHistoLego(TEveWindowSlot *slot)
{
   // Viewer initialization, tab creation
   if (fHisto2dAllEventsv0 == 0) {

      slot->MakeCurrent();
      fHisto2dAllEventsv0 = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      fHisto2dAllEventss0 = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      fHisto2dAllEventsv0->AddScene(fHisto2dAllEventss0);
      fHisto2dAllEventsv0->SetElementName("2D Lego Viewer");
      fHisto2dAllEventss0->SetElementName("2D Lego Scene");

      TGLViewer* glv = fHisto2dAllEventsv0->GetGLViewer();
      fHisto2dAllEventsLegoOverlay = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(fHisto2dAllEventsLegoOverlay);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);

      // Plotting histogram lego
      fLegoAllEvents = new TEveCaloLego(fDataAllEvents);
      fHisto2dAllEventss0->AddElement(fLegoAllEvents);

      // Move to real world coordinates
      fLegoAllEvents->InitMainTrans();
      Float_t sc = TMath::Min(fLegoAllEvents->GetEtaRng(), fLegoAllEvents->GetPhiRng());
      fLegoAllEvents->RefMainTrans().SetScale(sc, sc, sc);

      // Set event handler to move from perspective to orthographic view.
      glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, fLegoAllEvents));

      fHisto2dAllEventsLegoOverlay->SetCaloLego(fLegoAllEvents);
   }

   return fLegoAllEvents;
}

//______________________________________________________________________________
TEveCalo3D* AliEveLego::Create3DView()
{
   //initialization
   if (fHisto2ds2 == 0) {
      fHisto2ds2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      gEve->GetDefaultViewer()->AddScene(fHisto2ds2);
      fHisto2ds2->SetElementName("3D Histogram Scene");
      fHisto2ds2->AddElement(fCalo3d);
   }

   return fCalo3d;
}

//______________________________________________________________________________
TEveCalo3D* AliEveLego::Create3DView(TEveWindowSlot *slot)
{
   // creates a 3d view for the 3d histogram
   if ( fHisto2dAllEventsv1 == 0 ) {

      slot->MakeCurrent();
      fHisto2dAllEventsv1 = gEve->SpawnNewViewer("3D Histogram", "3D Histogram");
      fHisto2dAllEventss1 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      fHisto2dAllEventsv1->AddScene(fHisto2dAllEventss1);
      fHisto2dAllEventsv1->SetElementName("3D Histogram Viewer");
      fHisto2dAllEventss1->SetElementName("3D Histogram Scene");

      fCalo3dAllEvents = new TEveCalo3D(fDataAllEvents);

      fCalo3dAllEvents->SetBarrelRadius(550);
      fCalo3dAllEvents->SetEndCapPos(550);
      fHisto2dAllEventss1->AddElement(fCalo3dAllEvents);
   }

   return fCalo3dAllEvents;
}

//______________________________________________________________________________
void AliEveLego::CreateProjections(TEveWindowSlot* slot1, TEveWindowSlot* slot2)
{
   // create projections
   if (fHisto2dAllEventsv2 == 0) {

      slot1->MakeCurrent();
      fHisto2dAllEventsv2 = gEve->SpawnNewViewer("RPhi projection", "RPhi projection");
      fHisto2dAllEventss2 = gEve->SpawnNewScene("RPhi projection", "RPhi projection");
      fHisto2dAllEventsv2->AddScene(fHisto2dAllEventss2);
      fHisto2dAllEventsv2->SetElementName("RPhi Projection Viewer");
      fHisto2dAllEventss2->SetElementName("RPhi Projection Scene");

      TEveProjectionManager* mng1 = new TEveProjectionManager();
      mng1->SetProjection(TEveProjection::kPT_RPhi);

      TEveProjectionAxes* axeghisto2dAllEventss1 = new TEveProjectionAxes(mng1);
      fHisto2dAllEventss2->AddElement(axeghisto2dAllEventss1);
      TEveCalo2D* fcalo2d1 = (TEveCalo2D*) mng1->ImportElements(fCalo3dAllEvents);
      fHisto2dAllEventss2->AddElement(fcalo2d1);

      fHisto2dAllEventsv2->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   }

   if (fHisto2dAllEventsv3 == 0) {

      slot2->MakeCurrent();
      fHisto2dAllEventsv3 = gEve->SpawnNewViewer("RhoZ projection", "RhoZ projection");
      fHisto2dAllEventss3 = gEve->SpawnNewScene("RhoZ projection", "RhoZ projection");
      fHisto2dAllEventsv3->AddScene(fHisto2dAllEventss3);
      fHisto2dAllEventsv3->SetElementName("RhoZ Projection Viewer");
      fHisto2dAllEventss3->SetElementName("RhoZ Projection Viewer");

      TEveProjectionManager* mng2 = new TEveProjectionManager();
      mng2->SetProjection(TEveProjection::kPT_RhoZ);

      TEveProjectionAxes* axeghisto2dAllEventss2 = new TEveProjectionAxes(mng2);
      fHisto2dAllEventss3->AddElement(axeghisto2dAllEventss2);
      TEveCalo2D* fcalo2d2 = (TEveCalo2D*) mng2->ImportElements(fCalo3dAllEvents);
      fHisto2dAllEventss3->AddElement(fcalo2d2);

      fHisto2dAllEventsv3->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   }

   return;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadAllEvents()
{
   // load all events data from ESD
   if ( fHisto2dAllEventsSlot == 0 ) {

      printf("Filling histogram...\n");
      TStopwatch timer;
      timer.Start();

      // Creating 2D histograms
      fHistoposAllEvents = new TH2F("fHistoposAllEvents","Histo 2d positive",
                                 100,-1.5,1.5,80,-kPi,kPi);
      fHistonegAllEvents = new TH2F("fHistonegAllEvents","Histo 2d negative",
                                 100,-1.5,1.5,80,-kPi,kPi);

      fHistoposAllEvents->SetDirectory(0);
      fHistonegAllEvents->SetDirectory(0);

      fDataAllEvents = new TEveCaloDataHist();
      fDataAllEvents->AddHistogram(fHistonegAllEvents);
      fDataAllEvents->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
      fDataAllEvents->AddHistogram(fHistoposAllEvents);
      fDataAllEvents->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
      fDataAllEvents->GetEtaBins()->SetTitleFont(120);
      fDataAllEvents->GetEtaBins()->SetTitle("h");
      fDataAllEvents->GetPhiBins()->SetTitleFont(120);
      fDataAllEvents->GetPhiBins()->SetTitle("f");
      fDataAllEvents->IncDenyDestroy();

      // Creating frames
      fHisto2dAllEventsSlot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
      TEveWindowPack* packH = fHisto2dAllEventsSlot->MakePack();
      packH->SetElementName("Projections");
      packH->SetHorizontal();
      packH->SetShowTitleBar(kFALSE);

      fHisto2dAllEventsSlot = packH->NewSlot();
      TEveWindowPack* pack0 = fHisto2dAllEventsSlot->MakePack();
      pack0->SetShowTitleBar(kFALSE);
      TEveWindowSlot*  slotLeftTop   = pack0->NewSlot();
      TEveWindowSlot* slotLeftBottom = pack0->NewSlot();

      fHisto2dAllEventsSlot = packH->NewSlot();
      TEveWindowPack* pack1 = fHisto2dAllEventsSlot->MakePack();
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
      timer.Stop();
      timer.Print();

   }

   return fDataAllEvents;
}

//______________________________________________________________________________
Float_t AliEveLego::GetPtMax()
{
   return fData->GetMaxVal(fLego->GetPlotEt());
}

//______________________________________________________________________________
Float_t AliEveLego::GetPtMaxAE()
{
   return fDataAllEvents->GetMaxVal(fLegoAllEvents->GetPlotEt());
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
   fDataAllEvents->SetSliceThreshold(0,val);
   fDataAllEvents->SetSliceThreshold(1,val);
   fDataAllEvents->DataChanged();

   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
void AliEveLego::SetEventSelection()
{
   // activate/deactivate info box
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
   // activate/deactivate info box
   if (show == 0)
   {
      gEve->GetDefaultGLViewer()->RemoveOverlayElement(fGButton);
      fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fGButton);
      fHisto2dv->GetGLViewer()->RemoveOverlayElement(fGButton);

      gEve->GetDefaultGLViewer()->RemoveOverlayElement(fB1);
      fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fB1);
      fHisto2dv->GetGLViewer()->RemoveOverlayElement(fB1);

      gEve->GetDefaultGLViewer()->RemoveOverlayElement(fB2);
      fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fB2);
      fHisto2dv->GetGLViewer()->RemoveOverlayElement(fB2);

   } else {

      //Collision candidate
      if (updateonly == kFALSE) {
         gEve->GetDefaultGLViewer()->AddOverlayElement(fGButton);
         fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fGButton);
         fHisto2dv->GetGLViewer()->AddOverlayElement(fGButton);
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
         fHisto2dv->GetGLViewer()->AddOverlayElement(fB1);

         gEve->GetDefaultGLViewer()->AddOverlayElement(fB2);
         fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fB2);
         fHisto2dv->GetGLViewer()->AddOverlayElement(fB2);
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
   // show trigger information
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


