// $Id$
// Author: Stefano Carrazza 2010, CERN, stefano.carrazza@cern.ch

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
#include "TStopwatch.h"

//______________________________________________________________________________
// This class provides the following features:
//
// 1) 2D and 3D calorimeter like histograms of tracks distribution,
// providing particle selection by charge (positive and negative),
// and by specie: electrons, muons, pions, kaons and protons.
//
// 2) Histograms are plotted around detectors, allowing track selection: all tracks,
// only primary tracks. It allows pT maximum and minimum threashold.
//
// 3) 2D and 3D calorimeter like histograms of particles distribution obtained
// from all events.
//
// 4) Possibility to use AliPhysicsSelection during the all events histograms creation.
// It is also possible to switch between real data and simulation data (MC).
//

ClassImp(AliEveLego)
Double_t kPi = TMath::Pi();

//______________________________________________________________________________
AliEveLego::AliEveLego(const char* name) :
  TEveElementList(name),
  fIsMC(kFALSE),
  fCollisionCandidatesOnly(kFALSE),
  fParticleTypeId(0),
  fParticleTypeIdAE(0),
  fTracksId(1),
  fMaxPt(10000),
  fTracksIdAE(0),
  fMaxPtAE(10000),
  fEsd(0),
  fPhysicsSelection(0),
  fHistopos(0),
  fHistoposAllEvents(0),
  fHistoneg(0),
  fHistonegAllEvents(0),
  fHistoElectrons(0),
  fHistoElectronsAllEvents(0),
  fHistoMuons(0),
  fHistoMuonsAllEvents(0),
  fHistoPions(0),
  fHistoPionsAllEvents(0),
  fHistoKaons(0),
  fHistoKaonsAllEvents(0),
  fHistoProtons(0),
  fHistoProtonsAllEvents(0),
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
  fHisto2dAllEventsSlot(0)
{
  // Constructor.
  gEve->AddToListTree(this,0);

  // Get Current ESD event
  fEsd = AliEveEventManager::AssertESD();

  // Particles types per default
  fParticleTypeId = new Bool_t[7];
  fParticleTypeIdAE = new Bool_t[7];

  for (Int_t s = 0; s < 7; s++)
  {
    if (s > 1)
    {
      fParticleTypeId[s] = kFALSE;
      fParticleTypeIdAE[s] = kFALSE;
    } else {
      fParticleTypeId[s] = kTRUE;
      fParticleTypeIdAE[s] = kTRUE;
    }
  }

  // Loading Physics Selection to determine the collision candidates
  fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection->Initialize(fEsd);

  fHistopos       = new TH2F("histopos","Histo 2d positive", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoneg       = new TH2F("histoneg","Histo 2d negative", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoElectrons = new TH2F("histoele","Histo 2d electron", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoMuons     = new TH2F("histomuo","Histo 2d muons   ", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoPions     = new TH2F("histopio","Histo 2d pions   ", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoKaons     = new TH2F("histokao","Histo 2d kaons   ", 100, -1.5, 1.5, 80, -kPi, kPi);
  fHistoProtons   = new TH2F("histopro","Histo 2d protons ", 100, -1.5, 1.5, 80, -kPi, kPi);

  fHistopos->SetDirectory(0);
  fHistoneg->SetDirectory(0);
  fHistoElectrons->SetDirectory(0);
  fHistoMuons->SetDirectory(0);
  fHistoPions->SetDirectory(0);
  fHistoKaons->SetDirectory(0);
  fHistoProtons->SetDirectory(0);

  // Colors from get_pdg_color() in /alice-macros/kine_tracks.C
  static Color_t fElectro = 5;
  static Color_t fMuon = 30;
  static Color_t fPion = 3;
  static Color_t fKaon = 38;
  static Color_t fProton = 10;

  // Adding data to TEveCaloDataHist
  fData = new TEveCaloDataHist();
  fData->AddHistogram(fHistoneg);
  fData->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
  fData->AddHistogram(fHistopos);
  fData->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
  fData->AddHistogram(fHistoElectrons);
  fData->RefSliceInfo(2).Setup("Elect.:", 0, fElectro);
  fData->AddHistogram(fHistoMuons);
  fData->RefSliceInfo(3).Setup("Muons:", 0, fMuon);
  fData->AddHistogram(fHistoPions);
  fData->RefSliceInfo(4).Setup("Pions:", 0, fPion);
  fData->AddHistogram(fHistoKaons);
  fData->RefSliceInfo(5).Setup("Kaons:", 0, fKaon);
  fData->AddHistogram(fHistoProtons);
  fData->RefSliceInfo(6).Setup("Proto.:", 0, fProton);

  fData->GetEtaBins()->SetTitleFont(120);
  fData->GetEtaBins()->SetTitle("h");
  fData->GetPhiBins()->SetTitleFont(120);
  fData->GetPhiBins()->SetTitle("f");
  fData->IncDenyDestroy();

  // Setting up the position of the 3D calorimeter histogram view.
  fCalo3d = new TEveCalo3D(fData);
  fCalo3d->SetBarrelRadius(550);
  fCalo3d->SetEndCapPos(550);

  // Adding data to the Lego object
  fLego = new TEveCaloLego(fData);

  // Creating projections
  fAl = AliEveMultiView::Instance();
  fAl->ImportEventRPhi(fCalo3d);
  fAl->ImportEventRhoZ(fCalo3d);

  // Update viewers and scenes
  Update();
}

//______________________________________________________________________________
AliEveLego::~AliEveLego()
{
   // Deleting variables
   delete fEsd;
   delete fPhysicsSelection;
   delete fHistopos;
   delete fHistoposAllEvents;
   delete fHistoneg;
   delete fHistonegAllEvents;
   delete fHistoElectrons;
   delete fHistoElectronsAllEvents;
   delete fHistoMuons;
   delete fHistoMuonsAllEvents;
   delete fHistoPions;
   delete fHistoPionsAllEvents;
   delete fHistoKaons;
   delete fHistoKaonsAllEvents;
   delete fHistoProtons;
   delete fHistoProtonsAllEvents;
   delete[] fParticleTypeId;
   delete[] fParticleTypeIdAE;
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
}

namespace
{
  //____________________________________________________________________________
  Double_t getphi(Double_t phi)
  {
    // Phi correction for alice

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
   fHistoElectrons->Reset();
   fHistoMuons->Reset();
   fHistoPions->Reset();
   fHistoKaons->Reset();
   fHistoProtons->Reset();

   // Getting current tracks, filling histograms
   for (int n = 0; n < fEsd->GetNumberOfTracks(); ++n) {

     AliESDtrack *track = fEsd->GetTrack(n);
     const Double_t sign = fEsd->GetTrack(n)->GetSign();
     const Int_t prob = GetParticleType(track);     

     // Filling histograms
     if (sign > 0)
       fHistopos->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

     if (sign < 0)
       fHistoneg->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

     if (prob == 0)
       fHistoElectrons->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

     if (prob == 1)
       fHistoMuons->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

     if (prob == 2)
       fHistoPions->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

     if (prob == 3)
       fHistoKaons->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

     if (prob == 4)
       fHistoProtons->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));
   }

   fData->DataChanged();

   FilterData();

   return fData;
}

//______________________________________________________________________________
TEveCaloDataHist* AliEveLego::LoadAllData()
{
   // Load data from all events ESD
   fHistoposAllEvents->Reset();
   fHistonegAllEvents->Reset();
   fHistoElectronsAllEvents->Reset();
   fHistoMuonsAllEvents->Reset();
   fHistoPionsAllEvents->Reset();
   fHistoKaonsAllEvents->Reset();
   fHistoProtonsAllEvents->Reset();

   TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

   // Getting current tracks for each event, filling histograms
   Int_t fAcceptedEvents = 0;
   for (int event = 0; event < t->GetEntries(); event++) {
      t->GetEntry(event);

      if (fCollisionCandidatesOnly == kTRUE)
        if (fPhysicsSelection->IsCollisionCandidate(fEsd) == kFALSE) continue;

      fAcceptedEvents++;

      for (int n = 0; n < fEsd->GetNumberOfTracks(); ++n) {

           AliESDtrack *track = fEsd->GetTrack(n);
           const Double_t sign = track->GetSign();
           const Int_t prob = GetParticleType(track);

           if (sign > 0)
             fHistoposAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

           if (sign < 0)
             fHistonegAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

           if (prob == 0)
             fHistoElectronsAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

           if (prob == 1)
             fHistoMuonsAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

           if (prob == 2)
             fHistoPionsAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

           if (prob == 3)
             fHistoKaonsAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));

           if (prob == 4)
             fHistoProtonsAllEvents->Fill(track->Eta(), getphi(track->Phi()), fabs(track->Pt()));
        }
   }

   // Setting the current view to the first event
   t->GetEntry(0);

   // Usefull information,
   // with this we can estimate the event efficiency
   printf("Number of events loaded: %i, with AliPhysicsSelection: %i\n",fAcceptedEvents,fCollisionCandidatesOnly);

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
      fHistoElectrons->Reset();
      fHistoMuons->Reset();
      fHistoPions->Reset();
      fHistoKaons->Reset();
      fHistoProtons->Reset();

      const AliESDVertex *pv  = fEsd->GetPrimaryVertex();

      for (Int_t n = 0; n < pv->GetNIndices(); n++ )
      {
         AliESDtrack *at = fEsd->GetTrack(pv->GetIndices()[n]);
         const Double_t sign = at->GetSign();
         const Int_t prob = GetParticleType(at);

         if (sign > 0)
           fHistopos->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

         if (sign < 0)
           fHistoneg->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

         if (prob == 0)
           fHistoElectrons->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

         if (prob == 1)
           fHistoMuons->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

         if (prob == 2)
           fHistoPions->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

         if (prob == 3)
           fHistoKaons->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

         if (prob == 4)
           fHistoProtons->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
      }
   }

   fData->DataChanged();

   // Max Pt threshold
   if (GetPtMax() >= fMaxPt){
      for (Int_t binx = 1; binx <= 100; binx++) {
         for (Int_t biny = 1; biny <= 80; biny++) {

           if (fHistopos->GetBinContent(binx, biny) >= fMaxPt)
             fHistopos->SetBinContent(binx, biny, fMaxPt);

           if (fHistoneg->GetBinContent(binx, biny) >= fMaxPt)
             fHistoneg->SetBinContent(binx, biny, fMaxPt);

           if (fHistoElectrons->GetBinContent(binx, biny) >= fMaxPt)
             fHistoElectrons->SetBinContent(binx, biny, fMaxPt);

           if (fHistoMuons->GetBinContent(binx, biny) >= fMaxPt)
             fHistoMuons->SetBinContent(binx, biny, fMaxPt);

           if (fHistoPions->GetBinContent(binx, biny) >= fMaxPt)
             fHistoPions->SetBinContent(binx, biny, fMaxPt);

           if (fHistoKaons->GetBinContent(binx, biny) >= fMaxPt)
             fHistoKaons->SetBinContent(binx, biny, fMaxPt);

           if (fHistoProtons->GetBinContent(binx, biny) >= fMaxPt)
             fHistoProtons->SetBinContent(binx, biny, fMaxPt);
         }
      }
   }

   // Particle type filter
   if ( fParticleTypeId[0] == kFALSE) fHistopos->Reset();
   if ( fParticleTypeId[1] == kFALSE) fHistoneg->Reset();
   if ( fParticleTypeId[2] == kFALSE) fHistoElectrons->Reset();
   if ( fParticleTypeId[3] == kFALSE) fHistoMuons->Reset();
   if ( fParticleTypeId[4] == kFALSE) fHistoPions->Reset();
   if ( fParticleTypeId[5] == kFALSE) fHistoKaons->Reset();
   if ( fParticleTypeId[6] == kFALSE) fHistoProtons->Reset();

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
      fHistoElectronsAllEvents->Reset();
      fHistoMuonsAllEvents->Reset();
      fHistoPionsAllEvents->Reset();
      fHistoKaonsAllEvents->Reset();
      fHistoProtonsAllEvents->Reset();

      TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

      // Getting current tracks for each event, filling histograms
      Int_t fAcceptedEvents = 0;
      for (int event = 0; event < t->GetEntries(); event++) {

        t->GetEntry(event);
        const AliESDVertex *pv  = fEsd->GetPrimaryVertex();

        if (fCollisionCandidatesOnly == kTRUE)
          if (fPhysicsSelection->IsCollisionCandidate(fEsd) == kFALSE) continue;

        fAcceptedEvents++;

        for (Int_t n = 0; n < pv->GetNIndices(); n++ )
        {
          AliESDtrack *at = fEsd->GetTrack(pv->GetIndices()[n]);
          const Double_t sign = at->GetSign();
          const Int_t prob = GetParticleType(at);

          if (sign > 0)
            fHistoposAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

          if (sign < 0)
            fHistonegAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

          if (prob == 0)
            fHistoElectronsAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

          if (prob == 1)
            fHistoMuonsAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

          if (prob == 2)
            fHistoPionsAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

          if (prob == 3)
            fHistoKaonsAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));

          if (prob == 4)
            fHistoProtonsAllEvents->Fill(at->Eta(), getphi(at->Phi()), fabs(at->Pt()));
        }
      }
      t->GetEntry(0);
      printf("Number of events loaded: %i, with AliPhysicsSelection: %i\n",fAcceptedEvents,fCollisionCandidatesOnly);

   } else {

      LoadAllData();

    }
   
   fDataAllEvents->DataChanged();

   // Max Pt threshold
   if (GetPtMaxAE() >= fMaxPtAE){
      for (Int_t binx = 1; binx <= 100; binx++) {
         for (Int_t biny = 1; biny <= 80; biny++) {

            if (fHistoposAllEvents->GetBinContent(binx, biny) >= fMaxPtAE)            
              fHistoposAllEvents->SetBinContent(binx, biny, fMaxPtAE);

            if (fHistonegAllEvents->GetBinContent(binx, biny) >= fMaxPtAE)            
              fHistonegAllEvents->SetBinContent(binx, biny, fMaxPtAE);

            if (fHistoElectronsAllEvents->GetBinContent(binx, biny) >= fMaxPt)
              fHistoElectronsAllEvents->SetBinContent(binx, biny, fMaxPt);

            if (fHistoMuonsAllEvents->GetBinContent(binx, biny) >= fMaxPt)
              fHistoMuonsAllEvents->SetBinContent(binx, biny, fMaxPt);

            if (fHistoPionsAllEvents->GetBinContent(binx, biny) >= fMaxPt)
              fHistoPionsAllEvents->SetBinContent(binx, biny, fMaxPt);

            if (fHistoKaonsAllEvents->GetBinContent(binx, biny) >= fMaxPt)
              fHistoKaonsAllEvents->SetBinContent(binx, biny, fMaxPt);

            if (fHistoProtonsAllEvents->GetBinContent(binx, biny) >= fMaxPt)
              fHistoProtonsAllEvents->SetBinContent(binx, biny, fMaxPt);
         }
      }
   }

   // Particles species and charges filter
   if ( fParticleTypeIdAE[0] == kFALSE ) fHistoposAllEvents->Reset();
   if ( fParticleTypeIdAE[1] == kFALSE ) fHistonegAllEvents->Reset();
   if ( fParticleTypeIdAE[2] == kFALSE ) fHistoElectronsAllEvents->Reset();
   if ( fParticleTypeIdAE[3] == kFALSE ) fHistoMuonsAllEvents->Reset();
   if ( fParticleTypeIdAE[4] == kFALSE ) fHistoPionsAllEvents->Reset();
   if ( fParticleTypeIdAE[5] == kFALSE ) fHistoKaonsAllEvents->Reset();
   if ( fParticleTypeIdAE[6] == kFALSE ) fHistoProtonsAllEvents->Reset();

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
   // Initialization
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
   // Creates a 3d view for the 3d histogram
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
   // Create projections
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

   if (fHisto2dAllEventsv3 == 0)
   {
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
      fHistoElectronsAllEvents = new TH2F("fHistoElectronsAllEvents","Histo 2d electrons",
                                       100,-1.5,1.5,80,-kPi,kPi);
      fHistoMuonsAllEvents = new TH2F("fHistoMuonsAllEvents","Histo 2d muons",
                                       100,-1.5,1.5,80,-kPi,kPi);
      fHistoPionsAllEvents = new TH2F("fHistoPionsAllEvents","Histo 2d pions",
                                       100,-1.5,1.5,80,-kPi,kPi);
      fHistoKaonsAllEvents = new TH2F("fHistoKaonsAllEvents","Histo 2d kaons",
                                       100,-1.5,1.5,80,-kPi,kPi);
      fHistoProtonsAllEvents = new TH2F("fHistoProtonsAllEvents","Histo 2d protons",
                                       100,-1.5,1.5,80,-kPi,kPi);

      fHistoposAllEvents->SetDirectory(0);
      fHistonegAllEvents->SetDirectory(0);
      fHistoElectronsAllEvents->SetDirectory(0);
      fHistoMuonsAllEvents->SetDirectory(0);
      fHistoPionsAllEvents->SetDirectory(0);
      fHistoKaonsAllEvents->SetDirectory(0);
      fHistoProtonsAllEvents->SetDirectory(0);

      // colors from get_pdg_color() in /alice-macros/kine_tracks.C
      static Color_t fElectro = 5;
      static Color_t fMuon = 30;
      static Color_t fPion = 3;
      static Color_t fKaon = 38;
      static Color_t fProton = 10;

      fDataAllEvents = new TEveCaloDataHist();
      fDataAllEvents->AddHistogram(fHistonegAllEvents);
      fDataAllEvents->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
      fDataAllEvents->AddHistogram(fHistoposAllEvents);
      fDataAllEvents->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
      fDataAllEvents->AddHistogram(fHistoElectronsAllEvents);
      fDataAllEvents->RefSliceInfo(2).Setup("Electrons:", 0, fElectro);
      fDataAllEvents->AddHistogram(fHistoMuonsAllEvents);
      fDataAllEvents->RefSliceInfo(3).Setup("Muons:", 0, fMuon);
      fDataAllEvents->AddHistogram(fHistoPionsAllEvents);
      fDataAllEvents->RefSliceInfo(4).Setup("Pions:", 0, fPion);
      fDataAllEvents->AddHistogram(fHistoKaonsAllEvents);
      fDataAllEvents->RefSliceInfo(5).Setup("Kaons:", 0, fKaon);
      fDataAllEvents->AddHistogram(fHistoProtonsAllEvents);
      fDataAllEvents->RefSliceInfo(6).Setup("Protons:", 0, fProton);

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

      FilterAllData();

      gEve->Redraw3D(kTRUE);

      printf("Filling histogram... Finished\n");
      timer.Stop();
      timer.Print();
   }
   return fDataAllEvents;
}

//______________________________________________________________________________
void AliEveLego::ApplyParticleTypeSelectionAE()
{
  // Reload all events applying particle type selection
  FilterAllData();
}

//______________________________________________________________________________
Float_t AliEveLego::GetPtMax()
{
   // Return pT maximum
   return fData->GetMaxVal(fLego->GetPlotEt());
}

//______________________________________________________________________________
Float_t AliEveLego::GetPtMaxAE()
{
   // Return pT max from all events
   return fDataAllEvents->GetMaxVal(fLegoAllEvents->GetPlotEt());
}

//______________________________________________________________________________
Int_t AliEveLego::GetParticleType(AliESDtrack *track)
{
  // Determine the particle type
	Double_t prob[5]={0.};
  track->GetESDpid(prob);

  Double_t max = prob[0];
  Int_t index = 0;

  for (Int_t i = 1 ; i < 5; i++)
  {
    if (prob[i] > max){
      max = prob[i];
      index = i;
    }
  }
  return index;
}

//______________________________________________________________________________
void AliEveLego::SetParticleType(Int_t id, Bool_t status)
{
  // Activate/deactivate particles types
  fParticleTypeId[id] = status;

  Update();
}

//______________________________________________________________________________
void AliEveLego::SetParticleTypeAE(Int_t id, Bool_t status)
{
  // Activate/deactivate particles types
  fParticleTypeIdAE[id] = status;
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
   fData->SetSliceThreshold(2,val);
   fData->SetSliceThreshold(3,val);
   fData->SetSliceThreshold(4,val);
   fData->SetSliceThreshold(5,val);
   fData->SetSliceThreshold(6,val);
   fData->DataChanged();

   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
void AliEveLego::SetThresholdAE(Double_t val)
{
   // Setting up the new threshold for all histograms
   fDataAllEvents->SetSliceThreshold(0,val);
   fDataAllEvents->SetSliceThreshold(1,val);
   fDataAllEvents->SetSliceThreshold(2,val);
   fDataAllEvents->SetSliceThreshold(3,val);
   fDataAllEvents->SetSliceThreshold(4,val);
   fDataAllEvents->SetSliceThreshold(5,val);
   fDataAllEvents->SetSliceThreshold(6,val);
   fDataAllEvents->DataChanged();

   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
void AliEveLego::SwitchDataType(Bool_t status)
{
  // Activate/deactivate MC / real data type
  fIsMC = status;

  // Removing defaul physics selection
  delete fPhysicsSelection;
  fPhysicsSelection = NULL;

  // Re-initialization of physics selection
  fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection->SetAnalyzeMC(fIsMC);
  fPhysicsSelection->Initialize(fEsd);
  FilterAllData();
}

//______________________________________________________________________________
void AliEveLego::SetCollisionCandidatesOnly()
{
  // Activate/deactivate MC / real data type
  if (fCollisionCandidatesOnly == 0)
  {
    fCollisionCandidatesOnly = 1;
  } else {
    fCollisionCandidatesOnly = 0;
  }
  FilterAllData();
}

/******************************************************************************/


