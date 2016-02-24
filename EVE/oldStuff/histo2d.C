// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGLViewer.h>
#include <TGLWidget.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <TEveBrowser.h>
#include <TEveCalo.h>
#include <TEveCaloData.h>
#include <TEveCaloLegoOverlay.h>
#include <TEveLegoEventHandler.h>
#include <TEveManager.h>
#include <TEveScene.h>
#include <TEveTrans.h>
#include <TEveViewer.h>
#include <TEveWindow.h>

#include <AliESDEvent.h>
#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

double pi = TMath::Pi();
TEveViewer *g_histo2d_v   = 0;
TEveScene  *g_histo2d_s   = 0;
TEveScene  *g_histo2d_s2  = 0;
TEveCaloLegoOverlay* g_histo2d_lego_overlay = 0;

Double_t GetPhi(Double_t phi);
TEveCaloLego* CreateHistoLego(TEveCaloData* data);
TEveCalo3D* Create3DView(TEveCaloData* data);
AliEveMultiView* CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d);

TEveCaloDataHist* histo2d()
{ 

   // Access to esdTree
   AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();

   // Creating 2D histograms
   TH2F *histopos = new TH2F("histopos","Histo 2d positive",100,-1.5,1.5,80,-pi,pi);
   TH2F *histoneg = new TH2F("histoneg","Histo 2d negative",100,-1.5,1.5,80,-pi,pi);

   Info("histo2d", "Event: %d, Number of tracks: %d\n", AliEveEventManager::Instance()->GetEventId(), esd->GetNumberOfTracks() );

   // Getting current tracks, filling histograms 
   for ( int n = 0; n < esd->GetNumberOfTracks(); ++n ) {    
  
      if (esd->GetTrack(n)->GetSign() > 0) {
         histopos->Fill(esd->GetTrack(n)->Eta(),
	      	        GetPhi(esd->GetTrack(n)->Phi()),
                        fabs(esd->GetTrack(n)->Pt()));
      } else {
         histoneg->Fill(esd->GetTrack(n)->Eta(),
                        GetPhi(esd->GetTrack(n)->Phi()),
                        fabs(esd->GetTrack(n)->Pt()));
      }
   }

   TEveCaloDataHist* data = new TEveCaloDataHist();
   AliEveEventManager::Instance()->RegisterTransient(data);
   
   data->AddHistogram(histoneg);
   data->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
   data->AddHistogram(histopos);
   data->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
   data->GetEtaBins()->SetTitleFont(120);
   data->GetEtaBins()->SetTitle("h");
   data->GetPhiBins()->SetTitleFont(120);
   data->GetPhiBins()->SetTitle("f");
   data->IncDenyDestroy();

   // Plotting the lego histogram in a new tab
   CreateHistoLego(data);
   
   // Plotting the 3D histogram
   TEveCalo3D *calo3d = Create3DView(data);

   // Plotting projections RPhi and RhoZ
   CreateProjections(data, calo3d);
   
   gEve->Redraw3D();

   return data;
}

//______________________________________________________________________________
Double_t GetPhi(Double_t phi)
{
   if (phi > pi) {
      phi -= 2*pi;
   }
   return phi;
}

//______________________________________________________________________________
TEveCaloLego* CreateHistoLego(TEveCaloData* data){

   TGLViewer* glv;

   // Viewer initialization, tab creation
   if ( g_histo2d_v == 0 ) {
      TEveWindowSlot *slot    = 0;
      TEveBrowser    *browser = gEve->GetBrowser();

      slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
      slot->MakeCurrent();
      g_histo2d_v = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      g_histo2d_s = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      g_histo2d_v->AddScene(g_histo2d_s);
      g_histo2d_v->SetElementName("2D Lego Viewer");
      g_histo2d_s->SetElementName("2D Lego Scene");

      glv = g_histo2d_v->GetGLViewer();
      g_histo2d_lego_overlay = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(g_histo2d_lego_overlay);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
   } else {
      glv = g_histo2d_v->GetGLViewer(); 
   }
   
   //plotting histo
   TEveCaloLego* lego = new TEveCaloLego(data);
   g_histo2d_s->AddElement(lego);
   AliEveEventManager::Instance()->RegisterTransient(lego);

   // move to real world coordinates
   lego->InitMainTrans();
   Float_t sc = TMath::Min(lego->GetEtaRng(), lego->GetPhiRng());
   lego->RefMainTrans().SetScale(sc, sc, sc);

   // set event handler to move from perspective to orthographic view.
   glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));

   g_histo2d_lego_overlay->SetCaloLego(lego);

   return lego;
}

//______________________________________________________________________________
TEveCalo3D* Create3DView(TEveCaloData* data){
  
   //initialization
   if ( g_histo2d_s2 == 0 ) {
      g_histo2d_s2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      gEve->GetDefaultViewer()->AddScene(g_histo2d_s2);
      g_histo2d_s2->SetElementName("3D Histogram Scene");
   }
 
   TEveCalo3D* calo3d = new TEveCalo3D(data);
   AliEveEventManager::Instance()->RegisterTransient(calo3d);
   
   calo3d->SetBarrelRadius(550);
   calo3d->SetEndCapPos(550);
   g_histo2d_s2->AddElement(calo3d);
 
   return calo3d;
}

//______________________________________________________________________________
AliEveMultiView* CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d){

   AliEveMultiView *al = AliEveMultiView::Instance();
   al->ImportEventRPhi(calo3d);
   al->ImportEventRhoZ(calo3d);

   return al;
}
