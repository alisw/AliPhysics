/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

double pi = TMath::Pi();
TEveViewer *v   = 0;
TEveScene  *s   = 0;
TEveScene  *s2  = 0;
TEveCaloLegoOverlay* lego_overlay = 0;


TEveCaloDataHist* histo2d()
{

   // Access to esdTree
   AliESDEvent* esd = AliEveEventManager::AssertESD();

   // Creating 2D histograms
   TH2F *histopos = new TH2F("histopos","Histo 2d positive",100,-1.5,1.5,80,-pi,pi);
   TH2F *histoneg = new TH2F("histoneg","Histo 2d negative",100,-1.5,1.5,80,-pi,pi);

   cout<<"Event: "<<AliEveEventManager::GetMaster()->GetEventId()
       <<", Number of tracks: "<<esd->GetNumberOfTracks()<<endl;

   // Getting current tracks, filling histograms 
   for (int n = 0; n < esd->GetNumberOfTracks(); ++n) {    
  
      if (esd->GetTrack(n)->GetSign() > 0) {
         histopos->Fill(esd->GetTrack(n)->Eta(),
	      	        getphi(esd->GetTrack(n)->Phi()),
                        fabs(esd->GetTrack(n)->Pt()));
      } else {
         histoneg->Fill(esd->GetTrack(n)->Eta(),
                        getphi(esd->GetTrack(n)->Phi()),
                        fabs(esd->GetTrack(n)->Pt()));
      }
   }

   TEveCaloDataHist* data = new TEveCaloDataHist();
   AliEveEventManager::RegisterTransient(data);
   
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
   Create_histo_lego(data);
   
   // Plotting the 3D histogram
   TEveCalo3D *calo3d = Create_3D_view(data);

   // Plotting projections RPhi and RhoZ
   Create_projections(data, calo3d);
   
   gEve->Redraw3D(kTRUE);

   return data;
}

//______________________________________________________________________________
Double_t getphi(Double_t phi)
{
   if (phi > pi) {
      phi -= 2*pi;
   }
   return phi;
}

//______________________________________________________________________________
TEveCaloLego* Create_histo_lego(TEveCaloData* data){

   TGLViewer* glv;

   // Viewer initialization, tab creation
   if (v == 0) {
      TEveWindowSlot *slot    = 0;
      TEveBrowser    *browser = gEve->GetBrowser();

      slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
      slot->MakeCurrent();
      v = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      s = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      v->AddScene(s);
      v->SetElementName("2D Lego Viewer");
      s->SetElementName("2D Lego Scene");

      glv = v->GetGLViewer();
      lego_overlay = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(lego_overlay);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
   } else {
      glv = v->GetGLViewer(); 
   }
   
   //plotting histo
   TEveCaloLego* lego = new TEveCaloLego(data);
   s->AddElement(lego);
   AliEveEventManager::RegisterTransient(lego);

   // move to real world coordinates
   lego->InitMainTrans();
   Float_t sc = TMath::Min(lego->GetEtaRng(), lego->GetPhiRng());
   lego->RefMainTrans().SetScale(sc, sc, sc);

   // set event handler to move from perspective to orthographic view.
   glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));

   lego_overlay->SetCaloLego(lego);

   return lego;
}

//______________________________________________________________________________
TEveCalo3D* Create_3D_view(TEveCaloData* data){
  
   //initialization
   if (s2 == 0) {
      s2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      gEve->GetDefaultViewer()->AddScene(s2);
      s2->SetElementName("3D Histogram Scene");
   }
 
   TEveCalo3D* calo3d = new TEveCalo3D(data);
   AliEveEventManager::RegisterTransient(calo3d);
   
   calo3d->SetBarrelRadius(550);
   calo3d->SetEndCapPos(550);
   s2->AddElement(calo3d);
 
   return calo3d;
}

//______________________________________________________________________________
AliEveMultiView* Create_projections(TEveCaloData* data, TEveCalo3D *calo3d){

   AliEveMultiView *al = AliEveMultiView::Instance();
   al->ImportEventRPhi(calo3d);
   al->ImportEventRhoZ(calo3d);

   return al;
}
