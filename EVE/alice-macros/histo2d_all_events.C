/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

double pi = TMath::Pi();
TEveViewer *v_0   = 0;
TEveViewer *v_1   = 0;
TEveViewer *v_2   = 0;
TEveViewer *v_3   = 0;
TEveScene  *s_0   = 0;
TEveScene  *s_1   = 0;
TEveScene  *s_2   = 0;
TEveScene  *s_3   = 0;
TEveCaloLegoOverlay* lego_overlay_t = 0;
TEveWindowSlot* slot_t = 0;


TEveCaloDataHist* histo2d_all_events()
{

   TEveCaloDataHist* data_t;
   
   if(slot_t == 0){
      cout<<"Filling histogram..."<<endl;
   
      // Access to esdTree
      AliESDEvent* esd = AliEveEventManager::AssertESD();
      TTree* t = AliEveEventManager::GetMaster()->GetESDTree();

      // Creating 2D histograms
      TH2F *histopos_t = new TH2F("histopos_t","Histo 2d positive",
                                 100,-1.5,1.5,80,-pi,pi);
      TH2F *histoneg_t = new TH2F("histoneg_t","Histo 2d negative",
                                 100,-1.5,1.5,80,-pi,pi);

         // Getting current tracks for each event, filling histograms 
         for (int event = 0; event < t->GetEntries(); event++) {	
            t->GetEntry(event);
               for (int n = 0; n < esd->GetNumberOfTracks(); ++n) {    
  
                  if (esd->GetTrack(n)->GetSign() > 0) {
                     histopos_t->Fill(esd->GetTrack(n)->Eta(),
	      	                      getphi(esd->GetTrack(n)->Phi()),
                                      fabs(esd->GetTrack(n)->Pt()));
                  } else {
                     histoneg_t->Fill(esd->GetTrack(n)->Eta(),
                                      getphi(esd->GetTrack(n)->Phi()),
                                      fabs(esd->GetTrack(n)->Pt()));
                  }
               }
         }

      data_t = new TEveCaloDataHist();
      data_t->AddHistogram(histoneg_t);
      data_t->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
      data_t->AddHistogram(histopos_t);
      data_t->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
      data_t->GetEtaBins()->SetTitleFont(120);
      data_t->GetEtaBins()->SetTitle("h");
      data_t->GetPhiBins()->SetTitleFont(120);
      data_t->GetPhiBins()->SetTitle("f");
      data_t->IncDenyDestroy();

      // Creating frames
      slot_t = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
      TEveWindowPack* packH = slot_t->MakePack();
      packH->SetElementName("Projections");
      packH->SetHorizontal();
      packH->SetShowTitleBar(kFALSE);

      slot_t = packH->NewSlot();
      TEveWindowPack* pack0 = slot_t->MakePack();
      pack0->SetShowTitleBar(kFALSE);
      TEveWindowSlot*  slotLeftTop   = pack0->NewSlot();
      TEveWindowSlot* slotLeftBottom = pack0->NewSlot();

      slot_t = packH->NewSlot();
      TEveWindowPack* pack1 = slot_t->MakePack();
      pack1->SetShowTitleBar(kFALSE);
      TEveWindowSlot* slotRightTop    = pack1->NewSlot();
      TEveWindowSlot* slotRightBottom = pack1->NewSlot();

      // Creating viewers and scenes   
      TEveCalo3D* calo3d = Create_3D_view(data_t, slotLeftTop);
      Create_histo_lego(data_t, slotLeftBottom);
      Create_projections(data_t, calo3d, slotRightTop, slotRightBottom);

      gEve->Redraw3D(kTRUE);

      cout<<"Filling histogram... Finished"<<endl;
   }

   return data_t;
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
TEveCaloLego* Create_histo_lego(TEveCaloData* data, TEveWindowSlot* slot){

   TEveCaloLego* lego;

   // Viewer initialization, tab creation
   if (v_0 == 0) {

      TEveBrowser *browser = gEve->GetBrowser();
      slot->MakeCurrent();
      v_0 = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      s_0 = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      v_0->AddScene(s_0);
      v_0->SetElementName("2D Lego Viewer");
      s_0->SetElementName("2D Lego Scene");

      TGLViewer* glv = v_0->GetGLViewer();
      lego_overlay_t = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(lego_overlay_t);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);

      // Plotting histogram lego
      lego = new TEveCaloLego(data);
      s_0->AddElement(lego);

      // Move to real world coordinates
      lego->InitMainTrans();
      Float_t sc = TMath::Min(lego->GetEtaRng(), lego->GetPhiRng());
      lego->RefMainTrans().SetScale(sc, sc, sc);

      // Set event handler to move from perspective to orthographic view.
      glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));

      lego_overlay_t->SetCaloLego(lego);
   }

   return lego;
}

//______________________________________________________________________________
TEveCalo3D* Create_3D_view(TEveCaloData* data, TEveWindowSlot* slot){

   TEveCalo3D* calo3d;

   if (v_1 == 0) {
      
      TEveBrowser *browser = gEve->GetBrowser();
      slot->MakeCurrent();
      v_1 = gEve->SpawnNewViewer("3D Histogram", "3D Histogram");
      s_1 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      v_1->AddScene(s_1);
      v_1->SetElementName("3D Histogram Viewer");
      s_1->SetElementName("3D Histogram Scene");

      calo3d = new TEveCalo3D(data);
   
      calo3d->SetBarrelRadius(550);
      calo3d->SetEndCapPos(550);
      s_1->AddElement(calo3d);
   } 

   return calo3d;
}

//______________________________________________________________________________
void Create_projections(TEveCaloData* data, TEveCalo3D *calo3d, TEveWindowSlot* slot1, TEveWindowSlot* slot2){

   if (v_2 == 0) {
      
      TEveBrowser *browser = gEve->GetBrowser();
      slot1->MakeCurrent();
      v_2 = gEve->SpawnNewViewer("RPhi projection", "RPhi projection");
      s_2 = gEve->SpawnNewScene("RPhi projection", "RPhi projection");
      v_2->AddScene(s_2);
      v_2->SetElementName("RPhi Projection Viewer");
      s_2->SetElementName("RPhi Projection Scene");

      TEveProjectionManager* mng1 = new TEveProjectionManager();
      mng1->SetProjection(TEveProjection::kPT_RPhi);

      TEveProjectionAxes* axes_1 = new TEveProjectionAxes(mng1);
      s_2->AddElement(axes_1);
      TEveCalo2D* calo2d1 = (TEveCalo2D*) mng1->ImportElements(calo3d);
      s_2->AddElement(calo2d1);

      v_2->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   } 

   if (v_3 == 0) {
      
      TEveBrowser *browser = gEve->GetBrowser();
      slot2->MakeCurrent();
      v_3 = gEve->SpawnNewViewer("RhoZ projection", "RhoZ projection");
      s_3 = gEve->SpawnNewScene("RhoZ projection", "RhoZ projection");
      v_3->AddScene(s_3);
      v_3->SetElementName("RhoZ Projection Viewer");
      s_3->SetElementName("RhoZ Projection Viewer");

      TEveProjectionManager* mng2 = new TEveProjectionManager();
      mng2->SetProjection(TEveProjection::kPT_RhoZ);

      TEveProjectionAxes* axes_2 = new TEveProjectionAxes(mng2);
      s_3->AddElement(axes_2);
      TEveCalo2D* calo2d2 = (TEveCalo2D*) mng2->ImportElements(calo3d);
      s_3->AddElement(calo2d2);

      v_3->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   } 

   return;
}
