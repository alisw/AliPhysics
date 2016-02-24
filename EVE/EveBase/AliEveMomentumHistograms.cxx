// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <AliEveMomentumHistograms.h>

#include <AliESDEvent.h>
#include <AliEveEventManager.h>

#include <TGLWidget.h>
#include <TH2.h>
#include <TEveLegoEventHandler.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TEveTrans.h>
#include <TEveProjectionAxes.h>

AliEveMomentumHistograms::AliEveMomentumHistograms()
{
    pi = TMath::Pi();
    fHistoViewer   = 0;
    g_histo2d_s   = 0;
    g_histo2d_s2  = 0;
    g_histo2d_lego_overlay = 0;
    
    g_histo2d_all_events_v0 = 0;
    g_histo2d_all_events_v1 = 0;
    g_histo2d_all_events_v2 = 0;
    g_histo2d_all_events_v3 = 0;
    g_histo2d_all_events_s0 = 0;
    g_histo2d_all_events_s1 = 0;
    g_histo2d_all_events_s2 = 0;
    g_histo2d_all_events_s3 = 0;
    g_histo2d_all_events_lego_overlay = 0;
    g_histo2d_all_events_slot = 0;
}

AliEveMomentumHistograms::~AliEveMomentumHistograms()
{
    if(fHistoViewer){delete fHistoViewer; fHistoViewer = 0;}
    if(g_histo2d_s){delete g_histo2d_s; g_histo2d_s = 0;}
    if(g_histo2d_s2){delete g_histo2d_s2; g_histo2d_s2 = 0;}
    if(g_histo2d_lego_overlay){delete g_histo2d_lego_overlay; g_histo2d_lego_overlay = 0;}
    
    if(g_histo2d_all_events_v0){delete g_histo2d_all_events_v0;g_histo2d_all_events_v0 = 0;}
    if(g_histo2d_all_events_v1){delete g_histo2d_all_events_v1;g_histo2d_all_events_v1 = 0;}
    if(g_histo2d_all_events_v2){delete g_histo2d_all_events_v2;g_histo2d_all_events_v2 = 0;}
    if(g_histo2d_all_events_v3){delete g_histo2d_all_events_v3;g_histo2d_all_events_v3 = 0;}
    if(g_histo2d_all_events_s0){delete g_histo2d_all_events_s0;g_histo2d_all_events_s0 = 0;}
    if(g_histo2d_all_events_s1){delete g_histo2d_all_events_s1;g_histo2d_all_events_s1 = 0;}
    if(g_histo2d_all_events_s2){delete g_histo2d_all_events_s2;g_histo2d_all_events_s2 = 0;}
    if(g_histo2d_all_events_s3){delete g_histo2d_all_events_s3;g_histo2d_all_events_s3 = 0;}
    
    if(g_histo2d_all_events_lego_overlay){delete g_histo2d_all_events_lego_overlay;g_histo2d_all_events_lego_overlay =0;}
    if(g_histo2d_all_events_slot){delete g_histo2d_all_events_slot;g_histo2d_all_events_slot = 0;}
}

TEveCaloDataHist* AliEveMomentumHistograms::Draw()
{
    
    // Access to esdTree
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    
    // Creating 2D histograms
    TH2F *histopos = new TH2F("histopos","Histo 2d positive",100,-1.5,1.5,80,-pi,pi);
    TH2F *histoneg = new TH2F("histoneg","Histo 2d negative",100,-1.5,1.5,80,-pi,pi);
    
    Info("AliEveMomentumHistograms::Draw", "Event: %d, Number of tracks: %d\n", AliEveEventManager::Instance()->GetEventId(), esd->GetNumberOfTracks() );
    
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
    AliEveMultiView *al = AliEveMultiView::Instance();
    al->ImportEventRPhi(calo3d);
    al->ImportEventRhoZ(calo3d);
    
    gEve->Redraw3D();
    
    return data;
}

Double_t AliEveMomentumHistograms::GetPhi(Double_t phi)
{
    if (phi > pi) {phi -= 2*pi;}
    return phi;
}

void AliEveMomentumHistograms::CreateHistoLego(TEveCaloData* data)
{
    TGLViewer* glv;
    
    // Viewer initialization, tab creation
    if ( fHistoViewer == 0 )
    {
        TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
        slot->MakeCurrent();
        
        fHistoViewer = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
        g_histo2d_s = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
        
        fHistoViewer->AddScene(g_histo2d_s);
        fHistoViewer->SetElementName("2D Lego Viewer");
        g_histo2d_s->SetElementName("2D Lego Scene");
        
        glv = fHistoViewer->GetGLViewer();
        g_histo2d_lego_overlay = new TEveCaloLegoOverlay();
        glv->AddOverlayElement(g_histo2d_lego_overlay);
        glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
    }
    else
    {
        glv = fHistoViewer->GetGLViewer();
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
}

TEveCalo3D* AliEveMomentumHistograms::Create3DView(TEveCaloData* data)
{
    
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

TEveCaloDataHist* AliEveMomentumHistograms::DrawAllEvents()
{
    
    TEveCaloDataHist* data_t;
    
    if ( g_histo2d_all_events_slot == 0 ) {
        Info("AliEveMomentumHistograms::DrawAllEvents", "Filling histogram...");
        
        // Access to esdTree
        AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
        TTree* t = AliEveEventManager::Instance()->GetESDTree();
        
        // Creating 2D histograms
        TH2F *histopos_t = new TH2F("histopos_t","Histo 2d positive",
                                    100,-1.5,1.5,80,-pi,pi);
        TH2F *histoneg_t = new TH2F("histoneg_t","Histo 2d negative",
                                    100,-1.5,1.5,80,-pi,pi);
        
        // Getting current tracks for each event, filling histograms
        for ( int event = 0; event < t->GetEntries(); event++ ) {
            t->GetEntry(event);
            for ( int n = 0; n < esd->GetNumberOfTracks(); ++n ) {
                
                if ( esd->GetTrack(n)->GetSign() > 0 ) {
                    histopos_t->Fill(esd->GetTrack(n)->Eta(),
                                     GetPhi(esd->GetTrack(n)->Phi()),
                                     fabs(esd->GetTrack(n)->Pt()));
                } else {
                    histoneg_t->Fill(esd->GetTrack(n)->Eta(),
                                     GetPhi(esd->GetTrack(n)->Phi()),
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
        g_histo2d_all_events_slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
        TEveWindowPack* packH = g_histo2d_all_events_slot->MakePack();
        packH->SetElementName("Projections");
        packH->SetHorizontal();
        packH->SetShowTitleBar(kFALSE);
        
        g_histo2d_all_events_slot = packH->NewSlot();
        TEveWindowPack* pack0 = g_histo2d_all_events_slot->MakePack();
        pack0->SetShowTitleBar(kFALSE);
        TEveWindowSlot*  slotLeftTop   = pack0->NewSlot();
        TEveWindowSlot* slotLeftBottom = pack0->NewSlot();
        
        g_histo2d_all_events_slot = packH->NewSlot();
        TEveWindowPack* pack1 = g_histo2d_all_events_slot->MakePack();
        pack1->SetShowTitleBar(kFALSE);
        TEveWindowSlot* slotRightTop    = pack1->NewSlot();
        TEveWindowSlot* slotRightBottom = pack1->NewSlot();
        
        // Creating viewers and scenes
        TEveCalo3D* calo3d = Create3DView(data_t, slotLeftTop);
        CreateHistoLego(data_t, slotLeftBottom);
        CreateProjections(data_t, calo3d, slotRightTop, slotRightBottom);
        
        gEve->Redraw3D();
        
        Info("AliEveMomentumHistograms::DrawAllEvents", "...Finished");
    }
    
    return data_t;
}

TEveCaloLego* AliEveMomentumHistograms::CreateHistoLego(TEveCaloData* data, TEveWindowSlot* slot)
{
    TEveCaloLego* lego;
    
    // Viewer initialization, tab creation
    if ( g_histo2d_all_events_v0 == 0 ) {
        
        TEveBrowser *browser = gEve->GetBrowser();
        slot->MakeCurrent();
        g_histo2d_all_events_v0 = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
        g_histo2d_all_events_s0 = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
        g_histo2d_all_events_v0->AddScene(g_histo2d_all_events_s0);
        g_histo2d_all_events_v0->SetElementName("2D Lego Viewer");
        g_histo2d_all_events_s0->SetElementName("2D Lego Scene");
        
        TGLViewer* glv = g_histo2d_all_events_v0->GetGLViewer();
        g_histo2d_all_events_lego_overlay = new TEveCaloLegoOverlay();
        glv->AddOverlayElement(g_histo2d_all_events_lego_overlay);
        glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
        
        // Plotting histogram lego
        lego = new TEveCaloLego(data);
        g_histo2d_all_events_s0->AddElement(lego);
        
        // Move to real world coordinates
        lego->InitMainTrans();
        Float_t sc = TMath::Min(lego->GetEtaRng(), lego->GetPhiRng());
        lego->RefMainTrans().SetScale(sc, sc, sc);
        
        // Set event handler to move from perspective to orthographic view.
        glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));
        
        g_histo2d_all_events_lego_overlay->SetCaloLego(lego);
    }
    
    return lego;
}

TEveCalo3D* AliEveMomentumHistograms::Create3DView(TEveCaloData* data, TEveWindowSlot* slot)
{
    TEveCalo3D* calo3d;
    
    if ( g_histo2d_all_events_v1 == 0 ) {
        
        TEveBrowser *browser = gEve->GetBrowser();
        slot->MakeCurrent();
        g_histo2d_all_events_v1 = gEve->SpawnNewViewer("3D Histogram", "3D Histogram");
        g_histo2d_all_events_s1 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
        g_histo2d_all_events_v1->AddScene(g_histo2d_all_events_s1);
        g_histo2d_all_events_v1->SetElementName("3D Histogram Viewer");
        g_histo2d_all_events_s1->SetElementName("3D Histogram Scene");
        
        calo3d = new TEveCalo3D(data);
        
        calo3d->SetBarrelRadius(550);
        calo3d->SetEndCapPos(550);
        g_histo2d_all_events_s1->AddElement(calo3d);
    }
    
    return calo3d;
}

void AliEveMomentumHistograms::CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d, TEveWindowSlot* slot1, TEveWindowSlot* slot2)
{
    
    if ( g_histo2d_all_events_v2 == 0 ) {
        
        TEveBrowser *browser = gEve->GetBrowser();
        slot1->MakeCurrent();
        g_histo2d_all_events_v2 = gEve->SpawnNewViewer("RPhi projection", "RPhi projection");
        g_histo2d_all_events_s2 = gEve->SpawnNewScene("RPhi projection", "RPhi projection");
        g_histo2d_all_events_v2->AddScene(g_histo2d_all_events_s2);
        g_histo2d_all_events_v2->SetElementName("RPhi Projection Viewer");
        g_histo2d_all_events_s2->SetElementName("RPhi Projection Scene");
        
        TEveProjectionManager* mng1 = new TEveProjectionManager();
        mng1->SetProjection(TEveProjection::kPT_RPhi);
        
        TEveProjectionAxes* axeg_histo2d_all_events_s1 = new TEveProjectionAxes(mng1);
        g_histo2d_all_events_s2->AddElement(axeg_histo2d_all_events_s1);
        TEveCalo2D* calo2d1 = (TEveCalo2D*) mng1->ImportElements(calo3d);
        g_histo2d_all_events_s2->AddElement(calo2d1);
        
        g_histo2d_all_events_v2->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    }
    
    if ( g_histo2d_all_events_v3 == 0 ) {
        
        TEveBrowser *browser = gEve->GetBrowser();
        slot2->MakeCurrent();
        g_histo2d_all_events_v3 = gEve->SpawnNewViewer("RhoZ projection", "RhoZ projection");
        g_histo2d_all_events_s3 = gEve->SpawnNewScene("RhoZ projection", "RhoZ projection");
        g_histo2d_all_events_v3->AddScene(g_histo2d_all_events_s3);
        g_histo2d_all_events_v3->SetElementName("RhoZ Projection Viewer");
        g_histo2d_all_events_s3->SetElementName("RhoZ Projection Viewer");
        
        TEveProjectionManager* mng2 = new TEveProjectionManager();
        mng2->SetProjection(TEveProjection::kPT_RhoZ);
        
        TEveProjectionAxes* axeg_histo2d_all_events_s2 = new TEveProjectionAxes(mng2);
        g_histo2d_all_events_s3->AddElement(axeg_histo2d_all_events_s2);
        TEveCalo2D* calo2d2 = (TEveCalo2D*) mng2->ImportElements(calo3d);
        g_histo2d_all_events_s3->AddElement(calo2d2);
        
        g_histo2d_all_events_v3->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    } 
    
    return;
}



