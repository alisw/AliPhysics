//
//  AliEveMomentumHistograms.h
//
//  Created by Jeremi Niedziela on 1/12/15
//  main author: Stefano Carrazza 2010
//

#ifndef __AliEveMomentumHistograms__
#define __AliEveMomentumHistograms__

#include <AliEveMultiView.h>

#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveCaloLegoOverlay.h>
#include <TEveWindow.h>
#include <TEveCalo.h>

class AliEveMomentumHistograms
{
public:
    AliEveMomentumHistograms();
    ~AliEveMomentumHistograms();
    
    TEveCaloDataHist* Draw();
    TEveCaloDataHist* DrawAllEvents();
    
private:
    double pi;
    TEveViewer *fHistoViewer;
    TEveScene  *g_histo2d_s;
    TEveScene  *g_histo2d_s2;
    TEveCaloLegoOverlay* g_histo2d_lego_overlay;
    
    TEveViewer *g_histo2d_all_events_v0;
    TEveViewer *g_histo2d_all_events_v1;
    TEveViewer *g_histo2d_all_events_v2;
    TEveViewer *g_histo2d_all_events_v3;
    TEveScene  *g_histo2d_all_events_s0;
    TEveScene  *g_histo2d_all_events_s1;
    TEveScene  *g_histo2d_all_events_s2;
    TEveScene  *g_histo2d_all_events_s3;
    TEveCaloLegoOverlay* g_histo2d_all_events_lego_overlay;
    TEveWindowSlot* g_histo2d_all_events_slot;
    
    Double_t GetPhi(Double_t phi);
    void CreateHistoLego(TEveCaloData* data);
    TEveCaloLego* CreateHistoLego(TEveCaloData* data, TEveWindowSlot* slot);
    TEveCalo3D* Create3DView(TEveCaloData* data);
    TEveCalo3D* Create3DView(TEveCaloData* data, TEveWindowSlot* slot);
    void CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d, TEveWindowSlot* slot1, TEveWindowSlot* slot2);
    
    AliEveMomentumHistograms(const AliEveMomentumHistograms&);
    AliEveMomentumHistograms& operator=(const AliEveMomentumHistograms&);
};

#endif