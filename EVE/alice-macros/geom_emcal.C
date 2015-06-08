// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>

#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

void geom_emcal()
{
    AliEveEventManager::GetMaster()->AssertGeometry();
    AliEveMultiView *mv = AliEveMultiView::Instance();
    
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
    if (!node) {
        Warning("geom_emcal()", "Node XEN1_1 not found.");
        return;
    }
    
    TEveGeoTopNode* emcal_re = new TEveGeoTopNode(gGeoManager, node);
    emcal_re->SetVisLevel(1);
    emcal_re->SetMainTransparency(70);
    
    // 3D view:
    gEve->AddGlobalElement(emcal_re);
    
    // RPhi
    mv->ImportGeomRPhi(emcal_re);
    emcal_re->ProjectAllChildren();
    emcal_re->PropagateRnrStateToProjecteds();
    g_proj = gEve->SpawnNewScene("emcal proj", "emcal proj");
    AliEveMultiView::Instance()->GetRPhiView()->AddScene(g_proj);
    g_proj->SetElementName("EMCal geom RPhi");
    g_proj->AddElement(emcal_re);
    
    gSystem->ProcessEvents();
    gEve->Redraw3D();
}
