// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMultiView.h"
#include "AliEveInit.h"

#include <TEveProjectionAxes.h>
#include <TEveBrowser.h>
#include <TBrowser.h>

#include <iostream>

using namespace std;

ClassImp(AliEveMultiView)

AliEveMultiView* AliEveMultiView::fgInstance = 0;

AliEveMultiView::AliEveMultiView() :
fRPhiMgr(0), fRhoZMgr(0),
f3DView(0), fRPhiView(0), fRhoZView(0),
f3DGeomScene(0), fRPhiGeomScene(0), fRhoZGeomScene(0),
f3DEventScene(0),fRPhiEventScene(0), fRhoZEventScene(0)
{
    // Constructor --- creates required scenes, projection managers and GL viewers
    
    if (fgInstance){
        throw TEveException("AliEveMultiView::AliEveMultiView already instantiated.");
    }
    fgInstance = this;
    
    // Scenes
    f3DGeomScene  = gEve->SpawnNewScene("3D Geometry","Scene holding 3D geometry.");
    fRPhiGeomScene  = gEve->SpawnNewScene("RPhi Geometry",
                                          "Scene holding projected geometry for the RPhi view.");
    fRhoZGeomScene  = gEve->SpawnNewScene("RhoZ Geometry",
                                          "Scene holding projected geometry for the RhoZ view.");
    
    f3DEventScene = gEve->SpawnNewScene("3D Event Data","Scene holding 3D event-data.");
    fRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data",
                                          "Scene holding projected event-data for the RPhi view.");
    fRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
                                          "Scene holding projected event-data for the RhoZ view.");
    
    // Projection managers
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    bool showAxes = settings.GetValue("axes.show", false);
    
    fRPhiMgr = new TEveProjectionManager();
    fRPhiMgr->SetProjection(TEveProjection::kPT_RPhi);
    gEve->AddToListTree(fRPhiMgr, kFALSE);
    
    if(showAxes)
    {
        TEveProjectionAxes* a = new TEveProjectionAxes(fRPhiMgr);
        a->SetMainColor(kWhite);
        a->SetTitle("R-Phi");
        a->SetTitleSize(0.05);
        a->SetTitleFont(102);
        a->SetLabelSize(0.025);
        a->SetLabelFont(102);
        fRPhiGeomScene->AddElement(a);
    }
    
    fRhoZMgr = new TEveProjectionManager();
    fRhoZMgr->SetProjection(TEveProjection::kPT_RhoZ);
    gEve->AddToListTree(fRhoZMgr, kFALSE);
    
    if(showAxes)
    {
        TEveProjectionAxes* a = new TEveProjectionAxes(fRhoZMgr);
        a->SetMainColor(kWhite);
        a->SetTitle("Rho-Z");
        a->SetTitleSize(0.05);
        a->SetTitleFont(102);
        a->SetLabelSize(0.025);
        a->SetLabelFont(102);
        fRhoZGeomScene->AddElement(a);
    }
    
    // Viewers
    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;
    
    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack(); // slot is destroyed here
    pack->SetElementName("Multi View");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);
    
    pack->NewSlotWithWeight(2)->MakeCurrent(); // new slot is created from pack
    f3DView = gEve->SpawnNewViewer("3D View MV", "");
    f3DView->AddScene(f3DGeomScene);
    f3DView->AddScene(f3DEventScene);
    
    pack = pack->NewSlot()->MakePack(); // new slot created from pack, then slot is destroyed and new pack returned
    pack->SetShowTitleBar(kFALSE);
    pack->NewSlot()->MakeCurrent(); // new slot from pack
    fRPhiView = gEve->SpawnNewViewer("RPhi View", "");
    fRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    fRPhiView->AddScene(fRPhiGeomScene);
    fRPhiView->AddScene(fRPhiEventScene);
    
    pack->NewSlot()->MakeCurrent(); // new slot from pack
    fRhoZView = gEve->SpawnNewViewer("RhoZ View", "");
    fRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    fRhoZView->AddScene(fRhoZGeomScene);
    fRhoZView->AddScene(fRhoZEventScene);
}

AliEveMultiView::~AliEveMultiView()
{
    DestroyAllGeometries();
    delete fRPhiMgr;
    delete fRhoZMgr;
}

void AliEveMultiView::InitSimpleGeom(TEveGeoShape* geom, bool threeD, bool rPhi, bool rhoZ)
{
    if(!geom)
    {
        cout<<"AliEveMultiView::InitSimpleGeom -- geometry is NULL!"<<endl;
        return;
    }
    
    fGeomVector.push_back(geom);
    
    if(threeD){
        gEve->AddElement(geom,f3DGeomScene);
    }
    if(rPhi){
        fRPhiMgr->SetCurrentDepth(-10);
        fRPhiMgr->ImportElements(geom, fRPhiGeomScene);
        fRPhiMgr->SetCurrentDepth(0);
    }
    if(rhoZ){
        fRhoZMgr->SetCurrentDepth(-10);
        fRhoZMgr->ImportElements(geom, fRhoZGeomScene);
        fRhoZMgr->SetCurrentDepth(0);
    }
}

void AliEveMultiView::ImportEvent(TEveElement* el)
{
    gEve->AddElement(el,f3DEventScene);
    fRPhiMgr->ImportElements(el, fRPhiEventScene);
    fRhoZMgr->ImportElements(el, fRhoZEventScene);
}

void AliEveMultiView::ImportEvent3D(TEveElement* el)
{
    gEve->AddElement(el, f3DEventScene);
}

void AliEveMultiView::ImportEventRPhi(TEveElement* el)
{
    fRPhiMgr->ImportElements(el, fRPhiEventScene);
}

void AliEveMultiView::ImportEventRhoZ(TEveElement* el)
{
    fRhoZMgr->ImportElements(el, fRhoZEventScene);
}

void AliEveMultiView::DestroyEvent3D()
{
    f3DEventScene->DestroyElements();
}

void AliEveMultiView::DestroyEventRPhi()
{
    fRPhiEventScene->DestroyElements();
}

void AliEveMultiView::DestroyEventRhoZ()
{
    fRhoZEventScene->DestroyElements();
}

void AliEveMultiView::DestroyAllEvents()
{
    f3DEventScene->DestroyElements();
    fRPhiEventScene->DestroyElements();
    fRhoZEventScene->DestroyElements();
}

void AliEveMultiView::DestroyAllGeometries()
{
    for(int i=0;i<fGeomVector.size();i++)
    {
        if(fGeomVector[i])
        {
            fGeomVector[i]->DestroyElements();
            gEve->RemoveElement(fGeomVector[i],f3DGeomScene);
            fGeomVector[i] = 0;
        }
    }
}

