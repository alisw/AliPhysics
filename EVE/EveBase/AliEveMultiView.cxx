// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMultiView.h"
#include "AliEveInit.h"
#include <TGPack.h>
#include <TBrowser.h>

#include <iostream>

using namespace std;

//______________________________________________________________________________
// Full description of AliEveMultiView
//

ClassImp(AliEveMultiView)

AliEveMultiView* AliEveMultiView::fgInstance = 0;

AliEveMultiView* AliEveMultiView::Instance()
{
    // Return static instance.
    
    return fgInstance;
}

AliEveMultiView::AliEveMultiView() :
fRPhiMgr(0), fRhoZMgr(0),
f3DView(0), fRPhiView(0), fRhoZView(0),
fRPhiGeomScene(0), fRhoZGeomScene(0),
fRPhiEventScene(0), fRhoZEventScene(0),
fGeomGentle(0), fGeomGentleRPhi(0), fGeomGentleRhoZ(0),
fGeomGentleTrd(0),
fPack(0)
{
    // Constructor --- creates required scenes, projection managers
    // and GL viewers.
    
    if (fgInstance)
        throw TEveException("AliEveMultiView::AliEveMultiView already instantiated.");
    fgInstance = this;
    
    // Scenes
    //========
    
    fRPhiGeomScene  = gEve->SpawnNewScene("RPhi Geometry",
                                          "Scene holding projected geometry for the RPhi view.");
    fRhoZGeomScene  = gEve->SpawnNewScene("RhoZ Geometry",
                                          "Scene holding projected geometry for the RhoZ view.");
    fRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data",
                                          "Scene holding projected event-data for the RPhi view.");
    fRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
                                          "Scene holding projected event-data for the RhoZ view.");
    
    // Projection managers
    //=====================
    
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
    //=========
    
    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;
    
    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack(); // slot is destroyed here
    pack->SetElementName("Multi View");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);
    
    pack->NewSlotWithWeight(2)->MakeCurrent(); // new slot is created from pack
    f3DView = gEve->SpawnNewViewer("3D View MV", "");
    f3DView->AddScene(gEve->GetGlobalScene());
    f3DView->AddScene(gEve->GetEventScene());
    
    pack = pack->NewSlot()->MakePack(); // new slot created from pack, then slot is destroyed and new pack returned
    pack->SetShowTitleBar(kFALSE);
    pack->NewSlot()->MakeCurrent(); // new slot from pack
    fRPhiView = gEve->SpawnNewViewer("RPhi View", "");
    fRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    fRPhiView->AddScene(fRPhiGeomScene);
    fRPhiView->AddScene(fRPhiEventScene);
    
    fPack = pack;
    
    pack->NewSlot()->MakeCurrent(); // new slot from pack
    fRhoZView = gEve->SpawnNewViewer("RhoZ View", "");
    fRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    fRhoZView->AddScene(fRhoZGeomScene);
    fRhoZView->AddScene(fRhoZEventScene);
}

AliEveMultiView::~AliEveMultiView()
{
    DestroyAllGeometries();
    
    delete fGeomGentle;
    delete fGeomGentleRPhi;
    delete fGeomGentleRhoZ;
    
    delete	fRPhiMgr;
    delete fRhoZMgr;
}

//-------------------------------------------------------------------------

void AliEveMultiView::InitGeomGentle(TEveGeoShape* g3d, TEveGeoShape* grphi, TEveGeoShape* grhoz)
{
    // Initialize gentle geometry.
    
    fGeomGentle     = g3d;
    fGeomGentleRPhi = grphi; fGeomGentleRPhi->IncDenyDestroy();
    fGeomGentleRhoZ = grhoz; fGeomGentleRhoZ->IncDenyDestroy();
    
    ImportGeomRPhi(fGeomGentleRPhi);
    ImportGeomRhoZ(fGeomGentleRhoZ);
}

void AliEveMultiView::InitGeomGentleTrd(TEveGeoShape* gtrd)
{
    // Initialize gentle geometry TRD.
    
    fGeomGentleTrd = gtrd;
    ImportGeomRPhi(fGeomGentleTrd);
    ImportGeomRhoZ(fGeomGentleTrd);
}

void AliEveMultiView::InitSimpleGeom(TEveGeoShape* geom, bool rPhi, bool rhoZ)
{
    fGeomVector.push_back(geom);
    if(rPhi) ImportGeomRPhi(geom);
    if(rhoZ) ImportGeomRhoZ(geom);
}

//-------------------------------------------------------------------------

void AliEveMultiView::SetDepth(Float_t d)
{
    // Set current depth on all projection managers.
    
    fRPhiMgr->SetCurrentDepth(d);
    fRhoZMgr->SetCurrentDepth(d);
}

//-------------------------------------------------------------------------

void AliEveMultiView::ImportGeomRPhi(TEveElement* el)
{
    // Import el into r-phi geometry scene.
    
    fRPhiMgr->ImportElements(el, fRPhiGeomScene);
}

void AliEveMultiView::ImportGeomRhoZ(TEveElement* el)
{
    // Import el into rho-z geometry scene.
    
    fRhoZMgr->ImportElements(el, fRhoZGeomScene);
}

void AliEveMultiView::ImportEventRPhi(TEveElement* el)
{
    // Import el into r-phi event scene.
    
    fRPhiMgr->ImportElements(el, fRPhiEventScene);
}

void AliEveMultiView::ImportEventRhoZ(TEveElement* el)
{
    // Import el into rho-z event scene.
    
    fRhoZMgr->ImportElements(el, fRhoZEventScene);
}

void AliEveMultiView::DestroyEventRPhi()
{
    // Destroy all elements in r-phi event scene.
    
    fRPhiEventScene->DestroyElements();
}

void AliEveMultiView::DestroyEventRhoZ()
{
    // Destroy all elements in rho-z event scene.
    
    fRhoZEventScene->DestroyElements();
}


//-------------------------------------------------------------------------

void AliEveMultiView::SetCenterRPhi(Double_t x, Double_t y, Double_t z)
{
    // Set center of r-phi manager.
    
    fRPhiMgr->SetCenter(x, y, z);
}

void AliEveMultiView::SetCenterRhoZ(Double_t x, Double_t y, Double_t z)
{
    // Set center of rho-z manager.
    
    fRhoZMgr->SetCenter(x, y, z);
}

void AliEveMultiView::DestroyAllGeometries()
{
    // Destroy 3d, r-phi and rho-z geometries.
    if(fGeomGentle){
        fGeomGentle->DestroyElements();
        gEve->RemoveElement(fGeomGentle,gEve->GetGlobalScene());
        fGeomGentle = 0;
    }
    if(fGeomGentleRPhi){
        fGeomGentleRPhi->DestroyElements();
        gEve->RemoveElement(fGeomGentleRPhi,gEve->GetGlobalScene());
        fGeomGentleRPhi = 0;
    }
    if(fGeomGentleRhoZ){
        fGeomGentleRhoZ->DestroyElements();
        gEve->RemoveElement(fGeomGentleRhoZ,gEve->GetGlobalScene());
        fGeomGentleRhoZ = 0;
    }
    if(fGeomGentleTrd){
        fGeomGentleTrd->DestroyElements();
        gEve->RemoveElement(fGeomGentleTrd,gEve->GetGlobalScene());
        fGeomGentleTrd = 0;
    }
        
    for(int i=0;i<fGeomVector.size();i++)
    {
        if(fGeomVector[i])
        {
            fGeomVector[i]->DestroyElements();
            gEve->RemoveElement(fGeomVector[i],gEve->GetGlobalScene());
            fGeomVector[i] = 0;
        }
    }
}

