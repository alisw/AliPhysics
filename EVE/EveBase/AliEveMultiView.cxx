// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMultiView.h"
#include <TGPack.h>
#include <TBrowser.h>

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

AliEveMultiView::AliEveMultiView(Bool_t setMuonView) :
fRPhiMgr(0), fRhoZMgr(0), fMuonMgr(0),
f3DView(0), fRPhiView(0), fRhoZView(0), fMuonView(0),
fRPhiGeomScene(0), fRhoZGeomScene(0), fMuonGeomScene(0),
fRPhiEventScene(0), fRhoZEventScene(0), fMuonEventScene(0),
fGeomGentle(0), fGeomGentleRPhi(0), fGeomGentleRhoZ(0),
fGeomGentleTrd(0),fGeomGentleEmcal(0),fGeomGentleZdc(0), fGeomGentleMuon(0), fIsMuonView(kFALSE),
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
    fMuonGeomScene  = gEve->SpawnNewScene("Muon Geometry",
                                          "Scene holding projected geometry for the Muon view.");
    fRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data",
                                          "Scene holding projected event-data for the RPhi view.");
    fRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
                                          "Scene holding projected event-data for the RhoZ view.");
    fMuonEventScene = gEve->SpawnNewScene("Muon Event Data",
                                          "Scene holding projected event-data for the Muon view.");
    
    fIsMuonView = setMuonView;
    
    // Projection managers
    //=====================
    
    fRPhiMgr = new TEveProjectionManager();
    fRPhiMgr->SetProjection(TEveProjection::kPT_RPhi);
    gEve->AddToListTree(fRPhiMgr, kFALSE);
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
    
    if(fIsMuonView)
    {
        fMuonMgr = new TEveProjectionManager();
        fMuonMgr->SetProjection(TEveProjection::kPT_RhoZ);
        gEve->AddToListTree(fMuonMgr, kFALSE);
        {
            TEveProjectionAxes* a = new TEveProjectionAxes(fMuonMgr);
            a->SetMainColor(kWhite);
            a->SetTitle("Rho-Z Muon");
            a->SetTitleSize(0.05);
            a->SetTitleFont(102);
            a->SetLabelSize(0.025);
            a->SetLabelFont(102);
            fMuonGeomScene->AddElement(a);
        }
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
    f3DView = gEve->SpawnNewViewer("3D View", "");
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
    
    if(fIsMuonView)
    {
        pack->NewSlot()->MakeCurrent(); // new slot from pack
        fMuonView = gEve->SpawnNewViewer("RhoZ View Muon", "");
        fMuonView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
        fMuonView->AddScene(fMuonGeomScene);
        fMuonView->AddScene(fMuonEventScene);
    }
}

AliEveMultiView::~AliEveMultiView()
{
    DestroyAllGeometries();
    
    delete fGeomGentle;
    delete fGeomGentleRPhi;
    delete fGeomGentleRhoZ;
    
    delete	fRPhiMgr;
    delete fRhoZMgr;
    delete fMuonMgr;
    
}

//-------------------------------------------------------------------------

void AliEveMultiView::InitGeomGentle(TEveGeoShape* g3d, TEveGeoShape* grphi, TEveGeoShape* grhoz, TEveGeoShape* gmuon)
{
    // Initialize gentle geometry.
    
    fGeomGentle     = g3d;
    fGeomGentleRPhi = grphi; fGeomGentleRPhi->IncDenyDestroy();
    fGeomGentleRhoZ = grhoz; fGeomGentleRhoZ->IncDenyDestroy();
    if(fIsMuonView) { fGeomGentleMuon = gmuon; fGeomGentleMuon->IncDenyDestroy(); }
    
    ImportGeomRPhi(fGeomGentleRPhi);
    ImportGeomRhoZ(fGeomGentleRhoZ);
    if(fIsMuonView) ImportGeomMuon(fGeomGentleMuon);
}

void AliEveMultiView::InitGeomGentleTrd(TEveGeoShape* gtrd)
{
    // Initialize gentle geometry TRD.
    
    fGeomGentleTrd = gtrd;
    ImportGeomRPhi(fGeomGentleTrd);
    ImportGeomRhoZ(fGeomGentleTrd);
    if(fIsMuonView) ImportGeomMuon(fGeomGentleTrd);
}

void AliEveMultiView::InitGeomGentleEmcal(TEveGeoShape* gemcal)
{
    // Initialize gentle geometry EMCal.
    
    fGeomGentleEmcal = gemcal;
    ImportGeomRPhi(fGeomGentleEmcal);
    ImportGeomRhoZ(fGeomGentleEmcal);
    if(fIsMuonView) ImportGeomMuon(fGeomGentleEmcal);
}

void AliEveMultiView::InitGeomGentleZdc(TEveGeoShape* gzdc)
{
    // Initialize gentle geometry ZDC.
    
    fGeomGentleZdc = gzdc;
    ImportGeomRPhi(fGeomGentleZdc);
    ImportGeomRhoZ(fGeomGentleZdc);
    if(fIsMuonView) ImportGeomMuon(fGeomGentleZdc);
}

void AliEveMultiView::InitGeomGentleMuon(TEveGeoShape* gmuon, Bool_t showRPhi, Bool_t showRhoZ, Bool_t showMuon)
{
    // Initialize gentle geometry for MUON.
    
    fGeomGentleMuon = gmuon;
    if (showRPhi) ImportGeomRPhi(fGeomGentleMuon);
    if (showRhoZ) ImportGeomRhoZ(fGeomGentleMuon);
    if (showMuon && fIsMuonView) ImportGeomMuon(fGeomGentleMuon);
    
}

//-------------------------------------------------------------------------

void AliEveMultiView::SetDepth(Float_t d)
{
    // Set current depth on all projection managers.
    
    fRPhiMgr->SetCurrentDepth(d);
    fRhoZMgr->SetCurrentDepth(d);
    if(fIsMuonView) fMuonMgr->SetCurrentDepth(d);
    
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

void AliEveMultiView::ImportGeomMuon(TEveElement* el)
{
    // Import el into muon geometry scene.
    
    if(fIsMuonView) fMuonMgr->ImportElements(el, fMuonGeomScene);
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

void AliEveMultiView::ImportEventMuon(TEveElement* el)
{
    // Import el into muon event scene.
    
    if(fIsMuonView) fMuonMgr->ImportElements(el, fMuonEventScene);
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

void AliEveMultiView::DestroyEventMuon()
{
    // Destroy all elements in rho-z event scene.
    
    if(fIsMuonView) fMuonEventScene->DestroyElements();
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

void AliEveMultiView::SetCenterMuon(Double_t x, Double_t y, Double_t z)
{
    // Set center of rho-z manager.
    
    if(fIsMuonView) fMuonMgr->SetCenter(x, y, z);
}

void AliEveMultiView::DestroyAllGeometries()
{
    // Destroy 3d, r-phi and rho-z geometries.
    
    fGeomGentle->DestroyElements();
    gEve->RemoveElement(fGeomGentle,gEve->GetGlobalScene());
    
    fGeomGentleRPhi->DestroyElements();
    gEve->RemoveElement(fGeomGentleRPhi,gEve->GetGlobalScene());
    
    fGeomGentleRhoZ->DestroyElements();
    gEve->RemoveElement(fGeomGentleRhoZ,gEve->GetGlobalScene());
    
    fGeomGentleMuon->DestroyElements();
    gEve->RemoveElement(fGeomGentleMuon,gEve->GetGlobalScene());
    
    fGeomGentleTrd->DestroyElements();
    gEve->RemoveElement(fGeomGentleTrd,gEve->GetGlobalScene());
    
    fGeomGentleEmcal->DestroyElements();
    gEve->RemoveElement(fGeomGentleEmcal,gEve->GetGlobalScene());

    if(fGeomGentleZdc){
        fGeomGentleZdc->DestroyElements();
        gEve->RemoveElement(fGeomGentleZdc,gEve->GetGlobalScene());
    }
//    if(fIsMuonView) fGeomGentleMuon->DestroyElements();
    
}

