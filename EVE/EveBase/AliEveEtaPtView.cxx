// Author: Jeremi Niedziela 2015

#include "AliEveEtaPtView.h"
#include <TBrowser.h>
#include <TMath.h>
#include <TEveStraightLineSet.h>
#include <TRandom.h>
#include <TEveBoxSet.h>
#include <TGeoTube.h>
#include <TEveGeoShape.h>

using namespace TMath;

ClassImp(AliEveEtaPtView)

AliEveEtaPtView* AliEveEtaPtView::fInstance = 0;

AliEveEtaPtView* AliEveEtaPtView::Instance()
{
    return fInstance;
}

AliEveEtaPtView::AliEveEtaPtView() :
    f3DView(0),
    fGeom(0)
{
    if (fInstance){
        throw TEveException("AliEveEtaPtView::AliEveEtaPtView already instantiated.");
    }
    fInstance = this;

//    f3DView = gEve->SpawnNewViewer("3D Eta-Pt View", "");
//    f3DView->AddScene(gEve->GetGlobalScene());
//    f3DView->AddScene(gEve->GetEventScene());
    
    f3DGeomScene = new TEveScene("3D Eta-pt scene");
    fRPhiGeomScene = new TEveScene("RPhi Eta-pt scene");
    fRhoZGeomScene = new TEveScene("PhoZ Eta-pt scene");
    
    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;
    
    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack(); // slot is destroyed here
    pack->SetElementName("Eta-pt view");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);
    
    pack->NewSlotWithWeight(2)->MakeCurrent(); // new slot is created from pack
    f3DView = gEve->SpawnNewViewer("3D Eta-pt View", "");
    f3DView->AddScene(f3DGeomScene);
    f3DView->AddScene(gEve->GetEventScene());
    
    pack = pack->NewSlot()->MakePack(); // new slot created from pack, then slot is destroyed and new pack returned
    pack->SetShowTitleBar(kFALSE);
    pack->NewSlot()->MakeCurrent(); // new slot from pack
    fRPhiView = gEve->SpawnNewViewer("RPhi Eta-pt View", "");
    fRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    fRPhiView->AddScene(fRPhiGeomScene);
//    fRPhiView->AddScene(fRPhiEventScene);
    
    fPack = pack;
    
    pack->NewSlot()->MakeCurrent(); // new slot from pack
    fRhoZView = gEve->SpawnNewViewer("RhoZ Eta-pt View", "");
    fRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
    fRhoZView->AddScene(fRhoZGeomScene);
//    fRhoZView->AddScene(fRhoZEventScene);

}

AliEveEtaPtView::~AliEveEtaPtView()
{
    fGeom->DestroyElements();
    delete fGeom;
}

void AliEveEtaPtView::InitGeom()
{
    Float_t x=0, y=0, z=0;

    TGeoTube *tube = new TGeoTube(0,500,100);
    TGeoTube *tube2 = new TGeoTube(0,400,100);
    TGeoTube *tube3 = new TGeoTube(0,400,100);
    TGeoTube *tube4 = new TGeoTube(0,300,100);
    TGeoTube *tube5 = new TGeoTube(0,300,100);

    TEveGeoShape *eveTube = new TEveGeoShape();
    TEveGeoShape *eveTube2 = new TEveGeoShape();
    TEveGeoShape *eveTube3 = new TEveGeoShape();
    TEveGeoShape *eveTube4 = new TEveGeoShape();
    TEveGeoShape *eveTube5 = new TEveGeoShape();
    
    eveTube->SetShape(tube);
    eveTube2->SetShape(tube2);
    eveTube3->SetShape(tube3);
    eveTube4->SetShape(tube4);
    eveTube5->SetShape(tube5);
    
    eveTube->SetMainColor(kGreen);
    eveTube2->SetMainColor(kRed);
    eveTube3->SetMainColor(kRed);
    eveTube4->SetMainColor(kYellow);
    eveTube5->SetMainColor(kYellow);
    
    eveTube->SetMainTransparency(70);
    eveTube2->SetMainTransparency(70);
    eveTube3->SetMainTransparency(70);
    eveTube4->SetMainTransparency(70);
    eveTube5->SetMainTransparency(70);
    
    
    TEveTrans& t = eveTube->RefMainTrans();
    TEveTrans& t2 = eveTube2->RefMainTrans();
    TEveTrans& t3 = eveTube3->RefMainTrans();
    TEveTrans& t4 = eveTube4->RefMainTrans();
    TEveTrans& t5 = eveTube5->RefMainTrans();

    t.SetPos(x, y, z);
    t2.SetPos(x, y, z-200);
    t3.SetPos(x, y, z+200);
    t4.SetPos(x, y, z-400);
    t5.SetPos(x, y, z+400);

    
    f3DGeomScene->AddElement(eveTube);
    f3DGeomScene->AddElement(eveTube2);
    f3DGeomScene->AddElement(eveTube3);
    f3DGeomScene->AddElement(eveTube4);
    f3DGeomScene->AddElement(eveTube5);

    ImportShapeToRPhi(eveTube);
    ImportShapeToRPhi(eveTube2);
    ImportShapeToRPhi(eveTube3);
    ImportShapeToRPhi(eveTube4);
    ImportShapeToRPhi(eveTube5);
    
    ImportShapeToRhoZ(eveTube);
    ImportShapeToRhoZ(eveTube2);
    ImportShapeToRhoZ(eveTube3);
    ImportShapeToRhoZ(eveTube4);
    ImportShapeToRhoZ(eveTube5);
    
    f3DView->GetGLViewer()->UpdateScene();
    gEve->Redraw3D(kTRUE);
    
}

void AliEveEtaPtView::ImportShapeToRhoZ(TEveGeoShape *shape)
{
    fRhoZGeomScene->AddElement(shape);
    fRhoZView->GetGLViewer()->UpdateScene();
}

void AliEveEtaPtView::ImportShapeToRPhi(TEveGeoShape *shape)
{
    fRPhiGeomScene->AddElement(shape);
    fRPhiView->GetGLViewer()->UpdateScene();
}

